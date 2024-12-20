//!!-----------------------------------------------------------------------
//!!     Copyright (C) 2019-2024 CEREA (ENPC) - INERIS
//!!     SSH-aerosol is distributed under the GNU General Public License v3
//!!-----------------------------------------------------------------------

using namespace ssh_soap;

void compute_kp_org_bins_ssh(model_config &config, vector<species>& surrogate,
                             Array <double, 3> &MOinit, double &Temperature, Array<double, 3> &MOW, int b)
{
  //compute partitioning constant of the organic phase by taking into account activity
  //coefficients and the kelvin effect
  double temp1,temp2,maxi,MOWsurf;
  double kelvin_effect=1.e0;
  int ilayer,i,iphase,jphase;
  int n=surrogate.size();

  if (config.fixed_density==false)
    {
      double tmp=0.;
      double tmp2=0.;
      ilayer=config.nlayer-1;
      for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
	for (i=0;i<n;++i)
	  if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophobic)
	    {
	      tmp+=surrogate[i].Ap_layer_init(b,ilayer,iphase)/surrogate[i].rho;
	      tmp2+=surrogate[i].Ap_layer_init(b,ilayer,iphase);
	    }

      if (tmp>0)
	config.rho_organic=tmp2/tmp;
      else
	config.rho_organic=1300.;
    }
  
  if (config.compute_kelvin_effect and config.diameters(b) > 0.0) //compute the kelvin_effect
    {
      //compute the mean molar mass of organic phases
      //it is used to have the same mean molar mass for all organic phases and to prevent
      //two organic phase from having different kelvin effect
      //(with different kelvin effect numerical problems could arised) 
      temp1=0.0;
      temp2=0.0;
      for (iphase=0;iphase<config.nphase(b,config.nlayer-1);++iphase)
        {
          temp1+=MOinit(b,config.nlayer-1,iphase);
          if(MOW(b,config.nlayer-1,iphase) <= 1.) 
            MOW(b,config.nlayer-1,iphase) = 200.0;
          temp2+=MOinit(b,config.nlayer-1,iphase)/MOW(b,config.nlayer-1,iphase);
        }
      if (temp1>0.0)
        MOWsurf=temp1/temp2;
      else
        MOWsurf=200.0;
      kelvin_effect = 2.0*config.surface_tension_org*
        MOWsurf/(8.314*Temperature*config.rho_organic*
                 0.5*config.diameters(b));
      if(kelvin_effect > 50.0)
        kelvin_effect = 50.0;
      kelvin_effect=exp(kelvin_effect);
    }
  for (i=0;i<n;++i)
    if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophobic)
      {		  
        //compute the characteristic time
        for (ilayer=0;ilayer<config.nlayer;++ilayer)
          {
            maxi=0.0;
            jphase=0;
            if (surrogate[i].kp_from_experiment)
              for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                if (surrogate[i].Ap_layer_init(b,ilayer,iphase)>0.0)
                  jphase=iphase;
		
            if (i!=config.iH2O)
             for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                {
                  if (surrogate[i].nonvolatile)
                    surrogate[i].Kp(b,ilayer,iphase)=1.e9*200.0/
                      (surrogate[i].gamma_org_layer(b,ilayer,iphase)*MOW(b,ilayer,iphase));
                  else if (surrogate[i].kp_from_experiment)
                    if(iphase==jphase)
                      surrogate[i].Kp(b,ilayer,iphase)=surrogate[i].kpi; //Kp_exp_org(Temperature);
                    else
                      surrogate[i].Kp(b,ilayer,iphase)=0.0;
                  else
                    {
                      if(MOW(b,ilayer,iphase) <= 1.) 
                        MOW(b,ilayer,iphase) = 200.0;
                      if(surrogate[i].gamma_org_layer(b,ilayer,iphase) <=1.e-20)   //KS Why 0 there?
                        surrogate[i].gamma_org_layer(b,ilayer,iphase) = 1.0;
                      surrogate[i].Kp(b,ilayer,iphase)=surrogate[i].kpi/MOW(b,ilayer,iphase)/
                        surrogate[i].gamma_org_layer(b,ilayer,iphase);
                    }		
                  surrogate[i].Kp(b,ilayer,iphase)/=kelvin_effect;                        
                }
            else
              for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                surrogate[i].Kp(b,ilayer,iphase)=
                  surrogate[i].kpi/MOW(b,ilayer,iphase)/
                  surrogate[i].gamma_org_layer(b,ilayer,iphase)/kelvin_effect;
	    
          }
      }
}

void compute_kp_org_ssh(model_config &config, vector<species>& surrogate,
			Array <double, 3> &MOinit, double &Temperature, Array<double, 3> &MOW)
{
  //compute partitioning constant of the organic phase by taking into account activity
  //coefficients and the kelvin effect
  int b;
  for (b=0;b<config.nbins;++b)
    compute_kp_org_bins_ssh(config, surrogate, MOinit, Temperature, MOW, b);
}

void compute_kp_aq_ssh(model_config &config, vector<species>& surrogate,
		       double &Temperature, Array <double, 1> &ionic,
		       Array <double, 1> &chp,Array<double, 1> &AQinit,Array<double, 1> &MMaq, int b1, int b2)
{
  //compute partitioning constant of the aqueous phase by taking into account activity
  //coefficients and the kelvin effect
  double kelvin_effect=1.;
  int b,i;
  int n=surrogate.size();
  double fion1,fion2;

  for (b=b1;b<b2;++b)      
    if (MMaq(b) > 0.0 and config.diameters(b)>0)
      {
	if (config.compute_kelvin_effect) //compute the kelvin effect
	  {
	    kelvin_effect=2.0*config.surface_tension_aq*MMaq(b)/
	      (8.314*Temperature*config.AQrho(b)*
	       0.5*config.diameters(b));
	    
	    if(kelvin_effect > 50.0)
	      kelvin_effect = 50.0;
	    kelvin_effect=exp(kelvin_effect);
	  }
	for (i=0;i<n;++i)
	  if(surrogate[i].is_organic and surrogate[i].hydrophilic and config.compute_organic)
	    {	  
	      surrogate[i].gamma_LR=surrogate[i].LR(b);
	      if (surrogate[i].nonvolatile)
		surrogate[i].Kaq(b)=1.0e9/surrogate[i].gamma_aq_bins(b)*18.0/MMaq(b);
	      else 
		{
		  surrogate[i].Kaq(b)=min(1.e20,
                                          surrogate[i].Kp_eff_aqrealdyn_ssh(config, Temperature, ionic(b), chp(b),surrogate[config.iHp].LR(b),
                                                                            surrogate[config.iHp].SRMR(b),MMaq(b),fion1,fion2,b))
		    /surrogate[i].gamma_aq_bins(b);
		
		  surrogate[i].fion1=fion1;
		  surrogate[i].fion2=fion2;
		}
	     
	      surrogate[i].Kaq(b)/=kelvin_effect;
	    }
	  else if (surrogate[i].is_inorganic_precursor and config.compute_inorganic and surrogate[i].is_solid==false)
	    {
	      if (i==config.iH2SO4)
		surrogate[i].Kaq(b)=1.0e10;          
	      else if (i==config.iNH3)
                {		
		  surrogate[i].Kaq(b)=surrogate[i].kaqi
		    *(1.0+surrogate[i].keqi*chp(b)*surrogate[config.iHp].gamma_aq_bins(b)/surrogate[config.iNH4p].gamma_aq_bins(b));
		  /*surrogate[i].dKaq(b)=
		    surrogate[i].Henry*exp(-deltaH_over_RT0*(T0/Temperature-1.0)-deltaCp0_over_R*(1+log(T0/Temperature)-T0/Temperature))
		    *R*Temperature/(1000.*1.0e6*1.013e5)
		    *(surrogate[i].Kequilibrium_ssh(Temperature)*surrogate[config.iHp].gamma_aq_bins(b)/surrogate[config.iNH4p].gamma_aq_bins(b));*/
		  if (config.compute_kelvin_effect) //compute the kelvin effect
		    {
		      surrogate[i].Kaq(b)/=kelvin_effect;
		      //surrogate[i].dKaq(b)/=kelvin_effect;
		    }
		}
	      else if (i==config.iHNO3)
		{
		  surrogate[i].Kaq(b)=surrogate[i].kaqi
		    *(1.0+surrogate[i].keqi/(chp(b)*surrogate[config.iHp].gamma_aq_bins(b)*surrogate[config.iNO3m].gamma_aq_bins(b)));
		  /*surrogate[i].dKaq(b)=-surrogate[i].Henry*exp(-deltaH_over_RT0*(T0/Temperature-1.0)-deltaCp0_over_R*(1+log(T0/Temperature)-T0/Temperature))
		    *R*Temperature/(1000.*1.0e6*1.013e5)
		    *(surrogate[i].Kequilibrium_ssh(Temperature)/(pow(chp(b),2)*surrogate[config.iHp].gamma_aq_bins(b)*surrogate[config.iClm].gamma_aq_bins(b)));*/
		  if (config.compute_kelvin_effect) //compute the kelvin effect
		    {
		      surrogate[i].Kaq(b)/=kelvin_effect;
		      //surrogate[i].dKaq(b)/=kelvin_effect;
		    }
		}
	      else if (i==config.iCO2)
		{
		  /*
		  surrogate[i].Kaq(b)=surrogate[i].kaqi
		  *(1.0+surrogate[i].keqi/(chp(b)*surrogate[config.iHp].gamma_aq_bins(b)*surrogate[config.iNO3m].gamma_aq_bins(b)));*/
		  surrogate[i].Kaq(b)=surrogate[config.iCO2].kaqi*surrogate[config.iH2O].Aaq_bins_init(b)/surrogate[config.iH2O].MM*MMaq(b)/AQinit(b)*surrogate[config.iH2O].gamma_aq_bins(b)*(surrogate[config.iCO2].keq/(chp(b)*surrogate[config.iHp].gamma_aq_bins(b)*surrogate[config.iHCO3m].gamma_aq_bins(b)))*
		    (1.+(surrogate[config.iCO2].keq2/(chp(b)*surrogate[config.iHp].gamma_aq_bins(b)*surrogate[config.iCO3mm].gamma_aq_bins(b)/surrogate[config.iHCO3m].gamma_aq_bins(b))));		  
		  //surrogate[i].Kaq(b)=surrogate[config.iCO2].kaqi*0.79*(surrogate[config.iCO2].keq/(chp(b)*surrogate[config.iHp].gamma_aq_bins(b)*surrogate[config.iHCO3m].gamma_aq_bins(b)))*
		  //  (1.+(surrogate[config.iCO2].keq2/(chp(b)*surrogate[config.iHp].gamma_aq_bins(b)*surrogate[config.iCO3mm].gamma_aq_bins(b)/surrogate[config.iHCO3m].gamma_aq_bins(b))));
		  //surrogate[i].Kaq(b)=surrogate[config.iCO2].kaqi*(surrogate[config.iCO2].keq/(chp(b)*surrogate[config.iHp].gamma_aq_bins(b)*surrogate[config.iHCO3m].gamma_aq_bins(b)))*
		  //  (1.+(surrogate[config.iCO2].keq2/(chp(b)*surrogate[config.iHp].gamma_aq_bins(b)*surrogate[config.iCO3mm].gamma_aq_bins(b)/surrogate[config.iHCO3m].gamma_aq_bins(b))));
		  //surrogate[i].Kaq(b)=surrogate[config.iCO2].kaqi*surrogate[config.iCO2].keq/chp(b)*(1.+surrogate[config.iCO2].keq2/chp(b));
		  //cout <<  "HCO2 " << surrogate[i].Kaq(b) << " " << surrogate[config.iCO2].kaqi << " " << surrogate[config.iCO2].keq << " " << surrogate[config.iCO2].keq2 << " " << chp(b) << endl;
		  /*surrogate[i].dKaq(b)=-surrogate[i].Henry*exp(-deltaH_over_RT0*(T0/Temperature-1.0)-deltaCp0_over_R*(1+log(T0/Temperature)-T0/Temperature))
		    *R*Temperature/(1000.*1.0e6*1.013e5)
		    *(surrogate[i].Kequilibrium_ssh(Temperature)/(pow(chp(b),2)*surrogate[config.iHp].gamma_aq_bins(b)*surrogate[config.iClm].gamma_aq_bins(b)));*/
		  if (config.compute_kelvin_effect) //compute the kelvin effect
		    {
		      surrogate[i].Kaq(b)/=kelvin_effect;
		      //surrogate[i].dKaq(b)/=kelvin_effect;
		    }
		}
	      else if (i==config.iHCl)
		{
		  surrogate[i].Kaq(b)=surrogate[i].kaqi
		    *(1.0+surrogate[i].keqi/(chp(b)*surrogate[config.iHp].gamma_aq_bins(b)*surrogate[config.iClm].gamma_aq_bins(b)));
                  /*
		  surrogate[i].dKaq(b)=-surrogate[i].Henry*exp(-deltaH_over_RT0*(T0/Temperature-1.0)-deltaCp0_over_R*(1+log(T0/Temperature)-T0/Temperature))
		    *R*Temperature/(1000.*1.0e6*1.013e5)
		    *(surrogate[i].Kequilibrium_ssh(Temperature)/(pow(chp(b),2)*surrogate[config.iHp].gamma_aq_bins(b)*surrogate[config.iClm].gamma_aq_bins(b)));*/
		  if (config.compute_kelvin_effect) //compute the kelvin effect
		    {
		      surrogate[i].Kaq(b)/=kelvin_effect;
		      //surrogate[i].dKaq(b)/=kelvin_effect;
		    }
		}
	    }
	  else if (i==config.iH2O)	    	    
	    surrogate[i].Kaq(b)=surrogate[i].kaqi/(MMaq(b)*surrogate[i].gamma_aq_bins(b))/kelvin_effect;			    	    
      }
}



void hydratation_dyn_ssh(model_config &config, vector<species>& surrogate, double &RH,Array<double, 1> &AQinit)  
{
  int n=surrogate.size();
  int i,b,ilayer,iphase;  
  Array <int, 1> ifound;
  ifound.resize(n);
  ifound=0;
    
  for (i=0;i<n;i++)    
    {
      int j=surrogate[i].iHyd;
      if (j>-1 and ifound(i)==0)
	{
	  ifound(i)=1;
	  ifound(j)=1;	 
	  int k=surrogate[j].iHyd;          
	  if (k>-1)
	    {
	      //cout << "before: " << surrogate[i].kprod_gas/surrogate[i].kloss_gas << endl;
	      surrogate[i].kprod_gas=0.;
	      surrogate[j].kprod_gas=0.;
	      surrogate[k].kprod_gas=0.;
	      if (surrogate[i].hydrophilic)
		for (b=0;b<config.nbins;++b)
		  {
		    double Keq1=RH*surrogate[i].Khyd*surrogate[i].gamma_aq_bins(b)/surrogate[j].gamma_aq_bins(b); //*surrogate[j].MM/surrogate[i].MM;
		    double Keq2=RH*surrogate[j].Khyd*surrogate[j].gamma_aq_bins(b)/surrogate[k].gamma_aq_bins(b); //*surrogate[k].MM/surrogate[j].MM;
		    //cout << Keq1 << " " << Keq2 << endl;
		    
		    double atot=surrogate[i].Aaq_bins_init(b)/surrogate[i].MM+surrogate[j].Aaq_bins_init(b)/surrogate[j].MM+surrogate[k].Aaq_bins_init(b)/surrogate[k].MM;
		    surrogate[i].Aaq_bins_init(b)=atot*surrogate[i].MM/(1+Keq1*(1+Keq2));
		    surrogate[j].Aaq_bins_init(b)=atot*surrogate[j].MM/(1+Keq1*(1+Keq2))*Keq1;
		    surrogate[k].Aaq_bins_init(b)=atot*surrogate[k].MM/(1+Keq1*(1+Keq2))*Keq1*Keq2;
		    
		    atot=surrogate[i].Aaq_bins_init0(b)/surrogate[i].MM+surrogate[j].Aaq_bins_init0(b)/surrogate[j].MM+surrogate[k].Aaq_bins_init0(b)/surrogate[k].MM;
		    surrogate[i].Aaq_bins_init0(b)=atot*surrogate[i].MM/(1+Keq1*(1+Keq2));
		    surrogate[j].Aaq_bins_init0(b)=atot*surrogate[j].MM/(1+Keq1*(1+Keq2))*Keq1;
		    surrogate[k].Aaq_bins_init0(b)=atot*surrogate[k].MM/(1+Keq1*(1+Keq2))*Keq1*Keq2;

		    double prod=surrogate[i].kprod_aq(b)/surrogate[i].MM+surrogate[j].kprod_aq(b)/surrogate[j].MM+surrogate[k].kprod_aq(b)/surrogate[k].MM;
		    if (surrogate[i].kprod_aq(b)>0.)
		      surrogate[i].kloss_aq(b)=prod/(1+Keq1*(1+Keq2))*surrogate[i].MM*surrogate[i].kloss_aq(b)/surrogate[i].kprod_aq(b);
		    //cout << surrogate[i].kprod_aq(b) << " " << surrogate[j].kprod_aq(b) << " " << surrogate[k].kprod_aq(b) << endl;
		    if (surrogate[j].kprod_aq(b)>0.)
		      surrogate[j].kloss_aq(b)=prod/(1+Keq1*(1+Keq2))*Keq1*surrogate[j].MM*surrogate[j].kloss_aq(b)/surrogate[j].kprod_aq(b);
		    if (surrogate[k].kprod_aq(b)>0.)
		      surrogate[k].kloss_aq(b)=prod/(1+Keq1*(1+Keq2))*Keq1*Keq2*surrogate[k].MM*surrogate[k].kloss_aq(b)/surrogate[k].kprod_aq(b);
		   
		    surrogate[i].kprod_aq(b)=prod/(1+Keq1*(1+Keq2))*surrogate[i].MM;
		    surrogate[j].kprod_aq(b)=prod/(1+Keq1*(1+Keq2))*Keq1*surrogate[j].MM;
		    surrogate[k].kprod_aq(b)=prod/(1+Keq1*(1+Keq2))*Keq1*Keq2*surrogate[k].MM;

		    surrogate[i].kprod_gas+=surrogate[i].kloss_aq(b)*surrogate[i].Aaq_bins_init(b);
		    surrogate[j].kprod_gas+=surrogate[j].kloss_aq(b)*surrogate[j].Aaq_bins_init(b);
		    surrogate[k].kprod_gas+=surrogate[k].kloss_aq(b)*surrogate[k].Aaq_bins_init(b);
		    
		  }

	      //cout << "after: " << surrogate[i].kprod_gas/surrogate[i].kloss_gas << endl;

	      if (surrogate[i].hydrophobic)
		for (b=0;b<config.nbins;++b)
		  for (ilayer=0;ilayer<config.nlayer;++ilayer)
		    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		      {
			double Keq1=RH*surrogate[i].Khyd*surrogate[i].gamma_org_layer(b,ilayer,iphase)/surrogate[i].GAMMAinf/surrogate[j].gamma_org_layer(b,ilayer,iphase)*surrogate[j].GAMMAinf;
			double Keq2=RH*surrogate[j].Khyd*surrogate[j].gamma_org_layer(b,ilayer,iphase)/surrogate[j].GAMMAinf/surrogate[k].gamma_org_layer(b,ilayer,iphase)*surrogate[k].GAMMAinf;

			double atot=surrogate[i].Ap_layer_init(b,ilayer,iphase)/surrogate[i].MM+surrogate[j].Ap_layer_init(b,ilayer,iphase)/surrogate[j].MM+surrogate[k].Ap_layer_init(b,ilayer,iphase)/surrogate[k].MM;

			surrogate[i].Ap_layer_init(b,ilayer,iphase)=atot*surrogate[i].MM/(1+Keq1*(1+Keq2));
			surrogate[j].Ap_layer_init(b,ilayer,iphase)=atot*surrogate[j].MM/(1+Keq1*(1+Keq2))*Keq1;
			surrogate[k].Ap_layer_init(b,ilayer,iphase)=atot*surrogate[k].MM/(1+Keq1*(1+Keq2))*Keq1*Keq2;

			atot=surrogate[i].Ap_layer_init0(b,ilayer,iphase)/surrogate[i].MM+surrogate[j].Ap_layer_init0(b,ilayer,iphase)/surrogate[j].MM+surrogate[k].Ap_layer_init0(b,ilayer,iphase)/surrogate[k].MM;
			surrogate[i].Ap_layer_init0(b,ilayer,iphase)=atot*surrogate[i].MM/(1+Keq1*(1+Keq2));
			surrogate[j].Ap_layer_init0(b,ilayer,iphase)=atot*surrogate[j].MM/(1+Keq1*(1+Keq2))*Keq1;
			surrogate[k].Ap_layer_init0(b,ilayer,iphase)=atot*surrogate[k].MM/(1+Keq1*(1+Keq2))*Keq1*Keq2;		 

			double prod=surrogate[i].kprod(b,ilayer,iphase)/surrogate[i].MM+surrogate[j].kprod(b,ilayer,iphase)/surrogate[j].MM+surrogate[k].kprod(b,ilayer,iphase)/surrogate[k].MM;
			if (surrogate[i].kprod(b,ilayer,iphase)>0.)
			  surrogate[i].kloss(b,ilayer,iphase)=prod/(1+Keq1*(1+Keq2))*surrogate[i].MM*surrogate[i].kloss(b,ilayer,iphase)/surrogate[i].kprod(b,ilayer,iphase);
			if (surrogate[j].kprod(b,ilayer,iphase)>0.)
			  surrogate[j].kloss(b,ilayer,iphase)=prod/(1+Keq1*(1+Keq2))*Keq1*surrogate[j].MM*surrogate[j].kloss(b,ilayer,iphase)/surrogate[j].kprod(b,ilayer,iphase);
			if (surrogate[k].kprod(b,ilayer,iphase)>0.)
			  surrogate[k].kloss(b,ilayer,iphase)=prod/(1+Keq1*(1+Keq2))*Keq1*Keq2*surrogate[k].MM*surrogate[k].kloss(b,ilayer,iphase)/surrogate[k].kprod(b,ilayer,iphase);
		   
			surrogate[i].kprod(b,ilayer,iphase)=prod/(1+Keq1*(1+Keq2))*surrogate[i].MM;
			surrogate[j].kprod(b,ilayer,iphase)=prod/(1+Keq1*(1+Keq2))*Keq1*surrogate[j].MM;
			surrogate[k].kprod(b,ilayer,iphase)=prod/(1+Keq1*(1+Keq2))*Keq1*Keq2*surrogate[k].MM;

			surrogate[i].kprod_gas+=surrogate[i].kloss(b,ilayer,iphase)*surrogate[i].Ap_layer_init(b,ilayer,iphase);
			surrogate[j].kprod_gas+=surrogate[j].kloss(b,ilayer,iphase)*surrogate[j].Ap_layer_init(b,ilayer,iphase);
			surrogate[k].kprod_gas+=surrogate[k].kloss(b,ilayer,iphase)*surrogate[k].Ap_layer_init(b,ilayer,iphase);


			//cout << surrogate[i].Ag << " " << surrogate[j].Ag << " " << surrogate[k].Ag << endl;
			//cout << b << " " << ilayer << " " << iphase << " " << surrogate[i].Ap_layer_init(b,ilayer,iphase) << " " << surrogate[j].Ap_layer_init(b,ilayer,iphase) << " " << surrogate[k].Ap_layer_init(b,ilayer,iphase) << endl;

		      }

	      //cout << surrogate[i].kprod_gas << " " << surrogate[j].kprod_gas << " " << surrogate[k].kprod_gas << endl;
		
	    }
	  else
	    {
	      //cout << "nooon..." << endl;
	      /*
	      surrogate[i].kprod_gas=0.;
	      surrogate[j].kprod_gas=0.;
	      if (surrogate[i].hydrophilic)
		for (b=0;b<config.nbins;++b)
		  {
		    double Keq1=RH*surrogate[i].Khyd*surrogate[i].gamma_aq_bins(b)/surrogate[j].gamma_aq_bins(b); 
		    //cout << Keq1 << " " << Keq2 << endl;
		    
		    double atot=surrogate[i].Aaq_bins_init(b)/surrogate[i].MM+surrogate[j].Aaq_bins_init(b)/surrogate[j].MM;
		    surrogate[i].Aaq_bins_init(b)=atot*surrogate[i].MM/(1+Keq1);
		    surrogate[j].Aaq_bins_init(b)=atot*surrogate[j].MM/(1+Keq1)*Keq1;
		    
		    atot=surrogate[i].Aaq_bins_init0(b)/surrogate[i].MM+surrogate[j].Aaq_bins_init0(b)/surrogate[j].MM;
		    surrogate[i].Aaq_bins_init0(b)=atot*surrogate[i].MM/(1+Keq1);
		    surrogate[j].Aaq_bins_init0(b)=atot*surrogate[j].MM/(1+Keq1)*Keq1;

		    double prod=surrogate[i].kprod_aq(b)/surrogate[i].MM+surrogate[j].kprod_aq(b)/surrogate[j].MM;
		    if (surrogate[i].kprod_aq(b)>0.)
		      surrogate[i].kloss_aq(b)=prod/(1+Keq1)*surrogate[i].MM*surrogate[i].kloss_aq(b)/surrogate[i].kprod_aq(b);
		    //cout << surrogate[i].kprod_aq(b) << " " << surrogate[j].kprod_aq(b) << " " << surrogate[k].kprod_aq(b) << endl;
		    if (surrogate[j].kprod_aq(b)>0.)
		      surrogate[j].kloss_aq(b)=prod/(1+Keq1)*Keq1*surrogate[j].MM*surrogate[j].kloss_aq(b)/surrogate[j].kprod_aq(b);
		   
		    surrogate[i].kprod_aq(b)=prod/(1+Keq1)*surrogate[i].MM;
		    surrogate[j].kprod_aq(b)=prod/(1+Keq1)*Keq1*surrogate[j].MM;

		    surrogate[i].kprod_gas+=surrogate[i].kloss_aq(b)*surrogate[i].Aaq_bins_init(b);
		    surrogate[j].kprod_gas+=surrogate[j].kloss_aq(b)*surrogate[j].Aaq_bins_init(b); 
		  }
*/

	      
	      surrogate[i].kprod_gas=0.;
	      surrogate[j].kprod_gas=0.;
	      if (surrogate[i].hydrophilic)
		for (b=0;b<config.nbins;++b)
		  {
		    double Keq1=RH*surrogate[i].Khyd*surrogate[i].gamma_aq_bins(b)/surrogate[j].gamma_aq_bins(b);
		    double atot=surrogate[i].Aaq_bins_init(b)/surrogate[i].MM+surrogate[j].Aaq_bins_init(b)/surrogate[j].MM;
		    surrogate[i].Aaq_bins_init(b)=atot*surrogate[i].MM/(1+Keq1);
		    surrogate[j].Aaq_bins_init(b)=atot*surrogate[j].MM/(1+Keq1)*Keq1;

		   
		    atot=surrogate[i].Aaq_bins_init0(b)/surrogate[i].MM+surrogate[j].Aaq_bins_init0(b)/surrogate[j].MM;
		    surrogate[i].Aaq_bins_init0(b)=atot*surrogate[i].MM/(1+Keq1);
		    surrogate[j].Aaq_bins_init0(b)=atot*surrogate[j].MM/(1+Keq1)*Keq1;

		    double prod=surrogate[i].kprod_aq(b)/surrogate[i].MM+surrogate[j].kprod_aq(b)/surrogate[j].MM;
		    if (surrogate[i].kprod_aq(b)>0.)
		      surrogate[i].kloss_aq(b)=prod/(1+Keq1)*surrogate[i].MM*surrogate[i].kloss_aq(b)/surrogate[i].kprod_aq(b);
		    if (surrogate[j].kprod_aq(b)>0.)
		      surrogate[j].kloss_aq(b)=prod/(1+Keq1)*Keq1*surrogate[j].MM*surrogate[j].kloss_aq(b)/surrogate[j].kprod_aq(b);
		   
		    surrogate[i].kprod_aq(b)=prod/(1+Keq1)*surrogate[i].MM;
		    surrogate[j].kprod_aq(b)=prod/(1+Keq1)*Keq1*surrogate[j].MM;		   

		    surrogate[i].kprod_gas+=surrogate[i].kloss_aq(b)*surrogate[i].Aaq_bins_init(b);
		    surrogate[j].kprod_gas+=surrogate[j].kloss_aq(b)*surrogate[j].Aaq_bins_init(b);
		    
		  }

	      

	      if (surrogate[i].hydrophobic)
		for (b=0;b<config.nbins;++b)
		  for (ilayer=0;ilayer<config.nlayer;++ilayer)
		    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		      {
			double Keq1=RH*surrogate[i].Khyd*surrogate[i].gamma_org_layer(b,ilayer,iphase)/surrogate[i].GAMMAinf/surrogate[j].gamma_org_layer(b,ilayer,iphase)*surrogate[j].GAMMAinf;
			
			double atot=surrogate[i].Ap_layer_init(b,ilayer,iphase)/surrogate[i].MM+surrogate[j].Ap_layer_init(b,ilayer,iphase)/surrogate[j].MM;
			
			surrogate[i].Ap_layer_init(b,ilayer,iphase)=atot*surrogate[i].MM/(1+Keq1);
			surrogate[j].Ap_layer_init(b,ilayer,iphase)=atot*surrogate[j].MM/(1+Keq1)*Keq1;

			atot=surrogate[i].Ap_layer_init0(b,ilayer,iphase)/surrogate[i].MM+surrogate[j].Ap_layer_init0(b,ilayer,iphase)/surrogate[j].MM;
			surrogate[i].Ap_layer_init0(b,ilayer,iphase)=atot*surrogate[i].MM/(1+Keq1);
			surrogate[j].Ap_layer_init0(b,ilayer,iphase)=atot*surrogate[j].MM/(1+Keq1)*Keq1;

			double prod=surrogate[i].kprod(b,ilayer,iphase)/surrogate[i].MM+surrogate[j].kprod(b,ilayer,iphase)/surrogate[j].MM;
			if (surrogate[i].kprod(b,ilayer,iphase)>0.)
			  surrogate[i].kloss(b,ilayer,iphase)=prod/(1+Keq1)*surrogate[i].MM*surrogate[i].kloss(b,ilayer,iphase)/surrogate[i].kprod(b,ilayer,iphase);
			if (surrogate[j].kprod(b,ilayer,iphase)>0.)
			  surrogate[j].kloss(b,ilayer,iphase)=prod/(1+Keq1)*Keq1*surrogate[j].MM*surrogate[j].kloss(b,ilayer,iphase)/surrogate[j].kprod(b,ilayer,iphase);
		   
			surrogate[i].kprod(b,ilayer,iphase)=prod/(1+Keq1)*surrogate[i].MM;
			surrogate[j].kprod(b,ilayer,iphase)=prod/(1+Keq1)*Keq1*surrogate[j].MM;

			surrogate[i].kprod_gas+=surrogate[i].kloss(b,ilayer,iphase)*surrogate[i].Ap_layer_init(b,ilayer,iphase);
			surrogate[j].kprod_gas+=surrogate[j].kloss(b,ilayer,iphase)*surrogate[j].Ap_layer_init(b,ilayer,iphase);
		    
		      }

	      
	    }
	}
    }

  
  
}

void characteristic_time_ssh(model_config &config, vector<species>& surrogate,
			     Array <double, 3> &MOinit, Array <double, 1> &AQinit,
			     double LWCtot)
{
  //compute characteristic time for each compound to reach an equilibrium between the gas phase and
  //the organic phase for each bin and each layer
  int i,b,ilayer,iphase;
  int n=surrogate.size();
  double sum1,sum2,sum3;
  double time_dif;
  double Vlayer_dif=0.0;

  //compute total particulate concentrations
  if (config.explicit_representation)
    {      
      for (i=0;i<n;++i)
	if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophobic)
	  for (b=0;b<config.nbins;++b)
	    {	      
	      time_dif=0.0;
	      Vlayer_dif=0.0;
	
	      for (ilayer=config.nlayer-1;ilayer>=0;--ilayer)
		{
		  Vlayer_dif+=config.Vlayer(ilayer);
	  
		  sum1=0.0;
		  for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		    sum1+=surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase);

		  sum2=0.0;
		  for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		    sum2+=surrogate[i].Ap_layer_init(b,ilayer,iphase);

		  sum3=0.0;
		  for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		    sum3+=MOinit(b,ilayer,iphase);

		  time_dif=pow(config.dbound(b,config.nlayer)-config.dbound(b,ilayer),2)/(surrogate[i].dif_org(b,ilayer));			   	      
		  for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		    surrogate[i].tau_diffusion(b,ilayer,iphase)=time_dif; 

		  if (sum2>0.0)
		    if(sum1/sum3<config.kp_low_volatility)	     
		      for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)			
			{			    
			  surrogate[i].time(b,ilayer,iphase)=
			    (surrogate[i].tau_diffusion(b,ilayer,iphase)+sum1/config.Vlayer(ilayer)*Vlayer_dif*surrogate[i].tau_air(b)*
			     (sum3/config.Vlayer(ilayer)+AQinit(b))/sum3*config.Vlayer(ilayer))
			    /(1.0+sum1/sum2*surrogate[i].Ap);
			}
		    else 
		      for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
			surrogate[i].time(b,ilayer,iphase)=2.0*config.tequilibrium;
		  else
		    if(sum1/sum3<config.kp_low_volatility)	     
		      for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)			
			{			    
			  surrogate[i].time(b,ilayer,iphase)=
			    (surrogate[i].tau_diffusion(b,ilayer,iphase)+sum1/config.Vlayer(ilayer)*Vlayer_dif*surrogate[i].tau_air(b)*
			     (sum3/config.Vlayer(ilayer)+AQinit(b))/sum3*config.Vlayer(ilayer))/(1.+sum1);			  
			}
		    else 
		      for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
			surrogate[i].time(b,ilayer,iphase)=2.0*config.tequilibrium;
		  /*
		    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		      surrogate[i].time(b,ilayer,iphase)=0.0; */      		  
		}
	    }

	else
	  for (b=0;b<config.nbins;++b)
	    for (ilayer=0;ilayer<config.nlayer;++ilayer)
	      for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		surrogate[i].time(b,ilayer,iphase)=2.0*config.tequilibrium;

    }
  else
    {
      for (i=0;i<n;++i)
	if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophobic)
	  {
	    surrogate[i].Ap=sum(surrogate[i].Ap_layer_init);
	    if (surrogate[i].hydrophilic and LWCtot>config.LWClimit)
	      surrogate[i].Ap+=sum(surrogate[i].Aaq_bins_init);
	  }

      //compute the characteristic time
      for (i=0;i<n;++i)
	if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophobic)
	  for (b=0;b<config.nbins;++b)	  
	    for (ilayer=0;ilayer<config.nlayer;++ilayer)
	      {
		sum1=0.0;
		for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		  sum1+=surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase);

		sum2=0.0;
		for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		  sum2+=surrogate[i].Ap_layer_init(b,ilayer,iphase);

		sum3=0.0;
		for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		  sum3+=MOinit(b,ilayer,iphase);

		if (sum2>0.0)
		  if(sum1/sum3<config.kp_low_volatility)
		    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		      {
			surrogate[i].time(b,ilayer,iphase)=
			  (surrogate[i].tau_diffusion(b,ilayer,iphase)+
			   sum1/config.Vlayer(ilayer)*surrogate[i].tau_air(b)*
			   (sum3/config.Vlayer(ilayer)+AQinit(b))/sum3*config.Vlayer(ilayer))
			  /(1.0+sum1/sum2*surrogate[i].Ap);
		      }

		  else 
		    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		      surrogate[i].time(b,ilayer,iphase)=2.0*config.tequilibrium;
		else
		  if(sum1/sum3<config.kp_low_volatility)
		    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		      surrogate[i].time(b,ilayer,iphase)=
			(surrogate[i].tau_diffusion(b,ilayer,iphase)+
			 sum1/config.Vlayer(ilayer)*surrogate[i].tau_air(b)*
			 (sum3/config.Vlayer(ilayer)+AQinit(b))/sum3*config.Vlayer(ilayer))
			/(1.0+sum1);
		  else 
		    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		      surrogate[i].time(b,ilayer,iphase)=2.0*config.tequilibrium;
		//surrogate[i].time(b,ilayer,iphase)=0.0;
	      }
    }  

  if (config.imethod==0)
    for (b=0;b<config.nbins;++b)	  
      for (ilayer=config.nlayer-1;ilayer>=0;--ilayer)
	for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
	  {
	    for (i=0;i<n;i++)
	      if (ilayer<config.nlayer-1 and surrogate[i].tau_diffusion(b,ilayer,iphase)<config.tequilibrium and surrogate[i].time(b,ilayer+1,iphase)<config.tequilibrium)
		surrogate[i].time(b,ilayer,iphase)=0.0;
	    
	    surrogate[config.iH2O].time(b,ilayer,iphase)=0.0; //H2O is force to be at equilibrium
	  }
}

void characteristic_time_aq_ssh(model_config &config, vector<species>& surrogate, double &Temperature,
				Array <double, 1> &chp,
				Array <double, 1> &LWC, Array <double, 1> &AQinit,
				Array<double, 3> &MOinit, Array <double, 1> &MMaq, Array <double, 1> &ionic)
{
  //compute characteristic time for each compound to reach an equilibrium between the gas phase and
  //the aqueous phase in a bin
  int i,b,ilayer,iphase;
  int n=surrogate.size();
  double Kaq;
  double LWCtot=0.0;
  Array <double, 1> max_time_inorg;
  max_time_inorg.resize(config.nbins);
  double conc_aq,sumMO;
  double AQ2;
  double a1,f1;
  double b1;
  double c1;
  double delta;

  for (b=0;b<config.nbins;++b)
    max_time_inorg(b)=0.0;
  
  for (b=0;b<config.nbins;++b)
    LWCtot+=LWC(b);

  Array <double, 1> conc_org;
  conc_org.resize(config.nbins);
  conc_org=0.;
  for (b=0;b<config.nbins;++b)
    compute_conc_org_bins_ssh(config, surrogate, Temperature,  MMaq,
			      AQinit, chp, ionic,LWC, b, conc_org(b));


  if(LWCtot>config.LWClimit) //if there is enough water to have an aqueous phase
    {     
      //compute the total paticulate concentration
      for (i=0;i<n;++i)
        if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophilic)
          {
            surrogate[i].Aaq=0.0;
            for (b=0;b<config.nbins;++b)
              surrogate[i].Aaq+=surrogate[i].Aaq_bins_init(b);
            
            if (surrogate[i].hydrophobic)
              for (b=0;b<config.nbins;++b)
                for (ilayer=0;ilayer<config.nlayer;++ilayer)
                  for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                    surrogate[i].Aaq+=surrogate[i].Ap_layer_init(b,ilayer,iphase);
          }
        else if(surrogate[i].is_inorganic_precursor and config.compute_inorganic and surrogate[i].is_solid==false)
          {
            surrogate[i].Aaq=0.0;
            if (i==config.iH2SO4)
              for (b=0;b<config.nbins;++b)
                surrogate[i].Aaq+=surrogate[config.iHSO4m].Aaq_bins_init(b)/surrogate[config.iHSO4m].MM*surrogate[i].MM
                  +surrogate[config.iSO4mm].Aaq_bins_init(b)/surrogate[config.iSO4mm].MM*surrogate[i].MM;
            else if (i==config.iNH3)
              {
                for (b=0;b<config.nbins;++b)
                  surrogate[i].Aaq+=surrogate[config.iNH4p].Aaq_bins_init(b)/surrogate[config.iNH4p].MM*surrogate[i].MM;
              }
            else if (i==config.iHNO3)
              {
                for (b=0;b<config.nbins;++b)
                  surrogate[i].Aaq+=surrogate[config.iNO3m].Aaq_bins_init(b)/surrogate[config.iNO3m].MM*surrogate[i].MM;		
              }
	    else if (i==config.iCO2)
              {
                for (b=0;b<config.nbins;++b)
                  surrogate[i].Aaq+=surrogate[config.iHCO3m].Aaq_bins_init(b)/surrogate[config.iHCO3m].MM*surrogate[i].MM+surrogate[config.iCO3mm].Aaq_bins_init(b)/surrogate[config.iCO3mm].MM*surrogate[i].MM;		
              }
            else if (i==config.iHCl)
              for (b=0;b<config.nbins;++b)
                surrogate[i].Aaq+=surrogate[config.iClm].Aaq_bins_init(b)/surrogate[config.iClm].MM*surrogate[i].MM;
          }
      
      //compute the characteristic time
      conc_aq,sumMO;
      for (i=0;i<n;++i)
        if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophilic)
          for (b=0;b<config.nbins;++b)
            {
              sumMO=0.0;
              for (ilayer=0;ilayer<config.nlayer;++ilayer)
                for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                  sumMO+=MOinit(b,ilayer,iphase);
              
              Kaq=surrogate[i].Kaq(b);
              if(Kaq<config.kp_low_volatility)
                if (surrogate[i].Aaq_bins_init(b)>1.e-10)
		  {
		    AQ2=AQinit(b);
		    surrogate[i].time_aq(b)=
		      (AQ2+sumMO)/AQ2*
		      (Kaq*AQ2*surrogate[i].tau_air(b))
		      /(1.0+Kaq*AQ2/surrogate[i].Aaq_bins_init(b)*surrogate[i].Aaq);

		    if (i==config.iH2O)
		      {		       
		        a1=Kaq;
			b1=1.0-Kaq*AQ2-Kaq*surrogate[i].Atot;
			c1=-(AQinit(b)-surrogate[i].Aaq_bins_init(b));
			delta=pow(b1,2.0)-4.0*a1*c1;

			if (delta>=0)
			  {
			    AQ2=max(AQinit(b),(-b1+pow(delta,0.5))/(2.0*a1));
			    f1=pow((AQ2+sumMO)/(AQinit(b)+sumMO),1.0/3.0);			
			    surrogate[i].time_aq(b)=
			      (AQ2+sumMO)/AQ2*
			      (Kaq*AQ2*surrogate[i].tau_air(b)/f1)
			      /(1.0+Kaq*AQ2/surrogate[i].Aaq_bins_init(b)*surrogate[i].Aaq);
			  }
		      }

		  }
                else
		  {
		    AQ2=AQinit(b); 
		    surrogate[i].time_aq(b)=
		      (AQ2+sumMO)/AQ2*
		      (Kaq*AQ2*surrogate[i].tau_air(b))
		      /(1.0+Kaq*AQ2);

		    if (i==config.iH2O)
		      {		       
		        a1=Kaq;
			b1=1.0-Kaq*AQ2-Kaq*surrogate[i].Atot;
			c1=-(AQinit(b)-surrogate[i].Aaq_bins_init(b));
			delta=pow(b1,2.0)-4.0*a1*c1;

			if (delta>=0)
			  {
			    AQ2=max(AQinit(b),(-b1+pow(delta,0.5))/(2.0*a1));
			    f1=pow((AQ2+sumMO)/(AQinit(b)+sumMO),1.0/3.0);			
			    surrogate[i].time_aq(b)=
			      (AQ2+sumMO)/AQ2*
			      (Kaq*AQ2*surrogate[i].tau_air(b)/f1)
			      /(1.0+Kaq*AQ2);
			  }
		      }
		  }                  
              else
                surrogate[i].time_aq(b)=2.0*config.tequilibrium;
            }
        else if (surrogate[i].is_inorganic_precursor and config.compute_inorganic and surrogate[i].is_solid==false)
          for (b=0;b<config.nbins;++b)
            {
              conc_aq=0.0;
              if (i==config.iH2SO4)
                conc_aq+=surrogate[config.iHSO4m].Aaq_bins_init(b)/surrogate[config.iHSO4m].MM*surrogate[i].MM
                  +surrogate[config.iSO4mm].Aaq_bins_init(b)/surrogate[config.iSO4mm].MM*surrogate[i].MM;
              else if (i==config.iNH3)
                conc_aq+=surrogate[config.iNH4p].Aaq_bins_init(b)/surrogate[config.iNH4p].MM*surrogate[i].MM;
              else if (i==config.iHNO3)
                conc_aq+=surrogate[config.iNO3m].Aaq_bins_init(b)/surrogate[config.iNO3m].MM*surrogate[i].MM;
	      else if (i==config.iCO2)
		conc_aq+=surrogate[config.iHCO3m].Aaq_bins_init(b)/surrogate[config.iHCO3m].MM*surrogate[i].MM+surrogate[config.iCO3mm].Aaq_bins_init(b)/surrogate[config.iCO3mm].MM*surrogate[i].MM;
              else if (i==config.iHCl)
                conc_aq+=surrogate[config.iClm].Aaq_bins_init(b)/surrogate[config.iClm].MM*surrogate[i].MM;
              Kaq=surrogate[i].Kaq(b);
	      
              if(Kaq<config.kp_low_volatility)		
                if (conc_aq>0.0)
                  surrogate[i].time_aq(b)=(Kaq*conc_org(b)*surrogate[i].tau_air(b))
                    /(1.0+Kaq*conc_org(b)/conc_aq*surrogate[i].Aaq);
                else 
                  surrogate[i].time_aq(b)=Kaq*conc_org(b)*surrogate[i].tau_air(b); //0.0*config.tequilibrium;
              else
                surrogate[i].time_aq(b)=2.0*config.tequilibrium;

	      /*if (i==config.iCO2)
		{
		  cout << Kaq << " " <<  config.kp_low_volatility << " " << surrogate[i].time_aq(b) << " " << Kaq*AQinit(b)*surrogate[i].tau_air(b) << endl;
		  exit(0);
		  }*/

	      /*
	      if (surrogate[i].time_aq(b)==0 and i==config.iCO2)
		{
		  cout << "time is zero" << " " << conc_aq << endl;
		  cout << surrogate[config.iHCO3m].Aaq_bins_init << endl;
		  cout << surrogate[config.iCO3mm].Aaq_bins_init << endl;
		  exit(0);
		}
	      if (Kaq>=config.kp_low_volatility and i==config.iCO2)
		{
		  cout << "time is above" << " " << conc_aq << endl;
		  cout << surrogate[config.iHCO3m].Aaq_bins_init << endl;
		  cout << surrogate[config.iCO3mm].Aaq_bins_init << endl;
		  cout << surrogate[i].Kaq(b) << endl;
		  cout << chp << endl;
		  cout << surrogate[config.iHp].gamma_aq_bins << endl;
		  cout << surrogate[config.iCO3mm].gamma_aq_bins << endl;
		  cout << surrogate[config.iHCO3m].gamma_aq_bins << endl;
		  exit(0);
		}*/

	      	      
              if (i==config.iHCl or i==config.iHNO3 or i==config.iNH3 or i==config.iCO2)		
		max_time_inorg(b)=max(max_time_inorg(b),surrogate[i].time_aq(b));	        		

              if (i==config.iH2SO4)
                surrogate[i].time_aq(b)=2.0*config.tequilibrium;
            }
    }
  else
    for (i=0;i<n;++i)
      if(surrogate[i].hydrophilic)
        for (b=0;b<config.nbins;++b)
          surrogate[i].time_aq(b)=config.tequilibrium;

  /*  
      if (config.compute_inorganic)
      {
      double sensicl=0.;
      double sensino3=0.;
      double sensinh4=0.;      
      for (b=0;b<config.nbins;++b)
      {
      double Jhp=0.0;
      double asulf=0.;
      double a1=0.;
      double b1=0.;
      double c1=0.;
      double conc_org=LWC(b);
      double faq=0.;
      double conc_aq;
      double organion=0.;
      double organion2=0.;
      for (i=0;i<n;++i)
      if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophilic)
      conc_org+=surrogate[i].Aaq_bins_init(b);
      
      double sum_mass=AQinit(b);
      for (ilayer=0;ilayer<config.nlayer;ilayer++)
      for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
      sum_mass+=MOinit(b,ilayer,iphase);

      if (sum_mass>0.0)
      faq=AQinit(b)/sum_mass;
      else
      faq=0.0;
          
      i=config.iH2SO4;          
      if (i>=0)
      {
      double K=surrogate[i].Kequilibrium_ssh(Temperature)*surrogate[config.iHSO4m].gamma_aq_bins(b)
      /(surrogate[config.iHp].gamma_aq_bins(b)*surrogate[config.iSO4mm].gamma_aq_bins(b));
              
      a1+=(2.0-1.0/(1.0+K/chp(b)))*surrogate[i].Ag/surrogate[i].tau_air(b)/surrogate[i].MM; 		
		      
      double total=surrogate[config.iSO4mm].Aaq_bins_init(b)/surrogate[config.iSO4mm].MM+
      surrogate[config.iHSO4m].Aaq_bins_init(b)/surrogate[config.iHSO4m].MM;	       		     
      asulf=config.AQrho(b)*total/AQinit(b)*1.0/pow(1.0+K/chp(b),2.0)*K/(chp(b)*chp(b));
              
      }
      i=config.iNH3;
      if (i>=0)
      {		
      Kaq=surrogate[i].Kaq(b)/chp(b);	       
      conc_aq=surrogate[config.iNH4p].Aaq_bins_init(b)/surrogate[config.iNH4p].MM;
      c1+=conc_aq/(Kaq*AQinit(b)*surrogate[i].tau_air(b))*faq;
      a1-=surrogate[i].Ag/surrogate[i].tau_air(b)/surrogate[i].MM*faq;
      }
      i=config.iHNO3;
      if (i>=0)
      {
      Kaq=surrogate[i].Kaq(b)*chp(b);	       
      conc_aq=surrogate[config.iNO3m].Aaq_bins_init(b)/surrogate[config.iNO3m].MM;
      b1-=conc_aq/(Kaq*AQinit(b)*surrogate[i].tau_air(b))*faq;
      a1+=surrogate[i].Ag/surrogate[i].tau_air(b)/surrogate[i].MM*faq;
      }
      i=config.iHCl;
      if (i>=0)            
      {
      Kaq=surrogate[i].Kaq(b)*chp(b);	       
      conc_aq=surrogate[config.iClm].Aaq_bins_init(b)/surrogate[config.iClm].MM;
      b1-=conc_aq/(Kaq*AQinit(b)*surrogate[i].tau_air(b))*faq;
      a1+=surrogate[i].Ag/surrogate[i].tau_air(b)/surrogate[i].MM*faq;
      }

      a1=a1/conc_org*1000.0;
      b1=b1/conc_org*1000.0;
      c1=c1/conc_org*1000.0;
          
      Jhp=(a1*chp(b)+b1*pow(chp(b),2)+c1+organion*chp(b))/(config.Ke/chp(b)+chp(b)+organion2*chp(b)+chp(b)*asulf);         
      double clconc=(surrogate[config.iClm].Aaq_bins_init(b)/surrogate[config.iClm].MM+surrogate[config.iHCl].Ag/surrogate[config.iHCl].MM);
      double nh4conc=(surrogate[config.iNH4p].Aaq_bins_init(b)/surrogate[config.iNH4p].MM+surrogate[config.iNH3].Ag/surrogate[config.iNH3].MM);
      double no3conc=(surrogate[config.iNO3m].Aaq_bins_init(b)/surrogate[config.iNO3m].MM+surrogate[config.iHNO3].Ag/surrogate[config.iHNO3].MM);
      double sensimax=0.;
      if (config.iHCl>=0)
      sensicl=(surrogate[config.iClm].Aaq_bins_init(b)/surrogate[config.iClm].MM+surrogate[config.iHCl].Ag/surrogate[config.iHCl].MM)*
      (surrogate[config.iHCl].dKaq(b)*AQinit(b)/(1+surrogate[config.iHCl].Kaq(b)*AQinit(b))-surrogate[config.iHCl].Kaq(b)*AQinit(b)/(pow(1+surrogate[config.iHCl].Kaq(b)*AQinit(b),2))*AQinit(b)*surrogate[config.iHCl].dKaq(b));
      if (config.iHNO3>=0)
      sensino3=(surrogate[config.iNO3m].Aaq_bins_init(b)/surrogate[config.iNO3m].MM+surrogate[config.iHNO3].Ag/surrogate[config.iHNO3].MM)*
      (surrogate[config.iHNO3].dKaq(b)*AQinit(b)/(1+surrogate[config.iHNO3].Kaq(b)*AQinit(b))-surrogate[config.iHNO3].Kaq(b)*AQinit(b)/(pow(1+surrogate[config.iHNO3].Kaq(b)*AQinit(b),2))*AQinit(b)*surrogate[config.iHNO3].dKaq(b));
      if (config.iNH3>=0)
      sensinh4=(surrogate[config.iNH4p].Aaq_bins_init(b)/surrogate[config.iNH4p].MM+surrogate[config.iNH3].Ag/surrogate[config.iNH3].MM)*
      (surrogate[config.iNH3].dKaq(b)*AQinit(b)/(1+surrogate[config.iNH3].Kaq(b)*AQinit(b))-surrogate[config.iNH3].Kaq(b)*AQinit(b)/(pow(1+surrogate[config.iNH3].Kaq(b)*AQinit(b),2))*AQinit(b)*surrogate[config.iNH3].dKaq(b)); 
          
      int isensi=config.iHCl;
      sensimax=abs(sensicl);
          
            
      if (abs(sensino3)>sensimax)
      {
      isensi=config.iHNO3;
      sensimax=abs(sensino3);
              
      }
          
      if (abs(sensinh4)>sensimax)
      {
      isensi=config.iNH3;
      sensimax=abs(sensinh4);             
      }

      // The species the more sensitive to pH is forced at equilibrium
      surrogate[isensi].time_aq(b)=0.;
      }
      }*/
  
  for (i=0;i<n;i++)
    if (surrogate[i].aqt==1 or surrogate[i].aqt==2)
      for (b=0;b<config.nbins;b++)         	
	max_time_inorg(b)=max(max_time_inorg(b),(surrogate[i].fion1+surrogate[i].fion2)*surrogate[i].time_aq(b));       

  for (b=0;b<config.nbins;++b)
    surrogate[config.iH2O].time_aq(b)=-10.0; //H2O is forced to be at equilibrium
    
  for (i=0;i<n;i++)
    if (surrogate[i].aqt==1 or surrogate[i].aqt==2)
      for (b=0;b<config.nbins;b++)         
	{	  
	  surrogate[i].time_aq(b)=(1.0-surrogate[i].fion1-surrogate[i].fion2)*surrogate[i].time_aq(b)
	    +(surrogate[i].fion1+surrogate[i].fion2)*max_time_inorg(b); //max(max_time_inorg(b),surrogate[i].time_aq(b));
	}

  if (config.nbins>1)
    for (i=0;i<n;i++)
      {
	int imin=-1;
	for (b=0;b<config.nbins;b++)
	  if (surrogate[i].time_aq(b)<config.tequilibrium)
	    imin=b;

	/*if (i==config.iNH3 or i==config.iHNO3)
	  {
	    cout << surrogate[i].name << " bef " << surrogate[i].time_aq << " " << surrogate[i].Kaq << " " << endl;
	    for (b=0;b<config.nbins;b++)
	      cout << b << " " << 1.0/(surrogate[i].Kaq(b)*conc_org(b)*surrogate[i].tau_air(b));
	      }*/

	if (i==config.iHNO3 or i==config.iNH3 or i==config.iHCl)
	  if (imin>0 and i!=config.iCO2)
	    for (b=imin-1;b>=0;b--)
	      //if (surrogate[i].time_aq(b)>surrogate[i].time_aq(b+1))
	      if (surrogate[i].Kaq(b)*AQinit(b)>0)
		surrogate[i].time_aq(b)=min(surrogate[i].time_aq(b),surrogate[i].time_aq(b+1)*surrogate[i].Kaq(b)*AQinit(b)/surrogate[i].Kaq(b+1)*AQinit(b+1));
	/*
	  for (b=config.nbins-2;b>=0;b--)
	    if (surrogate[i].time_aq(b)>surrogate[i].time_aq(b+1))
	    surrogate[i].time_aq(b)=surrogate[i].time_aq(b+1)*surrogate[i].Kaq(b)*AQinit(b)/(1.+surrogate[i].Kaq(b)*AQinit(b));*/
	//surrogate[i].time_aq(b)=min(surrogate[i].time_aq(b),surrogate[i].time_aq(b+1));

	//if (i==config.iNH3 or i==config.iHNO3)
	  //  cout << surrogate[i].name << " " << surrogate[i].time_aq << " " << surrogate[i].Kaq << endl;
      }
  /*
    for (b=0;b<config.nbins;++b)
    {
    surrogate[config.iNH3].time_aq(b)=min(surrogate[config.iNH3].time_aq(b),surrogate[config.iHNO3].time_aq(b));
    surrogate[config.iHNO3].time_aq(b)=surrogate[config.iNH3].time_aq(b);
    }

    surrogate[config.iNH3].time_aq=surrogate[config.iHNO3].time_aq;*/
  /*
    surrogate[config.iNH3].time_aq=2.*config.tequilibrium;
    surrogate[config.iHNO3].time_aq=2.*config.tequilibrium;*/
}

void compute_flux_chem_ssh(model_config &config, vector<species>& surrogate,
			   Array<double, 3> &MOinit, Array<double, 3> &MOW, 
			   Array<double, 1> &AQinit, Array<double, 1> &LWC, Array<double, 1> &MMaq, 
			   Array<double, 1> &chp, Array<double, 1> &ionic, double &deltat, double &tiny, double &Temperature, int index, double &RH)
{    
  int n=surrogate.size();
  int i,j,b,ilayer,iphase,jion,jmol;
  double gamma=1.+0.5*pow(2,0.5);  
  double sum_flux_gas=0.0;
  double Xmono=0.0;          
  double XH2O;
  double xmin=1.0e-10;
  double conc_org;

  for (i=0;i<n;++i)
    {      
      sum_flux_gas=0.0;
      for (b=0;b<config.nbins;++b)
        {
          for (ilayer=0;ilayer<config.nlayer;ilayer++)
            for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)       
              {
                surrogate[i].flux_chem(b,ilayer,iphase,index)=0.0;
                if (surrogate[i].Ag>tiny)
                  sum_flux_gas-=surrogate[i].k1(b,ilayer,iphase,index)/surrogate[i].Ag;
              }
          surrogate[i].flux_chem_aq(b,index)=0.0;
          if (surrogate[i].Ag>tiny)
            sum_flux_gas-=surrogate[i].k1_aq(b,index)/surrogate[i].Ag;          
        }

      surrogate[i].Jdn_gas(index)=min(sum_flux_gas,0.0);
      surrogate[i].flux_chem_gas(index)=0.0;      
    }

  if (config.chemistry)
    {
      for (b=0;b<config.nbins;++b)
        for (ilayer=0;ilayer<config.nlayer;ilayer++)
          for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
            {        
              Xmono=0.0;          
              for (i=0;i<n;i++)
                if (surrogate[i].is_monomer and surrogate[i].is_organic)              
                  Xmono+=surrogate[i].gamma_org_layer(b,ilayer,iphase)*surrogate[i].Ap_layer_init(b,ilayer,iphase)*MOW(b,ilayer,iphase)/surrogate[i].MM/max(MOinit(b,ilayer,iphase),1.0e-10);              

              XH2O=surrogate[config.iH2O].gamma_org_layer(b,ilayer,iphase)*surrogate[config.iH2O].Ap_layer_init(b,ilayer,iphase)*MOW(b,ilayer,iphase)/surrogate[config.iH2O].MM/max(MOinit(b,ilayer,iphase),1.0e-10);
              xmin=1.0e-10;

              for (i=0;i<n;++i)
                if (surrogate[i].is_organic and surrogate[i].hydrophobic)
		  {		   
		    if (surrogate[i].is_monomer)
		      {
			j=surrogate[i].ioligo;
			double Keq2=surrogate[j].MM/surrogate[i].MM*pow(surrogate[i].Keq_oligo,surrogate[i].moligo-1)*pow(max(Xmono,xmin),surrogate[i].moligo-2)/pow(max(XH2O,xmin),surrogate[i].moligo-1);                      
			double flux=surrogate[i].koligo*(surrogate[i].gamma_org_layer(b,ilayer,iphase)*surrogate[i].Ap_layer_init(b,ilayer,iphase)*Xmono-surrogate[j].gamma_org_layer(b,ilayer,iphase)*surrogate[j].Ap_layer_init(b,ilayer,iphase)/Keq2)*deltat;
			double fac=1.0;
			double fac2=1.0;                      
			if (surrogate[i].time(b,ilayer,iphase)<config.tequilibrium)                                                  
			  if (surrogate[i].Ag+surrogate[i].Ap_layer_init(b,ilayer,iphase)>0.0)                        
			    fac=surrogate[i].Ap_layer_init(b,ilayer,iphase)/(surrogate[i].Ag+surrogate[i].Ap_layer_init(b,ilayer,iphase));                                                                              
                        
			if (surrogate[j].time(b,ilayer,iphase)<config.tequilibrium)                        
			  if (surrogate[j].Ag+surrogate[j].Ap_layer_init(b,ilayer,iphase)>0.0)                                                                                       
			    fac2=surrogate[j].Ap_layer_init(b,ilayer,iphase)/(surrogate[j].Ag+surrogate[j].Ap_layer_init(b,ilayer,iphase));  
			surrogate[i].flux_chem(b,ilayer,iphase,index)+=-fac*flux;
			surrogate[i].flux_chem_gas(index)+=-(1.0-fac)*flux;
			surrogate[j].flux_chem(b,ilayer,iphase,index)+=fac2*flux;			                      
			surrogate[j].flux_chem_gas(index)+=(1.0-fac2)*flux;	                    
		      }

		    for (jmol=0;jmol<surrogate[i].n_irreversible;jmol++)
		      if (surrogate[i].i_irreversible[jmol]>=0)
			{
			  double flux=surrogate[i].k_irreversible[jmol]*deltat*surrogate[i].gamma_org_layer(b,ilayer,iphase)*surrogate[i].Ap_layer_init(b,ilayer,iphase);
			  if (surrogate[i].irr_catalyzed_water[jmol])
			    flux=flux*RH;
			  if (surrogate[i].irr_catalyzed_pH[jmol])
			    flux=flux*config.chp_org_ref;	    
			
			  double fac=1.0;
			  if (surrogate[i].time(b,ilayer,iphase)<config.tequilibrium)                                                  
			    if (surrogate[i].Ag+surrogate[i].Ap_layer_init(b,ilayer,iphase)>0.0)                        
			      fac=surrogate[i].Ap_layer_init(b,ilayer,iphase)/(surrogate[i].Ag+surrogate[i].Ap_layer_init(b,ilayer,iphase)); 
			
			  surrogate[i].flux_chem(b,ilayer,iphase,index)+=-flux*fac;
			  surrogate[i].flux_chem_gas(index)+=-flux*(1.0-fac);
			  //cout << "aaaa111... " << endl;
			  if (surrogate[surrogate[i].i_irreversible[jmol]].hydrophobic)
			    if (surrogate[i].irr_mass_conserving[jmol])
			      surrogate[surrogate[i].i_irreversible[jmol]].flux_chem(b,ilayer,iphase,index)+=flux;
			    else
			      surrogate[surrogate[i].i_irreversible[jmol]].flux_chem(b,ilayer,iphase,index)+=flux*surrogate[surrogate[i].i_irreversible[jmol]].MM/surrogate[i].MM;
			  else
			    if (surrogate[i].irr_mass_conserving[jmol])
			      surrogate[surrogate[i].i_irreversible[jmol]].flux_chem_aq(b,index)+=flux;
			    else
			      surrogate[surrogate[i].i_irreversible[jmol]].flux_chem_aq(b,index)+=flux*surrogate[surrogate[i].i_irreversible[jmol]].MM/surrogate[i].MM;

			  if (surrogate[i].iother_irreversible[jmol]>=0)
			    if (surrogate[surrogate[i].iother_irreversible[jmol]].hydrophobic)
			      surrogate[surrogate[i].iother_irreversible[jmol]].flux_chem(b,ilayer,iphase,index)+=flux*surrogate[surrogate[i].iother_irreversible[jmol]].MM/surrogate[i].MM;
			  else
			    surrogate[surrogate[i].iother_irreversible[jmol]].flux_chem_aq(b,index)+=flux*surrogate[surrogate[i].iother_irreversible[jmol]].MM/surrogate[i].MM;
			}

		  }
            }
  
      for (b=0;b<config.nbins;++b)
        {        
          double Xmonoaq=0.0;
          for (i=0;i<n;i++)
            if (surrogate[i].is_monomer and surrogate[i].is_organic and surrogate[i].hydrophilic)  
              Xmonoaq+=surrogate[i].gamma_aq_bins(b)*surrogate[i].GAMMAinf*surrogate[i].Aaq_bins_init(b)*MMaq(b)/surrogate[i].MM/max(AQinit(b),1.0e-10);

          compute_conc_org_bins_ssh(config, surrogate, Temperature,  MMaq, AQinit, chp, ionic, LWC, b, conc_org);
      
          double XH2Oaq=surrogate[config.iH2O].gamma_aq_bins(b)*surrogate[config.iH2O].Aaq_bins_init(b)*MMaq(b)/surrogate[config.iH2O].MM/max(AQinit(b),1.0e-10);
          double xmin=1.0e-10;

          for (i=0;i<n;++i)
            {
              if (surrogate[i].is_organic)
                if (surrogate[i].is_organic and surrogate[i].hydrophilic) 
                  if (surrogate[i].is_monomer)
                    {
                      j=surrogate[i].ioligo;                                        
                      double Kaq2=surrogate[j].MM/surrogate[i].MM*pow(surrogate[i].Keq_oligo,surrogate[i].moligo-1)*pow(max(Xmonoaq,xmin),surrogate[i].moligo-2)/pow(max(XH2Oaq,xmin),surrogate[i].moligo-1);                                            
                      double flux=surrogate[i].koligo*(surrogate[i].gamma_aq_bins(b)*surrogate[i].GAMMAinf*surrogate[i].Aaq_bins_init(b)*Xmonoaq-surrogate[j].gamma_aq_bins(b)*surrogate[j].GAMMAinf*surrogate[j].Aaq_bins_init(b)/Kaq2)*deltat;
		      if (surrogate[i].catalyzed_ph and config.isorropia_ph==false)			
			flux=flux*chp(b)*surrogate[config.iHp].gamma_aq_bins(b)/config.chp_org_ref;			 			
		      else if (surrogate[i].catalyzed_ph and config.isorropia_ph==true)			
			flux=flux*chp(b)/config.chp_org_ref;	
                      double fac=1.0;
                      double fac2=1.0;
                      
                      if (surrogate[i].time_aq(b)<config.tequilibrium)                                                  
                        if (surrogate[i].Ag+surrogate[i].Aaq_bins_init(b)>0.0)
                          fac=surrogate[i].Aaq_bins_init(b)/(surrogate[i].Ag+surrogate[i].Aaq_bins_init(b));                                                    
                        
                      if (surrogate[j].time_aq(b)<config.tequilibrium)                        
                        if (surrogate[j].Ag+surrogate[j].Aaq_bins_init(b)>0.0)
                          fac2=surrogate[j].Aaq_bins_init(b)/(surrogate[j].Ag+surrogate[j].Aaq_bins_init(b));
                      
                      surrogate[i].flux_chem_aq(b,index)+=-fac*flux;
                      surrogate[i].flux_chem_gas(index)+=-(1.0-fac)*flux;
                      surrogate[j].flux_chem_aq(b,index)+=fac2*flux;		
                      surrogate[j].flux_chem_gas(index)+=(1.0-fac2)*flux;
                    }
	      
              if (surrogate[i].rion and surrogate[i].Aaq_bins_init(b)>0.0 and config.compute_inorganic)
                for (jion=0;jion<surrogate[i].nion_chem;jion++)
                  {            
                    double molality=surrogate[surrogate[i].iion(jion)].gamma_aq_bins(b)*surrogate[surrogate[i].iion(jion)].Aaq_bins_init(b)/surrogate[i].MM/conc_org*1000.0;
                    double flux=surrogate[i].kion[jion]*molality*deltat*surrogate[i].GAMMAinf*surrogate[i].gamma_aq_bins(b)*surrogate[i].Aaq_bins_init(b);
		    if (surrogate[i].rion_ph_catalyzed[jion] and config.isorropia_ph==false)
		      flux=flux*chp(b)*surrogate[config.iHp].gamma_aq_bins(b);
		    else if (surrogate[i].rion_ph_catalyzed[jion] and config.isorropia_ph==true)
		      flux=flux*chp(b);
		    if (surrogate[i].rion_water_catalyzed[jion])
		      flux=flux*RH;
		    
                    double fac=1.0;
                    double fac2=1.0;
                    j=surrogate[i].iproduct(jion);
                    if (surrogate[i].time_aq(b)<config.tequilibrium)                                                  
                      if (surrogate[i].Ag+surrogate[i].Aaq_bins_init(b)>0.0)
                        fac=surrogate[i].Aaq_bins_init(b)/(surrogate[i].Ag+surrogate[i].Aaq_bins_init(b)); 

                    if (surrogate[j].time_aq(b)<config.tequilibrium)                        
                      if (surrogate[j].Ag+surrogate[j].Aaq_bins_init(b)>0.0)
                        fac2=surrogate[j].Aaq_bins_init(b)/(surrogate[j].Ag+surrogate[j].Aaq_bins_init(b));
                                 
                    surrogate[i].flux_chem_aq(b,index)+=-fac*flux;
                    surrogate[i].flux_chem_gas(index)+=-(1.0-fac)*flux;

                    if (surrogate[i].rion_catalyzed[jion]==false and molality>0.)
                      {
                        if (surrogate[i].iion(jion)==config.iHSO4m or surrogate[i].iion(jion)==config.iSO4mm)
                          {
                            double totsulf=surrogate[config.iHSO4m].Aaq_bins_init(b)/surrogate[config.iHSO4m].MM+surrogate[config.iSO4mm].Aaq_bins_init(b)/surrogate[config.iSO4mm].MM;
                            surrogate[config.iHSO4m].flux_chem_aq(b,index)-=surrogate[config.iHSO4m].Aaq_bins_init(b)/totsulf*flux/surrogate[i].MM;
                            surrogate[config.iSO4mm].flux_chem_aq(b,index)-=surrogate[config.iSO4mm].Aaq_bins_init(b)/totsulf*flux/surrogate[i].MM;
                          }
                        if (surrogate[i].iion(jion)==config.iNO3m)
                          {
                            if (surrogate[j].time_aq(b)<config.tequilibrium)
                              {
                                double totno3=surrogate[config.iNO3m].Aaq_bins_init(b)+surrogate[config.iHNO3].Ag/surrogate[config.iHNO3].MM*surrogate[config.iNO3m].MM;
                                surrogate[config.iNO3m].flux_chem_aq(b,index)-=surrogate[config.iNO3m].Aaq_bins_init(b)/totno3*flux/surrogate[i].MM*surrogate[config.iNO3m].MM;
                                surrogate[config.iHNO3].flux_chem_gas(index)-=surrogate[config.iHNO3].Ag/totno3*flux/surrogate[i].MM*surrogate[config.iHNO3].MM;                            
                              }
                            else                          
                              surrogate[config.iNO3m].flux_chem_aq(b,index)-=flux/surrogate[i].MM*surrogate[config.iNO3m].MM;                                                      
                          }

                        if (surrogate[i].iion(jion)==config.iClm)
                          {
                            if (surrogate[j].time_aq(b)<config.tequilibrium)
                              {
                                double totcl=surrogate[config.iClm].Aaq_bins_init(b)+surrogate[config.iHCl].Ag/surrogate[config.iHCl].MM*surrogate[config.iClm].MM;
                                surrogate[config.iClm].flux_chem_aq(b,index)-=surrogate[config.iClm].Aaq_bins_init(b)/totcl*flux/surrogate[i].MM*surrogate[config.iClm].MM;
                                surrogate[config.iHCl].flux_chem_gas(index)-=surrogate[config.iHCl].Ag/totcl*flux/surrogate[i].MM*surrogate[config.iHCl].MM;
                              }
                            else
                              surrogate[config.iClm].flux_chem_aq(b,index)-=flux/surrogate[i].MM*surrogate[config.iClm].MM;
                          }
                        if (surrogate[i].iion(jion)==config.iNH4p)
                          {
                            if (surrogate[j].time_aq(b)<config.tequilibrium)
                              {
                                double totnh4=surrogate[config.iNH4p].Aaq_bins_init(b)+surrogate[config.iNH3].Ag/surrogate[config.iNH3].MM*surrogate[config.iNH4p].MM;
                                surrogate[config.iNH4p].flux_chem_aq(b,index)-=surrogate[config.iNH4p].Aaq_bins_init(n)/totnh4*flux/surrogate[i].MM*surrogate[config.iNH4p].MM;
                                surrogate[config.iNH3].flux_chem_gas(index)-=surrogate[config.iNH3].Ag/totnh4*flux/surrogate[i].MM*surrogate[config.iNH3].MM;
                              }
                            else
                              surrogate[config.iNH4p].flux_chem_aq(b,index)-=flux/surrogate[i].MM*surrogate[config.iNH4p].MM;
                          }
                      }                
                    surrogate[j].flux_chem_aq(b,index)+=fac2*flux*surrogate[j].MM/surrogate[i].MM;		
                    surrogate[j].flux_chem_gas(index)+=(1.0-fac2)*flux*surrogate[j].MM/surrogate[i].MM;
                  }

	      if (surrogate[i].is_organic)
		for (jmol=0;jmol<surrogate[i].n_irreversible;jmol++)
		  if (surrogate[i].i_irreversible[jmol]>=0)
		    {
		      double flux=surrogate[i].k_irreversible[jmol]*deltat*surrogate[i].GAMMAinf*surrogate[i].gamma_aq_bins(b)*surrogate[i].Aaq_bins_init(b);
		      if (surrogate[i].irr_catalyzed_water[jmol])
			flux=flux*RH;
		      if (surrogate[i].irr_catalyzed_pH[jmol] and config.isorropia_ph==false)
			flux=flux*chp(b)*surrogate[config.iHp].gamma_aq_bins(b);
		      else if (surrogate[i].irr_catalyzed_pH[jmol] and config.isorropia_ph==true)
			flux=flux*chp(b);
		    
		      double fac=1.0;
		      if (surrogate[i].time_aq(b)<config.tequilibrium)                                                  
			if (surrogate[i].Ag+surrogate[i].Aaq_bins_init(b)>0.0)
			  fac=surrogate[i].Aaq_bins_init(b)/(surrogate[i].Ag+surrogate[i].Aaq_bins_init(b)); 

		      surrogate[i].flux_chem_aq(b,index)+=-flux*fac;
		      surrogate[i].flux_chem_gas(index)+=-flux*(1.0-fac);

		      if (surrogate[i].irr_mass_conserving[jmol]==false)
			flux=flux*surrogate[surrogate[i].i_irreversible[jmol]].MM/surrogate[i].MM;
		      if (surrogate[surrogate[i].i_irreversible[jmol]].hydrophilic)
			surrogate[surrogate[i].i_irreversible[jmol]].flux_chem_aq(b,index)+=flux;
		      else
			for (ilayer=0;ilayer<config.nlayer;ilayer++)
			  surrogate[surrogate[i].i_irreversible[jmol]].flux_chem(b,ilayer,0,index)+=flux*config.Vlayer(ilayer);

		      if (surrogate[i].iother_irreversible[jmol]>=0)
			{
			  flux=flux*surrogate[surrogate[i].iother_irreversible[jmol]].MM/surrogate[surrogate[i].i_irreversible[jmol]].MM;
			  if (surrogate[surrogate[i].iother_irreversible[jmol]].hydrophilic)
			    surrogate[surrogate[i].iother_irreversible[jmol]].flux_chem_aq(b,index)+=flux;
			  else
			    for (ilayer=0;ilayer<config.nlayer;ilayer++)
			      surrogate[surrogate[i].iother_irreversible[jmol]].flux_chem(b,ilayer,0,index)+=flux*config.Vlayer(ilayer);
			}
		      
		    }
	      
            }    
        }

  
      for (i=0;i<n;++i)
	if (surrogate[i].is_organic or i==config.iH2O)
        {
      
          double sum_flux_gas=0.;
          if (surrogate[i].Ag>tiny)
            sum_flux_gas=surrogate[i].flux_chem_gas(index)/surrogate[i].Ag;

          for (b=0;b<config.nbins;++b)
            {              
		if (surrogate[i].hydrophobic)
		  for (ilayer=0;ilayer<config.nlayer;ilayer++)
		    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		      {
			//ut << surrogate[i].name << " " << b << " " << ilayer << " " << iphase << " " << surrogate[i].time(b,ilayer,iphase) << " " << config.tequilibrium << endl;
			if (surrogate[i].time(b,ilayer,iphase)>=config.tequilibrium)
			  {
			    //ut << "A" << endl;
			    if (surrogate[i].flux_chem(b,ilayer,iphase,index)/deltat+surrogate[i].k1(b,ilayer,iphase,index)<0.0 and
				surrogate[i].Ap_layer_init(b,ilayer,iphase)>tiny)
			      surrogate[i].Jdn(b,ilayer,iphase)=(surrogate[i].flux_chem(b,ilayer,iphase,index)/deltat+surrogate[i].k1(b,ilayer,iphase,index))/
				surrogate[i].Ap_layer_init(b,ilayer,iphase); 
			    else
			      surrogate[i].Jdn(b,ilayer,iphase)=0.0;
			    
			    if (surrogate[i].Ag>tiny)
			      sum_flux_gas-=surrogate[i].k1(b,ilayer,iphase,index)/surrogate[i].Ag;
			  }       
			else
			  {
			    //ut << "B" << endl;
			    if (surrogate[i].flux_chem(b,ilayer,iphase,index)/deltat<0.0 and surrogate[i].Ap_layer_init(b,ilayer,iphase)>tiny)
			      surrogate[i].Jdn(b,ilayer,iphase)=surrogate[i].flux_chem(b,ilayer,iphase,index)/deltat/
				surrogate[i].Ap_layer_init(b,ilayer,iphase); 
			    else
			      surrogate[i].Jdn(b,ilayer,iphase)=0.0;                                    
			  }
		      }

		if (surrogate[i].hydrophilic)
              if (surrogate[i].time_aq(b)>=config.tequilibrium)
                {
                  if (surrogate[i].flux_chem_aq(b,index)/deltat+surrogate[i].k1_aq(b,index)<0.0 and surrogate[i].Aaq_bins_init(b)>tiny)
                    surrogate[i].Jdn_aq(b,index)=(surrogate[i].flux_chem_aq(b,index)/deltat+surrogate[i].k1_aq(b,index))/surrogate[i].Aaq_bins_init(b);
                  else
                    surrogate[i].Jdn_aq(b,index)=0.0;
            
                  if (surrogate[i].Ag>tiny)
                    sum_flux_gas-=surrogate[i].k1_aq(b,index)/surrogate[i].Ag;          
                }
              else
                {
                  if (surrogate[i].flux_chem_aq(b,index)/deltat<0.0 and surrogate[i].Aaq_bins_init(b)>tiny)
                    surrogate[i].Jdn_aq(b,index)=surrogate[i].flux_chem_aq(b,index)/deltat/surrogate[i].Aaq_bins_init(b);
                  else
                    surrogate[i].Jdn_aq(b,index)=0.0;                 
                }

            }
     
          surrogate[i].Jdn_gas(index)=min(sum_flux_gas,0.0);           
        }

      for (i=0;i<n;++i)    
        {      
        
          surrogate[i].Jdn_tot=surrogate[i].flux_chem_gas(index);
          for (b=0;b<config.nbins;++b)
            {
              for (ilayer=0;ilayer<config.nlayer;ilayer++)
                for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)       
                  surrogate[i].Jdn_tot+=surrogate[i].flux_chem(b,ilayer,iphase,index);
              surrogate[i].Jdn_tot+=surrogate[i].flux_chem_aq(b,index);
            }     
          if (surrogate[i].Atot>tiny and surrogate[i].Jdn_tot<0.0)
            surrogate[i].Jdn_tot=surrogate[i].Jdn_tot/deltat/surrogate[i].Atot;
          else
            surrogate[i].Jdn_tot=0.0;
        }    
  

      for (i=0;i<n;++i)
        {
          surrogate[i].flux_chem_gas(index)=0.0;
          for (b=0;b<config.nbins;++b)
            {
              for (ilayer=0;ilayer<config.nlayer;ilayer++)
                for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)       
                  surrogate[i].flux_chem(b,ilayer,iphase,index)=0.0;
              surrogate[i].flux_chem_aq(b,index)=0.0;
            }     
        }

      for (b=0;b<config.nbins;++b)
        for (ilayer=0;ilayer<config.nlayer;ilayer++)
          for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
            {        
              double Xmono=0.0;          
              for (i=0;i<n;i++)
                if (surrogate[i].is_monomer and surrogate[i].is_organic)              
                  Xmono+=surrogate[i].gamma_org_layer(b,ilayer,iphase)*surrogate[i].Ap_layer_init(b,ilayer,iphase)*MOW(b,ilayer,iphase)/surrogate[i].MM/max(MOinit(b,ilayer,iphase),1.0e-10);              

              double XH2O=surrogate[config.iH2O].gamma_org_layer(b,ilayer,iphase)*surrogate[config.iH2O].Ap_layer_init(b,ilayer,iphase)*MOW(b,ilayer,iphase)/surrogate[config.iH2O].MM/max(MOinit(b,ilayer,iphase),1.0e-10);
              double xmin=1.0e-10;

              for (i=0;i<n;++i)
                if (surrogate[i].is_organic and surrogate[i].hydrophobic)
		  {
		    if (surrogate[i].is_monomer)
		      {
			j=surrogate[i].ioligo;                       
			double Keq2=surrogate[j].MM/surrogate[i].MM*pow(surrogate[i].Keq_oligo,surrogate[i].moligo-1)*pow(max(Xmono,xmin),surrogate[i].moligo-2)/pow(max(XH2O,xmin),surrogate[i].moligo-1);                      
			double flux=surrogate[i].koligo*(surrogate[i].gamma_org_layer(b,ilayer,iphase)*surrogate[i].Ap_layer_init(b,ilayer,iphase)*Xmono-surrogate[j].gamma_org_layer(b,ilayer,iphase)*surrogate[j].Ap_layer_init(b,ilayer,iphase)/Keq2)*deltat;                                                       
			double fac=1.0;
			double fac2=1.0;                      
			if (surrogate[i].time(b,ilayer,iphase)<config.tequilibrium)                                                  
			  if (surrogate[i].Ag+surrogate[i].Ap_layer_init(b,ilayer,iphase)>0.0)                        
			    fac=surrogate[i].Ap_layer_init(b,ilayer,iphase)/(surrogate[i].Ag+surrogate[i].Ap_layer_init(b,ilayer,iphase));                                                                              
                        
			if (surrogate[j].time(b,ilayer,iphase)<config.tequilibrium)                        
			  if (surrogate[j].Ag+surrogate[j].Ap_layer_init(b,ilayer,iphase)>0.0)                                                                                      
			    fac2=surrogate[j].Ap_layer_init(b,ilayer,iphase)/(surrogate[j].Ag+surrogate[j].Ap_layer_init(b,ilayer,iphase));         
			surrogate[i].flux_chem(b,ilayer,iphase,index)+=-fac*flux;
			surrogate[i].flux_chem_gas(index)+=-(1.0-fac)*flux;
			surrogate[j].flux_chem(b,ilayer,iphase,index)+=fac2*flux;			                      
			surrogate[j].flux_chem_gas(index)+=(1.0-fac2)*flux;	                    
		      }

		    for (jmol=0;jmol<surrogate[i].n_irreversible;jmol++)
		      if (surrogate[i].i_irreversible[jmol]>=0)
			{
			  double flux=surrogate[i].k_irreversible[jmol]*deltat*surrogate[i].gamma_org_layer(b,ilayer,iphase)*surrogate[i].Ap_layer_init(b,ilayer,iphase);
			  if (surrogate[i].irr_catalyzed_water[jmol])
			    flux=flux*RH;
			  if (surrogate[i].irr_catalyzed_pH[jmol])
			    flux=flux*config.chp_org_ref;
			  double fac=1.0;
			  if (surrogate[i].time(b,ilayer,iphase)<config.tequilibrium)                                                  
			    if (surrogate[i].Ag+surrogate[i].Ap_layer_init(b,ilayer,iphase)>0.0)                        
			      fac=surrogate[i].Ap_layer_init(b,ilayer,iphase)/(surrogate[i].Ag+surrogate[i].Ap_layer_init(b,ilayer,iphase)); 
			
			  surrogate[i].flux_chem(b,ilayer,iphase,index)+=-flux*fac;
			  surrogate[i].flux_chem_gas(index)+=-flux*(1.0-fac);

			  if (surrogate[surrogate[i].i_irreversible[jmol]].hydrophobic)
			    if (surrogate[i].irr_mass_conserving[jmol])
			      surrogate[surrogate[i].i_irreversible[jmol]].flux_chem(b,ilayer,iphase,index)+=flux;
			    else
			      surrogate[surrogate[i].i_irreversible[jmol]].flux_chem(b,ilayer,iphase,index)+=flux*surrogate[surrogate[i].i_irreversible[jmol]].MM/surrogate[i].MM;
			  else
			    if (surrogate[i].irr_mass_conserving[jmol])
			      surrogate[surrogate[i].i_irreversible[jmol]].flux_chem_aq(b,index)+=flux;
			    else
			      surrogate[surrogate[i].i_irreversible[jmol]].flux_chem_aq(b,index)+=flux*surrogate[surrogate[i].i_irreversible[jmol]].MM/surrogate[i].MM;
			  //ut << "ok2 " << endl;
			  if (surrogate[i].iother_irreversible[jmol]>=0)
			    if (surrogate[surrogate[i].iother_irreversible[jmol]].hydrophobic)
			      surrogate[surrogate[i].iother_irreversible[jmol]].flux_chem(b,ilayer,iphase,index)+=flux*surrogate[surrogate[i].iother_irreversible[jmol]].MM/surrogate[i].MM;
			    else
			      surrogate[surrogate[i].iother_irreversible[jmol]].flux_chem_aq(b,index)+=flux*surrogate[surrogate[i].iother_irreversible[jmol]].MM/surrogate[i].MM;
			  
			}
		  
		  }
            }
  
      for (b=0;b<config.nbins;++b)
        {        
          double Xmonoaq=0.0;
          for (i=0;i<n;i++)
            if (surrogate[i].is_monomer and surrogate[i].is_organic and surrogate[i].hydrophilic)  
              Xmonoaq+=surrogate[i].gamma_aq_bins(b)*surrogate[i].GAMMAinf*surrogate[i].Aaq_bins_init(b)*MMaq(b)/surrogate[i].MM/max(AQinit(b),1.0e-10);          

          compute_conc_org_bins_ssh(config, surrogate, Temperature,  MMaq, AQinit, chp, ionic, LWC, b, conc_org);
          /*
          double conc_org=0.0;//LWC(b);
          for (i=0;i<n;++i)
            if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophilic)
              conc_org+=surrogate[i].Aaq_bins_init(b);

          conc_org=max(conc_org,1.e-5*AQinit(b)); //config.MOmin);
          //conc_org=max(conc_org,config.MOmin);*/
      
          double XH2Oaq=surrogate[config.iH2O].gamma_aq_bins(b)*surrogate[config.iH2O].Aaq_bins_init(b)*MMaq(b)/surrogate[config.iH2O].MM/max(AQinit(b),1.0e-10);
          double xmin=1.0e-10;

          for (i=0;i<n;++i)
            if (surrogate[i].is_organic)
              {
                if (surrogate[i].is_organic and surrogate[i].hydrophilic) 
                  if (surrogate[i].is_monomer)
                    {
                      j=surrogate[i].ioligo;                                        
                      double Kaq2=surrogate[j].MM/surrogate[i].MM*pow(surrogate[i].Keq_oligo,surrogate[i].moligo-1)*pow(max(Xmonoaq,xmin),surrogate[i].moligo-2)/pow(max(XH2Oaq,xmin),surrogate[i].moligo-1);                                            
                      double flux=surrogate[i].koligo*(surrogate[i].gamma_aq_bins(b)*surrogate[i].GAMMAinf*surrogate[i].Aaq_bins_init(b)*Xmonoaq-surrogate[j].gamma_aq_bins(b)*surrogate[j].GAMMAinf*surrogate[j].Aaq_bins_init(b)/Kaq2)*deltat;
		      if (surrogate[i].catalyzed_ph and config.isorropia_ph==false)			
			flux=flux*chp(b)*surrogate[config.iHp].gamma_aq_bins(b)/config.chp_org_ref;
		      else if (surrogate[i].catalyzed_ph and config.isorropia_ph==true)			
			flux=flux*chp(b)/config.chp_org_ref;
                      double fac=1.0;
                      double fac2=1.0;
                      
                      if (surrogate[i].time_aq(b)<config.tequilibrium)                                                  
                        if (surrogate[i].Ag+surrogate[i].Aaq_bins_init(b)>0.0)
                          fac=surrogate[i].Aaq_bins_init(b)/(surrogate[i].Ag+surrogate[i].Aaq_bins_init(b));                                                    
                        
                      if (surrogate[j].time_aq(b)<config.tequilibrium)                        
                        if (surrogate[j].Ag+surrogate[j].Aaq_bins_init(b)>0.0)
                          fac2=surrogate[j].Aaq_bins_init(b)/(surrogate[j].Ag+surrogate[j].Aaq_bins_init(b));

		      /*
			if (flux>0.0)                        
                        flux=flux/(1.0-gamma*surrogate[i].Jdn(b,ilayer,iphase,index)*deltat-gamma*surrogate[i].Jdn_tot*deltat);                                                  
			else
			{
			cout << 
			flux=flux/(1.0-gamma*surrogate[j].Jdn(b,ilayer,iphase,index)*deltat-gamma*surrogate[j].Jdn_tot*deltat);
			}*/

                      if (flux>0.0)                        
                        flux=flux/(1.0-gamma*surrogate[i].Jdn_aq(b,index)*deltat-gamma*surrogate[i].Jdn_tot*deltat);                                                  
                      else
                        flux=flux/(1.0-gamma*surrogate[j].Jdn_aq(b,index)*deltat-gamma*surrogate[j].Jdn_tot*deltat);   
                      
                      surrogate[i].flux_chem_aq(b,index)+=-fac*flux;
                      surrogate[i].flux_chem_gas(index)+=-(1.0-fac)*flux;
                      surrogate[j].flux_chem_aq(b,index)+=fac2*flux;		
                      surrogate[j].flux_chem_gas(index)+=(1.0-fac2)*flux;
                    }

                if (surrogate[i].rion and surrogate[i].Aaq_bins_init(b)>0.0 and config.compute_inorganic)
                  for (jion=0;jion<surrogate[i].nion_chem;jion++)
                    {            
                      double molality=surrogate[surrogate[i].iion(jion)].gamma_aq_bins(b)*surrogate[surrogate[i].iion(jion)].Aaq_bins_init(b)/surrogate[i].MM/conc_org*1000.0;
                      double flux=surrogate[i].kion[jion]*molality*deltat*surrogate[i].GAMMAinf*surrogate[i].gamma_aq_bins(b)*surrogate[i].Aaq_bins_init(b);
		      if (surrogate[i].rion_ph_catalyzed[jion] and config.isorropia_ph==false)
			flux=flux*chp(b)*surrogate[config.iHp].gamma_aq_bins(b);
		      else if (surrogate[i].rion_ph_catalyzed[jion] and config.isorropia_ph==true)
			flux=flux*chp(b);
		      if (surrogate[i].rion_water_catalyzed[jion])
			flux=flux*RH;
		      
                      double fac=1.0;
                      double fac2=1.0;
                      j=surrogate[i].iproduct(jion);
                      if (surrogate[i].time_aq(b)<config.tequilibrium)                                                  
                        if (surrogate[i].Ag+surrogate[i].Aaq_bins_init(b)>0.0)
                          fac=surrogate[i].Aaq_bins_init(b)/(surrogate[i].Ag+surrogate[i].Aaq_bins_init(b)); 

                      if (surrogate[j].time_aq(b)<config.tequilibrium)                        
                        if (surrogate[j].Ag+surrogate[j].Aaq_bins_init(b)>0.0)
                          fac2=surrogate[j].Aaq_bins_init(b)/(surrogate[j].Ag+surrogate[j].Aaq_bins_init(b));                                 

                      if (surrogate[i].rion_catalyzed[jion]==false and molality>0.)
                        {
                          if (surrogate[i].iion(jion)==config.iHSO4m or surrogate[i].iion(jion)==config.iSO4mm)
                            {
                              flux=flux/(1.0-gamma*surrogate[i].Jdn_aq(b,index)*deltat-gamma*(1.0-fac)*surrogate[i].Jdn_gas(index)*deltat-gamma*surrogate[i].Jdn_tot*deltat
                                         -gamma*surrogate[config.iHSO4m].Jdn_aq(b,index)/surrogate[config.iHSO4m].MM*surrogate[i].MM*deltat
                                         -gamma*surrogate[config.iSO4mm].Jdn_aq(b,index)/surrogate[config.iSO4mm].MM*surrogate[i].MM*deltat);
                              double totsulf=surrogate[config.iHSO4m].Aaq_bins_init(b)/surrogate[config.iHSO4m].MM+surrogate[config.iSO4mm].Aaq_bins_init(b)/surrogate[config.iSO4mm].MM;
                              surrogate[config.iHSO4m].flux_chem_aq(b,index)-=surrogate[config.iHSO4m].Aaq_bins_init(b)/totsulf*flux/surrogate[i].MM;
                              surrogate[config.iSO4mm].flux_chem_aq(b,index)-=surrogate[config.iSO4mm].Aaq_bins_init(b)/totsulf*flux/surrogate[i].MM;                               
                            }
                          if (surrogate[i].iion(jion)==config.iNO3m)
                            {
                              if (surrogate[j].time_aq(b)<config.tequilibrium)
                                {
                                  double totno3=surrogate[config.iNO3m].Aaq_bins_init(b)+surrogate[config.iHNO3].Ag/surrogate[config.iHNO3].MM*surrogate[config.iNO3m].MM;
                                  double fac3=surrogate[config.iHNO3].Ag/surrogate[config.iHNO3].MM*surrogate[config.iNO3m].MM/totno3;
                                  flux=flux/(1.0-gamma*surrogate[i].Jdn_aq(b,index)*deltat-gamma*(1.0-fac)*surrogate[i].Jdn_gas(index)*deltat-gamma*surrogate[i].Jdn_tot*deltat
                                             -gamma*surrogate[config.iNO3m].Jdn_aq(b,index)/surrogate[config.iNO3m].MM*surrogate[i].MM*deltat
                                             -gamma*fac3*surrogate[config.iHNO3].Jdn_gas(index)/surrogate[config.iHNO3].MM*surrogate[i].MM*deltat);                             
                                  surrogate[config.iNO3m].flux_chem_aq(b,index)-=surrogate[config.iNO3m].Aaq_bins_init(b)/totno3*flux/surrogate[i].MM*surrogate[config.iNO3m].MM;
                                  surrogate[config.iHNO3].flux_chem_gas(index)-=surrogate[config.iHNO3].Ag/totno3*flux/surrogate[i].MM*surrogate[config.iHNO3].MM;                            
                                }
                              else                          
                                {
                                  flux=flux/(1.0-gamma*surrogate[i].Jdn_aq(b,index)*deltat-gamma*(1.0-fac)*surrogate[i].Jdn_gas(index)*deltat-gamma*surrogate[i].Jdn_tot*deltat
                                             -gamma*surrogate[config.iNO3m].Jdn_aq(b,index)/surrogate[config.iNO3m].MM*surrogate[i].MM*deltat); 
                                  surrogate[config.iNO3m].flux_chem_aq(b,index)-=flux/surrogate[i].MM*surrogate[config.iNO3m].MM;                                                      
                                }
                            }

                          if (surrogate[i].iion(jion)==config.iClm)
                            {
                              if (surrogate[j].time_aq(b)<config.tequilibrium)
                                {
                                  double totcl=surrogate[config.iClm].Aaq_bins_init(b)+surrogate[config.iHCl].Ag/surrogate[config.iHCl].MM*surrogate[config.iClm].MM;
                                  double fac3=surrogate[config.iHCl].Ag/surrogate[config.iHCl].MM*surrogate[config.iClm].MM/totcl;
                                  flux=flux/(1.0-gamma*surrogate[i].Jdn_aq(b,index)*deltat-gamma*(1.0-fac)*surrogate[i].Jdn_gas(index)*deltat-gamma*surrogate[i].Jdn_tot*deltat
                                             -gamma*surrogate[config.iClm].Jdn_aq(b,index)/surrogate[config.iClm].MM*surrogate[i].MM*deltat
                                             -gamma*fac3*surrogate[config.iHCl].Jdn_gas(index)/surrogate[config.iHCl].MM*surrogate[i].MM*deltat);
                                  surrogate[config.iClm].flux_chem_aq(b,index)-=surrogate[config.iClm].Aaq_bins_init(b)/totcl*flux/surrogate[i].MM*surrogate[config.iClm].MM;
                                  surrogate[config.iHCl].flux_chem_gas(index)-=surrogate[config.iHCl].Ag/totcl*flux/surrogate[i].MM*surrogate[config.iHCl].MM;
                                }
                              else
                                {
                                  flux=flux/(1.0-gamma*surrogate[i].Jdn_aq(b,index)*deltat-gamma*(1.0-fac)*surrogate[i].Jdn_gas(index)*deltat-gamma*surrogate[i].Jdn_tot*deltat
                                             -gamma*surrogate[config.iClm].Jdn_aq(b,index)/surrogate[config.iClm].MM*surrogate[i].MM*deltat); 
                                  surrogate[config.iClm].flux_chem_aq(b,index)-=flux/surrogate[i].MM*surrogate[config.iClm].MM;
                                }
                            }
                          if (surrogate[i].iion(jion)==config.iNH4p)
                            {
                              if (surrogate[j].time_aq(b)<config.tequilibrium)
                                {
                                  double totnh4=surrogate[config.iNH4p].Aaq_bins_init(b)+surrogate[config.iNH3].Ag/surrogate[config.iNH3].MM*surrogate[config.iNH4p].MM;
                                  double fac3=surrogate[config.iNH3].Ag/surrogate[config.iNH3].MM*surrogate[config.iNH4p].MM/totnh4;
                                  flux=flux/(1.0-gamma*surrogate[i].Jdn_aq(b,index)*deltat-gamma*(1.0-fac)*surrogate[i].Jdn_gas(index)*deltat-gamma*surrogate[i].Jdn_tot*deltat
                                             -gamma*surrogate[config.iNH4p].Jdn_aq(b,index)/surrogate[config.iNH4p].MM*surrogate[i].MM*deltat
                                             -gamma*fac3*surrogate[config.iNH3].Jdn_gas(index)/surrogate[config.iNH3].MM*surrogate[i].MM*deltat);
                                  surrogate[config.iNH4p].flux_chem_aq(b,index)-=surrogate[config.iNH4p].Aaq_bins_init(n)/totnh4*flux/surrogate[i].MM*surrogate[config.iNH4p].MM;
                                  surrogate[config.iNH3].flux_chem_gas(index)-=surrogate[config.iNH3].Ag/totnh4*flux/surrogate[i].MM*surrogate[config.iNH3].MM;
                                }
                              else
                                {
                                  flux=flux/(1.0-gamma*surrogate[i].Jdn_aq(b,index)*deltat-gamma*(1.0-fac)*surrogate[i].Jdn_gas(index)*deltat-gamma*surrogate[i].Jdn_tot*deltat
                                             -gamma*surrogate[config.iNH4p].Jdn_aq(b,index)/surrogate[config.iNH4p].MM*surrogate[i].MM*deltat); 
                                  surrogate[config.iNH4p].flux_chem_aq(b,index)-=flux/surrogate[i].MM*surrogate[config.iNH4p].MM;
                                }
                            }
                        }                

                      surrogate[i].flux_chem_aq(b,index)+=-fac*flux;
                      surrogate[i].flux_chem_gas(index)+=-(1.0-fac)*flux;
                      surrogate[j].flux_chem_aq(b,index)+=fac2*flux*surrogate[j].MM/surrogate[i].MM;		
                      surrogate[j].flux_chem_gas(index)+=(1.0-fac2)*flux*surrogate[j].MM/surrogate[i].MM;
                    }

		
		if (surrogate[i].is_organic)
		  for (jmol=0;jmol<surrogate[i].n_irreversible;jmol++)
		    if (surrogate[i].i_irreversible[jmol]>=0)
		      {
			double flux=surrogate[i].k_irreversible[jmol]*deltat*surrogate[i].GAMMAinf*surrogate[i].gamma_aq_bins(b)*surrogate[i].Aaq_bins_init(b);
			if (surrogate[i].irr_catalyzed_water[jmol])
			  flux=flux*RH;
			if (surrogate[i].irr_catalyzed_pH[jmol] and config.isorropia_ph==false)
			  flux=flux*chp(b)*surrogate[config.iHp].gamma_aq_bins(b);
			else if (surrogate[i].irr_catalyzed_pH[jmol] and config.isorropia_ph==true)
			  flux=flux*chp(b);
		      
			double fac=1.0;
			if (surrogate[i].time_aq(b)<config.tequilibrium)                                                  
			  if (surrogate[i].Ag+surrogate[i].Aaq_bins_init(b)>0.0)
			    fac=surrogate[i].Aaq_bins_init(b)/(surrogate[i].Ag+surrogate[i].Aaq_bins_init(b)); 
		  
			surrogate[i].flux_chem_aq(b,index)+=-flux*fac;
			surrogate[i].flux_chem_gas(index)+=-flux*(1.0-fac);
			if (surrogate[i].irr_mass_conserving[jmol]==false)
			  flux=flux*surrogate[surrogate[i].i_irreversible[jmol]].MM/surrogate[i].MM;
			if (surrogate[surrogate[i].i_irreversible[jmol]].hydrophilic)
			  surrogate[surrogate[i].i_irreversible[jmol]].flux_chem_aq(b,index)+=flux;
			else
			  for (ilayer=0;ilayer<config.nlayer;ilayer++)
			    surrogate[surrogate[i].i_irreversible[jmol]].flux_chem(b,ilayer,0,index)+=flux*config.Vlayer(ilayer);

			//cout << "ok3 " << endl;
			if (surrogate[i].iother_irreversible[jmol]>=0)
			  {
			    flux=flux*surrogate[surrogate[i].iother_irreversible[jmol]].MM/surrogate[surrogate[i].i_irreversible[jmol]].MM;
			    if (surrogate[surrogate[i].iother_irreversible[jmol]].hydrophilic)
			      surrogate[surrogate[i].iother_irreversible[jmol]].flux_chem_aq(b,index)+=flux;
			    else
			      for (ilayer=0;ilayer<config.nlayer;ilayer++)
				surrogate[surrogate[i].iother_irreversible[jmol]].flux_chem(b,ilayer,0,index)+=flux*config.Vlayer(ilayer);
			  }
		      }
    
              }
        } 
    }    
}



void prodloss_chem_ssh(model_config &config, vector<species>& surrogate,
                       Array<double, 3> &MOinit, Array<double, 3> &MOW, 
                       Array<double, 1> &AQinit, Array<double, 1> &LWC, Array<double, 1> &MMaq, 
                       Array<double, 1> &chp, Array<double, 1> &ionic, double &tiny, double &Temperature, int index, double &RH)
{    
  int n=surrogate.size();
  int i,j,b,ilayer,iphase,jion,jmol; 
  double Xmono=0.0;          
  double XH2O;
  double xmin=1.0e-10;
  double conc_org;

  if (config.chemistry and config.compute_organic)
    {
      for (b=0;b<config.nbins;++b)
        for (ilayer=0;ilayer<config.nlayer;ilayer++)
          for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
            {        
              Xmono=0.0;          
              for (i=0;i<n;i++)
                if (surrogate[i].is_monomer and surrogate[i].is_organic)              
                  Xmono+=surrogate[i].gamma_org_layer(b,ilayer,iphase)*surrogate[i].Ap_layer_init(b,ilayer,iphase)*MOW(b,ilayer,iphase)/surrogate[i].MM/max(MOinit(b,ilayer,iphase),1.0e-10);              
              XH2O=surrogate[config.iH2O].gamma_org_layer(b,ilayer,iphase)*surrogate[config.iH2O].Ap_layer_init(b,ilayer,iphase)*MOW(b,ilayer,iphase)/surrogate[config.iH2O].MM/max(MOinit(b,ilayer,iphase),1.0e-10);
              xmin=1.0e-10;

              for (i=0;i<n;++i)
                if (surrogate[i].is_organic and surrogate[i].hydrophobic)
		  {
		    if (surrogate[i].is_monomer)
		      {
			j=surrogate[i].ioligo;                       
			double Keq2=surrogate[j].MM/surrogate[i].MM*pow(surrogate[i].Keq_oligo,surrogate[i].moligo-1)*pow(max(Xmono,xmin),surrogate[i].moligo-2)/pow(max(XH2O,xmin),surrogate[i].moligo-1);                      
			double flux1=surrogate[i].koligo*surrogate[i].gamma_org_layer(b,ilayer,iphase)*surrogate[i].Ap_layer_init(b,ilayer,iphase)*Xmono;
			double flux2=surrogate[i].koligo*surrogate[j].gamma_org_layer(b,ilayer,iphase)*surrogate[j].Ap_layer_init(b,ilayer,iphase)/Keq2;
			
			surrogate[i].kprod(b,ilayer,iphase)+=flux2;
			surrogate[j].kprod(b,ilayer,iphase)+=flux1;
			surrogate[i].kloss(b,ilayer,iphase)+=surrogate[i].koligo*surrogate[i].gamma_org_layer(b,ilayer,iphase)*Xmono;
			surrogate[j].kloss(b,ilayer,iphase)+=surrogate[i].koligo*surrogate[j].gamma_org_layer(b,ilayer,iphase)/Keq2;                      
		      }

		    
		    for (jmol=0;jmol<surrogate[i].n_irreversible;jmol++)
		      if (surrogate[i].i_irreversible[jmol]>=0)
			{
			  double flux=surrogate[i].k_irreversible[jmol]*surrogate[i].gamma_org_layer(b,ilayer,iphase);
			  if (surrogate[i].irr_catalyzed_water[jmol])
			    flux=flux*RH;
			  if (surrogate[i].irr_catalyzed_pH[jmol])
			    flux=flux*config.chp_org_ref;
			
			  surrogate[i].kloss(b,ilayer,iphase)+=flux;
			  if (surrogate[i].irr_mass_conserving[jmol]==false)
			    {
			      flux=flux*surrogate[surrogate[i].i_irreversible[jmol]].MM/surrogate[i].MM;
			    }
			  if (surrogate[surrogate[i].i_irreversible[jmol]].hydrophobic)
			    surrogate[surrogate[i].i_irreversible[jmol]].kprod(b,ilayer,iphase)+=flux*surrogate[i].Ap_layer_init(b,ilayer,iphase);
			  else
			    surrogate[surrogate[i].i_irreversible[jmol]].kprod_aq(b)+=flux*surrogate[i].Ap_layer_init(b,ilayer,iphase);

			  if (surrogate[i].iother_irreversible[jmol]>=0)
			    {
			      flux=flux*surrogate[surrogate[i].iother_irreversible[jmol]].MM/surrogate[surrogate[i].i_irreversible[jmol]].MM;
			      if (surrogate[surrogate[i].iother_irreversible[jmol]].hydrophobic)
				surrogate[surrogate[i].iother_irreversible[jmol]].kprod(b,ilayer,iphase)+=flux*surrogate[i].Ap_layer_init(b,ilayer,iphase);
			      else
				{
				  surrogate[surrogate[i].iother_irreversible[jmol]].kprod_aq(b)+=flux*surrogate[i].Ap_layer_init(b,ilayer,iphase);
				}
			    }
			}
		  
		  }
            }
  
      for (b=0;b<config.nbins;++b)
        {        
          double Xmonoaq=0.0;
          for (i=0;i<n;i++)
            if (surrogate[i].is_monomer and surrogate[i].is_organic and surrogate[i].hydrophilic)  
              Xmonoaq+=surrogate[i].gamma_aq_bins(b)*surrogate[i].GAMMAinf*surrogate[i].Aaq_bins_init(b)*MMaq(b)/surrogate[i].MM/max(AQinit(b),1.0e-10);          

          compute_conc_org_bins_ssh(config, surrogate, Temperature,  MMaq, AQinit, chp, ionic, LWC, b, conc_org);
          /*
          double conc_org=0.0;//LWC(b);
          for (i=0;i<n;++i)
            if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophilic)
              conc_org+=surrogate[i].Aaq_bins_init(b);

          conc_org=max(conc_org,1.e-5*AQinit(b)); //config.MOmin);
          //conc_org=max(conc_org,config.MOmin);*/
      
          double XH2Oaq=surrogate[config.iH2O].gamma_aq_bins(b)*surrogate[config.iH2O].Aaq_bins_init(b)*MMaq(b)/surrogate[config.iH2O].MM/max(AQinit(b),1.0e-10);
          double xmin=1.0e-10;

          for (i=0;i<n;++i)
            {
              if (surrogate[i].is_organic)
                if (surrogate[i].is_organic and surrogate[i].hydrophilic)
		  {
		    if (surrogate[i].is_monomer)
		      {
			j=surrogate[i].ioligo;                                        
			double Kaq2=surrogate[j].MM/surrogate[i].MM*pow(surrogate[i].Keq_oligo,surrogate[i].moligo-1)*pow(max(Xmonoaq,xmin),surrogate[i].moligo-2)/pow(max(XH2Oaq,xmin),surrogate[i].moligo-1);                                            
			double flux1=surrogate[i].koligo*surrogate[i].gamma_aq_bins(b)*surrogate[i].GAMMAinf*surrogate[i].Aaq_bins_init(b)*Xmonoaq;
			double flux2=surrogate[i].koligo*surrogate[j].gamma_aq_bins(b)*surrogate[j].GAMMAinf*surrogate[j].Aaq_bins_init(b)/Kaq2;			
			if (surrogate[i].catalyzed_ph and config.isorropia_ph==false)
			  {
			    flux1=flux1*chp(b)*surrogate[config.iHp].gamma_aq_bins(b)/config.chp_org_ref;
			    flux2=flux2*chp(b)*surrogate[config.iHp].gamma_aq_bins(b)/config.chp_org_ref;
			  }
			else if (surrogate[i].catalyzed_ph and config.isorropia_ph==true)
			  {
			    flux1=flux1*chp(b)/config.chp_org_ref;
			    flux2=flux2*chp(b)/config.chp_org_ref;
			  }
			surrogate[i].kprod_aq(b)+=flux2;
			surrogate[j].kprod_aq(b)+=flux1;
			surrogate[i].kloss_aq(b)+=surrogate[i].koligo*surrogate[i].gamma_aq_bins(b)*surrogate[i].GAMMAinf*Xmonoaq;
			surrogate[j].kloss_aq(b)+=surrogate[i].koligo*surrogate[j].gamma_aq_bins(b)*surrogate[j].GAMMAinf/Kaq2;                                            
		      }

		    for (jmol=0;jmol<surrogate[i].n_irreversible;jmol++)
		      if (surrogate[i].i_irreversible[jmol]>=0)
			{
			  double flux=surrogate[i].k_irreversible[jmol]*surrogate[i].gamma_aq_bins(b)*surrogate[i].GAMMAinf;
			  if (surrogate[i].irr_catalyzed_water[jmol])
			    flux=flux*RH;
			  if (surrogate[i].irr_catalyzed_pH[jmol] and config.isorropia_ph==false)
			    flux=flux*chp(b)*surrogate[config.iHp].gamma_aq_bins(b);
			  else if (surrogate[i].irr_catalyzed_pH[jmol] and config.isorropia_ph==true)
			    flux=flux*chp(b);
			
			  surrogate[i].kloss_aq(b)+=flux;
			  if (surrogate[i].irr_mass_conserving[jmol]==false)
			    flux=flux*surrogate[surrogate[i].i_irreversible[jmol]].MM/surrogate[i].MM;
			  if (surrogate[surrogate[i].i_irreversible[jmol]].hydrophilic)
			    surrogate[surrogate[i].i_irreversible[jmol]].kprod_aq(b)+=flux*surrogate[i].Aaq_bins_init(b);
			  else
			    for (ilayer=0;ilayer<config.nlayer;ilayer++)
			      surrogate[surrogate[i].i_irreversible[jmol]].kprod(b,ilayer,0)+=flux*surrogate[i].Aaq_bins_init(b)*config.Vlayer(ilayer);
			
			  if (surrogate[i].iother_irreversible[jmol]>=0)
			    {
			      ///
			        
			      if (surrogate[i].irr_mass_conserving[jmol]==false)
				flux=flux*surrogate[surrogate[i].iother_irreversible[jmol]].MM/surrogate[surrogate[i].i_irreversible[jmol]].MM;
			      else
				flux=flux*surrogate[surrogate[i].iother_irreversible[jmol]].MM/surrogate[i].MM;
			      if (surrogate[surrogate[i].iother_irreversible[jmol]].hydrophilic)
				surrogate[surrogate[i].iother_irreversible[jmol]].kprod_aq(b)+=flux*surrogate[i].Aaq_bins_init(b);
			      else
				for (ilayer=0;ilayer<config.nlayer;ilayer++)
				  surrogate[surrogate[i].iother_irreversible[jmol]].kprod(b,ilayer,0)+=flux*surrogate[i].Aaq_bins_init(b)*config.Vlayer(ilayer);
			    }
			}
		  }

	      
	      //cout << surrogate[i].name << " " << surrogate[i].rion << " " << surrogate[i].Aaq_bins_init(b) << " " << surrogate[i].nion_chem << " " << config.compute_inorganic << endl;
              if (surrogate[i].rion and surrogate[i].Aaq_bins_init(b)>0.0) // and config.compute_inorganic)
                for (jion=0;jion<surrogate[i].nion_chem;jion++)
                  {            
                    double molality=surrogate[surrogate[i].iion(jion)].gamma_aq_bins(b)*surrogate[surrogate[i].iion(jion)].Aaq_bins_init(b)/surrogate[i].MM/conc_org*1000.0;
                    double flux=surrogate[i].kion[jion]*molality*surrogate[i].GAMMAinf*surrogate[i].gamma_aq_bins(b);
		    if (surrogate[i].rion_ph_catalyzed[jion] and config.isorropia_ph==false)
		      flux=flux*chp(b)*surrogate[config.iHp].gamma_aq_bins(b);
		    else if (surrogate[i].rion_ph_catalyzed[jion] and config.isorropia_ph==true)
		      flux=flux*chp(b);
		      
		    if (surrogate[i].rion_water_catalyzed[jion])
		      flux=flux*RH;
		    j=surrogate[i].iproduct(jion);
		    surrogate[i].kloss_aq(b)+=flux;
                    if (surrogate[i].rion_catalyzed[jion]==false and molality>0.)
                      {                                                  
			surrogate[surrogate[i].iion(jion)].kloss_aq(b)+=flux*surrogate[i].Aaq_bins_init(b)/surrogate[i].MM*surrogate[surrogate[i].iion(jion)].MM/surrogate[surrogate[i].iion(jion)].Aaq_bins_init(b);
                      }                
                    surrogate[j].kprod_aq(b)+=flux*surrogate[i].Aaq_bins_init(b)*surrogate[j].MM/surrogate[i].MM;
                  }
            }    
        }
    }
}



void prodloss_aq_ssh(model_config &config, vector<species>& surrogate, Array<double, 1> &AQinit,  Array<double, 1> &LWC, Array<double, 1> &MMaq, Array<double, 1> &chp, Array<double, 1> &ionic, Array<double, 3> &MOinit, int index, double &Temperature, double deltat)
{
  //compute kinetic rates for the absorption of a compound in the aqueous phase in a bin
  //index = 0 : first evaluation of rates
  //index = 1 : second evaluation of rates
  double Kaq;
  int n=surrogate.size();
  int i,b,ilayer,iphase;
  double sum_mass;  
  double conc_aq;
  double conc_org;
  double sumkpositive=0.0;
  double sumknegative=0.0;

  for (i=0;i<n;++i)
    {
      surrogate[i].kloss_gas=0.0;
      surrogate[i].kprod_gas=0.0;
      for (b=0;b<config.nbins;b++)
	{
	  surrogate[i].kprod_aq(b)=0.0;
	  surrogate[i].kloss_aq(b)=0.0;
	}	
    }
      

  if (sum(LWC)>config.LWClimit)
    for (b=0;b<config.nbins;++b)
      {
        compute_conc_org_bins_ssh(config, surrogate, Temperature,  MMaq, AQinit, chp, ionic, LWC, b, conc_org);

        sum_mass=AQinit(b);
        for (ilayer=0;ilayer<config.nlayer;ilayer++)
          for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
            sum_mass+=MOinit(b,ilayer,iphase);

        for (i=0;i<n;++i)
          if(surrogate[i].is_organic and surrogate[i].hydrophilic and config.compute_organic)
            {
	      double tau_airloc=surrogate[i].tau_air(b);
	      if (surrogate[i].time_aq(b)<config.tequilibrium)
		tau_airloc=surrogate[i].tau_air(b)*config.tequilibrium/surrogate[i].time_aq(b);
	      
              Kaq=surrogate[i].Kaq(b);          
              if (Kaq > 1.e-20 and AQinit(b)> 1e-20 and sum_mass > 1e-20 and tau_airloc > 0.)
                {
                  surrogate[i].kprod_aq(b)=surrogate[i].Ag*AQinit(b)/sum_mass/tau_airloc;
                  surrogate[i].kloss_aq(b)=AQinit(b)/sum_mass/(Kaq*AQinit(b)*tau_airloc);
                  surrogate[i].kprod_gas+=AQinit(b)/sum_mass*surrogate[i].Aaq_bins_init(b)/(Kaq*AQinit(b)*tau_airloc);
                  surrogate[i].kloss_gas+=AQinit(b)/sum_mass/tau_airloc;
                }
            }
          else if (surrogate[i].is_inorganic_precursor and config.compute_inorganic and surrogate[i].is_solid==false) //for inorganic compounds      
            {	     
              if (i==config.iH2SO4)
                {		
                  surrogate[config.iSO4mm].kprod_aq(b)=surrogate[i].Ag/surrogate[i].tau_air(b)*surrogate[config.iSO4mm].MM/surrogate[i].MM;
                  surrogate[config.iSO4mm].kloss_aq(b)=0.0;
                  surrogate[i].kprod_gas+=0.0;
                  surrogate[i].kloss_gas+=1.0/surrogate[i].tau_air(b);		
                }
              else if (i==config.iNH3)
                {
		  double tau_airloc=surrogate[i].tau_air(b);
		  if (surrogate[i].time_aq(b)<config.tequilibrium)
		    tau_airloc=surrogate[i].tau_air(b)*config.tequilibrium/surrogate[i].time_aq(b);
                  Kaq=surrogate[i].Kaq(b);
                  //compute kinetic rate of absorption
                  conc_aq=surrogate[config.iNH4p].Aaq_bins_init(b)/surrogate[config.iNH4p].MM*surrogate[i].MM;
		  double fcorr;
		  /*if (1.0/(Kaq*conc_org*tau_airloc)*AQinit(b)/sum_mass>config.tequilibrium)
		    {
		      fcorr=(Kaq*conc_org*tau_airloc)/AQinit(b)*sum_mass;
		    }
		    else*/
		    fcorr=1.;

                  surrogate[config.iNH4p].kprod_aq(b)=surrogate[i].Ag/tau_airloc*AQinit(b)/sum_mass*surrogate[config.iNH4p].MM/surrogate[i].MM*fcorr;
                  surrogate[config.iNH4p].kloss_aq(b)=1.0/(Kaq*conc_org*tau_airloc)*AQinit(b)/sum_mass*fcorr;
                  surrogate[i].kprod_gas+=conc_aq/(Kaq*conc_org*tau_airloc)*AQinit(b)/sum_mass*fcorr;
                  surrogate[i].kloss_gas+=1.0/tau_airloc*AQinit(b)/sum_mass*fcorr;

                }
              else if (i==config.iHNO3)
                {
		  double tau_airloc=surrogate[i].tau_air(b);
		  if (surrogate[i].time_aq(b)<config.tequilibrium)
		    tau_airloc=surrogate[i].tau_air(b)*config.tequilibrium/surrogate[i].time_aq(b);
                  Kaq=surrogate[i].Kaq(b);
		  conc_aq=surrogate[config.iNO3m].Aaq_bins_init(b)/surrogate[config.iNO3m].MM*surrogate[i].MM;
                  //compute kinetic rate of absorption
		  double fcorr;
		  /*if (1.0/(Kaq*conc_org*tau_airloc)*AQinit(b)/sum_mass>config.tequilibrium)
		    {
		      fcorr=(Kaq*conc_org*tau_airloc)/AQinit(b)*sum_mass;
		    }
		  else*/
		    fcorr=1.;
		  surrogate[config.iNO3m].kprod_aq(b)=surrogate[i].Ag/tau_airloc*AQinit(b)/sum_mass*surrogate[config.iNO3m].MM/surrogate[i].MM*fcorr;
		  surrogate[config.iNO3m].kloss_aq(b)=1.0/(Kaq*conc_org*tau_airloc)*AQinit(b)/sum_mass*fcorr;
		  surrogate[i].kprod_gas+=conc_aq/(Kaq*conc_org*tau_airloc)*AQinit(b)/sum_mass*fcorr;
		  surrogate[i].kloss_gas+=1.0/tau_airloc*AQinit(b)/sum_mass*fcorr;
                }
	      else if (i==config.iCO2)
                {
		  double tau_airloc=surrogate[i].tau_air(b);
		  if (surrogate[i].time_aq(b)<config.tequilibrium)
		    tau_airloc=surrogate[i].tau_air(b)*config.tequilibrium/surrogate[i].time_aq(b);
		  //cout << b << " 1 " << surrogate[i].tau_air(b) << " " << tau_airloc << " " << surrogate[i].time_aq(b) << endl;
                  Kaq=surrogate[i].Kaq(b);
		  //cout << "Kaq " << surrogate[i].Kaq(b) << endl;
                  //compute kinetic rate of absorption
                  //conc_aq=surrogate[config.iNO3m].Aaq_bins_init(b)/surrogate[config.iNO3m].MM*surrogate[i].MM;
		  double frac1=1.0/(1.0+(surrogate[config.iCO2].keq2/(chp(b)*surrogate[config.iHp].gamma_aq_bins(b)*surrogate[config.iCO3mm].gamma_aq_bins(b)/surrogate[config.iHCO3m].gamma_aq_bins(b))));
                  surrogate[config.iHCO3m].kprod_aq(b)=frac1*surrogate[i].Ag/tau_airloc*AQinit(b)/sum_mass*surrogate[config.iHCO3m].MM/surrogate[i].MM;
		  surrogate[config.iCO3mm].kprod_aq(b)=(1.0-frac1)*surrogate[i].Ag/tau_airloc*AQinit(b)/sum_mass*surrogate[config.iCO3mm].MM/surrogate[i].MM;
                  surrogate[config.iHCO3m].kloss_aq(b)=1.0/(Kaq*conc_org*tau_airloc)*AQinit(b)/sum_mass;
		  surrogate[config.iCO3mm].kloss_aq(b)=1.0/(Kaq*conc_org*tau_airloc)*AQinit(b)/sum_mass;
                  //surrogate[i].kprod_gas+=(surrogate[config.iHCO3m].Aaq_bins_init(b)/surrogate[config.iHCO3m].MM+surrogate[config.iCO3mm].Aaq_bins_init(b)/surrogate[config.iCO3mm].MM)*surrogate[i].MM/(Kaq*conc_org*surrogate[i].tau_air(b))*AQinit(b)/sum_mass;
                  //surrogate[i].kloss_gas+=1.0/surrogate[i].tau_air(b)*AQinit(b)/sum_mass;
                }
              else if (i==config.iHCl)
                {
		  double tau_airloc=surrogate[i].tau_air(b);
		  if (surrogate[i].time_aq(b)<config.tequilibrium)
		    tau_airloc=surrogate[i].tau_air(b)*config.tequilibrium/surrogate[i].time_aq(b);
                  Kaq=surrogate[i].Kaq(b);
                  //compute kinetic rate of absorption
                  conc_aq=surrogate[config.iClm].Aaq_bins_init(b)/surrogate[config.iClm].MM*surrogate[i].MM;

                  surrogate[config.iClm].kprod_aq(b)=surrogate[i].Ag/tau_airloc*AQinit(b)/sum_mass*surrogate[i].fac_corr_ph(b)*surrogate[config.iClm].MM/surrogate[i].MM;
                  surrogate[config.iClm].kloss_aq(b)=1.0/(Kaq*conc_org*tau_airloc)*AQinit(b)/sum_mass;
                  surrogate[i].kprod_gas+=conc_aq/(Kaq*conc_org*tau_airloc)*AQinit(b)/sum_mass;
                  surrogate[i].kloss_gas+=1.0/tau_airloc*AQinit(b)/sum_mass;		
                }


          
            }
      }
}




void prodloss_aq_bins_ssh(model_config &config, vector<species>& surrogate, Array<double, 1> &AQinit,  Array<double, 1> &LWC, Array<double, 1> &MMaq, Array<double, 1> &chp, Array<double, 1> &ionic, Array<double, 3> &MOinit, double &conc_org, int index, double &Temperature, double deltat, int b)
{
  //compute kinetic rates for the absorption of a compound in the aqueous phase in a bin
  //index = 0 : first evaluation of rates
  //index = 1 : second evaluation of rates
  double Kaq;
  int n=surrogate.size();
  int i,ilayer,iphase;
  double sum_mass;  
  double conc_aq;
  double sumkpositive=0.0;
  double sumknegative=0.0;

  //cout << surrogate[config.iCO2].time_aq << endl;
  if (sum(LWC)>config.LWClimit)   
    {
      //compute_conc_org_bins_ssh(config, surrogate, Temperature,  MMaq, AQinit, chp, ionic, LWC, b, conc_org);
      
      sum_mass=AQinit(b);
      for (ilayer=0;ilayer<config.nlayer;ilayer++)
        for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
          sum_mass+=MOinit(b,ilayer,iphase);
      
      for (i=0;i<n;++i)
        if (surrogate[i].is_inorganic_precursor and config.compute_inorganic and surrogate[i].is_solid==false) //for inorganic compounds      
          {                  
            if (i==config.iNH3)
              {
		double tau_airloc=surrogate[i].tau_air(b);
		if (surrogate[i].time_aq(b)<config.tequilibrium)
		  tau_airloc=surrogate[i].tau_air(b)*config.tequilibrium/surrogate[i].time_aq(b);
		
                Kaq=surrogate[i].Kaq(b);
                //compute kinetic rate of absorption
                conc_aq=surrogate[config.iNH4p].Aaq_bins_init(b)/surrogate[config.iNH4p].MM*surrogate[i].MM;
                
                surrogate[config.iNH4p].kprod_aq(b)=surrogate[i].Ag/tau_airloc*AQinit(b)/sum_mass*surrogate[config.iNH4p].MM/surrogate[i].MM;
                surrogate[config.iNH4p].kloss_aq(b)=1.0/(Kaq*conc_org*tau_airloc)*AQinit(b)/sum_mass;
              }
            else if (i==config.iHNO3)
              {
		double tau_airloc=surrogate[i].tau_air(b);
		if (surrogate[i].time_aq(b)<config.tequilibrium)
		  tau_airloc=surrogate[i].tau_air(b)*config.tequilibrium/surrogate[i].time_aq(b);
		
                Kaq=surrogate[i].Kaq(b);
                //compute kinetic rate of absorption
                conc_aq=surrogate[config.iNO3m].Aaq_bins_init(b)/surrogate[config.iNO3m].MM*surrogate[i].MM;
                surrogate[config.iNO3m].kprod_aq(b)=surrogate[i].Ag/tau_airloc*AQinit(b)/sum_mass*surrogate[config.iNO3m].MM/surrogate[i].MM;
                surrogate[config.iNO3m].kloss_aq(b)=1.0/(Kaq*conc_org*tau_airloc)*AQinit(b)/sum_mass;
              }
	    else if (i==config.iCO2)
              {
		double tau_airloc=surrogate[i].tau_air(b);
		if (surrogate[i].time_aq(b)<config.tequilibrium)
		  tau_airloc=surrogate[i].tau_air(b)*config.tequilibrium/surrogate[i].time_aq(b);
		if (tau_airloc==0.)
		  cout << tau_airloc << " is zero" << endl;
		  //cout << b << " 2 " << surrogate[i].tau_air(b) << " " << tau_airloc << " " << surrogate[i].time_aq(b) << endl;
		
		//conc_aq=surrogate[config.iNO3m].Aaq_bins_init(b)/surrogate[config.iNO3m].MM*surrogate[i].MM;
		Kaq=surrogate[i].Kaq(b);
		double frac1=1.0/(1.0+(surrogate[config.iCO2].keq2/(chp(b)*surrogate[config.iHp].gamma_aq_bins(b)*surrogate[config.iCO3mm].gamma_aq_bins(b)/surrogate[config.iHCO3m].gamma_aq_bins(b))));
		surrogate[config.iHCO3m].kprod_aq(b)=frac1*surrogate[i].Ag/tau_airloc*AQinit(b)/sum_mass*surrogate[config.iHCO3m].MM/surrogate[i].MM;
		surrogate[config.iCO3mm].kprod_aq(b)=(1.0-frac1)*surrogate[i].Ag/tau_airloc*AQinit(b)/sum_mass*surrogate[config.iCO3mm].MM/surrogate[i].MM;
		surrogate[config.iHCO3m].kloss_aq(b)=1.0/(Kaq*conc_org*tau_airloc)*AQinit(b)/sum_mass;
		surrogate[config.iCO3mm].kloss_aq(b)=1.0/(Kaq*conc_org*tau_airloc)*AQinit(b)/sum_mass;
		surrogate[i].kprod_gas+=(surrogate[config.iHCO3m].Aaq_bins_init(b)/surrogate[config.iHCO3m].MM+surrogate[config.iCO3mm].Aaq_bins_init(b)/surrogate[config.iCO3mm].MM)*surrogate[i].MM/(Kaq*conc_org*surrogate[i].tau_air(b))*AQinit(b)/sum_mass;
		surrogate[i].kloss_gas+=1.0/surrogate[i].tau_air(b)*AQinit(b)/sum_mass;
	      }
            else if (i==config.iHCl)
              {
		double tau_airloc=surrogate[i].tau_air(b);
		if (surrogate[i].time_aq(b)<config.tequilibrium)
		  tau_airloc=surrogate[i].tau_air(b)*config.tequilibrium/surrogate[i].time_aq(b);
		
                Kaq=surrogate[i].Kaq(b);
                //compute kinetic rate of absorption
                conc_aq=surrogate[config.iClm].Aaq_bins_init(b)/surrogate[config.iClm].MM*surrogate[i].MM;
                
                surrogate[config.iClm].kprod_aq(b)=surrogate[i].Ag/tau_airloc*AQinit(b)/sum_mass*surrogate[i].fac_corr_ph(b)*surrogate[config.iClm].MM/surrogate[i].MM;
                surrogate[config.iClm].kloss_aq(b)=1.0/(Kaq*conc_org*tau_airloc)*AQinit(b)/sum_mass;		
              }
          
          }
    }
}


void flux_aq_ssh(model_config &config, vector<species>& surrogate, Array<double, 1> &AQinit,  Array<double, 1> &LWC, Array<double, 1> &MMaq, Array<double, 1> &chp, Array<double, 1> &ionic, Array<double, 3> &MOinit, double &tiny, double &Temperature, int index)
{
  //compute kinetic rates for the absorption of a compound in the aqueous phase in a bin
  //index = 0 : first evaluation of rates
  //index = 1 : second evaluation of rates
  double Kaq;
  int n=surrogate.size();
  int i,b,ilayer,iphase;
  double sum_mass;  
  double conc_aq;
  double conc_org;
  double sumkpositive=0.0;
  double sumknegative=0.0;

  for (i=0;i<n;++i)      		
    for (b=0;b<config.nbins;b++)
      {
	surrogate[i].k1_aq(b,index)=0.0;
	surrogate[i].Jdn_aq(b,index)=0.0;
      }

  for (b=0;b<config.nbins;++b)
    {
      compute_conc_org_bins_ssh(config, surrogate, Temperature,  MMaq, AQinit, chp, ionic, LWC, b, conc_org);

      for (i=0;i<n;++i)
        if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophilic)
          {
            sum_mass=AQinit(b);
            for (ilayer=0;ilayer<config.nlayer;ilayer++)
              for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                sum_mass+=MOinit(b,ilayer,iphase);
          
            Kaq=surrogate[i].Kaq(b);          
	    if (Kaq > 1.e-20 and AQinit(b)> 1e-20 and sum_mass > 1e-20 and surrogate[i].tau_air(b) > 0.)
	      {
		surrogate[i].k1_aq(b,index)=AQinit(b)/sum_mass*
		  (surrogate[i].Ag*Kaq*AQinit(b)-surrogate[i].Aaq_bins_init(b))/
		  (Kaq*AQinit(b)*surrogate[i].tau_air(b));
		  
		surrogate[i].Jdn_aq(b,index)=0.0;
		if (surrogate[i].k1_aq(b,index)<0.0 and surrogate[i].Aaq_bins_init(b)>tiny)
		  surrogate[i].Jdn_aq(b,index)=surrogate[i].k1_aq(b,index)/
		    surrogate[i].Aaq_bins_init(b);
	      }
          }
        else if (surrogate[i].is_inorganic_precursor and config.compute_inorganic and surrogate[i].is_solid==false) //for inorganic compounds      
          {
            sum_mass=AQinit(b);
            for (ilayer=0;ilayer<config.nlayer;ilayer++)
              for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		sum_mass+=MOinit(b,ilayer,iphase);       

            conc_aq;
            if (i==config.iH2SO4)
              {	      
                surrogate[i].k1_aq(b,index)=surrogate[i].Ag/surrogate[i].tau_air(b);      
                surrogate[i].Jdn_aq(b,index)=0.0;
              }
            else if (i==config.iNH3)
              {
                Kaq=surrogate[i].Kaq(b);
                //compute kinetic rate of absorption
                conc_aq=surrogate[config.iNH4p].Aaq_bins_init(b)/surrogate[config.iNH4p].MM*surrogate[i].MM;
	      
                surrogate[i].k1_aq(b,index)=(surrogate[i].Ag*Kaq*conc_org-conc_aq)/(Kaq*conc_org*surrogate[i].tau_air(b))*AQinit(b)/sum_mass*surrogate[i].fac_corr_ph(b);	   
 
                surrogate[i].Jdn_aq(b,index)=0.0;
                if (surrogate[i].k1_aq(b,index)<0.0 and surrogate[config.iNH4p].Aaq_bins_init(b)>tiny)
                  surrogate[i].Jdn_aq(b,index)=surrogate[i].k1_aq(b,index)/surrogate[config.iNH4p].Aaq_bins_init(b)*surrogate[config.iNH4p].MM/surrogate[i].MM;
              }
            else if (i==config.iHNO3)
              {
                Kaq=surrogate[i].Kaq(b);
                //compute kinetic rate of absorption
                conc_aq=surrogate[config.iNO3m].Aaq_bins_init(b)/surrogate[config.iNO3m].MM*surrogate[i].MM;
	      
                surrogate[i].k1_aq(b,index)=(surrogate[i].Ag*Kaq*conc_org-conc_aq)/(Kaq*conc_org*surrogate[i].tau_air(b))*AQinit(b)/sum_mass*surrogate[i].fac_corr_ph(b);     
              
                surrogate[i].Jdn_aq(b,index)=0.0;
                if (surrogate[i].k1_aq(b,index)<0.0 and surrogate[config.iNO3m].Aaq_bins_init(b)>tiny)
                  surrogate[i].Jdn_aq(b,index)=surrogate[i].k1_aq(b,index)/surrogate[config.iNO3m].Aaq_bins_init(b)*surrogate[config.iNO3m].MM/surrogate[i].MM;
              }
	    else if (i==config.iCO2)
              {
		//conc_aq=surrogate[config.iNO3m].Aaq_bins_init(b)/surrogate[config.iNO3m].MM*surrogate[i].MM;
		/*double frac1=1.0/(1.0+(surrogate[config.iCO2].keq2/(chp(b)*surrogate[config.iHp].gamma_aq_bins(b)*surrogate[config.iCO3mm].gamma_aq_bins(b)/surrogate[config.iHCO3m].gamma_aq_bins(b))));
		surrogate[config.iHCO3m].kprod_aq(b)=frac1*surrogate[i].Ag/surrogate[i].tau_air(b)*AQinit(b)/sum_mass*surrogate[config.iHCO3m].MM/surrogate[i].MM;
		surrogate[config.iCO3mm].kprod_aq(b)=(1.0-frac1)*surrogate[i].Ag/surrogate[i].tau_air(b)*AQinit(b)/sum_mass*surrogate[config.iCO3mm].MM/surrogate[i].MM;
		surrogate[config.iHCO3m].kloss_aq(b)=1.0/(Kaq*conc_org*surrogate[i].tau_air(b))*AQinit(b)/sum_mass;
		surrogate[config.iCO3mm].kloss_aq(b)=1.0/(Kaq*conc_org*surrogate[i].tau_air(b))*AQinit(b)/sum_mass;
		surrogate[i].kprod_gas+=(surrogate[config.iHCO3m].Aaq_bins_init(b)/surrogate[config.iHCO3m].MM+surrogate[config.iCO3mm].Aaq_bins_init(b)/surrogate[config.iCO3mm].MM)*surrogate[i].MM/(Kaq*conc_org*surrogate[i].tau_air(b))*AQinit(b)/sum_mass;
		surrogate[i].kloss_gas+=1.0/surrogate[i].tau_air(b)*AQinit(b)/sum_mass;*/
		cout << "Carbonates Inorganics in dynamic with imethod=0. Should not happen" << endl;
		exit(1);
	      }
            else if (i==config.iHCl)
              {
                Kaq=surrogate[i].Kaq(b);
                //compute kinetic rate of absorption
                conc_aq=surrogate[config.iClm].Aaq_bins_init(b)/surrogate[config.iClm].MM*surrogate[i].MM;
	      
                surrogate[i].k1_aq(b,index)=(surrogate[i].Ag*Kaq*conc_org-conc_aq)/(Kaq*conc_org*surrogate[i].tau_air(b))*AQinit(b)/sum_mass*surrogate[i].fac_corr_ph(b);       

                surrogate[i].Jdn_aq(b,index)=0.0;
                if (surrogate[i].k1_aq(b,index)<0.0 and surrogate[config.iClm].Aaq_bins_init(b)>tiny)
                  surrogate[i].Jdn_aq(b,index)=surrogate[i].k1_aq(b,index)/surrogate[config.iClm].Aaq_bins_init(b)*surrogate[config.iClm].MM/surrogate[i].MM;
              }

          
          }
    }
  
  for (i=0;i<n;++i)      
    if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophilic)
      {
        sumkpositive=0.0;
        sumknegative=0.0;
        for (b=0;b<config.nbins;b++)
          if (surrogate[i].time_aq(b)>=config.tequilibrium)			
            {
              if (surrogate[i].k1_aq(b,index)>0.0)  
                sumkpositive+=surrogate[i].k1_aq(b,index);
              else if (surrogate[i].k1_aq(b,index)<0.0) 
                sumknegative+=surrogate[i].k1_aq(b,index);
            }
		
        for (b=0;b<config.nbins;b++)
          {
            //Flux limitations				
            if (surrogate[i].k1_aq(b,index)>0.0 and surrogate[i].Ag*surrogate[i].Kaq(b)*AQinit(b)-surrogate[i].Aaq_bins_init(b)>0.0)
              surrogate[i].k1_aq(b,index)=min(surrogate[i].k1_aq(b,index),1.0/config.deltatmin*
                                              ((surrogate[i].Ag+max(-sumknegative-sumkpositive,0.0)*config.deltatmin+surrogate[i].Aaq_bins_init(b))*surrogate[i].Kaq(b)*AQinit(b)/
                                               (1.0+surrogate[i].Kaq(b)*AQinit(b))-surrogate[i].Aaq_bins_init(b)));
            else if (surrogate[i].k1_aq(b,index)<0.0 and surrogate[i].Ag*surrogate[i].Kaq(b)*AQinit(b)-surrogate[i].Aaq_bins_init(b)<0.0)
              surrogate[i].k1_aq(b,index)=max(surrogate[i].k1_aq(b,index),1.0/config.deltatmin*
                                              ((surrogate[i].Ag+min(sumknegative+sumkpositive,0.0)*config.deltatmin+surrogate[i].Aaq_bins_init(b))*surrogate[i].Kaq(b)*AQinit(b)/
                                               (1.0+surrogate[i].Kaq(b)*AQinit(b))-surrogate[i].Aaq_bins_init(b)));
		
            //compute fluxes
            surrogate[i].Jdn_aq(b,index)=0.0; 
            if (surrogate[i].k1_aq(b,index)<0.0 and
                surrogate[i].Aaq_bins_init(b)>tiny)
              surrogate[i].Jdn_aq(b,index)=surrogate[i].k1_aq(b,index)/
						 surrogate[i].Aaq_bins_init(b);
          }
      }



}


void solidification_bins_ssh(model_config &config, vector<species>& surrogate, double &Temperature, Array<double, 1> &MMaq, Array<double,1> &AQinit, Array<double,1> &chp, Array<double,1> &ionic, Array<double,1> &LWC, double& factor)
{
  int n=surrogate.size();
  int i,b;
  Array <double, 1> excess,ratio;
  double Ke,prod_conc,xmol,total;     
  excess.resize(n);
  ratio.resize(n);
  /*
  for (i=0;i<n;i++)
    surrogate[i].Ag2=surrogate[i].Ag; */

  for (b=0;b<config.nbins;b++)
    {

      excess=0.0; 
      double sum_excess=0.0;
      ratio=1.0;
      double sum_error_ap=10;
      double factor2=1.;
      
      double conc_org=0.;               

      /*
      for (i=0;i<n;i++)
	if (i==config.iNa or i==config.iCa or i==config.iK or i==config.iMg)
	  surrogate[i].Aaq_bins_init2(b)=surrogate[i].Aaq_bins_init0(b);
	else
	  {
	    //surrogate[i].Ag2=surrogate[i].Ag;
	    //surrogate[i].Asol_bins_init2(b)=0.;
	    surrogate[i].Aaq_bins_init2(b)=surrogate[i].Aaq_bins_init(b);
	    //surrogate[i].molality2=surrogate[i].molality;
	    //if (surrogate[i].nion>2 and surrogate[i].Ap>config.MOmin)
	    //  factor2=1.0; //0.5;
	  }

      if (config.solids)
	for (int j=0;j<n;j++)
	  if (surrogate[j].is_solid)
	    {			
	      for (int j2=0;j2<surrogate[j].nion;j2++)
		{
		  //cout << j2 << " " << surrogate[j].iion1 << " " << surrogate[j].iion2 << " " << config.iHCO3m << " " << condifendl;
		  if (j2==0 and (surrogate[j].iion1==config.iNa or surrogate[j].iion1==config.iCa or surrogate[j].iion1==config.iK or surrogate[j].iion1==config.iMg))
		    {
		      surrogate[surrogate[j].pion1].Aaq_bins_init2(b)+=(surrogate[j].Asol_bins_init0(b))/surrogate[j].MM*surrogate[j].pion1;		      
		    }
		  else if (j2==2 and (surrogate[j].iion2==config.iNa or surrogate[j].iion2==config.iCa or surrogate[j].iion2==config.iK or surrogate[j].iion2==config.iMg))
		    {
		      surrogate[surrogate[j].pion1].Aaq_bins_init2(b)+=(surrogate[j].Asol_bins_init0(b))/surrogate[j].MM*surrogate[j].pion2;		      
		    }
		  else if (j2==2 and (surrogate[j].iion3==config.iNa or surrogate[j].iion3==config.iCa or surrogate[j].iion3==config.iK or surrogate[j].iion3==config.iMg))
		    {
		      surrogate[surrogate[j].pion1].Aaq_bins_init2(b)+=(surrogate[j].Asol_bins_init0(b))/surrogate[j].MM*surrogate[j].pion3;		      
		    }
		}
	      surrogate[j].Asol_bins_init2(b)=0.;
	    }

      compute_conc_org_bins_ssh(config, surrogate, Temperature,  MMaq, AQinit, chp, ionic, LWC, b, conc_org);
      for (i=0;i<n;++i)
	if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)      
	surrogate[i].molality2=surrogate[i].Aaq_bins_init2(b)/surrogate[i].MM/conc_org*1000.0;     */

      compute_conc_org_bins_ssh(config, surrogate, Temperature,  MMaq, AQinit, chp, ionic, LWC, b, conc_org);
  
      for (i=0;i<n;++i)
	if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)      
	  surrogate[i].molality=surrogate[i].Aaq_bins_init(b)/surrogate[i].MM/conc_org*1000.0; 
      
      for (i=0;i<n;i++)	
	{
	  //surrogate[i].Ag2=surrogate[i].Ag;
	  surrogate[i].Asol_bins_init2(b)=surrogate[i].Asol_bins_init(b);
	  surrogate[i].Aaq_bins_init2(b)=surrogate[i].Aaq_bins_init(b);
	  surrogate[i].molality2=surrogate[i].molality;	  
	  
	  //if (surrogate[i].nion>2 and surrogate[i].Ap>config.MOmin)
	  //  factor2=1.0; //0.5;	  
	}


      
      
      /*
      for (i=0;i<n;i++)       
	{
	  //surrogate[i].Ag2=surrogate[i].Ag;
	  surrogate[i].Asol_bins_init2(b)=surrogate[i].Asol_bins_init(b);
	  surrogate[i].Aaq_bins_init2(b)=surrogate[i].Aaq_bins_init(b);
	  surrogate[i].molality2=surrogate[i].molality;
	  //if (surrogate[i].nion>2 and surrogate[i].Ap>config.MOmin)
	  //  factor2=1.0; //0.5;
	}*/
  
      int k=0;
      while (sum_error_ap/factor2>0.00001 and k<1)
	{
	  //double sum_error_ap_save=sum_error_ap;
	  sum_error_ap=0;
	  for (i=0;i<n;i++)
	    if (surrogate[i].is_solid)
	      {
		int iion1=surrogate[i].iion1;
		int iion2=surrogate[i].iion2;
		int pion1=surrogate[i].pion1;
		int pion2=surrogate[i].pion2;
		int iion3=surrogate[i].iion3;
		int pion3=surrogate[i].pion3;
		
		prod_conc=pow(surrogate[iion1].molality2,pion1)*pow(surrogate[iion2].molality2,pion2);
		xmol=-surrogate[i].Asol_bins_init2(b)/surrogate[i].MM/conc_org*1000.;            
		Ke=surrogate[i].keq/(pow(surrogate[iion1].gamma_aq,pion1)*pow(surrogate[iion2].gamma_aq,pion2));
				
		if (surrogate[i].nion>2)
		  {
		    Ke=Ke/pow(surrogate[iion3].gamma_aq,pion3);
		    prod_conc*=pow(surrogate[iion3].molality2,pion3);               
		  }
		else
		  pion3=0.0;
		excess(i)=prod_conc-Ke;
        
		//cout << "excess: " << excess(i) << endl; 
		if ((excess(i)>0. and prod_conc>0.) or (excess(i)<0. and xmol<0.)) 
		  {    
		    //cout << "in " << endl;
		    //ratio=1.;
		    /*
		    if (iion1==config.iNH4p or iion2==config.iNH4p or iion3==config.iNH4p)
		      {
			total=surrogate[config.iNH3].Ag2/surrogate[config.iNH3].MM+surrogate[config.iNH4p].Aaq_bins_init2(b)/surrogate[config.iNH4p].MM;
			if(total>0.) ratio(config.iNH4p)=surrogate[config.iNH4p].Aaq_bins_init2(b)/surrogate[config.iNH4p].MM/total;
		      }

		    if  (iion1==config.iNO3m or iion2==config.iNO3m or iion3==config.iNO3m)
		      {
			total=surrogate[config.iHNO3].Ag2/surrogate[config.iHNO3].MM+surrogate[config.iNO3m].Aaq_bins_init2(b)/surrogate[config.iNO3m].MM;
			if(total>0.) ratio(config.iNO3m)=surrogate[config.iNO3m].Aaq_bins_init2(b)/surrogate[config.iNO3m].MM/total;
		      }


		    if (iion1==config.iClm or iion2==config.iClm or iion3==config.iClm)
		      {
			total=surrogate[config.iHCl].Ag2/surrogate[config.iHCl].MM+surrogate[config.iClm].Aaq_bins_init2(b)/surrogate[config.iClm].MM;
			if(total>0.) ratio(config.iClm)=surrogate[config.iClm].Aaq_bins_init2(b)/surrogate[config.iClm].MM/total;
		      }*/

		    if (iion1==config.iSO4mm or iion1==config.iHSO4m or iion2==config.iSO4mm or iion2==config.iHSO4m or iion3==config.iSO4mm or iion3==config.iHSO4m)
		      {
			total=surrogate[config.iSO4mm].Aaq_bins_init2(b)/surrogate[config.iSO4mm].MM+surrogate[config.iHSO4m].Aaq_bins_init2(b)/surrogate[config.iHSO4m].MM;
			if(total>0.) ratio(config.iSO4mm)=surrogate[config.iSO4mm].Aaq_bins_init2(b)/surrogate[config.iSO4mm].MM/total;
			if(total>0.) ratio(config.iHSO4m)=1.0-ratio(config.iSO4mm);
		      }

		    
		    if (iion1==config.iCO3mm or iion1==config.iHCO3m or iion2==config.iCO3mm or iion2==config.iHCO3m or iion3==config.iCO3mm or iion3==config.iHCO3m)
		      {
			total=surrogate[config.iCO3mm].Aaq_bins_init2(b)/surrogate[config.iCO3mm].MM+surrogate[config.iHCO3m].Aaq_bins_init2(b)/surrogate[config.iHCO3m].MM;
			if(total>0.)
			  {
			    ratio(config.iCO3mm)=(surrogate[config.iCO3mm].Aaq_bins_init2(b)/surrogate[config.iCO3mm].MM/total);
			    ratio(config.iHCO3m)=1.-ratio(config.iCO3mm); //surrogate[config.iHCO3m].Aaq_bins_init2(b)/surrogate[config.iHCO3m].MM/total;
			  }

			/*
			cout << ratio(config.iCO3mm) << endl;
			cout << surrogate[config.iCO3mm].Aaq_bins_init2(b)/surrogate[config.iCO3mm].MM << " " << surrogate[config.iHCO3m].Aaq_bins_init2(b)/surrogate[config.iHCO3m].MM << endl;
			cout << surrogate[config.iCO3mm].Aaq_bins(b)/surrogate[config.iCO3mm].MM << " " << surrogate[config.iHCO3m].Aaq_bins(b)/surrogate[config.iHCO3m].MM << endl;
			cout << surrogate[config.iCO3mm].Aaq_bins_init(b)/surrogate[config.iCO3mm].MM << " " << surrogate[config.iHCO3m].Aaq_bins_init(b)/surrogate[config.iHCO3m].MM << endl;*/
		      }

		    

		    if (excess(i)<0.) // and iion3<1)
		      xmol=-surrogate[i].Asol_bins_init2(b)/surrogate[i].MM/conc_org*1000.; //0.; //-pow(Ke,1./3);
		    //xmol=max(xmol,-pow(-excess(i),1.0/(pion1+pion2+pion3)));
		    else
		      xmol=0.;
		    //xmol=max(xmol,pow(prod_conc-Ke,1.0/(pion1+pion2+pion3)));

		    double error=100.;
		    //double error_save=0.0;
		    double derror=0.;
		    int iter=0;
		    //double xmol_save=xmol;
		    while (abs(error)/max(Ke,1.0e-4)>1.0e-8 and iter<2000)
		      {
			double m1=surrogate[iion1].molality2-ratio(iion1)*xmol*pion1;
			double m2=surrogate[iion2].molality2-ratio(iion2)*xmol*pion2;                                                
			double m3=1.;                   
			if (iion3>0) 
			  m3=surrogate[iion3].molality2-ratio(iion3)*xmol*pion3;

			//cout << "mol " << m1 << " " << m2 << " " << endl;
			//cout << iter << endl;
                           
			error=pow(m1,pion1)*pow(m2,pion2);
			derror=-ratio(iion1)*pion1*pow(m1,pion1-1)*pion1*pow(m2,pion2)*pow(m3,pion3)-ratio(iion2)*pion2*pion2*pow(m2,pion2-1)*pow(m1,pion1)*pow(m3,pion3);
			if (iion3>0) 
			  {
			    derror-=ratio(iion3)*pion3*pion3*pow(m3,pion3-1)*pow(m2,pion2)*pow(m1,pion1);                   
			    error*=pow(m3,pion3);
			  }
			error=error-Ke;
			double xmol_s=xmol;
			if (abs(derror)>0.)                      
			  xmol=xmol-error/derror;

			xmol=min(xmol,surrogate[iion1].molality2/ratio(iion1)/pion1*0.999);
			xmol=min(xmol,surrogate[iion2].molality2/ratio(iion2)/pion2*0.999);
			if (iion3>0)
			  xmol=min(xmol,surrogate[iion3].molality2/ratio(iion3)/pion3*0.999);
			if (xmol==xmol_s)
			  error=0;
			iter++;
			//cout << "error: " << error << " " << error/Ke << endl;
			//xmol_save=xmol;
			//cout << "X: " << iter << " " << xmol << " " << error << " " << derror << " " << m1 << " " << m2 << " " << m3 << " " << surrogate[iion1].molality2 << endl; //Ke << " " << pow(m1,pion1)*pow(m2,pion2) << endl;
			
                             
		      }
		    if (iter==2000)
		      {
			xmol=0.99*xmol;
		      }

		    excess(iion1)-=xmol*ratio(iion1)*pion1;
		    excess(iion2)-=xmol*ratio(iion2)*pion2;	       
		    excess(i)=xmol;
		    //cout << xmol << " " << excess(iion1) << " " << excess(iion2) << endl;
		  }             

		if (prod_conc==0.0) excess(i)=0.0;             
		if (xmol>0.) sum_excess+=xmol;

		//if (xmol<0.) excess(i)=xmol_save;
		//cout << surrogate[i].name << " " << excess(i) << " " << xmol << " " << surrogate[i].Ap << " " << surrogate[i].keq << endl;

		if (iion1==config.iSO4mm or iion1==config.iHSO4m or iion2==config.iSO4mm or iion2==config.iHSO4m or iion3==config.iSO4mm or iion3==config.iHSO4m)
		  {
		    total=surrogate[config.iSO4mm].Aaq_bins_init(b)/surrogate[config.iSO4mm].MM+surrogate[config.iHSO4m].Aaq_bins_init(b)/surrogate[config.iHSO4m].MM;
		    if(total>0.) ratio(config.iSO4mm)=surrogate[config.iSO4mm].Aaq_bins_init(b)/surrogate[config.iSO4mm].MM/total;
		    if(total>0.) ratio(config.iHSO4m)=1.0-ratio(config.iSO4mm);       
		  }

		if (iion1==config.iCO3mm or iion1==config.iHCO3m or iion2==config.iCO3mm or iion2==config.iHCO3m or iion3==config.iCO3mm or iion3==config.iHCO3m)
		  {
		    total=surrogate[config.iCO3mm].Aaq_bins_init2(b)/surrogate[config.iCO3mm].MM+surrogate[config.iHCO3m].Aaq_bins_init2(b)/surrogate[config.iHCO3m].MM;
		    if(total>0.)
		      {
			ratio(config.iCO3mm)=(surrogate[config.iCO3mm].Aaq_bins_init2(b)/surrogate[config.iCO3mm].MM/total);
			ratio(config.iHCO3m)=surrogate[config.iHCO3m].Aaq_bins_init2(b)/surrogate[config.iHCO3m].MM/total;
		      }
		  }

		//excess*=factor2;
		if (excess(i)!=excess(i))
		  {
		    cout << "weird excess " << excess(i) << endl;
		    exit(0);
		  }
        
		if (excess(i)<0.) //dissolution
		  {		    
		    excess(i)=factor2*max(-surrogate[i].Asol_bins_init2(b)/surrogate[i].MM/conc_org*1000,excess(i));		    
		    //cout << "dissolution " << surrogate[i].Ap << " " << factor << " " << surrogate[config.iH2O].Aaq << " " << surrogate[i].gamma_aq << " "<< excess(i) << endl;
		    //surrogate[i].Ap=factor*excess(i)*surrogate[i].MM*conc_org/1000.+surrogate[i].Ap;
		    if (surrogate[i].Asol_bins_init2(b)>0.) sum_error_ap+=abs(excess(i)*surrogate[i].MM*conc_org/1000./surrogate[i].Asol_bins_init2(b));
		    surrogate[i].Asol_bins_init2(b)=max(excess(i)*surrogate[i].MM*conc_org/1000.+surrogate[i].Asol_bins_init2(b),0.);
		    //cout << surrogate[i].Ap << endl;            
		    //int iion1=surrogate[i].iion1;
		    //int pion1=surrogate[i].pion1;
		    surrogate[iion1].Aaq_bins_init2(b)=-excess(i)*surrogate[iion1].MM*conc_org/1000.*pion1+surrogate[iion1].Aaq_bins_init2(b);
		    //int iion2=surrogate[i].iion2;
		    //int pion2=surrogate[i].pion2;
		    surrogate[iion2].Aaq_bins_init2(b)=-excess(i)*surrogate[iion2].MM*conc_org/1000.*pion2+surrogate[iion2].Aaq_bins_init2(b);
		    if (surrogate[i].nion>2)
		      {
			//int iion3=surrogate[i].iion3;
			//int pion3=surrogate[i].pion3;
			surrogate[iion3].Aaq_bins_init2(b)=-excess(i)*surrogate[iion3].MM*conc_org/1000.*pion3+surrogate[iion3].Aaq_bins_init2(b);
		      }

		    /*
		    if (iion1==config.iNH4p or iion2==config.iNH4p or iion3==config.iNH4p)
		      {
			total=surrogate[config.iNH3].Ag2/surrogate[config.iNH3].MM+surrogate[config.iNH4p].Aaq_bins_init2(b)/surrogate[config.iNH4p].MM;
			surrogate[config.iNH3].Ag2=total*surrogate[config.iNH3].MM*(1.0-ratio(config.iNH4p));
			surrogate[config.iNH4p].Aaq_bins_init2(b)=total*surrogate[config.iNH4p].MM*ratio(config.iNH4p);
		      }

		    if (iion1==config.iNO3m or iion2==config.iNO3m or iion3==config.iNO3m)
		      {
			total=surrogate[config.iHNO3].Ag2/surrogate[config.iHNO3].MM+surrogate[config.iNO3m].Aaq_bins_init2(b)/surrogate[config.iNO3m].MM;
			surrogate[config.iHNO3].Ag2=total*surrogate[config.iHNO3].MM*(1.0-ratio(config.iNO3m));
			surrogate[config.iNO3m].Aaq_bins_init2(b)=total*surrogate[config.iNO3m].MM*ratio(config.iNO3m);
		      }


		    if (iion1==config.iClm or iion2==config.iClm or iion3==config.iClm)
		      {
			total=surrogate[config.iHCl].Ag2/surrogate[config.iHCl].MM+surrogate[config.iClm].Aaq_bins_init2(b)/surrogate[config.iClm].MM;
			surrogate[config.iHCl].Ag2=total*surrogate[config.iHCl].MM*(1.0-ratio(config.iClm));
			surrogate[config.iClm].Aaq_bins_init2(b)=total*surrogate[config.iClm].MM*ratio(config.iClm);
		      }*/

		    if (iion1==config.iSO4mm or iion1==config.iHSO4m or iion2==config.iSO4mm or iion2==config.iHSO4m or iion3==config.iSO4mm or iion3==config.iHSO4m)
		      {
			total=surrogate[config.iSO4mm].Aaq_bins_init2(b)/surrogate[config.iSO4mm].MM+surrogate[config.iHSO4m].Aaq_bins_init2(b)/surrogate[config.iHSO4m].MM;
			surrogate[config.iSO4mm].Aaq_bins_init2(b)=total*surrogate[config.iSO4mm].MM*ratio(config.iSO4mm);
			surrogate[config.iHSO4m].Aaq_bins_init2(b)=total*surrogate[config.iHSO4m].MM*(1.-ratio(config.iSO4mm));
		      }
		    
		    if (iion1==config.iCO3mm or iion1==config.iHCO3m or iion2==config.iCO3mm or iion2==config.iHCO3m or iion3==config.iCO3mm or iion3==config.iHCO3m)
		      {
			total=surrogate[config.iCO3mm].Aaq_bins_init2(b)/surrogate[config.iCO3mm].MM+surrogate[config.iHCO3m].Aaq_bins_init2(b)/surrogate[config.iHCO3m].MM;
			//cout << "diss: " << total << " " <<  ratio(config.iCO3mm) << endl;
			surrogate[config.iCO3mm].Aaq_bins_init2(b)=total*surrogate[config.iCO3mm].MM*ratio(config.iCO3mm);
			surrogate[config.iHCO3m].Aaq_bins_init2(b)=total*surrogate[config.iHCO3m].MM*ratio(config.iHCO3m);
			//cout << "diss2: " << 	surrogate[config.iCO3mm].Aaq_bins_init2(b) << " " << 	surrogate[config.iHCO3m].Aaq_bins_init2(b) << ' ' << total*surrogate[config.iHCO3m].MM*(1.-ratio(config.iCO3mm)) << " " << (1.-ratio(config.iCO3mm)) << " " << ratio(config.iHCO3m) << endl;
		      }
		  }

  
		//cout << "avant " << min_excess << " " << max_excess << " " << 0.225*surrogate[42].Ap << " " << 0.225*surrogate[42].Ap+0.27272727*surrogate[43].Ap+surrogate[config.iNH3].Ag*18./17.+surrogate[config.iNH4p].Aaq << endl;

		//cout << "prendant " << surrogate[config.iNa].Aaq << " " << surrogate[44].Ap*2*23./surrogate[44].MM+surrogate[config.iNa].Aaq << endl;
        
		if (excess(i)>0.) //solidification
		  {
		    //cout << "solid " << endl;
		    //cout << "rat2: " << ratio(config.iCO3mm) << " " << ratio(config.iHCO3m) << endl;
		    //cout << "solid " << surrogate[i].Ap << " " << factor << " " << surrogate[config.iH2O].Aaq << " " << surrogate[config.iH2O].gamma_aq << " " << config.iH2O << endl;
		    //surrogate[i].Ap=max(factor*excess(i)*surrogate[i].MM*conc_org/1000.+surrogate[i].Ap,0.)
		    
		    if (surrogate[i].Asol_bins_init2(b)>0.) sum_error_ap+=abs(excess(i)*surrogate[i].MM*conc_org/1000./surrogate[i].Asol_bins_init2(b));
		    excess(i)*=factor2;
		    surrogate[i].Asol_bins_init2(b)+=excess(i)*surrogate[i].MM*conc_org/1000.;
		    //int iion1=surrogate[i].iion1;
		    //int pion1=surrogate[i].pion1;		    
		    surrogate[iion1].Aaq_bins_init2(b)-=excess(i)*surrogate[iion1].MM*conc_org/1000.*pion1;
		    //cout << "Ca: " << surrogate[iion1].Aaq_bins_init2(b) << endl;
		    //int iion2=surrogate[i].iion2;
		    //int pion2=surrogate[i].pion2;
		    
		    surrogate[iion2].Aaq_bins_init2(b)-=excess(i)*surrogate[iion2].MM*conc_org/1000.*pion2;
		    //surrogate[config.iCO3mm].Aaq_bins_init2(b)-=excess(i)*surrogate[config.iCO3mm].MM*conc_org/1000.*pion2*ratio(config.iCO3mm);
		    //surrogate[config.iHCO3m].Aaq_bins_init2(b)-=excess(i)*surrogate[config.iHCO3m].MM*conc_org/1000.*pion2*(1.-ratio(config.iCO3mm));
		    if (surrogate[i].nion>2)
		      {
			//int iion3=surrogate[i].iion3;
			//int pion3=surrogate[i].pion3;
			surrogate[iion3].Aaq_bins_init2(b)=-excess(i)*surrogate[iion3].MM*conc_org/1000.*pion3+surrogate[iion3].Aaq_bins_init2(b);
		      }

		    //cout << "pendant " << min_excess << " " << max_excess << " " << 0.27272727*surrogate[43].Ap << " " << 0.27272727*surrogate[43].Ap+surrogate[config.iNH3].Ag*18./17.+surrogate[config.iNH4p].Aaq << endl;
		    double total=0.0;
		    //cout << surrogate[config.iNH3].Ag/surrogate[config.iNH3].MM+surrogate[config.iNH4p].Aaq/surrogate[config.iNH4p].MM << endl;
		    /*
		    if (iion1==config.iNH4p or iion2==config.iNH4p or iion3==config.iNH4p)
		      {
			total=surrogate[config.iNH3].Ag2/surrogate[config.iNH3].MM+surrogate[config.iNH4p].Aaq_bins_init2(b)/surrogate[config.iNH4p].MM;
			surrogate[config.iNH3].Ag2=max(total*surrogate[config.iNH3].MM*(1.0-ratio(config.iNH4p)),0.);
			surrogate[config.iNH4p].Aaq_bins_init2(b)=total*surrogate[config.iNH4p].MM*ratio(config.iNH4p);
		      }
		    //cout << surrogate[config.iNH3].Ag/surrogate[config.iNH3].MM << " " << surrogate[config.iNH4p].Aaq/surrogate[config.iNH4p].MM << endl;

		    //cout << surrogate[config.iNH3].Ag/surrogate[config.iNH3].MM+surrogate[config.iNH4p].Aaq/surrogate[config.iNH4p].MM << endl;
		    if (iion1==config.iNO3m or iion2==config.iNO3m or iion3==config.iNO3m)
		      {
			total=surrogate[config.iHNO3].Ag2/surrogate[config.iHNO3].MM+surrogate[config.iNO3m].Aaq_bins_init2(b)/surrogate[config.iNO3m].MM;
			surrogate[config.iHNO3].Ag2=max(total*surrogate[config.iHNO3].MM*(1.0-ratio(config.iNO3m)),0.);
			surrogate[config.iNO3m].Aaq_bins_init2(b)=total*surrogate[config.iNO3m].MM*ratio(config.iNO3m);
		      }


		    if (iion1==config.iClm or iion2==config.iClm or iion3==config.iClm)
		      {
			total=surrogate[config.iHCl].Ag2/surrogate[config.iHCl].MM+surrogate[config.iClm].Aaq_bins_init2(b)/surrogate[config.iClm].MM;
			surrogate[config.iHCl].Ag2=max(total*surrogate[config.iHCl].MM*(1.0-ratio(config.iClm)),0.);
			surrogate[config.iClm].Aaq_bins_init2(b)=total*surrogate[config.iClm].MM*ratio(config.iClm);
			}*/

		    if (iion1==config.iSO4mm or iion1==config.iHSO4m or iion2==config.iSO4mm or iion2==config.iHSO4m or iion3==config.iSO4mm or iion3==config.iHSO4m)
		      {
			total=surrogate[config.iSO4mm].Aaq_bins_init2(b)/surrogate[config.iSO4mm].MM+surrogate[config.iHSO4m].Aaq_bins_init2(b)/surrogate[config.iHSO4m].MM;
			surrogate[config.iSO4mm].Aaq_bins_init2(b)=total*surrogate[config.iSO4mm].MM*ratio(config.iSO4mm);
			surrogate[config.iHSO4m].Aaq_bins_init2(b)=total*surrogate[config.iHSO4m].MM*ratio(config.iHSO4m);
		      }
		    
		    if (iion1==config.iCO3mm or iion1==config.iHCO3m or iion2==config.iCO3mm or iion2==config.iHCO3m or iion3==config.iCO3mm or iion3==config.iHCO3m)
		      {
			total=surrogate[config.iCO3mm].Aaq_bins_init2(b)/surrogate[config.iCO3mm].MM+surrogate[config.iHCO3m].Aaq_bins_init2(b)/surrogate[config.iHCO3m].MM;
			//cout << "rat3: " << ratio(config.iCO3mm) << " " << ratio(config.iHCO3m) << " " << total << endl;
			surrogate[config.iCO3mm].Aaq_bins_init2(b)=total*surrogate[config.iCO3mm].MM*ratio(config.iCO3mm);
			surrogate[config.iHCO3m].Aaq_bins_init2(b)=total*surrogate[config.iHCO3m].MM*ratio(config.iHCO3m);

			//cout << "la " << surrogate[config.iCO3mm].Aaq_bins_init2(b) << " " << surrogate[config.iHCO3m].Aaq_bins_init2(b) << endl; 
		      }
		  }
   
		//cout << "fin " << min_excess << " " << max_excess << " " << 0.27272727*surrogate[43].Ap << " " << 0.27272727*surrogate[43].Ap+surrogate[config.iNH3].Ag*18./17.+surrogate[config.iNH4p].Aaq << endl;
		//cout << surrogate[i].Ap << endl;
		//cout << "la " << surrogate[i].Ap << " " << factor << " " << surrogate[config.iH2O].Aaq << " " << surrogate[config.iH2O].gamma_aq << " " << config.iH2O << endl;
        
		if (surrogate[i].Asol_bins_init2(b)<0. and surrogate[config.iH2O].Aaq_bins_init2(b)>1.0e4) exit(0);
		int j;
		for (j=0;j<n;++j)
		  if (surrogate[j].is_organic==false and j!=config.iH2O and surrogate[j].is_inorganic_precursor==false)
		    if (j!=config.iHp and j!=config.iOHm)
		      {

			//if (surrogate[j].Aaq>0.) cout << surrogate[j].name << " " << surrogate[j].Aaq << endl;
			surrogate[j].molality2=surrogate[j].Aaq_bins_init2(b)/surrogate[j].MM/conc_org*1000.;
			//if (surrogate[j].molality<0.) cout << surrogate[j].name << " " << surrogate[j].molality << " " << excess(i) << endl;
			//if (surrogate[j].molality<0.) exit(0);
		      }

		

		//cout << "apres " << surrogate[config.iNa].Aaq << " " << surrogate[44].Ap*2*23./surrogate[44].MM+surrogate[config.iNa].Aaq << endl;

	      }
	  //cout << sum_error_ap/factor2 << endl;
	  //if ((sum_error_ap_save-sum_error_ap)<1.e-4*sum_error_ap) factor2=max(factor2/2,0.5);
	  k++;
	  
	  for (i=0;i<n;i++)  
	    if (surrogate[i].Aaq_bins_init2(b)<0.) // and surrogate[i].Aaq_bins(b)>0.)
	      {
		cout << "exit wtf " << endl;
		cout << "error " << surrogate[i].name << " " << surrogate[i].Aaq_bins_init(b) << " " << surrogate[i].Aaq_bins(b) << " " << surrogate[i].Aaq_bins_init2(b) << " " << surrogate[i].molality << " " << surrogate[i].molality2 << endl;
		cout << "CO3mm " << surrogate[config.iCO3mm].Aaq_bins_init(b) << " " << surrogate[config.iCO3mm].Aaq_bins(b) << " " << surrogate[config.iCO3mm].Aaq_bins_init2(b) << " " << surrogate[config.iCO3mm].molality << " " << surrogate[config.iCO3mm].molality2 << endl;
		cout << ratio(config.iCO3mm) << " " << ratio(config.iHCO3m) << endl;
		exit(0);
	      }
	}      
            
      for (i=0;i<n;i++)       
	{	    
	  surrogate[i].Asol_bins_init(b)=surrogate[i].Asol_bins_init2(b);
	  surrogate[i].Aaq_bins_init(b)=surrogate[i].Aaq_bins_init2(b);
	}
      
      for (i=0;i<n;i++)
	{
	  //surrogate[i].Ag=factor*surrogate[i].Ag2+(1.-factor)*surrogate[i].Ag;
	  //surrogate[i].Asol_bins_init(b)=factor*surrogate[i].Asol_bins_init2(b)+(1.-factor)*surrogate[i].Asol_bins_init(b);
	  //surrogate[i].Aaq_bins_init(b)=factor*surrogate[i].Aaq_bins_init2(b)+(1.-factor)*surrogate[i].Aaq_bins_init(b);
	  //surrogate[i].molality=factor*surrogate[i].molality2+(1.-factor)*surrogate[i].molality;
	  if (surrogate[i].Aaq_bins_init(b)<0.) // and surrogate[i].Aaq_bins(b)>0.)
	    {
	      cout << "error " << surrogate[i].name << " " << surrogate[i].Aaq_bins_init(b) << " " << surrogate[i].Aaq_bins(b) << " " << surrogate[i].Aaq_bins_init2(b) << " " << surrogate[i].molality << " " << surrogate[i].molality2 << endl;
	      cout << "CO3mm " << surrogate[config.iCO3mm].Aaq_bins_init(b) << " " << surrogate[config.iCO3mm].Aaq_bins(b) << " " << surrogate[config.iCO3mm].Aaq_bins_init2(b) << " " << surrogate[config.iCO3mm].molality << " " << surrogate[config.iCO3mm].molality2 << endl;
	      cout << ratio(config.iCO3mm) << " " << ratio(config.iHCO3m) << endl;
	      exit(0);
	    }
	  
	}
    }

	  /*
  for (b=0;b<config.nbins;b++)
    {
      if (surrogate[config.iHCO3m].Aaq_bins_init(b)<0.)
	{
	  cout << surrogate[config.iHCO3m].Aaq_bins_init(b) << " " << surrogate[config.iHCO3m].Aaq_bins(b) << endl;	  
	  cout << "neg conc" << endl;
	  exit(0);
	}
	}*/
    
}

void compute_morphology_ssh(model_config &config, Array<double,1> &vsol, Array<double, 1> & Number)
{
  int b,ilayer;
  double vorg;
  //compute morpholgy
  for (b=0;b<config.nbins;++b)
    if (Number(b) > 0.0)  
      {
        config.dbound(b,0)=pow(3.0/(4.0*3.14159)*vsol(b),1.0/3.0)/Number(b);
        vorg=(4.0/3.0*3.14159*pow(0.5e-6*config.diameters(b),3.0)-vsol(b)/Number(b));  
        for (ilayer=0;ilayer<config.nlayer;++ilayer)
          {
            config.dbound(b,ilayer+1)=pow(pow(config.dbound(b,ilayer),3)+3.0/(4.0*3.14159)*vorg*config.Vlayer(ilayer),1.0/3.0);
            config.Radius(b,ilayer)=0.5*(config.dbound(b,ilayer+1)+config.dbound(b,ilayer));	  
          }
      }
}

void flux_org_ssh(model_config &config, vector<species>& surrogate,
		  Array<double, 3> &MOinit, Array<double, 1> &AQinit,
		  double &tiny, int index)
{
  //compute kinetic rates for the absorption of a compound in the organic phase in a bin and a layer
  //index = 0 : first evaluation of rates)
  //index = 1 : second evaluation of rates
  //double Kp;
  int n=surrogate.size();
  int i,b,ilayer,iphase,jphase;
  double sum,sum_mass;
  Array<double, 1> kpmo_interface;
  kpmo_interface.resize(config.nlayer);
  Array <double,2> rJ_interface;
  rJ_interface.resize(config.nlayer,config.max_number_of_phases);
  double sumkpositive=0.0;
  double sumknegative=0.0;
  double sumkpmo_interface=0.0;
  double sumk=0.0;
  double sum1=0.0;
  double sum2=0.0;
  double sumap=0.0;
  double kpmo=0.0;
  double F1=0.0;
  double F2=0.0;
  double surf1=0.0;
  double surf2=0.0;
  double dorg=0.0;
  bool eq_last_layer=true;
  double tau_interface= 0.0;
  double ap_interface=0.0;		  
  double Vinterface=0.0;
  int ilayer_interface=0;		  
  double Jinterface=0.0;
  double ktot1=0.0;
  double kcond=0.0;	    	      
  double ktot=0.0;
  double b2=0.0;
  double c2=0.0;
  double a=0.0;
  int jlayer;
   
  for (i=0;i<n;++i)      		
    for (b=0;b<config.nbins;b++)
      for (ilayer=0;ilayer<config.nlayer;ilayer++)
	for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
	  {
	    surrogate[i].k1(b,ilayer,iphase,index)=0.0;
	    surrogate[i].Jdn(b,ilayer,iphase,index)=0.0;
	  }
  
  if (config.explicit_representation)
    {
      for (i=0;i<n;++i)
	if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophobic)
	  for (b=0;b<config.nbins;++b)
	    {	     
	      sum_mass=AQinit(b);
	      for (ilayer=0;ilayer<config.nlayer;ilayer++)
		for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		  {
		    sum_mass+=MOinit(b,ilayer,iphase);
		    surrogate[i].k1(b,ilayer,iphase,index)=0.0;
		  }
	   
	      //compute kinetic rate of absorption	    
	      for (ilayer=0;ilayer<config.nlayer;ilayer++)
		for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		  if(ilayer==config.nlayer-1 and surrogate[i].time(b,ilayer,iphase)>=config.tequilibrium)
		    {		      
		      sum=0.0;
		      for (jphase=0;jphase<config.nphase(b,ilayer);++jphase)
			sum+=surrogate[i].Kp(b,ilayer,jphase)*MOinit(b,ilayer,jphase);
		      
		      if (abs(surrogate[i].Ag*surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)-surrogate[i].Ap_layer_init(b,ilayer,iphase))/surrogate[i].Ap_layer_init(b,ilayer,iphase)<0.001)
                        surrogate[i].k1(b,ilayer,iphase,index)=
                          (surrogate[i].Ag*surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)
                           -surrogate[i].Ap_layer_init(b,ilayer,iphase))/
                          (0.0 //surrogate[i].tau_diffusion(b,ilayer,iphase)
                           +sum/(1.0-AQinit(b)/sum_mass)*surrogate[i].tau_air(b));		
		      else
                        surrogate[i].k1(b,ilayer,iphase,index)=
                          (surrogate[i].Ag*surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)
                           -surrogate[i].Ap_layer_init(b,ilayer,iphase))/
                          (0.0 //surrogate[i].tau_diffusion(b,ilayer,iphase)
                           +sum/(1.0-AQinit(b)/sum_mass)*surrogate[i].tau_air(b));					      		     
	
		      surrogate[i].Jdn(b,ilayer,iphase,index)=0.0;	  
		      if (surrogate[i].k1(b,ilayer,iphase,index)<0.0 and
			  surrogate[i].Ap_layer_init(b,ilayer,iphase)>tiny)
			surrogate[i].Jdn(b,ilayer,iphase,index)=surrogate[i].k1(b,ilayer,iphase,index)/
			  surrogate[i].Ap_layer_init(b,ilayer,iphase);
		    }

	    
	      for (ilayer=0;ilayer<config.nlayer;++ilayer)		    
		{
		  F1=0.0;
		  F2=0.0;
		  sum=0.0;
		  for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		    sum+=surrogate[i].Ap_layer_init(b,ilayer,iphase)/(surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase));

		  if (ilayer>0)
		    {
		      surf1=4.0*3.14159*pow(config.dbound(b,ilayer),2);
		      sum1=0.0;
		      for (iphase=0;iphase<config.nphase(b,ilayer-1);++iphase)
			sum1+=surrogate[i].Ap_layer_init(b,ilayer-1,iphase)/(surrogate[i].Kp(b,ilayer-1,iphase)*MOinit(b,ilayer-1,iphase));

		      dorg=(config.dbound(b,ilayer)-config.Radius(b,ilayer-1))/(config.Radius(b,ilayer)-config.Radius(b,ilayer-1))*surrogate[i].dif_org(b,ilayer-1)
			+(config.Radius(b,ilayer)-config.dbound(b,ilayer))/(config.Radius(b,ilayer)-config.Radius(b,ilayer-1))*surrogate[i].dif_org(b,ilayer);
			
		      F1=surf1*dorg*(sum1-sum)/(config.Radius(b,ilayer)-config.Radius(b,ilayer-1));
		    }

		  if (ilayer<config.nlayer-1)
		    {
		      surf2=4.0*3.14159*pow(config.dbound(b,ilayer+1),2);
		      sum2=0.0;
		      for (iphase=0;iphase<config.nphase(b,ilayer+1);++iphase)
			sum2+=surrogate[i].Ap_layer_init(b,ilayer+1,iphase)/(surrogate[i].Kp(b,ilayer+1,iphase)*MOinit(b,ilayer+1,iphase));

		      dorg=(config.dbound(b,ilayer+1)-config.Radius(b,ilayer))/(config.Radius(b,ilayer+1)-config.Radius(b,ilayer))*surrogate[i].dif_org(b,ilayer)
			+(config.Radius(b,ilayer+1)-config.dbound(b,ilayer+1))/(config.Radius(b,ilayer+1)-config.Radius(b,ilayer))*surrogate[i].dif_org(b,ilayer+1);
			
		      F2=surf2*dorg*(sum2-sum)/(config.Radius(b,ilayer+1)-config.Radius(b,ilayer));
		    }
		  
		  sum=0.0;
		  for (jphase=0;jphase<config.nphase(b,ilayer);++jphase)
		    sum+=surrogate[i].Kp(b,ilayer,jphase)*MOinit(b,ilayer,jphase);

		  for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		    if (surrogate[i].time(b,ilayer,iphase)>=config.tequilibrium)
		      {
			surrogate[i].k1(b,ilayer,iphase,index)+=pow(surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase),2.0)/sum*(F1+F2)/(4.0/3.0*3.14159*pow(config.diameters(b)*0.5e-6,3.0)*config.Vlayer(ilayer));

			//compute fluxes
			surrogate[i].Jdn(b,ilayer,iphase,index)=0.0;
			
			if (surrogate[i].k1(b,ilayer,iphase,index)<0.0 and
			    surrogate[i].Ap_layer_init(b,ilayer,iphase)>tiny)
			  surrogate[i].Jdn(b,ilayer,iphase,index)=surrogate[i].k1(b,ilayer,iphase,index)/
			    surrogate[i].Ap_layer_init(b,ilayer,iphase);
		      }
		    
		}
	    }
    }
  else // Use implicit representation (low number of layers).
    {
      if (config.nlayer>1)
	{
	  for (i=0;i<n;++i)
	    if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophobic)
	      for (b=0;b<config.nbins;++b)
		{
		  sum_mass=0.0;

		  if (AQinit(b)>config.MOmin)
		    sum_mass+=AQinit(b);
		  for (ilayer=0;ilayer<config.nlayer;ilayer++)
		    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		      sum_mass+=MOinit(b,ilayer,iphase);

		  eq_last_layer=true;
		  sumkpositive=0.0;
		  sumknegative=0.0;
		  tau_interface=max(1.0,surrogate[i].tau_diffusion(b,0,0));
		  ap_interface=0.0;		  
		  Vinterface=0.0;
		  sumkpmo_interface=0.0;
		  ilayer_interface=0;		  
		  kpmo_interface=0.0;
		  rJ_interface=0.0;
		  ilayer_interface=config.nlayer-1;
		  for (ilayer=config.nlayer-1;ilayer>=0;--ilayer)
		    {
		      for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
			if (surrogate[i].tau_diffusion(b,ilayer,iphase)<config.tequilibrium)
			  {
			    tau_interface=surrogate[i].tau_diffusion(b,ilayer,iphase);
			    ilayer_interface=ilayer;
			  }
		    }
		
		  //ilayer_interface=1;

		  for (ilayer=ilayer_interface;ilayer<config.nlayer;++ilayer)		   
		    {
		      for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
			{
			  ap_interface+=surrogate[i].Ap_layer_init(b,ilayer,iphase);
			  kpmo_interface(ilayer)+=surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase);
			  sumkpmo_interface+=kpmo_interface(ilayer);
			}
		      Vinterface+=config.Vlayer(ilayer);
		    }

		  Jinterface=(surrogate[i].Ag*sumkpmo_interface-ap_interface)/tau_interface;
		  
		  for (ilayer=0;ilayer<config.nlayer;++ilayer)
		    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		      if (surrogate[i].time(b,ilayer,iphase)>=config.tequilibrium)
			{		      
			  //compute kinetic rate of absorption
			  if (ilayer<ilayer_interface)
			    {			      
			      surrogate[i].k1(b,ilayer,iphase,index)=
				(surrogate[i].Ag*surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)
				 -surrogate[i].Ap_layer_init(b,ilayer,iphase))/surrogate[i].tau_diffusion(b,ilayer,iphase);
			    }
			  else
			    {
			      //mix concentrations over the interface
			      surrogate[i].Ap_layer_init(b,ilayer,iphase)=config.Vlayer(ilayer)/Vinterface*surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)/kpmo_interface(ilayer)*ap_interface;
			      surrogate[i].k1(b,ilayer,iphase,index)=config.Vlayer(ilayer)/Vinterface*surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)/kpmo_interface(ilayer)*Jinterface;
			    }
			  surrogate[i].Jdn(b,ilayer,iphase,index)=0.0;
			}
		      else
			surrogate[i].k1(b,ilayer,iphase,index)=0.0;
	      
		  
		  ktot1=0.0;
		  for (ilayer=0;ilayer<ilayer_interface;++ilayer)
		    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		      if (surrogate[i].time(b,ilayer,iphase)>=config.tequilibrium)		
			ktot1+=surrogate[i].k1(b,ilayer,iphase,index);
		     		  
		  for (ilayer=0;ilayer<config.nlayer;++ilayer)
		    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		      if (surrogate[i].time(b,ilayer,iphase)>=config.tequilibrium)			
			if (ilayer<ilayer_interface or
			    (Jinterface*ktot1<=0.0 and sumkpmo_interface*surrogate[i].tau_air(b)<surrogate[i].tau_diffusion(b,max(ilayer_interface-1,0),iphase)))
			  {
			    if (surrogate[i].k1(b,ilayer,iphase,index)>0.0)  
			      sumkpositive+=surrogate[i].k1(b,ilayer,iphase,index);
			    else if (surrogate[i].k1(b,ilayer,iphase,index)<0.0) 
			      sumknegative+=surrogate[i].k1(b,ilayer,iphase,index);
			  }
			else
			  {
			    sum=0.0;
			    for (jphase=0;jphase<config.nphase(b,ilayer);++jphase)
			      sum+=surrogate[i].Kp(b,ilayer,jphase)*MOinit(b,ilayer,jphase);
			
			    eq_last_layer=false;
			    if (surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)>0.0)
			      {
				surrogate[i].k1(b,ilayer,iphase,index)=
				  (1.0-AQinit(b)/sum_mass)*
				  config.Vlayer(ilayer)/Vinterface*surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)/kpmo_interface(ilayer)*
				  (surrogate[i].Ag-ap_interface/sumkpmo_interface)/surrogate[i].tau_air(b);			    
			      }
			    else
			      surrogate[i].k1(b,ilayer,iphase,index)=0.0;

			    if (surrogate[i].k1(b,ilayer,iphase,index)>0.0) 
			      sumkpositive+=surrogate[i].k1(b,ilayer,iphase,index);
			    else if (surrogate[i].k1(b,ilayer,iphase,index)<0.0)  
			      sumknegative+=surrogate[i].k1(b,ilayer,iphase,index);
			  }
      	      
		  kcond=0.0;	    	      
		  ktot=sumkpositive+sumknegative;
		  if (ktot>0.0)			
		    {
		      sum1=0.0;
		      sum2=0.0;
		      for (ilayer=0;ilayer<config.nlayer;++ilayer)		   
			if (surrogate[i].time(b,ilayer,0)>=config.tequilibrium)
			  {					  
			    sumk=0.0;
			    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
			      sumk+=surrogate[i].k1(b,ilayer,iphase,index);
			
			    sumap=0.0;
			    kpmo=0.0;
	
			    if (sumk>0.0)
			      {
				for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
				  {
				    sumap+=surrogate[i].Ap_layer_init(b,ilayer,iphase);
				    kpmo+=surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase);
				  }
				sum1+=sumk*sumap/kpmo;
				sum2+=sumk;
			      }		       	  
			  }
		      kcond=surrogate[i].Ag-sum1/sum2;
		    }
		  else if (ktot<0.0)
		    {
		      sum1=0.0;
		      sum2=0.0;
		      for (ilayer=0;ilayer<config.nlayer;++ilayer)
			for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
			  if (surrogate[i].time(b,ilayer,0)>=config.tequilibrium)
			    {
			      sumk=0.0;
			      for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
				sumk+=surrogate[i].k1(b,ilayer,iphase,index);
			  			 
			      sumap=0.0;
			      kpmo=0.0;
			  
			      if (sumk<0.0)
				{
				  for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
				    {
				      sumap+=surrogate[i].Ap_layer_init(b,ilayer,iphase);
				      kpmo+=surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase);
				    }
				  sum1+=sumk*sumap/kpmo;
				  sum2+=sumk;
				}


			    }    
		      kcond=surrogate[i].Ag-sum1/sum2;
		    }

		  if (sum_mass>0.)
		    kcond=(1.0-AQinit(b)/sum_mass)*(kcond)/surrogate[i].tau_air(b);
		  else
		    kcond=(kcond)/surrogate[i].tau_air(b);
		  b2=0.0;
		  c2=0.0;
		  if (eq_last_layer)
		    {
		      for (ilayer=0;ilayer<config.nlayer;++ilayer)
			for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)		
			  if (surrogate[i].time(b,ilayer,iphase)>=config.tequilibrium)		  	  
			    if (surrogate[i].k1(b,ilayer,iphase,index)>0.0 and ktot>0.0 and kcond>0.0)			  
			      surrogate[i].k1(b,ilayer,iphase,index)=1.0/(1.0/(surrogate[i].k1(b,ilayer,iphase,index)/sumkpositive*kcond)+1.0/surrogate[i].k1(b,ilayer,iphase,index))		    
				-sumknegative*surrogate[i].k1(b,ilayer,iphase,index)/sumkpositive;			  
			    else if (surrogate[i].k1(b,ilayer,iphase,index)<0.0 and ktot<0.0 and kcond<0.0)			  
			      surrogate[i].k1(b,ilayer,iphase,index)=1.0/
											 (1.0/(surrogate[i].k1(b,ilayer,iphase,index)/sumknegative*kcond)
											  +1.0/surrogate[i].k1(b,ilayer,iphase,index))-sumkpositive*surrogate[i].k1(b,ilayer,iphase,index)/sumknegative;	      
		    }
		  else
		    {
		      jlayer;
		      sumk=0.0;
		      for (jlayer=ilayer_interface;jlayer<config.nlayer;jlayer++)
			for (jphase=0;jphase<config.nphase(b,jlayer);++jphase)
			  sumk+=surrogate[i].k1(b,jlayer,jphase,index);

		      for (ilayer=0;ilayer<config.nlayer;++ilayer)
			for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)		
			  if (surrogate[i].time(b,ilayer,iphase)>=config.tequilibrium)		  	  
			    if (surrogate[i].k1(b,ilayer,iphase,index)>0.0 and ktot>0.0 and kcond>0.0)
			      {			   			 
                                if (kcond-max(sumk,0.0)<=0.0)
                                  {
                                    if (ilayer>=ilayer_interface)
                                      a=config.Vlayer(ilayer)/Vinterface*surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)/kpmo_interface(ilayer)*(sumk-sumknegative);				  			      	
                                    else
                                      a=0.0;
                                  }
				else				  
				  if (ilayer>=ilayer_interface) // interface 
				    if (sumkpositive>kcond) 
				      a=config.Vlayer(ilayer)/Vinterface*surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)/kpmo_interface(ilayer)*(kcond-b2-sumknegative*sumk/sumkpositive);
				    else
				      a=surrogate[i].k1(b,ilayer,iphase,index)-sumknegative*surrogate[i].k1(b,ilayer,iphase,index)/sumkpositive;			      
				  else // inner layers
				    {                                  
                                      if(sumkpositive > 0.0) {
		  		        a=1.0/(1.0/(surrogate[i].k1(b,ilayer,iphase,index)/(sumkpositive-max(sumk,0.0))*(kcond-max(sumk,0.0)))+1.0/surrogate[i].k1(b,ilayer,iphase,index))-sumknegative*surrogate[i].k1(b,ilayer,iphase,index)/sumkpositive;
				        b2+=a;
				      }
                                      else
                                        a = 0;
				    }
				     	    
			    
				surrogate[i].k1(b,ilayer,iphase,index)=a;
			      }		    
			    else if (surrogate[i].k1(b,ilayer,iphase,index)<0.0 and ktot<0.0 and kcond<0.0)
			      {					
				if (kcond-min(sumk,0.0)>=0.0)
                                  {
                                    if (ilayer>=ilayer_interface)
				      surrogate[i].k1(b,ilayer,iphase,index)=config.Vlayer(ilayer)/Vinterface*surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)/kpmo_interface(ilayer)*(sumk-sumkpositive);				  
				    else 
				      surrogate[i].k1(b,ilayer,iphase,index)=0.0;
				  }
				else
                                  if (ilayer>=ilayer_interface)
                                    {
                                      if (sumknegative<kcond)
                                        surrogate[i].k1(b,ilayer,iphase,index)=config.Vlayer(ilayer)/Vinterface*surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)/kpmo_interface(ilayer)*(kcond-c2-sumkpositive*sumk/sumknegative);
                                      else 
                                        surrogate[i].k1(b,ilayer,iphase,index)=surrogate[i].k1(b,ilayer,iphase,index)-sumkpositive*surrogate[i].k1(b,ilayer,iphase,index)/sumknegative;
                                    }
                                  else
                                    {
				      surrogate[i].k1(b,ilayer,iphase,index)=1.0/
					(1.0/(surrogate[i].k1(b,ilayer,iphase,index)/(sumknegative-min(sumk,0.0))*(kcond-min(sumk,0.0)))
					 +1.0/surrogate[i].k1(b,ilayer,iphase,index))-sumkpositive*surrogate[i].k1(b,ilayer,iphase,index)/sumknegative;
                                      c2+=surrogate[i].k1(b,ilayer,iphase,index);
                                    }
			      }	  
		    }
		}

	}
      else // There is a single layer.
	{      
	  for (b=0;b<config.nbins;++b)
	    {	     
	      sum_mass=AQinit(b);
	      for (ilayer=0;ilayer<config.nlayer;ilayer++)
		for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		  sum_mass+=MOinit(b,ilayer,iphase);
	      sum_mass=max(sum_mass,config.MOmin);
	   
	      //compute kinetic rate of absorption	    
	      ilayer = 0;
	      if (AQinit(b) != sum_mass) 
		for (i=0;i<n;++i)
		  if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophobic)
		    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		      if(surrogate[i].time(b,ilayer,iphase)>=config.tequilibrium)
			{		      
			  sum=0.0;
			  for (jphase=0;jphase<config.nphase(b,ilayer);++jphase)
			    sum+=surrogate[i].Kp(b,ilayer,jphase)*MOinit(b,ilayer,jphase);
			  if (sum > 0.0 and surrogate[i].tau_air(b) > 0.0)
			    surrogate[i].k1(b,ilayer,iphase,index)=			
			      (surrogate[i].Ag*surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)
			       -surrogate[i].Ap_layer_init(b,ilayer,iphase))/
			      (sum/(1.0-AQinit(b)/sum_mass)*surrogate[i].tau_air(b));
			}
	    }
	}
    }

  for (i=0;i<n;++i)      
    if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophobic)
      {
	sumkpositive=0.0;
	sumknegative=0.0;
	for (b=0;b<config.nbins;b++)
	  for (ilayer=0;ilayer<config.nlayer;++ilayer)
	    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
	      if (surrogate[i].time(b,ilayer,iphase)>=config.tequilibrium)			
		{
		  if (surrogate[i].k1(b,ilayer,iphase,index)>0.0)  
		    sumkpositive+=surrogate[i].k1(b,ilayer,iphase,index);
		  else if (surrogate[i].k1(b,ilayer,iphase,index)<0.0) 
		    sumknegative+=surrogate[i].k1(b,ilayer,iphase,index);
		}
		
	for (b=0;b<config.nbins;b++)
	  for (ilayer=0;ilayer<config.nlayer;ilayer++)
	    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
	      {
		//Flux limitations				
		if (surrogate[i].k1(b,ilayer,iphase,index)>0.0 and surrogate[i].Ag*surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)-surrogate[i].Ap_layer_init(b,ilayer,iphase)>0.0)
		  surrogate[i].k1(b,ilayer,iphase,index)=min(surrogate[i].k1(b,ilayer,iphase,index),1.0/config.deltatmin*
							     ((surrogate[i].Ag+max(-sumknegative-sumkpositive,0.0)*config.deltatmin+surrogate[i].Ap_layer_init(b,ilayer,iphase))*surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)/
							      (1.0+surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase))-surrogate[i].Ap_layer_init(b,ilayer,iphase)));
		else if (surrogate[i].k1(b,ilayer,iphase,index)<0.0 and surrogate[i].Ag*surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)-surrogate[i].Ap_layer_init(b,ilayer,iphase)<0.0)
		  surrogate[i].k1(b,ilayer,iphase,index)=max(surrogate[i].k1(b,ilayer,iphase,index),1.0/config.deltatmin*
							     ((surrogate[i].Ag+min(sumknegative+sumkpositive,0.0)*config.deltatmin+surrogate[i].Ap_layer_init(b,ilayer,iphase))*surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)/
							      (1.0+surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase))-surrogate[i].Ap_layer_init(b,ilayer,iphase)));
	      }
      }

  for (i=0;i<n;++i)      
    if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophobic)
      {
	for (b=0;b<config.nbins;b++)
	  for (ilayer=0;ilayer<config.nlayer;ilayer++)
	    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
	      {		
		//compute fluxes
		surrogate[i].Jdn(b,ilayer,iphase,index)=0.0;
		if (surrogate[i].k1(b,ilayer,iphase,index)<0.0 and
		    surrogate[i].Ap_layer_init(b,ilayer,iphase)>tiny)
		  surrogate[i].Jdn(b,ilayer,iphase,index)=surrogate[i].k1(b,ilayer,iphase,index)/
		    surrogate[i].Ap_layer_init(b,ilayer,iphase);
	      }
      }
}


void prodloss_org_ssh(model_config &config, vector<species>& surrogate,
		      Array<double, 3> &MOinit, Array<double, 1> &AQinit,
		      double &tiny, int index, double deltat)
{
  //compute kinetic rates for the absorption of a compound in the organic phase in a bin and a layer
  //index = 0 : first evaluation of rates)
  //index = 1 : second evaluation of rates
  //double Kp;
  int n=surrogate.size();
  int i,b,ilayer,iphase,jphase;
  double sum,sum_mass;
  Array<double, 1> kpmo_interface;
  kpmo_interface.resize(config.nlayer);
  Array <double,2> rJ_interface;
  rJ_interface.resize(config.nlayer,config.max_number_of_phases);
  double sumkpositive=0.0;
  double sumknegative=0.0;
  double sumkpmo_interface=0.0;
  double sumk=0.0;
  double sum1=0.0;
  double sum2=0.0;
  double sumap=0.0;
  double kpmo=0.0;
  double F1=0.0;
  double F2=0.0;
  double F1b,F2b;
  double surf1=0.0;
  double surf2=0.0;
  double dorg=0.0;
  bool eq_last_layer=true;
  double tau_interface= 0.0;
  double ap_interface=0.0;		  
  double Vinterface=0.0;
  int ilayer_interface=0;		  
  double Jinterface=0.0;
  double ktot1=0.0;
  double kcond=0.0;	    	      
  double ktot=0.0;
  double b2=0.0;
  double c2=0.0;
  double a=0.0;
  int jlayer;

  
  for (i=0;i<n;++i)
    {
      //surrogate[i].kloss_gas=0.0;
      //surrogate[i].kprod_gas=0.0;    
      surrogate[i].kprod=0.0;
      surrogate[i].kloss=0.0;       
    }
  
  if (config.explicit_representation)
    {
      for (i=0;i<n;++i)
	if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophobic)
	  for (b=0;b<config.nbins;++b)
	    {	     
	      sum_mass=AQinit(b);
	      for (ilayer=0;ilayer<config.nlayer;ilayer++)
		for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		  {
		    sum_mass+=MOinit(b,ilayer,iphase);
		    surrogate[i].k1(b,ilayer,iphase,index)=0.0;
		  }
	   
	      //compute kinetic rate of absorption	    
	      for (ilayer=0;ilayer<config.nlayer;ilayer++)
		for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		  if(ilayer==config.nlayer-1 and surrogate[i].time(b,ilayer,iphase)>=config.tequilibrium)
		    {		      
		      sum=0.0;
		      for (jphase=0;jphase<config.nphase(b,ilayer);++jphase)
			sum+=surrogate[i].Kp(b,ilayer,jphase)*MOinit(b,ilayer,jphase);

		      double tau_airloc=surrogate[i].tau_air(b);
		      if (surrogate[i].time(b,ilayer,iphase)<config.tequilibrium)
			tau_airloc=surrogate[i].tau_air(b)*config.tequilibrium/surrogate[i].time(b,ilayer,iphase);

		      surrogate[i].kprod(b,ilayer,iphase)+=surrogate[i].Ag*surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)/
			(sum/(1.0-AQinit(b)/sum_mass)*tau_airloc);
		      surrogate[i].kloss(b,ilayer,iphase)+=1./(sum/(1.0-AQinit(b)/sum_mass)*tau_airloc);
		      surrogate[i].kprod_gas+=surrogate[i].Ap_layer_init(b,ilayer,iphase)/
			(sum/(1.0-AQinit(b)/sum_mass)*tau_airloc);
		      surrogate[i].kprod_gas+=surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)/
                        (sum/(1.0-AQinit(b)/sum_mass)*tau_airloc);

		      if (index==0)
			surrogate[i].k1(b,ilayer,iphase,index)=
			  (surrogate[i].Ag*surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)
			   -surrogate[i].Ap_layer_init(b,ilayer,iphase))/
			  (sum/(1.0-AQinit(b)/sum_mass)*tau_airloc);					      		     
		    }

	    
	      for (ilayer=0;ilayer<config.nlayer;++ilayer)		    
		{
		  F1=0.0;
		  F2=0.0;
		  F1b=0.0;
		  F2b=0.0;
		  sum=0.0;
		  for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		    sum+=surrogate[i].Ap_layer_init(b,ilayer,iphase)/(surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase));

		  if (ilayer>0)
		    {
		      surf1=4.0*3.14159*pow(config.dbound(b,ilayer),2);
		      sum1=0.0;
		      for (iphase=0;iphase<config.nphase(b,ilayer-1);++iphase)
			sum1+=surrogate[i].Ap_layer_init(b,ilayer-1,iphase)/(surrogate[i].Kp(b,ilayer-1,iphase)*MOinit(b,ilayer-1,iphase));

		      dorg=(config.dbound(b,ilayer)-config.Radius(b,ilayer-1))/(config.Radius(b,ilayer)-config.Radius(b,ilayer-1))*surrogate[i].dif_org(b,ilayer-1)
			+(config.Radius(b,ilayer)-config.dbound(b,ilayer))/(config.Radius(b,ilayer)-config.Radius(b,ilayer-1))*surrogate[i].dif_org(b,ilayer);
			
		      F1=surf1*dorg*sum1/(config.Radius(b,ilayer)-config.Radius(b,ilayer-1));
		      F1b=surf1*dorg*sum/(config.Radius(b,ilayer)-config.Radius(b,ilayer-1));
		    }

		  if (ilayer<config.nlayer-1)
		    {
		      surf2=4.0*3.14159*pow(config.dbound(b,ilayer+1),2);
		      sum2=0.0;
		      for (iphase=0;iphase<config.nphase(b,ilayer+1);++iphase)
			sum2+=surrogate[i].Ap_layer_init(b,ilayer+1,iphase)/(surrogate[i].Kp(b,ilayer+1,iphase)*MOinit(b,ilayer+1,iphase));

		      dorg=(config.dbound(b,ilayer+1)-config.Radius(b,ilayer))/(config.Radius(b,ilayer+1)-config.Radius(b,ilayer))*surrogate[i].dif_org(b,ilayer)
			+(config.Radius(b,ilayer+1)-config.dbound(b,ilayer+1))/(config.Radius(b,ilayer+1)-config.Radius(b,ilayer))*surrogate[i].dif_org(b,ilayer+1);
			
		      F2=surf2*dorg*sum2/(config.Radius(b,ilayer+1)-config.Radius(b,ilayer));
		      F2b=surf2*dorg*sum/(config.Radius(b,ilayer+1)-config.Radius(b,ilayer));
		      
		    }
		  
		  sum=0.0;
		  for (jphase=0;jphase<config.nphase(b,ilayer);++jphase)
		    sum+=surrogate[i].Kp(b,ilayer,jphase)*MOinit(b,ilayer,jphase);

		  for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		    if (surrogate[i].time(b,ilayer,iphase)>=config.tequilibrium)
		      {
			surrogate[i].kprod(b,ilayer,iphase)+=pow(surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase),2.0)/sum*(F1+F2)/(4.0/3.0*3.14159*pow(config.diameters(b)*0.5e-6,3.0)*config.Vlayer(ilayer));
			if (surrogate[i].Ap_layer_init(b,ilayer,iphase)>0.0)
			  surrogate[i].kloss(b,ilayer,iphase)+=pow(surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase),2.0)/sum*(F1b+F2b)/(4.0/3.0*3.14159*pow(config.diameters(b)*0.5e-6,3.0)*config.Vlayer(ilayer))/surrogate[i].Ap_layer_init(b,ilayer,iphase);
			if (index==0)
			  surrogate[i].k1(b,ilayer,iphase,index)+=pow(surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase),2.0)/sum*(F1+F2-F1b-F2b)/(4.0/3.0*3.14159*pow(config.diameters(b)*0.5e-6,3.0)*config.Vlayer(ilayer));
		      }
		    
		}
	    }
    }
  else // Use implicit representation (low number of layers).
    {
      //cout << config.nphase << " " << config.nbins << endl;
      if (config.nlayer>1)
	{	  
	  for (i=0;i<n;i++)
	    if (surrogate[i].is_organic or i==config.iH2O)
	      if (surrogate[i].hydrophobic)
		for (b=0;b<config.nbins;b++)
		  for (ilayer=0;ilayer<config.nlayer-1;ilayer++)
		    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		      {
			surrogate[i].kprod(b,ilayer,iphase)+=surrogate[i].Ap_layer_init(b,config.nlayer-1,iphase)/(surrogate[i].Kp(b,config.nlayer-1,iphase)*MOinit(b,config.nlayer-1,iphase))
			  *surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)/max(surrogate[i].tau_diffusion(b,ilayer,iphase),config.tequilibrium);
			surrogate[i].kloss(b,ilayer,iphase)+=1./max(surrogate[i].tau_diffusion(b,ilayer,iphase),config.tequilibrium);
		    
			surrogate[i].kprod(b,config.nlayer-1,iphase)+=surrogate[i].Ap_layer_init(b,ilayer,iphase)/max(surrogate[i].tau_diffusion(b,ilayer,iphase),config.tequilibrium);
			surrogate[i].kloss(b,config.nlayer-1,iphase)+=1./(surrogate[i].Kp(b,config.nlayer-1,iphase)*MOinit(b,config.nlayer-1,iphase))
			  *surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)/max(surrogate[i].tau_diffusion(b,ilayer,iphase),config.tequilibrium);
			//cout << surrogate[i].kprod(b,config.nlayer-1,iphase) << endl;
		      }
	  
	  for (b=0;b<config.nbins;++b)
	    {	     
	      sum_mass=AQinit(b);
	      for (ilayer=0;ilayer<config.nlayer;ilayer++)
		for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		  sum_mass+=MOinit(b,ilayer,iphase);
	      sum_mass=max(sum_mass,config.MOmin);
	      
	      //compute kinetic rate of absorption	    
	      ilayer = config.nlayer-1;	     
	      if (AQinit(b) != sum_mass) 
		for (i=0;i<n;++i)
		  if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophobic)
		    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		      {
			double tau_airloc=surrogate[i].tau_air(b);
			if (surrogate[i].time(b,ilayer,iphase)<config.tequilibrium)
			  {
			    //cout << surrogate[i].name << " " << surrogate[i].time << endl;
			    tau_airloc=surrogate[i].tau_air(b)*config.tequilibrium/surrogate[i].time(b,ilayer,iphase);
			  }
			
			sum=0.0;
			for (jphase=0;jphase<config.nphase(b,ilayer);++jphase)
			  sum+=surrogate[i].Kp(b,ilayer,jphase)*MOinit(b,ilayer,jphase);
			
			surrogate[i].kprod(b,ilayer,iphase)+=surrogate[i].Ag*surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)
			  /(sum/(1.0-AQinit(b)/sum_mass)*tau_airloc);
			surrogate[i].kloss(b,ilayer,iphase)+=1./(sum/(1.0-AQinit(b)/sum_mass)*tau_airloc);
			surrogate[i].kloss_gas+=surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)
			  /(sum/(1.0-AQinit(b)/sum_mass)*tau_airloc);
			surrogate[i].kprod_gas+=surrogate[i].Ap_layer_init(b,ilayer,iphase)/(sum/(1.0-AQinit(b)/sum_mass)*tau_airloc);			
			if (sum > 0.0 and tau_airloc > 0.0 and index==0)
			  surrogate[i].k1(b,ilayer,iphase,index)=			
			    (surrogate[i].Ag*surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)
			     -surrogate[i].Ap_layer_init(b,ilayer,iphase))/
			    (sum/(1.0-AQinit(b)/sum_mass)*tau_airloc);			
		      }
	    }
	}
      else // There is a single layer.
	{      
	  for (b=0;b<config.nbins;++b)
	    {	     
	      sum_mass=AQinit(b);
	      for (ilayer=0;ilayer<config.nlayer;ilayer++)
		for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		  sum_mass+=MOinit(b,ilayer,iphase);
	      sum_mass=max(sum_mass,config.MOmin);
	   
	      //compute kinetic rate of absorption	    
	      ilayer = 0;
	      if (AQinit(b) != sum_mass) 
		for (i=0;i<n;++i)
		  if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophobic)
		    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		      {
			double tau_airloc=surrogate[i].tau_air(b);
			if (surrogate[i].time(b,ilayer,iphase)<config.tequilibrium)
			  tau_airloc=surrogate[i].tau_air(b)*config.tequilibrium/surrogate[i].time(b,ilayer,iphase);
	      
			sum=0.0;
			for (jphase=0;jphase<config.nphase(b,ilayer);++jphase)
			  sum+=surrogate[i].Kp(b,ilayer,jphase)*MOinit(b,ilayer,jphase);

			surrogate[i].kprod(b,ilayer,iphase)=surrogate[i].Ag*surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)
			  /(sum/(1.0-AQinit(b)/sum_mass)*tau_airloc);
			surrogate[i].kloss(b,ilayer,iphase)=1./(sum/(1.0-AQinit(b)/sum_mass)*tau_airloc);
			surrogate[i].kloss_gas+=surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)
			  /(sum/(1.0-AQinit(b)/sum_mass)*tau_airloc);
			surrogate[i].kprod_gas+=surrogate[i].Ap_layer_init(b,ilayer,iphase)/(sum/(1.0-AQinit(b)/sum_mass)*tau_airloc);			
			if (sum > 0.0 and surrogate[i].tau_air(b) > 0.0 and index==0)
			  surrogate[i].k1(b,ilayer,iphase,index)=			
			    (surrogate[i].Ag*surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)
			     -surrogate[i].Ap_layer_init(b,ilayer,iphase))/
			    (sum/(1.0-AQinit(b)/sum_mass)*tau_airloc);
			
			
		      }
	    }
	}
    }
}

void error_ph_bins_ssh(model_config &config, vector<species> &surrogate, int index_b, double Temperature, Array <double, 1> chp, 
		       double &error, double &derivative, Array <double, 1> AQinit, Array<double , 1> &ionic, Array<double, 1> &MMaq,
		       Array<double, 1> &LWC,
		       double &chp_new, double &conc_organion, double deltat)
{      
  int n=surrogate.size();
  int i;
  double inorganion=0.0;
  double organion=0.0;
  double total,dtotal;
  double sum_K_AQ;
  int b;
  //double kelvin_effect=1.0;
  double Kaq;
  double ratio_gamma1=0.0;
  double ratio_gamma2=0.0;
  double Kac1=0.0;
  double Kac2=0.0;
  double Kh=0.0;
  double dK=0.0;
  double fion1=0.0;
  double fion2=0.0;
  double K, ratio_gamma, Kac;
  /*
  if (config.imethod==0)
    if (config.compute_kelvin_effect) //compute the kelvin effect    
      {
        kelvin_effect=2.0*config.surface_tension_aq*MMaq(index_b)/
          (8.314*Temperature*config.AQrho(index_b)*
           0.5*config.diameters(index_b));
        if(kelvin_effect > 50.0)
          kelvin_effect = 50.0;
        kelvin_effect=exp(kelvin_effect);
      }*/

  //if (config.compute_organic)
  organion=conc_organion;
  //cout << organion << endl;
  
  derivative=0.0;
  //cout << "ici " << deltat << endl;
  if (config.imethod>0)
    {
      double conc_org;
      compute_conc_org_bins_ssh(config, surrogate, Temperature,  MMaq, AQinit, chp, ionic, LWC, index_b, conc_org);
      
      //cout << "laaaaa" << endl;
      for (i=0;i<n;++i)
	if (surrogate[i].is_inorganic_precursor and surrogate[i].is_solid==false)
	  {
	    if (i==config.iNH3)
	      {
                total=(surrogate[config.iNH4p].Aaq_bins_init0(index_b)
                       +deltat*surrogate[config.iNH4p].kprod_aq(index_b))
                  /(1.0+deltat*surrogate[config.iNH4p].kloss_aq(index_b))
                  /conc_org*1000.*surrogate[config.iNH4p].charge/surrogate[config.iNH4p].MM;

		dtotal=-total/(1.0+deltat*surrogate[config.iNH4p].kloss_aq(index_b))*
		  deltat*surrogate[config.iNH4p].kloss_aq(index_b)/chp(index_b);
		derivative+=dtotal;
		inorganion-=total;
		//cout << "NH3: " << total << endl;
	      }
		
	    //cout << "NH3: " << inorganion << endl;
          
	
	    else if (i==config.iHNO3)
	      {
                total=(surrogate[config.iNO3m].Aaq_bins_init0(index_b)
                       +deltat*surrogate[config.iNO3m].kprod_aq(index_b))
                  /(1.0+deltat*surrogate[config.iNO3m].kloss_aq(index_b))
                  /conc_org*1000.*surrogate[config.iNO3m].charge/surrogate[config.iNO3m].MM;

		dtotal=total/(1.0+deltat*surrogate[config.iNO3m].kloss_aq(index_b))*
		  deltat*surrogate[config.iNO3m].kloss_aq(index_b)/chp(index_b);
                
		derivative+=dtotal;
		inorganion-=total; 
		//cout << "NO3: " << total << endl;  
	      }
	    else if (i==config.iCO2)
              {
		//conc_aq=surrogate[config.iNO3m].Aaq_bins_init(b)/surrogate[config.iNO3m].MM*surrogate[i].MM;
		double A=surrogate[config.iCO2].keq2/(surrogate[config.iHp].gamma_aq_bins(index_b)*surrogate[config.iCO3mm].gamma_aq_bins(index_b)/surrogate[config.iHCO3m].gamma_aq_bins(index_b));

		//frac1: fraction of carbonates as HCO3-
		double frac1=1.0/(1.0+A/chp(index_b));

		//dfrac1: derivative of frac1. Set to as as it seems to lead to instabilities
		double dfrac1=0.; //1.0/pow(1.0+A/chp(index_b),2.)*A/pow(chp(index_b),2);


		
		//double dfrac1=surrogate[config.iCO2].keq2/(pow(chp(index_b),2)*surrogate[config.iHp].gamma_aq_bins(index_b)*surrogate[config.iCO3mm].gamma_aq_bins(index_b)/surrogate[config.iHCO3m].gamma_aq_bins(index_b))
		//  /pow(1.0+(surrogate[config.iCO2].keq2/(chp(index_b)*surrogate[config.iHp].gamma_aq_bins(index_b)*surrogate[config.iCO3mm].gamma_aq_bins(index_b)/surrogate[config.iHCO3m].gamma_aq_bins(index_b))),2);
		//frac1=0.03;
		//dfrac1=0;
		
		double temp=surrogate[config.iCO2].kaqi*surrogate[config.iH2O].Aaq_bins_init(index_b)/surrogate[config.iH2O].MM*MMaq(index_b)/AQinit(index_b)*surrogate[config.iH2O].gamma_aq_bins(index_b)*(surrogate[config.iCO2].keq/(chp(index_b)*surrogate[config.iHp].gamma_aq_bins(index_b)*surrogate[config.iHCO3m].gamma_aq_bins(index_b)))*
		    (1.+(surrogate[config.iCO2].keq2/(chp(index_b)*surrogate[config.iHp].gamma_aq_bins(index_b)*surrogate[config.iCO3mm].gamma_aq_bins(index_b)/surrogate[config.iHCO3m].gamma_aq_bins(index_b))));

		double Keffect=surrogate[i].Kaq(index_b)/temp;
		
		/*surrogate[i].Kaq(b)=surrogate[config.iCO2].kaqi*surrogate[config.iH2O].Aaq_bins_init(b)/surrogate[config.iH2O].MM*MMaq(b)/AQinit(b)*surrogate[config.iH2O].gamma_aq_bins(b)*(surrogate[config.iCO2].keq/(chp(b)*surrogate[config.iHp].gamma_aq_bins(b)*surrogate[config.iHCO3m].gamma_aq_bins(b)))*
		  (1.+(surrogate[config.iCO2].keq2/(chp(b)*surrogate[config.iHp].gamma_aq_bins(b)*surrogate[config.iCO3mm].gamma_aq_bins(b)/surrogate[config.iHCO3m].gamma_aq_bins(b))));*/
		
		Kaq=surrogate[i].Kaq(index_b);		
		double dKaq=-surrogate[i].Kaq(index_b)/chp(index_b)
		  -Keffect*surrogate[config.iCO2].kaqi*surrogate[config.iH2O].Aaq_bins_init(index_b)/surrogate[config.iH2O].MM*MMaq(index_b)/AQinit(index_b)*surrogate[config.iH2O].gamma_aq_bins(index_b)*(surrogate[config.iCO2].keq/(chp(index_b)*surrogate[config.iHp].gamma_aq_bins(index_b)*surrogate[config.iHCO3m].gamma_aq_bins(index_b))*(surrogate[config.iCO2].keq2/(pow(chp(index_b),2)*surrogate[config.iHp].gamma_aq_bins(index_b)*surrogate[config.iCO3mm].gamma_aq_bins(index_b)/surrogate[config.iHCO3m].gamma_aq_bins(index_b))));

		double total2b=0.;
		double total0=0.;
		
		if (config.solids)
		  for (int j=0;j<n;j++)
		    if (surrogate[j].is_solid)
		      {			
			for (int j2=0;j2<surrogate[j].nion;j2++)
			  {
			    //cout << j2 << " " << surrogate[j].iion1 << " " << surrogate[j].iion2 << " " << config.iHCO3m << " " << condifendl;
			    if (j2==0 and (surrogate[j].iion1==config.iHCO3m or surrogate[j].iion1==config.iCO3mm))
			      {
				total0+=(surrogate[j].Asol_bins_init0(index_b))/surrogate[j].MM*surrogate[j].pion1;
				total2b+=(surrogate[j].Asol_bins_init(index_b))/surrogate[j].MM*surrogate[j].pion1;
			      }
			    else if (j2==1 and (surrogate[j].iion2==config.iHCO3m or surrogate[j].iion2==config.iCO3mm))
			      {			      
				total0+=(surrogate[j].Asol_bins_init0(index_b))/surrogate[j].MM*surrogate[j].pion2;
				total2b+=(surrogate[j].Asol_bins_init(index_b))/surrogate[j].MM*surrogate[j].pion2;
			      }
			    else if (j2==2 and (surrogate[j].iion3==config.iHCO3m or surrogate[j].iion3==config.iCO3mm))
			      {
				total0+=(surrogate[j].Asol_bins_init0(index_b))/surrogate[j].MM*surrogate[j].pion3;
				total2b+=(surrogate[j].Asol_bins_init(index_b))/surrogate[j].MM*surrogate[j].pion3;
			      }
			  }
		      }
		
		double total1b=surrogate[config.iHCO3m].Aaq_bins_init(index_b)/surrogate[config.iHCO3m].MM+surrogate[config.iCO3mm].Aaq_bins_init(index_b)/surrogate[config.iHCO3m].MM;
		double frac2=(total1b)/(total2b+total1b);		
		total=(frac2*(surrogate[config.iHCO3m].Aaq_bins_init0(index_b)+frac1*total0*surrogate[config.iHCO3m].MM)
                       +deltat*surrogate[config.iHCO3m].kprod_aq(index_b))
                  /(1.0+deltat*surrogate[config.iHCO3m].kloss_aq(index_b))
                  /conc_org*1000./surrogate[config.iHCO3m].MM;


		
		double total1=total;
		dtotal=total/(1.0+deltat*surrogate[config.iHCO3m].kloss_aq(index_b))*deltat*surrogate[config.iHCO3m].kloss_aq(index_b)/Kaq*dKaq;
		
		double total2=(frac2*(surrogate[config.iCO3mm].Aaq_bins_init0(index_b)+(1-frac1)*total0*surrogate[config.iCO3mm].MM)
                       +deltat*surrogate[config.iCO3mm].kprod_aq(index_b))
                  /(1.0+deltat*surrogate[config.iCO3mm].kloss_aq(index_b))
                  /conc_org*1000./surrogate[config.iCO3mm].MM;
		total+=total2;
		dtotal+=total2/(1.0+deltat*surrogate[config.iCO3mm].kloss_aq(index_b))*deltat*surrogate[config.iCO3mm].kloss_aq(index_b)/Kaq*dKaq;
		//dtotal+=total2/(1.0+deltat*surrogate[config.iCO3mm].kloss_aq(index_b))*deltat*surrogate[config.iCO3mm].kloss_aq(index_b)/(Kaq*conc_org*surrogate[i].tau_air(index_b))*dKaq;		
		total=total*(1.-frac1)*surrogate[config.iCO3mm].charge+total*frac1*surrogate[config.iHCO3m].charge;
		dtotal=dtotal*(1.-frac1)*surrogate[config.iCO3mm].charge+dtotal*frac1*surrogate[config.iHCO3m].charge-total*dfrac1*surrogate[config.iCO3mm].charge+total*dfrac1*surrogate[config.iHCO3m].charge;

		derivative-=dtotal;
		inorganion-=total;
	      }	       	      
	    else if (i==config.iHCl)
	      {
                 total=(surrogate[config.iClm].Aaq_bins_init0(index_b)
                       +deltat*surrogate[config.iClm].kprod_aq(index_b))
                  /(1.0+deltat*surrogate[config.iClm].kloss_aq(index_b))
                  /conc_org*1000.*surrogate[config.iClm].charge/surrogate[config.iClm].MM;

		dtotal=total/(1.0+deltat*surrogate[config.iClm].kloss_aq(index_b))*
		  deltat*surrogate[config.iClm].kloss_aq(index_b)/chp(index_b);
                
		derivative+=dtotal;
		inorganion-=total; 
	      }
	    else if (i==config.iH2SO4)
	      {
                total=surrogate[config.iHSO4m].Aaq_bins_init(index_b)/surrogate[config.iHSO4m].MM
                  +surrogate[config.iSO4mm].Aaq_bins_init(index_b)/surrogate[config.iSO4mm].MM;
                
                K=surrogate[i].keqi*surrogate[config.iHSO4m].gamma_aq_bins(index_b)
                  /(surrogate[config.iHp].gamma_aq_bins(index_b)*surrogate[config.iSO4mm].gamma_aq_bins(index_b));
                
                derivative-=1000.*total/conc_org*1.0/pow(1.0+K/chp(index_b),2.0)*K/(chp(index_b)*chp(index_b)); //HSO4m+SO4mm
                inorganion+=1000.*total/conc_org*(2.0-1.0/(1.0+K/chp(index_b)));
		//cout << "H2SO4: " << inorganion << " " << total << " " << conc_org << " " << endl;
	      }
	  }
      /*
        else if (surrogate[i].is_organic and config.compute_organic and surrogate[i].hydrophilic)
          {
            if (surrogate[i].aqt==1) //_type=="monoacid")
              {
              
                total=surrogate[i].Aaq_bins_init(index_b)/surrogate[i].MM;
                fion1=0.0;
                fion2=0.0;
                surrogate[i].gamma_LR=surrogate[i].LR(index_b);
                Kaq=surrogate[i].Kp_eff_aqrealdyn_ssh(config, Temperature, ionic(index_b), chp(index_b),surrogate[config.iHp].LR(index_b),
                                                      surrogate[config.iHp].SRMR(index_b),MMaq(index_b), fion1, fion2,index_b)/kelvin_effect/surrogate[i].gamma_aq_bins(index_b);
                organion+=1000.*fion1*total/conc_org(index_b);
                derivative-=fion1*(1-fion1)/chp(index_b)*1000.*total/conc_org(index_b);
              }
            else if (surrogate[i].aqt==2) //_type=="diacid")
              {
                  
                total=surrogate[i].Aaq_bins_init(index_b)/surrogate[i].MM;
                fion1=0.0;
                fion2=0.0;
                surrogate[i].gamma_LR=surrogate[i].LR(index_b);	      
                Kaq=surrogate[i].Kp_eff_aqrealdyn_ssh(config, Temperature, ionic(index_b), chp(index_b),surrogate[config.iHp].LR(index_b),
                                                      surrogate[config.iHp].SRMR(index_b),MMaq(index_b), fion1, fion2,index_b)/kelvin_effect/surrogate[i].gamma_aq_bins(index_b);
                organion+=1000.*(fion1+2.0*fion2)*total/conc_org(index_b);
                derivative+=(fion1*(fion1-1)+4.0*fion1*fion2+4.0*fion2*(fion2-1))/chp(index_b)*1000.*total/conc_org(index_b);
                
              }
              }*/

      if (config.iNa>=0) inorganion-=surrogate[config.iNa].Aaq_bins_init(index_b)/surrogate[config.iNa].MM/conc_org*1000.*surrogate[config.iNa].charge;
      if (config.iMg>=0) inorganion-=surrogate[config.iMg].Aaq_bins_init(index_b)/surrogate[config.iMg].MM/conc_org*1000.*surrogate[config.iMg].charge;
      if (config.iCa>=0) inorganion-=surrogate[config.iCa].Aaq_bins_init(index_b)/surrogate[config.iCa].MM/conc_org*1000.*surrogate[config.iCa].charge;  
      if (config.iK>=0) inorganion-=surrogate[config.iK].Aaq_bins_init(index_b)/surrogate[config.iK].MM/conc_org*1000.*surrogate[config.iK].charge;

      //cout << "in: "<< inorganion << " " << derivative << " " << endl;

    }
  else
    {
      Array <double, 1> conc_org;
      conc_org.resize(config.nbins);
      
      for (b=0;b<config.nbins;b++)
        compute_conc_org_bins_ssh(config, surrogate, Temperature,  MMaq, AQinit, chp, ionic, LWC, index_b, conc_org(b));
    for (i=0;i<n;++i)
      if (surrogate[i].is_inorganic_precursor and surrogate[i].is_solid==false)
        {
          if (i==config.iNH3)
            {
              if (surrogate[i].time_aq(index_b)<config.tequilibrium)
                {
                  total=surrogate[i].Ag/surrogate[i].MM;
		  sum_K_AQ=1.0;
		  for (b=0;b<config.nbins;b++)
		    if (surrogate[i].time_aq(b)<config.tequilibrium)
		      {
			total+=surrogate[config.iNH4p].Aaq_bins_init(b)/surrogate[config.iNH4p].MM;
			sum_K_AQ+=surrogate[i].Kaq(b)*conc_org(b);
		      }
		  derivative-=1000.0*total*
		    (surrogate[i].Kaq(index_b)/chp(index_b)/sum_K_AQ
		     -surrogate[i].Kaq(index_b)/pow(sum_K_AQ,2.0)*surrogate[i].Kaq(index_b)/chp(index_b)*conc_org(index_b));
		  inorganion-=1000.0*total*surrogate[i].Kaq(index_b)/sum_K_AQ; 
		}
              else
		inorganion-=surrogate[config.iNH4p].Aaq_bins_init(index_b)/surrogate[config.iNH4p].MM/conc_org(index_b)*1000.0*surrogate[config.iNH4p].charge;
	      //cout << "NH3: " << inorganion << endl;
	    }
	
	  else if (i==config.iHNO3)
	    {
	      if (surrogate[i].time_aq(index_b)<config.tequilibrium)
		{
		  total=surrogate[i].Ag/surrogate[i].MM;
		  sum_K_AQ=1.0;
		  for (b=0;b<config.nbins;b++)
		    if (surrogate[i].time_aq(b)<config.tequilibrium)
		      {
			total+=surrogate[config.iNO3m].Aaq_bins_init(b)/surrogate[config.iNO3m].MM;
			sum_K_AQ+=surrogate[i].Kaq(b)*conc_org(b);
		      }
                
		  derivative+=1000.0*total*
		    (-surrogate[i].Kaq(index_b)/chp(index_b)/sum_K_AQ
		     +surrogate[i].Kaq(index_b)/pow(sum_K_AQ,2.0)*surrogate[i].Kaq(index_b)/chp(index_b)*conc_org(index_b));
		  inorganion+=1000.*total*surrogate[i].Kaq(index_b)/sum_K_AQ;
		}
	      else
		inorganion-=surrogate[config.iNO3m].Aaq_bins_init(index_b)/surrogate[config.iNO3m].MM/conc_org(index_b)*1000.*surrogate[config.iNO3m].charge;
	      //cout << "NO3: " << inorganion << endl;
	    }
	  else if (i==config.iCO2)
              {
		cout << "Carbonates Inorganics in dynamic with imethod=0. Should not happen" << endl;
		exit(1);
	      }
	  else if (i==config.iHCl)
	    {           
	      if (surrogate[i].time_aq(index_b)<config.tequilibrium)
		{
		  total=surrogate[i].Ag/surrogate[i].MM;
		  sum_K_AQ=1.0;
		  for (b=0;b<config.nbins;b++)
		    if (surrogate[i].time_aq(b)<config.tequilibrium)
		      {
			total+=surrogate[config.iClm].Aaq_bins_init(b)/surrogate[config.iClm].MM;
			sum_K_AQ+=surrogate[i].Kaq(b)*conc_org(b);
		      }
                
		  derivative+=1000.0*total*
		    (-surrogate[i].Kaq(index_b)/chp(index_b)/(sum_K_AQ)
		     +surrogate[i].Kaq(index_b)/pow(sum_K_AQ,2.0)*surrogate[i].Kaq(index_b)/chp(index_b)*conc_org(index_b));
		  inorganion+=1000.*total*surrogate[i].Kaq(index_b)/sum_K_AQ;
		}            
	      else
		{
                  inorganion-=surrogate[config.iClm].Aaq_bins_init(index_b)/surrogate[config.iClm].MM/conc_org(index_b)*1000.*surrogate[config.iClm].charge;
                }	
	    
	    }
	  else if (i==config.iH2SO4)
	    {
	      total=surrogate[config.iHSO4m].Aaq_bins_init(index_b)/surrogate[config.iHSO4m].MM
		+surrogate[config.iSO4mm].Aaq_bins_init(index_b)/surrogate[config.iSO4mm].MM;

	      K=surrogate[i].Kequilibrium_ssh(Temperature)*surrogate[config.iHSO4m].gamma_aq_bins(index_b)
		/(surrogate[config.iHp].gamma_aq_bins(index_b)*surrogate[config.iSO4mm].gamma_aq_bins(index_b));
	  
	      derivative-=1000.*total/conc_org(index_b)*1.0/pow(1.0+K/chp(index_b),2.0)*K/(chp(index_b)*chp(index_b)); //HSO4m+SO4mm
	      inorganion+=1000.*total/conc_org(index_b)*(2.0-1.0/(1.0+K/chp(index_b)));
	    }
	}
  /*
      else if (surrogate[i].is_organic and config.compute_organic and surrogate[i].hydrophilic)
	{
          
	  if (surrogate[i].aq_type=="monoacid")
	    {
	      if (surrogate[i].time_aq(index_b)<config.tequilibrium)
		{	       
		  total=surrogate[i].Ag/surrogate[i].MM;
		  sum_K_AQ=1.0;
		  for (b=0;b<config.nbins;b++)
		    if (surrogate[i].time_aq(b)<config.tequilibrium)
		      {
			total+=surrogate[i].Aaq_bins_init(b)/surrogate[i].MM;
			sum_K_AQ+=surrogate[i].Kaq(b)*AQinit(b);
		      }

		  fion1=0.0;
		  fion2=0.0;
		  surrogate[i].gamma_LR=surrogate[i].LR(index_b);
		  Kaq=surrogate[i].Kp_eff_aqrealdyn_ssh(config, Temperature, ionic(index_b), chp(index_b),surrogate[config.iHp].LR(index_b),
							surrogate[config.iHp].SRMR(index_b),MMaq(index_b), fion1, fion2,index_b)/kelvin_effect/surrogate[i].gamma_aq_bins(index_b);

		  ratio_gamma=pow(surrogate[config.iHp].LR(index_b),2.0)*surrogate[config.iHp].SRMR(index_b)/surrogate[i].gamma_LR;
		  Kac=surrogate[i].Kacidity1/ratio_gamma;
		  Kh=surrogate[i].Kpart_aq_ssh(Temperature, MMaq(index_b))/kelvin_effect/surrogate[i].gamma_aq_bins(index_b);
		  dK=-Kac*Kh/pow(chp(index_b),2.0);
		
		  derivative+=1000.*AQinit(index_b)/conc_org(index_b)*total*fion1*
		    (-(1.0-fion1)/chp(index_b)*surrogate[i].Kaq(index_b)/sum_K_AQ+dK/sum_K_AQ
		     -Kaq/pow(sum_K_AQ,2.0)*AQinit(index_b)*dK);                 
		  organion+=1000.*fion1*total*surrogate[i].Kaq(index_b)*AQinit(index_b)/sum_K_AQ/conc_org(index_b);		
		}
	      else
                {
              
		  total=surrogate[i].Aaq_bins_init(index_b)/surrogate[i].MM;
		  fion1=0.0;
		  fion2=0.0;
		  surrogate[i].gamma_LR=surrogate[i].LR(index_b);
		  Kaq=surrogate[i].Kp_eff_aqrealdyn_ssh(config, Temperature, ionic(index_b), chp(index_b),surrogate[config.iHp].LR(index_b),
							surrogate[config.iHp].SRMR(index_b),MMaq(index_b), fion1, fion2,index_b)/kelvin_effect/surrogate[i].gamma_aq_bins(index_b);
		  organion+=1000.*fion1*total/conc_org(index_b);
		  derivative-=fion1*(1-fion1)/chp(index_b)*1000.*total/conc_org(index_b);
                }

	    }
	  else if (surrogate[i].aq_type=="diacid")
	    {
	      if (surrogate[i].time_aq(index_b)<config.tequilibrium)
		{	       
		  total=surrogate[i].Ag/surrogate[i].MM;
		  sum_K_AQ=1.0;
		  for (b=0;b<config.nbins;b++)
		    if (surrogate[i].time_aq(b)<config.tequilibrium)
		      {
			total+=surrogate[i].Aaq_bins_init(b)/surrogate[i].MM;
			sum_K_AQ+=surrogate[i].Kaq(b)*AQinit(b);
		      }

		  fion1=0.0;
		  fion2=0.0;
		  surrogate[i].gamma_LR=surrogate[i].LR(index_b);
		  Kaq=surrogate[i].Kp_eff_aqrealdyn_ssh(config, Temperature, ionic(index_b), chp(index_b),surrogate[config.iHp].LR(index_b),
							surrogate[config.iHp].SRMR(index_b),MMaq(index_b), fion1, fion2,index_b)/kelvin_effect/surrogate[i].gamma_aq_bins(index_b);

		  ratio_gamma1=pow(surrogate[config.iHp].LR(index_b),2.0)*surrogate[config.iHp].SRMR(index_b)/surrogate[i].gamma_LR;
		  ratio_gamma2=pow(surrogate[config.iHp].LR(index_b),2.0)*surrogate[config.iHp].SRMR(index_b);
		  Kac1=surrogate[i].Kacidity1/ratio_gamma1;
		  Kac2=surrogate[i].Kacidity2/ratio_gamma2;
		  Kh=surrogate[i].Kpart_aq_ssh(Temperature, MMaq(index_b))/kelvin_effect/surrogate[i].gamma_aq_bins(index_b);
		  dK=-Kac1*Kh/pow(chp(index_b),2.0)-2.0*Kac1*Kac2*Kh/pow(chp(index_b),3.0);
		
		  derivative+=1000.*total*AQinit(index_b)/conc_org(index_b)*
		    ((fion1*(fion1-1.0)+4.0*fion1*fion2+4.0*fion2*(fion2-1))/chp(index_b)*surrogate[i].Kaq(index_b)/sum_K_AQ
		     +dK/sum_K_AQ-Kaq/pow(sum_K_AQ,2.0)*AQinit(index_b)*dK);                 
		  organion+=1000.0*(fion1+2.0*fion2)*total*surrogate[i].Kaq(index_b)*AQinit(index_b)/sum_K_AQ/conc_org(index_b);		
		}
	      else
                {
		  total=surrogate[i].Aaq_bins_init(index_b)/surrogate[i].MM;
		  fion1=0.0;
		  fion2=0.0;
		  surrogate[i].gamma_LR=surrogate[i].LR(index_b);	      
		  Kaq=surrogate[i].Kp_eff_aqrealdyn_ssh(config, Temperature, ionic(index_b), chp(index_b),surrogate[config.iHp].LR(index_b),
							surrogate[config.iHp].SRMR(index_b),MMaq(index_b), fion1, fion2,index_b)/kelvin_effect/surrogate[i].gamma_aq_bins(index_b);
		  organion+=1000.*(fion1+2.0*fion2)*total/conc_org(index_b);
		  derivative+=(fion1*(fion1-1)+4.0*fion1*fion2+4.0*fion2*(fion2-1))/chp(index_b)*1000.*total/conc_org(index_b);
                }
                }

                }*/

  if (config.iNa>=0) inorganion-=surrogate[config.iNa].Aaq_bins_init(index_b)/surrogate[config.iNa].MM/conc_org(index_b)*1000.*surrogate[config.iNa].charge;
  if (config.iMg>=0) inorganion-=surrogate[config.iMg].Aaq_bins_init(index_b)/surrogate[config.iMg].MM/conc_org(index_b)*1000.*surrogate[config.iMg].charge;
  if (config.iCa>=0) inorganion-=surrogate[config.iCa].Aaq_bins_init(index_b)/surrogate[config.iCa].MM/conc_org(index_b)*1000.*surrogate[config.iCa].charge;  
  if (config.iK>=0) inorganion-=surrogate[config.iK].Aaq_bins_init(index_b)/surrogate[config.iK].MM/conc_org(index_b)*1000.*surrogate[config.iK].charge;
    }

  if (abs(inorganion+organion)>=0.)
    {
      derivative=derivative*0.5*(1.0+(organion+inorganion)/pow(pow(organion+inorganion,2)+4*config.Ke/surrogate[config.iHp].gamma_aq_bins(index_b)/surrogate[config.iOHm].gamma_aq_bins(index_b),0.5))-1.0;
      chp_new=0.5*(organion+inorganion+pow(pow(organion+inorganion,2)+4*config.Ke/surrogate[config.iHp].gamma_aq_bins(index_b)/surrogate[config.iOHm].gamma_aq_bins(index_b),0.5));
    }
  else
    {
      chp_new=config.Ke/surrogate[config.iHp].gamma_aq_bins(index_b)/surrogate[config.iOHm].gamma_aq_bins(index_b)/
	(0.5*(-organion-inorganion+pow(pow(organion+inorganion,2)+4*config.Ke/surrogate[config.iHp].gamma_aq_bins(index_b)/surrogate[config.iOHm].gamma_aq_bins(index_b),0.5)));
      derivative=(chp_new/(0.5*(-organion-inorganion+pow(pow(organion+inorganion,2)+4*config.Ke/surrogate[config.iHp].gamma_aq_bins(index_b)/surrogate[config.iOHm].gamma_aq_bins(index_b),0.5))))*0.5*derivative*(1.-(organion+inorganion)/pow(pow(organion+inorganion,2)+4*config.Ke/surrogate[config.iHp].gamma_aq_bins(index_b)/surrogate[config.iOHm].gamma_aq_bins(index_b),0.5))-1.;
    }
  //cout << "in2: " << inorganion << " " << derivative << " " << chp_new <<endl;
  error=chp_new-chp(index_b);
  //exit(0);
}

void error_ph_dyn_ssh(model_config &config, vector<species> &surrogate, int index_b, double Temperature, Array <double, 1> chp, 
		      double &error, double &derivative, Array <double, 1> AQinit, Array<double , 1> &ionic, Array<double, 1> &MMaq,
		      Array<double, 1> &LWC,
		      double &chp_new)
{      
  int n=surrogate.size();
  int i;
  double inorganion=0.0;
  double organion=0.0;
  double total;
  double kelvin_effect=1.0;
  double Kaq,conc_org,K,fion1,fion2;
  if (config.compute_kelvin_effect) //compute the kelvin effect
    {
      kelvin_effect=2.0*config.surface_tension_aq*MMaq(index_b)/
	(8.314*Temperature*config.AQrho(index_b)*
	 0.5*config.diameters(index_b));
      if(kelvin_effect > 50.0)
	kelvin_effect = 50.0;
      kelvin_effect=exp(kelvin_effect);
    }

  compute_conc_org_bins_ssh(config, surrogate, Temperature,  MMaq, AQinit, chp, ionic, LWC, index_b, conc_org);

  derivative=0.0;
  for (i=0;i<n;++i)
    if (surrogate[i].is_inorganic_precursor and surrogate[i].is_solid==false)
      {
        if (i==config.iNH3)
	  {
	    inorganion-=surrogate[config.iNH4p].Aaq_bins_init(index_b)/surrogate[config.iNH4p].MM/conc_org*1000.*surrogate[config.iNH4p].charge;	   
	  }
        else if (i==config.iHNO3)
	  {
	    inorganion-=surrogate[config.iNO3m].Aaq_bins_init(index_b)/surrogate[config.iNO3m].MM/conc_org*1000.*surrogate[config.iNO3m].charge;	    
	  }
	else if (i==config.iCO2)
	  {
	    inorganion-=surrogate[config.iCO3mm].Aaq_bins_init(index_b)/surrogate[config.iCO3mm].MM/conc_org*1000.*surrogate[config.iCO3mm].charge+
	      surrogate[config.iHCO3m].Aaq_bins_init(index_b)/surrogate[config.iHCO3m].MM/conc_org*1000.*surrogate[config.iHCO3m].charge;
	  }
        else if (i==config.iHCl)
	  {
	    inorganion-=surrogate[config.iClm].Aaq_bins_init(index_b)/surrogate[config.iClm].MM/conc_org*1000.*surrogate[config.iClm].charge;	  
	  }
        else if (i==config.iH2SO4)
          {
            total=surrogate[config.iHSO4m].Aaq_bins_init(index_b)/surrogate[config.iHSO4m].MM
              +surrogate[config.iSO4mm].Aaq_bins_init(index_b)/surrogate[config.iSO4mm].MM;

	    K=surrogate[i].Kequilibrium_ssh(Temperature)*surrogate[config.iHSO4m].gamma_aq_bins(index_b)
              /(surrogate[config.iHp].gamma_aq_bins(index_b)*surrogate[config.iSO4mm].gamma_aq_bins(index_b));
            derivative-=1000.*total/conc_org*1.0/pow(1.0+K/chp(index_b),2.0)*K/(chp(index_b)*chp(index_b)); //HSO4m+SO4mm
	    inorganion+=1000.*total/conc_org*(2.0-1.0/(1.0+K/chp(index_b)));
          }
      }
    else if (surrogate[i].is_organic)
      {
	if (surrogate[i].aqt==1)
	  {
	    total=surrogate[i].Aaq_bins_init(index_b)/surrogate[i].MM;
	    fion1=0.0;
	    fion2=0.0;
	    surrogate[i].gamma_LR=surrogate[i].LR(index_b);
	    Kaq=surrogate[i].Kp_eff_aqrealdyn_ssh(config,Temperature, ionic(index_b), chp(index_b),surrogate[config.iHp].LR(index_b),
						  surrogate[config.iHp].SRMR(index_b),MMaq(index_b), fion1, fion2,index_b)/kelvin_effect/surrogate[i].gamma_aq_bins(index_b);
	    organion+=1000.*fion1*total/conc_org;
	    derivative-=fion1*(1-fion1)/chp(index_b)*1000.*total/conc_org;	   
	  }
	else if (surrogate[i].aqt==2)
	  {
	    total=surrogate[i].Aaq_bins_init(index_b)/surrogate[i].MM;
	    fion1=0.0;
	    fion2=0.0;
	    surrogate[i].gamma_LR=surrogate[i].LR(index_b);	    
	    Kaq=surrogate[i].Kp_eff_aqrealdyn_ssh(config, Temperature, ionic(index_b), chp(index_b),surrogate[config.iHp].LR(index_b),
						  surrogate[config.iHp].SRMR(index_b),MMaq(index_b), fion1, fion2, index_b)/kelvin_effect/surrogate[i].gamma_aq_bins(index_b);
	    organion+=1000.*(fion1+2.0*fion2)*total/conc_org;
	    derivative+=(fion1*(fion1-1.)+4.0*fion1*fion2+4.0*fion2*(fion2-1.))/chp(index_b)*1000.*total/conc_org;	   
	  }

      }

  if (config.iNa>=0) inorganion-=surrogate[config.iNa].Aaq_bins_init(index_b)/surrogate[config.iNa].MM/conc_org*1000.*surrogate[config.iNa].charge;
  if (config.iMg>=0) inorganion-=surrogate[config.iMg].Aaq_bins_init(index_b)/surrogate[config.iMg].MM/conc_org*1000.*surrogate[config.iMg].charge;
  if (config.iCa>=0) inorganion-=surrogate[config.iCa].Aaq_bins_init(index_b)/surrogate[config.iCa].MM/conc_org*1000.*surrogate[config.iCa].charge;  
  if (config.iK>=0) inorganion-=surrogate[config.iK].Aaq_bins_init(index_b)/surrogate[config.iK].MM/conc_org*1000.*surrogate[config.iK].charge;

  derivative=derivative*0.5*(1.0+(organion+inorganion)/pow(pow(organion+inorganion,2)+4*config.Ke/surrogate[config.iHp].gamma_aq_bins(index_b)/surrogate[config.iOHm].gamma_aq_bins(index_b),0.5))-1.0;
  chp_new=0.5*(organion+inorganion+pow(pow(organion+inorganion,2)+4*config.Ke/surrogate[config.iHp].gamma_aq_bins(index_b)/surrogate[config.iOHm].gamma_aq_bins(index_b),0.5));
  error=chp_new-chp(index_b);
}

void compute_ph_dyn_ssh(model_config &config, vector<species> &surrogate, double Temperature, Array <double, 1> &chp, 
			Array <double, 1> AQinit, Array<double , 1> &ionic, Array<double, 1> &MMaq, Array<double, 1> &LWC)
{
  int b,i;
  double error_tot=1000.0;
  double factor=1.0;
  double total;
  int nh=1;
  int n=surrogate.size(); 
  Array<double, 1> conc_org;
  conc_org.resize(config.nbins);
  conc_org=LWC;
  Array <double, 1> chp2;
  chp2.resize(config.nbins);

  for (b=0;b<config.nbins;b++)
    {
      compute_conc_org_bins_ssh(config, surrogate, Temperature,  MMaq, AQinit, chp, ionic, LWC, b, conc_org(b));
      /*
      for (i=0;i<n;i++)
        if (surrogate[i].is_organic or i==config.iH2O)
          conc_org(b)+=surrogate[i].Aaq_bins_init(b);
          conc_org(b)=max(conc_org(b),1.e-5*AQinit(b)); //config.MOmin);*/
    }

  while (error_tot>1.0e-3 and nh<=config.nh_max)
    {
      factor=1.0/nh;
      error_tot=0.0;
      for (b=0;b<config.nbins;++b)
        {
          if (config.compute_inorganic)
            {
              double error_h=1000.0;
              double derivative_h;
              int b2;
              int index=0;
              for (b2=0;b2<config.nbins;b2++)
                chp2(b2)=chp(b2);
              double chp_new;
              
	      while(abs(error_h/chp2(b))>1.0e-3 and index<1000)
		{
		  index++;
		  compute_kp_aq_ssh(config, surrogate, Temperature, ionic, chp2, AQinit, MMaq, b, b+1);
		  error_ph_dyn_ssh(config, surrogate, b, Temperature, chp2, error_h, derivative_h,AQinit,ionic,MMaq,LWC,chp_new);
		  if (chp2(b)==1.0e-14 and chp_new<=0.0)
		    error_h=0.0;
		  else if (chp_new <= 0.0 or chp_new!=chp_new)
		    {
		      error_h=1.0e-14-chp2(b);
		      chp2(b)=1.0e-14;
		    }
		  else
		    if (chp2(b)-error_h/derivative_h>0.0 and derivative_h!=0.0) 
		      chp2(b)=chp2(b)-error_h/derivative_h;
		    else
		      {
			if (chp2(b)+error_h<1.0e-14)
			  error_h=1.0e-14-chp2(b);
			chp2(b)=max(chp2(b)+error_h,1.0e-14);
		      }
		}	
	      error_tot=max(error_tot,abs(error_h)/chp2(b));
	      chp(b)=factor*min(chp2(b),100.)+(1.0-factor)*chp(b);	      
            }	  
          surrogate[config.iHp].Aaq_bins_init(b)=chp(b)*conc_org(b)/1000.0;
	  surrogate[config.iOHm].Aaq_bins_init(b)=config.Ke/chp(b)/surrogate[config.iHp].gamma_aq_bins(b)/surrogate[config.iOHm].gamma_aq_bins(b)*conc_org(b)/1000.0;
        }
      nh++;
    }

  for (i=0;i<n;++i)
    if (surrogate[i].is_inorganic_precursor and surrogate[i].is_solid==false)
      {
        if (i==config.iH2SO4)
	  for (b=0;b<config.nbins;++b)
	    {
	      total=surrogate[config.iHSO4m].Aaq_bins_init(b)/surrogate[config.iHSO4m].MM
		+surrogate[config.iSO4mm].Aaq_bins_init(b)/surrogate[config.iSO4mm].MM;
	      
	      double K=surrogate[i].Kequilibrium_ssh(Temperature)*surrogate[config.iHSO4m].gamma_aq_bins(b)
		/(surrogate[config.iHp].gamma_aq_bins(b)*surrogate[config.iSO4mm].gamma_aq_bins(b));
	      
	      surrogate[config.iHSO4m].Aaq_bins_init(b)=total*surrogate[config.iHSO4m].MM*1.0/(1.0+K/chp(b));
	      surrogate[config.iSO4mm].Aaq_bins_init(b)=total*surrogate[config.iSO4mm].MM*K/chp(b)/(1.0+K/chp(b));
	    }

      }
}	

void error_ph_dyn2_ssh(model_config &config, vector<species> &surrogate, int index_b, double Temperature, Array <double, 1> chp, 
		       double &error, double &derivative, Array <double, 1> AQ, Array<double , 1> &ionic, Array<double, 1> &MMaq,
		       Array<double, 1> &LWC, double &chp_new)
{      
  int n=surrogate.size();
  int i;
  double inorganion=0.0;
  double organion=0.0;
  double total;
  double kelvin_effect=1.0;
  double Kaq=0.;

  if (config.compute_kelvin_effect) //compute the kelvin effect
    {
      kelvin_effect=2.0*config.surface_tension_aq*MMaq(index_b)/
	(8.314*Temperature*config.AQrho(index_b)*
	 0.5*config.diameters(index_b));
      if(kelvin_effect > 50.0)
	kelvin_effect = 50.0;
      kelvin_effect=exp(kelvin_effect);
    }

  double conc_org=LWC(index_b);
  for (i=0;i<n;i++)
    if (surrogate[i].is_organic or i==config.iH2O)
      if (surrogate[i].hydrophilic)
        if (surrogate[i].aqt==2 or surrogate[i].aqt==1)
          {
            double fion1,fion2;
            double Kp=surrogate[i].Kp_eff_aqrealdyn_ssh(config, Temperature, ionic(index_b), chp(index_b),surrogate[config.iHp].LR(index_b),
                                                 surrogate[config.iHp].SRMR(index_b),MMaq(index_b), fion1, fion2,index_b);
            conc_org+=surrogate[i].Aaq_bins(index_b)*(1.-fion1-fion2);
          }
        else
          conc_org+=surrogate[i].Aaq_bins(index_b);

  derivative=0.0;
  for (i=0;i<n;++i)
    if (surrogate[i].is_inorganic_precursor and surrogate[i].is_solid==false)
      {
        if (i==config.iNH3)
	  inorganion-=surrogate[config.iNH4p].Aaq_bins(index_b)/surrogate[config.iNH4p].MM/conc_org*1000.0*surrogate[config.iNH4p].charge;         
        else if (i==config.iHNO3)
	  inorganion-=surrogate[config.iNO3m].Aaq_bins(index_b)/surrogate[config.iNO3m].MM/conc_org*1000.0*surrogate[config.iNO3m].charge;
	else if (i==config.iCO2)
	  inorganion-=surrogate[config.iCO3mm].Aaq_bins(index_b)/surrogate[config.iCO3mm].MM/conc_org*1000.0*surrogate[config.iCO3mm].charge+
	    surrogate[config.iHCO3m].Aaq_bins(index_b)/surrogate[config.iHCO3m].MM/conc_org*1000.0*surrogate[config.iHCO3m].charge;
        else if (i==config.iHCl)
	  inorganion-=surrogate[config.iClm].Aaq_bins(index_b)/surrogate[config.iClm].MM/conc_org*1000.0*surrogate[config.iClm].charge;
        else if (i==config.iH2SO4)
          {
            total=surrogate[config.iHSO4m].Aaq_bins(index_b)/surrogate[config.iHSO4m].MM
              +surrogate[config.iSO4mm].Aaq_bins(index_b)/surrogate[config.iSO4mm].MM;               
	    double K=surrogate[i].Kequilibrium_ssh(Temperature)*surrogate[config.iHSO4m].gamma_aq_bins(index_b)
              /(surrogate[config.iHp].gamma_aq_bins(index_b)*surrogate[config.iSO4mm].gamma_aq_bins(index_b));
            derivative-=1000.0*total/conc_org*1.0/pow(1.0+K/chp(index_b),2.0)*K/(chp(index_b)*chp(index_b)); //HSO4m+SO4mm1
	    inorganion+=1000.0*total/conc_org*(2.0-1.0/(1.0+K/chp(index_b)));


          }
      }
    else if (surrogate[i].is_organic)
      {
	if (surrogate[i].aqt==1)
	  {
	    total=surrogate[i].Aaq_bins(index_b)/surrogate[i].MM;
	    double fion1=0.0;
	    double fion2=0.0;
	    surrogate[i].gamma_LR=surrogate[i].LR(index_b);
	    Kaq=surrogate[i].Kp_eff_aqrealdyn_ssh(config, Temperature, ionic(index_b), chp(index_b),surrogate[config.iHp].LR(index_b),
						  surrogate[config.iHp].SRMR(index_b),MMaq(index_b), fion1, fion2, index_b)/kelvin_effect/surrogate[i].gamma_aq_bins(index_b);
	    organion+=1000.*fion1*total/conc_org;
	    derivative-=fion1*(1-fion1)/chp(index_b)*1000.0*total/conc_org;	   
	  }
	else if (surrogate[i].aqt==2)
	  {
	    total=surrogate[i].Aaq_bins(index_b)/surrogate[i].MM;
	    double fion1=0.0;
	    double fion2=0.0;
	    surrogate[i].gamma_LR=surrogate[i].LR(index_b);
	    Kaq=surrogate[i].Kp_eff_aqrealdyn_ssh(config, Temperature, ionic(index_b), chp(index_b),surrogate[config.iHp].LR(index_b),
						  surrogate[config.iHp].SRMR(index_b),MMaq(index_b), fion1, fion2, index_b)/kelvin_effect/surrogate[i].gamma_aq_bins(index_b);
	    organion+=1000.0*(fion1+2.0*fion2)*total/conc_org;
	    derivative+=(fion1*(fion1-1)+4.0*fion1*fion2+4.0*fion2*(fion2-1))/chp(index_b)*1000.*total/conc_org;	   
	  }

      }

  if (config.iNa>=0) inorganion-=surrogate[config.iNa].Aaq_bins(index_b)/surrogate[config.iNa].MM/conc_org*1000.*surrogate[config.iNa].charge;
  if (config.iMg>=0) inorganion-=surrogate[config.iMg].Aaq_bins(index_b)/surrogate[config.iMg].MM/conc_org*1000.*surrogate[config.iMg].charge;
  if (config.iCa>=0) inorganion-=surrogate[config.iCa].Aaq_bins(index_b)/surrogate[config.iCa].MM/conc_org*1000.*surrogate[config.iCa].charge;  
  if (config.iK>=0) inorganion-=surrogate[config.iK].Aaq_bins(index_b)/surrogate[config.iK].MM/conc_org*1000.*surrogate[config.iK].charge;


  derivative=derivative*0.5*(1.0+(organion+inorganion)/pow(pow(organion+inorganion,2)+4*config.Ke/surrogate[config.iHp].gamma_aq_bins(index_b)/surrogate[config.iOHm].gamma_aq_bins(index_b),0.5))-1.0;
  chp_new=0.5*(organion+inorganion+pow(pow(organion+inorganion,2)+4*config.Ke/surrogate[config.iHp].gamma_aq_bins(index_b)/surrogate[config.iOHm].gamma_aq_bins(index_b),0.5));
  error=chp_new-chp(index_b);
}
   
void compute_ph_dyn2_ssh(model_config &config, vector<species> &surrogate, double Temperature, Array <double, 1> &chp, 
			 Array <double, 1> AQ, Array<double , 1> &ionic, Array<double, 1> &MMaq, Array<double, 1> &LWC)
{
  int b,i;
  double error_tot=1000.0;
  double factor=1.0;
  double total;
  int nh=1;
  int n=surrogate.size();
  Array<double, 1> conc_org;
  conc_org.resize(config.nbins);
  conc_org=LWC;
  double error_h=1000.0;
  double derivative_h;
  int b2;
  int index=0;
  double chp_new;
  Array <double, 1> chp2;
  chp2.resize(config.nbins);

  for (b=0;b<config.nbins;b++)
    {
      compute_conc_org_bins_ssh(config, surrogate, Temperature,  MMaq, AQ, chp, ionic, LWC, b, conc_org(b));
      /*
      for (i=0;i<n;i++)
        if (surrogate[i].is_organic or i==config.iH2O)
          conc_org(b)+=surrogate[i].Aaq_bins(b);
          conc_org(b)=max(conc_org(b),1.e-5*AQ(b)); //config.MOmin);*/
    }

  while (error_tot>1.0e-3 and nh<=config.nh_max)
    {      
      factor=1.0/nh;
      error_tot=0.0;      
      for (b=0;b<config.nbins;++b)
        {
          if (config.compute_inorganic)
            {
              error_h=1000.0;
              index=0;
              for (b2=0;b2<config.nbins;b2++)
                chp2(b2)=chp(b2);
              
	      while(abs(error_h/chp2(b))>1.0e-3 and index<1000)
		{
		  index++;
		  compute_kp_aq_ssh(config, surrogate, Temperature, ionic, chp2, AQ, MMaq, b, b+1);
		  error_ph_dyn2_ssh(config, surrogate, b, Temperature, chp2, error_h, derivative_h,AQ,ionic,MMaq,LWC,chp_new);
		  if (chp2(b)==1.0e-14 and chp_new<=0.0)
		    error_h=0.0;
		  else if (chp_new <= 0.0 or chp_new!=chp_new)
		    {
		      error_h=1.0e-14-chp2(b);
		      chp2(b)=1.0e-14;
		    }
		  else
		    if (chp2(b)-error_h/derivative_h>0.0 and derivative_h!=0.0) 
		      chp2(b)=chp2(b)-error_h/derivative_h;
		    else
		      {
			if (chp2(b)+error_h<1.0e-14)
			  error_h=1.0e-14-chp2(b);
			chp2(b)=max(chp2(b)+error_h,1.0e-14);
		      }
		}	
	      error_tot=max(error_tot,abs(error_h)/chp2(b));
	      chp(b)=factor*min(chp2(b),100.)+(1.0-factor)*chp(b);	      
            }	  
          surrogate[config.iHp].Aaq_bins(b)=chp(b)*conc_org(b)/1000.;
	  surrogate[config.iOHm].Aaq_bins(b)=config.Ke/chp(b)/surrogate[config.iOHm].gamma_aq_bins(b)/surrogate[config.iHp].gamma_aq_bins(b)*conc_org(b)/1000.;
        }
      nh++;
    }

  for (i=0;i<n;++i)
    if (surrogate[i].is_inorganic_precursor and surrogate[i].is_solid==false)
      {
        if (i==config.iH2SO4)
	  for (b=0;b<config.nbins;++b)
	    {
	      total=surrogate[config.iHSO4m].Aaq_bins(b)/surrogate[config.iHSO4m].MM
		+surrogate[config.iSO4mm].Aaq_bins(b)/surrogate[config.iSO4mm].MM;
	      
	      double K=surrogate[i].Kequilibrium_ssh(Temperature)*surrogate[config.iHSO4m].gamma_aq_bins(b)
		/(surrogate[config.iHp].gamma_aq_bins(b)*surrogate[config.iSO4mm].gamma_aq_bins(b));
	      
	      surrogate[config.iHSO4m].Aaq_bins(b)=total*surrogate[config.iHSO4m].MM*1.0/(1.0+K/chp(b));
	      surrogate[config.iSO4mm].Aaq_bins(b)=total*surrogate[config.iSO4mm].MM*K/chp(b)/(1.0+K/chp(b));
	    }
      }

}


void activity_coefficients_dyn_aq_bins_ssh(model_config &config, vector<species>& surrogate,
                                           double &Temperature, Array<double, 1> &AQinit,
                                           Array<double, 3> &MOinit,
                                           Array<double, 1> &conc_inorganic,
                                           Array<double, 1> &ionic, Array<double, 1> &ionic_organic,
                                           Array<double, 1> &organion, Array<double, 1> &chp,
                                           Array<double, 1> &LWC, Array<double, 1> &MMaq, //Array<double, 1> &conc_org,
                                           double factor, double deltat, int index2, int b)
{
  //compute the activity coefficients with UNIFAC (short range interactions) for the organic phase
  //MOW: mean molar mass of the organic phase
  int i;
  int n=surrogate.size();
  double XH2O;
  Array<double, 1> X_unifac,gamma_unifac;
  double AQ=0.0;
  double error_h=1000.0;
  double derivative_h;
  Array <double, 1> chp2;
  chp2.resize(config.nbins);
  int b2;
  int index=0;
  double chp_new;
  double conc_org;
 
  conc_inorganic(b)=0.0;
  for (i=0;i<n;i++)
    if (surrogate[i].is_inorganic_precursor==false and surrogate[i].is_organic==false and i!=config.iH2O)
      conc_inorganic(b)+=surrogate[i].Aaq_bins_init(b);
 
  AQ=conc_inorganic(b);
	  
  for (i=0;i<n;++i)
    if (surrogate[i].hydrophilic)	      
      surrogate[i].Aaq=surrogate[i].Aaq_bins_init(b);	       	      
  surrogate[config.iH2O].Aaq+=LWC(b);

  for (i=0;i<n;++i)
    if (surrogate[i].hydrophilic)
      if(surrogate[i].is_organic or i==config.iH2O)
        AQ+=surrogate[i].Aaq;

  //compute_conc_org_bins_ssh(config, surrogate, Temperature,  MMaq, AQinit, chp, ionic, LWC, b, conc_org);

  compute_organion_bins_ssh(config, surrogate, Temperature,  MMaq, AQinit, ionic,
                            chp, ionic_organic, organion, LWC, conc_org, b);
  /*
  conc_org=LWC(b);	  
  for (i=0;i<n;i++)
    if (surrogate[i].is_organic or i==config.iH2O)
      conc_org+=surrogate[i].Aaq_bins_init(b);
      conc_org=max(conc_org,1.e-5*AQinit(b)); //config.MOmin);*/
	  	
  //config.rho_aqueous=config.AQrho(b);
  compute_ionic_strenght2_ssh(config, surrogate, Temperature, AQ, conc_inorganic(b), ionic(b), chp(b),
                              organion(b), ionic_organic(b), conc_org, factor, config.Ke/surrogate[config.iHp].gamma_aq_bins(b)/surrogate[config.iOHm].gamma_aq_bins(b));                    
          
  activity_coefficients_aq_ssh(config,surrogate,Temperature,0.0,MMaq(b),XH2O, conc_org);

  /*
    for (i=0;i<config.nion_aiomfac;i++)            
    config.gamma_MR_ions(b,i)=config.gamma_MR_ions_bins(b,i);*/
	  
  if (config.compute_long_and_medium_range_interactions)
    activity_coefficients_LR_MR_ssh(config, surrogate, Temperature, 0.0, ionic(b));

  /*
    for (i=0;i<config.nion_aiomfac;i++)            
    config.gamma_MR_ions_bins(b,i)=config.gamma_MR_ions(i);*/
	  
  for (i=0;i<n;++i)
    {
      surrogate[i].gamma_aq_bins(b)=surrogate[i].gamma_aq;              
      if (config.compute_long_and_medium_range_interactions)
        {
          surrogate[i].LR(b)=surrogate[i].gamma_LR;         
          surrogate[i].SRMR(b)=surrogate[i].gamma_SRMR;                  
        }
      else
        {
          surrogate[i].LR(b)=pow(10,-0.511*pow(298.0/Temperature,1.5)*
                                 pow(ionic(b),0.5)/(1.0+pow(ionic(b),0.5)));
          surrogate[i].SRMR(b)=1.0;
        }
    }        

  if (factor>0.)
    {	  
      //chp2(b)=max(chp(b),1.e-14);
      chp2=chp;     
      if (config.compute_inorganic)
        {
          error_h=1000.0;
          index=0;
          double chp_save=-1;
          double chp_save2=-1; 
          double factor2=1.;
          double error_h_save=10000;
          while(abs((chp_save-chp2(b))/chp2(b))/factor2>1.0e-3 and index<5 and abs(chp2(b)-chp_save)>0)
            {
              index++;
              chp_save2=chp_save;
              chp_save=chp2(b);	    
              compute_kp_aq_ssh(config, surrogate, Temperature, ionic, chp2, AQinit,MMaq, b, b+1);	      
              if (config.imethod>0)
                prodloss_aq_bins_ssh(config, surrogate, AQinit, LWC, MMaq, chp2, ionic, MOinit, conc_org, index2, Temperature, deltat, b);	      
              error_h_save=error_h;	   
              error_ph_bins_ssh(config, surrogate, b, Temperature, chp2, error_h, derivative_h,AQinit,ionic,MMaq,LWC,chp_new,organion(b),deltat);	     
              //cout << index << " " << chp2(b) << " " << chp_new << ' ' << derivative_h << " " << chp2(b)-error_h/derivative_h << endl;

              if (chp2(b)==1.0e-14 and chp_new<=0.0)
                error_h=0.0;
              else if (chp_new <= 0.0 or chp_new!=chp_new)
                {
                  error_h=1.0e-14-chp2(b);
                  chp2(b)=1.0e-14;
                }
              else
                {
                  //chp2(b)=chp_new;
                        
                  if (chp2(b)-error_h/derivative_h>0.0 and derivative_h!=0.0)
                    chp2(b)=chp2(b)-error_h/derivative_h;
                  else
                    {
                      /*
                        if (chp2(b)+error_h<1.0e-14)
                        error_h=1.0e-14-chp2(b);*/
                      chp2(b)=chp2(b)+error_h;
                    }
                       
                }

	      //chp2(b)=chp_new;
              /*
                if (abs(error_h)>abs(error_h_save))
                factor2=max(factor2/2,0.1);*/
              if (abs(chp2(b)-chp_save2)/max(chp2(b),1.e-14)<1.0e-4 and abs(chp2(b)-chp_save)/max(chp2(b),1.0e-14)>1.0e-4)
                factor2=max(factor2/2,0.1);

	      chp2(b)=min(max(chp2(b),0.5*chp_save),2*chp_save);
                    	     
              chp2(b)=min(max(factor2*chp2(b)+(1-factor2)*chp_save,1.e-14),10000.);
              //cout << index << " " << chp2(b) << " " << factor2 << " " << error_h << " " << chp_new << " " << derivative_h << endl;
            }
	  //cout << chp2(b) << " " << surrogate[config.iOHm].gamma_aq_bins(b) << " " << surrogate[config.iHp].gamma_aq_bins(b) << endl;
	  //cout << b << " " << chp2(b) << endl;

	  //if (index==50)
	  //exit(0);
          if (abs((chp_save-chp2(b))/chp2(b))/factor2>1.0e-3)
            {
              //cout << chp2(b) << " " << error_h << endl;
              //cout << "bug ph" << endl;
              //exit(0);
            }
          //cout << "chp_new: " << chp2(b) << " " << chp(b) << " " << b << endl;
          double varph=10.;
          chp2(b)=max(min(chp2(b),varph*chp(b)),chp(b)/varph);
                
          chp(b)=factor*min(chp2(b),100.)+(1.0-factor)*chp(b);
          chp(b)=max(chp(b),1.0e-14);

          //compute_conc_org_bins_ssh(config, surrogate, Temperature,  MMaq, AQinit, chp, ionic, LWC, b, conc_org);
          /*
          conc_org=LWC(b);	  
          for (i=0;i<n;i++)
            if (surrogate[i].is_organic or i==config.iH2O)
              conc_org+=surrogate[i].Aaq_bins_init(b);
              conc_org=max(conc_org,1.e-5*AQinit(b)); //config.MOmin);*/
	      
          // Concentrations of H+ is not calculated to avoid numerical issues when pH given by ISORROPIA is too low
          surrogate[config.iHp].Aaq_bins_init(b)=chp(b)*conc_org/1000.0;
	  surrogate[config.iOHm].Aaq_bins_init(b)=config.Ke/chp(b)/surrogate[config.iHp].gamma_aq_bins(b)/surrogate[config.iOHm].gamma_aq_bins(b)*conc_org/1000.0;
        }	  
    }
}









void activity_coefficients_dyn_aq_ssh(model_config &config, vector<species>& surrogate,
				      double &Temperature, Array<double, 1> &AQinit,
				      Array<double, 3> &MOinit,
				      Array<double, 1> &conc_inorganic,
				      Array<double, 1> &ionic, Array<double, 1> &ionic_organic,
				      Array<double, 1> &organion, Array<double, 1> &chp,
				      Array<double, 1> &LWC, Array<double, 1> &MMaq, double factor, double deltat, int index2)
{
  //compute the activity coefficients with UNIFAC (short range interactions) for the organic phase
  //MOW: mean molar mass of the organic phase
  int b;
 
  if (config.use_global_dynamic_parameters)
    {
      int i;
      int n=surrogate.size();
      double XH2O;
      Array<double, 1> X_unifac,gamma_unifac;
      double AQ=0.0;
      double error_h=1000.0;
      double derivative_h;
      Array <double, 1> chp2;
      chp2.resize(config.nbins);
      int b2;
      int index=0;
      double chp_new;
      double conc_org;
      double MMaqtemp=0.0;
      double sum_inorganic=0.0;
      double ionic_tmp,chp_tmp,organion_tmp,ionic_organic_tmp;
      organion_tmp=0.0;
      for (b=0;b<config.nbins;++b)
        organion_tmp+=organion(b)/config.nbins;

      ionic_organic_tmp=0.0;
      for (b=0;b<config.nbins;++b)
        ionic_organic_tmp+=ionic_organic(b)/config.nbins;

      for (b=0;b<config.nbins;++b)
        sum_inorganic+=conc_inorganic(b);

      chp_tmp=0.0;
      for (b=0;b<config.nbins;++b)
        chp_tmp+=chp(b)/config.nbins;

      AQ+=sum_inorganic;
	  
      for (i=0;i<n;++i)
        if (surrogate[i].hydrophilic)
          {
            surrogate[i].Aaq=0.0;
            for (b=0;b<config.nbins;++b)
              surrogate[i].Aaq+=surrogate[i].Aaq_bins_init(b);
            if(i==config.iH2O)
              for (b=0;b<config.nbins;++b)
                surrogate[i].Aaq+=LWC(b);
          }
	  
      for (i=0;i<n;++i)
        if (surrogate[i].hydrophilic)
          if(surrogate[i].is_organic or i==config.iH2O)
            AQ+=surrogate[i].Aaq;

      conc_org=0.;	 
      for (b=0;b<config.nbins;++b)
	compute_conc_org_bins_ssh(config, surrogate, Temperature,  MMaq, AQinit, chp, ionic, LWC, b, conc_org);
        
      conc_org=max(conc_org,1.0e-5*sum(AQinit)); //config.MOmin);

      //config.rho_aqueous=config.AQrho(b);
      compute_ionic_strenght_ssh(config, surrogate, Temperature, AQ, sum_inorganic, ionic_tmp, chp_tmp,
				 organion_tmp, ionic_organic_tmp, conc_org, factor);      

      activity_coefficients_aq_ssh(config,surrogate,Temperature,0.0,MMaqtemp,XH2O,conc_org);
      
      if (config.compute_long_and_medium_range_interactions)
        activity_coefficients_LR_MR_ssh(config, surrogate, Temperature, 0.0, ionic_tmp);
						  
      for (i=0;i<n;++i)
        for (b=0;b<config.nbins;++b)
          {
            surrogate[i].gamma_aq_bins(b)=surrogate[i].gamma_aq;
            if (config.compute_long_and_medium_range_interactions)
              {
                surrogate[i].LR(b)=surrogate[i].gamma_LR;                
                surrogate[i].SRMR(b)=surrogate[i].gamma_SRMR;
              }
            else
              {
                surrogate[i].LR(b)=1.0;
                surrogate[i].SRMR(b)=1.0;
              } 
          }
	  
      for (b=0;b<config.nbins;++b)
        {
          MMaq(b)=MMaqtemp;
          ionic(b)=ionic_tmp;
          chp(b)=chp_tmp;
          organion(b)=organion_tmp;
          ionic_organic(b)=ionic_organic_tmp;
	  i=config.iHp;
	  if (AQ>0.0)
	    surrogate[i].Aaq_bins_init(b)=AQinit(b)/AQ*surrogate[i].Aaq;
	  else
	    surrogate[i].Aaq_bins_init(b)=0.0;
        }
      
      compute_organion_ssh(config, surrogate, Temperature,  MMaq, AQinit, ionic,
                           chp, ionic_organic, organion, LWC);
    }
  else
    {
      for (b=0;b<config.nbins;++b)
        activity_coefficients_dyn_aq_bins_ssh(config, surrogate, Temperature, AQinit,
                                              MOinit,
                                              conc_inorganic,
                                              ionic, ionic_organic, organion,
                                              chp, LWC, MMaq, factor, deltat, index2, b);
      
      //compute_kp_aq_ssh(config, surrogate, Temperature, ionic, chp, MMaq, 0, config.nbins);
      if (config.imethod==0.)
        characteristic_time_aq_ssh(config, surrogate, Temperature, chp, LWC, AQinit, MOinit, MMaq, ionic);
    } 

}


void twostep_tot_ssh(model_config &config, vector<species>& surrogate, double &tequilibrium,
                     Array <double, 1> &AQinit,  Array <double, 1> &AQ,
                     Array <double, 1> &conc_inorganic, Array <double, 1> &ionic,
                     Array <double, 1> &ionic_organic, Array <double, 1> &organion,
                     Array <double, 1> &chp, Array <double, 1> &LWC,
                     Array <double, 3> &MOinit,  Array <double, 3> &MO, Array <double, 3> &MOW,
                     double &Temperature, double &RH, Array <double, 1> &MMaq,
                     bool compute_activity_coefficients, double factor, double t, double deltat,
		     int index)
{
  //This routine computes the aqueous-phase and organic-phase concentrations being at equilibrium
  //FOR THE COUPLED SYSTEM
  //factor: weight factor to decrease variations in case of strong variations of composition
  //                between two iterations
  //compute_activity_coefficients: should the activity coefficients be recomputed?
  //tequilibrium: criterium to distinguish cases under equilibrium (characteristic time lower
  //                  than tequilibrium)
  int n=surrogate.size();
  int i,b,ilayer,iphase;
  //double Kaq;
  Array <double, 1> Kp,chp0,chp_save;
  Kp.resize(config.max_number_of_phases);
  chp0.resize(config.nbins);
  chp_save.resize(config.nbins);
  
  double conc_equilibrium;
  double sum1=0.0;
  bool is_equilibrium=false;
  double LWCtot=0.0;
  int index_b; //,index_layer;
  double Mads,conceq;
  int iter,index2;
  double error_chp,conc_org,total,total2;
  double tiny=1.0e-26;
  
  for (b=0;b<config.nbins;++b)
    LWCtot+=LWC(b);

  /*
  if (index==0)
  prodloss_aq_ssh(config, surrogate, AQinit, LWC, MOinit, 0, deltat);*/
    
  
  //cout << chp << endl;
  chp0=chp;
  chp_save=chp;

  for (b=0;b<config.nbins;b++)
    if (chp(b)<1.0e-15)
      {
        cout << "la " << chp << endl;
        throw string("stop");
      }
  
  iter=0;
  error_chp=1000.;
  /*
    while (iter<50 and error_chp>1.0e-3)
    {
    error_chp=0.0;
    iter++;
    if (config.compute_inorganic)	
    {
    compute_ph_dyn_ssh(config, surrogate, Temperature, chp, AQinit, ionic, MMaq, LWC);
    for (b=0;b<config.nbins;++b)
    {	  
    conc_org=LWC(b);
    for (i=0;i<n;++i)
    if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophilic)
    conc_org+=surrogate[i].Aaq_bins_init(b);
    conc_org=max(conc_org,1.e-5*AQinit(b)); //config.MOmin);
              
              
    //config.rho_aqueous=config.AQrho(b);
    compute_ionic_strenght2_ssh(config, surrogate, Temperature, AQinit(b), conc_inorganic(b), ionic(b), chp(b),
    organion(b), ionic_organic(b), conc_org, 1.0);
    //cout << chp(b) << " " << chp0(b) << endl;
    error_chp=max(error_chp,abs(chp0(b)-chp(b))/chp0(b));
    chp0(b)=chp(b);	      
    }      
    }
	  
    if (compute_activity_coefficients)
    {
    if (config.compute_organic)
    {
    activity_coefficients_dyn_org_ssh(config, surrogate, Temperature, MOW);
    compute_kp_org_ssh(config, surrogate, MOinit, Temperature, MOW);
    }
	  
    if (LWCtot>config.LWClimit)
    {
    activity_coefficients_dyn_aq_ssh(config, surrogate, Temperature,AQinit,
    MOinit,conc_inorganic, ionic, ionic_organic,
    organion,chp,LWC,MMaq, 1.0);
    compute_kp_aq_ssh(config, surrogate, Temperature, ionic, chp, MMaq);   
    }
    }
    }

    for (b=0;b<config.nbins;b++)
    {
    chp(b)=factor*max(min(chp(b),2.*chp_save(b)),0.5*chp_save(b))+(1.-factor)*chp_save(b);
    compute_ionic_strenght2_ssh(config, surrogate, Temperature, AQinit(b), conc_inorganic(b), ionic(b), chp(b),
    organion(b), ionic_organic(b), conc_org, 1.0);
    }
    if (compute_activity_coefficients)
    {
    if (config.compute_organic)
    {
    activity_coefficients_dyn_org_ssh(config, surrogate, Temperature, MOW);
    compute_kp_org_ssh(config, surrogate, MOinit, Temperature, MOW);
    }
	  
    if (LWCtot>config.LWClimit)
    {
    activity_coefficients_dyn_aq_ssh(config, surrogate, Temperature,AQinit,
    MOinit,conc_inorganic, ionic, ionic_organic,
    organion,chp,LWC,MMaq, 1.0);
   compute_kp_aq_ssh(config, surrogate, Temperature, ionic, chp, MMaq);   
    }
    }
  */

  
  for (b=0;b<config.nbins;b++)
    {
      compute_conc_org_bins_ssh(config, surrogate, Temperature,  MMaq, AQinit, chp, ionic, LWC, b, conc_org);
      organion(b)=0.0;      
      for (i=0;i<n;++i)
        if (surrogate[i].hydrophilic and surrogate[i].is_organic)
          if(surrogate[i].nonvolatile==false and (surrogate[i].aqt==1 or surrogate[i].aqt==2))
            {
              double fion1=0.;
              double fion2=0.;
              double Kploc=surrogate[i].Kp_eff_aqrealdyn_ssh(config, Temperature, ionic(b), chp(b),surrogate[config.iHp].LR(b),
                                                             surrogate[config.iHp].SRMR(b),MMaq(b), fion1, fion2,b);
              //surrogate[i].Kp_eff_aq_ssh(config, Temperature, ionic(b), chp(b),surrogate[config.iHp].LR(b),
              //                                     surrogate[config.iHp].SRMR(b),MMaq(b), fion1, fion2);
              //molality1: molality of ions HA- or A-  
              double molality1=surrogate[i].Aaq_bins_init(b)*fion1/surrogate[i].MM/conc_org*1000.0;
              //molality2: molality of ions A2-
              double molality2=surrogate[i].Aaq_bins_init(b)*fion2/surrogate[i].MM/conc_org*1000.0;
              //compute ionic_organic and organion
              organion(b)+=molality1+2*molality2;
            }
          
      if (config.compute_inorganic==false and config.isorropia_ph==false)
        { 
          double inorganion=0.0;
          for (i=0;i<n;++i)
            if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)
              {
                surrogate[i].molality=surrogate[i].Aaq_bins_init(b)/surrogate[i].MM/conc_org*1000.0;           
                if (i!=config.iHp and i!=config.iOHm)
                  inorganion-=surrogate[i].molality*surrogate[i].charge; 
              }


          double chp_new=0.5*(organion(b)+inorganion+pow(pow(organion(b)+inorganion,2)+4*config.Ke/surrogate[config.iHp].gamma_aq_bins(b)/surrogate[config.iOHm].gamma_aq_bins(b),0.5));
          //cout << "new " << chp_new << endl;
          double change_ph=1.1;
          chp_new=max(max(min(chp_new,change_ph*chp(b)),1./change_ph*chp(b)),1.e-14);
          chp_new=min(chp_new,100.);
          //chp(b)=chp_new;
          chp(b)=factor*chp_new+(1.0-factor)*chp(b);
          //chp(b)=pow(chp_new,factor)*pow(chp(b),1-factor);
          if (chp(b)==0.0)
            chp(b)=pow(10.0,-5.6);
          //cout << b << " " << chp_new << " " << chp(b) << " " << AQinit(b) << " " << organion(b) << " " << inorganion << endl;
                
        }

      compute_ionic_strenght2_ssh(config, surrogate, Temperature, AQinit(b), conc_inorganic(b), ionic(b), chp(b),
                                  organion(b), ionic_organic(b), conc_org, 1.0, config.Ke/surrogate[config.iHp].gamma_aq_bins(b)/surrogate[config.iOHm].gamma_aq_bins(b));     
    }
    
  for (b=0;b<config.nbins;b++)
    if (chp(b)<1.0e-15)
      {
        cout << "la2 " << chp << endl;
        throw string("stop");
      }

  if (compute_activity_coefficients)
    {
      if (config.compute_organic)
        {
          activity_coefficients_dyn_org_ssh(config, surrogate, Temperature, MOW);
          compute_kp_org_ssh(config, surrogate, MOinit, Temperature, MOW);
	  characteristic_time_ssh(config, surrogate, MOinit, AQinit, LWCtot); 
        }
     
      if (LWCtot>config.LWClimit)
        {
          activity_coefficients_dyn_aq_ssh(config, surrogate, Temperature,AQinit,
					   MOinit, conc_inorganic, ionic, ionic_organic,
					   organion,chp,LWC,MMaq, factor, deltat, 1);
          compute_kp_aq_ssh(config, surrogate, Temperature, ionic, chp, AQinit, MMaq, 0, config.nbins);
	  characteristic_time_aq_ssh(config, surrogate, Temperature, chp, LWC, AQinit, MOinit, MMaq, ionic);  
        }    
    } 

  for (b=0;b<config.nbins;b++)
    if (chp(b)<1.0e-15)
      {
        cout << "la3 " << chp << endl;
        throw string("stop");
      }     

  for (i=0;i<n;i++)
    {
      surrogate[i].kprod_gas=0.0;
      surrogate[i].kloss_gas=0.0;
    }

  //compute kinetic rates
  if (LWCtot>config.LWClimit)
    prodloss_aq_ssh(config, surrogate, AQinit, LWC, MMaq, chp, ionic, MOinit, 1, Temperature, deltat);
  
  if (config.compute_organic)
    {
      if (index==0)
	prodloss_org_ssh(config, surrogate, MOinit, AQinit, tiny, 0, deltat);
      else
	prodloss_org_ssh(config, surrogate, MOinit, AQinit, tiny, 1, deltat);
    }

  if (config.compute_organic)
    hydratation_dyn_ssh(config, surrogate, RH, AQinit);

  if (config.chemistry)
    prodloss_chem_ssh(config, surrogate, MOinit, MOW, AQinit, LWC, MMaq, chp, ionic, tiny, Temperature, 0, RH);   

  double apnew;
  Array <int, 1> ifound;
  ifound.resize(n);
  ifound=0;
  
  for (i=0;i<n;++i)
    {  
        
      if(surrogate[i].is_organic and config.compute_organic)
        if(surrogate[i].hydrophobic)
          for (b=0;b<config.nbins;++b)
            for (ilayer=0;ilayer<config.nlayer;++ilayer)
              for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                {		    	  		    
                  apnew=(surrogate[i].Ap_layer_init0(b,ilayer,iphase)
                         +deltat*surrogate[i].kprod(b,ilayer,iphase))/(1.0+deltat*surrogate[i].kloss(b,ilayer,iphase));
                  surrogate[i].Ap_layer_init(b,ilayer,iphase)=apnew; //factor*apnew+(1.0-factor)*surrogate[i].Ap_layer_init(b,ilayer,iphase);
                }

      if (LWCtot>config.LWClimit)
	if((surrogate[i].is_organic and config.compute_organic) or ((i==config.iClm or i==config.iNO3m or i==config.iNH4p))) // and config.compute_inorganic))
	      if (surrogate[i].hydrophilic)
		for (b=0;b<config.nbins;++b)
		  {
		    apnew=(surrogate[i].Aaq_bins_init0(b)+deltat*surrogate[i].kprod_aq(b))
		      /(1.0+deltat*surrogate[i].kloss_aq(b));
		    surrogate[i].Aaq_bins_init(b)=apnew; //factor*apnew+(1.0-factor)*surrogate[i].Aaq_bins_init(b);
		  }

      if (LWCtot>config.LWClimit)
        if (i==config.iSO4mm) // and config.compute_inorganic)
          for (b=0;b<config.nbins;++b)
            {
	      double Keq;
	      if (config.isorropia_ph and surrogate[config.iHSO4m].Aaq_bins_init0(b)>0.)
		Keq=surrogate[config.iSO4mm].Aaq_bins_init0(b)/surrogate[config.iSO4mm].MM/surrogate[config.iHSO4m].Aaq_bins_init0(b)*surrogate[config.iHSO4m].MM;
	      else if (config.isorropia_ph)
		Keq=1.e15;
	      else
		Keq=surrogate[config.iH2SO4].keqi/chp(b)*surrogate[config.iHSO4m].gamma_aq_bins(b)/
		  (surrogate[config.iHp].gamma_aq_bins(b)*surrogate[config.iSO4mm].gamma_aq_bins(b));
	     
              total2=surrogate[config.iHSO4m].Aaq_bins_init(b)/surrogate[config.iHSO4m].MM+surrogate[config.iSO4mm].Aaq_bins_init(b)/surrogate[config.iSO4mm].MM;
              total=(surrogate[config.iHSO4m].Aaq_bins_init0(b)+deltat*surrogate[config.iHSO4m].kprod_aq(b))/
                (1.0+deltat*surrogate[config.iHSO4m].kloss_aq(b))/surrogate[config.iHSO4m].MM;
              total+=(surrogate[config.iSO4mm].Aaq_bins_init0(b)+deltat*surrogate[config.iSO4mm].kprod_aq(b))/
                (1.0+deltat*surrogate[config.iSO4mm].kloss_aq(b))/surrogate[config.iSO4mm].MM;

	      surrogate[config.iHSO4m].Aaq_bins_init(b)=total*surrogate[config.iHSO4m].MM*1.0/(1.0+Keq); //(1.0-factor)*surrogate[config.iHSO4m].Aaq_bins_init(b)+factor*total*surrogate[config.iHSO4m].MM*1.0/(1.0+Keq);
              surrogate[config.iSO4mm].Aaq_bins_init(b)=total*surrogate[config.iSO4mm].MM*Keq/(1.0+Keq); //(1.0-factor)*surrogate[config.iSO4mm].Aaq_bins_init(b)+factor*total*surrogate[config.iSO4mm].MM*Keq/(1.0+Keq);
            }

      
      if (LWCtot>config.LWClimit)
	if (config.compute_inorganic)
	  if (i==config.iCO3mm or i==config.iHCO3m)
	    for (b=0;b<config.nbins;++b)
	      {

		double frac1=1.0/(1.0+(surrogate[config.iCO2].keq2/(chp(b)*surrogate[config.iHp].gamma_aq_bins(b)*surrogate[config.iCO3mm].gamma_aq_bins(b)/surrogate[config.iHCO3m].gamma_aq_bins(b))));
		//cout << "frac1: " << b << " " << frac1 << " " << surrogate[config.iHCO3m].Aaq_bins_init(b) << " " << surrogate[config.iCO3mm].Aaq_bins_init(b) << endl;
		//cout << "g: " << surrogate[config.iHp].gamma_aq_bins(b) << " " << surrogate[config.iCO3mm].gamma_aq_bins(b) << " " << surrogate[config.iHCO3m].gamma_aq_bins(b) << endl;
		//cout << surrogate[config.iCO2].keq2 << " " << surrogate[config.iCO2].keq2/chp(b) << " " << chp(b) << endl;
		double Kaq=surrogate[config.iCO2].Kaq(b);
		//cout << Kaq << " " << Kaq*surrogate[config.iH2O].Aaq*surrogate[config.iCO2].Ag << " " << AQinit(b) << " " << surrogate[config.iH2O].Aaq_bins_init(b) << " " << surrogate[config.iHCO3m].Aaq_bins_init(b) <<  endl;

		double total2=0.;
		double total0=0.;
		
		if (config.solids)
		  for (int j=0;j<n;j++)
		    if (surrogate[j].is_solid)
		      {			
			for (int j2=0;j2<surrogate[j].nion;j2++)
			  {
			    //cout << j2 << " " << surrogate[j].iion1 << " " << surrogate[j].iion2 << " " << config.iHCO3m << " " << condifendl;
			    if (j2==0 and (surrogate[j].iion1==config.iHCO3m or surrogate[j].iion1==config.iCO3mm))
			      {
				total0+=(surrogate[j].Asol_bins_init0(b))/surrogate[j].MM*surrogate[j].pion1;
				total2+=(surrogate[j].Asol_bins_init(b))/surrogate[j].MM*surrogate[j].pion1;
			      }
			    else if (j2==1 and (surrogate[j].iion2==config.iHCO3m or surrogate[j].iion2==config.iCO3mm))
			      {			      
				total0+=(surrogate[j].Asol_bins_init0(b))/surrogate[j].MM*surrogate[j].pion2;
				total2+=(surrogate[j].Asol_bins_init(b))/surrogate[j].MM*surrogate[j].pion2;
			      }
			    else if (j2==2 and (surrogate[j].iion3==config.iHCO3m or surrogate[j].iion3==config.iCO3mm))
			      {
				total0+=(surrogate[j].Asol_bins_init0(b))/surrogate[j].MM*surrogate[j].pion3;
				total2+=(surrogate[j].Asol_bins_init(b))/surrogate[j].MM*surrogate[j].pion3;
			      }
			  }
		      }
		
		double total1=surrogate[config.iHCO3m].Aaq_bins_init(b)/surrogate[config.iHCO3m].MM+surrogate[config.iCO3mm].Aaq_bins_init(b)/surrogate[config.iHCO3m].MM;
		double frac2=(total1)/(total2+total1);
		total=(frac2*(surrogate[config.iHCO3m].Aaq_bins_init0(b)+frac1*total0*surrogate[config.iHCO3m].MM)+deltat*surrogate[config.iHCO3m].kprod_aq(b))
                  /(1.0+deltat*surrogate[config.iHCO3m].kloss_aq(b))/surrogate[config.iHCO3m].MM
		  +(frac2*(surrogate[config.iCO3mm].Aaq_bins_init0(b)+(1-frac1)*total0*surrogate[config.iCO3mm].MM)+deltat*surrogate[config.iCO3mm].kprod_aq(b))
                  /(1.0+deltat*surrogate[config.iCO3mm].kloss_aq(b))/surrogate[config.iCO3mm].MM;

		/*
		if (total>2*total1 and total1>0)
		  total=2*total1;
		if (total<0.5*total1 and total1>0)
		total=0.5*total1;*/
		
		/*
		double total2=0.;
		if (config.solids)
		  for (int j=0;j<n;j++)
		    if (surrogate[j].is_solid)
		      {			
			for (int j2=0;j2<surrogate[j].nion;j2++)
			  {
			    //cout << j2 << " " << surrogate[j].iion1 << " " << surrogate[j].iion2 << " " << config.iHCO3m << " " << condifendl;
			    if (j2==0 and (surrogate[j].iion1==config.iHCO3m or surrogate[j].iion1==config.iCO3mm))
			      total2+=(surrogate[j].Asol_bins_init0(b)-surrogate[j].Asol_bins_init(b))/surrogate[j].MM*surrogate[j].pion1;
			    else if (j2==1 and (surrogate[j].iion2==config.iHCO3m or surrogate[j].iion2==config.iCO3mm))
			      {
				//cout << "ok in: " << surrogate[j].Asol_bins_init0(b) << " " << surrogate[j].Asol_bins_init(b) << endl;
				total2+=(surrogate[j].Asol_bins_init0(b)-surrogate[j].Asol_bins_init(b))/surrogate[j].MM*surrogate[j].pion2;
			      }
			    else if (j2==2 and (surrogate[j].iion3==config.iHCO3m or surrogate[j].iion3==config.iCO3mm))
			      total2+=(surrogate[j].Asol_bins_init0(b)-surrogate[j].Asol_bins_init(b))/surrogate[j].MM*surrogate[j].pion3;
			  }
		      }

		total=(surrogate[config.iHCO3m].Aaq_bins_init0(b)+deltat*surrogate[config.iHCO3m].kprod_aq(b))
                  /(1.0+deltat*surrogate[config.iHCO3m].kloss_aq(b))/surrogate[config.iHCO3m].MM
		  +(total2*surrogate[config.iCO3mm].MM+surrogate[config.iCO3mm].Aaq_bins_init0(b)+deltat*surrogate[config.iCO3mm].kprod_aq(b))
                  /(1.0+deltat*surrogate[config.iCO3mm].kloss_aq(b))/surrogate[config.iCO3mm].MM;
		
		//total+=total2;
		cout << "rtota: " << total << " " << surrogate[config.iHCO3m].Aaq_bins_init(b)/61+surrogate[config.iCO3mm].Aaq_bins_init(b)/60 << " " << total2 << " " << total+total2 << endl;

		if (total<0)
		  {
		    cout<< surrogate[config.iHCO3m].Aaq_bins_init0(b)/61+surrogate[config.iCO3mm].Aaq_bins_init0(b)/60 << endl;
		    for (int j=0;j<n;j++)
		      if (surrogate[j].is_solid)
			cout << surrogate[j].Asol_bins_init0(b)/surrogate[j].MM << endl;
		    
		    exit(0);
		  }*/
		
		surrogate[config.iHCO3m].Aaq_bins_init(b)=total*frac1*surrogate[config.iHCO3m].MM;
		surrogate[config.iCO3mm].Aaq_bins_init(b)=total*(1.-frac1)*surrogate[config.iCO3mm].MM;		
		//cout << "toal: " << b << " " << frac1 << "  " << total << " " << surrogate[config.iHCO3m].kprod_aq(b)+surrogate[config.iCO3mm].kprod_aq(b) << " " << surrogate[config.iHCO3m].kloss_aq(b) << endl;
	      }      
	
      if ((surrogate[i].is_organic and config.compute_organic) or (surrogate[i].is_inorganic_precursor and config.compute_inorganic))
        {
          apnew=(surrogate[i].Ag0+deltat*surrogate[i].kprod_gas)/
            (1.0+deltat*surrogate[i].kloss_gas);
          //if (surrogate[i].Ag>tiny)
          //  apnew=max(min(apnew,10.*surrogate[i].Ag),0.1*surrogate[i].Ag);
          surrogate[i].Ag=apnew; //factor*apnew+(1.-factor)*surrogate[i].Ag;
        }

      /*
      if (surrogate[i].name=="GLY" or surrogate[i].name=="GLYOH" or surrogate[i].name=="GLYOHOH")
	cout << "update: " <<  surrogate[i].name << " " << surrogate[i].Ag << " " << surrogate[i].MM << " " << surrogate[i].Aaq_bins_init(0) << " " << surrogate[i].Aaq_bins_init(0)/surrogate[i].Kaq(0)/AQinit(0) << endl;*/
 
    }
  

  if (config.solids)
    {
      solidification_bins_ssh(config,surrogate, Temperature, MMaq, AQinit, chp, ionic, LWC, factor);
      //compute_conc_org_bins_ssh(config, surrogate, Temperature,  MMaq, AQinit, chp, ionic, LWC, b, conc_org);
    }
  
  //cout << "ici: " << surrogate[config.iCO3mm].Aaq_bins_init << endl;
  //cout << surrogate[config.iHCO3m].Aaq_bins_init << endl;


  if (config.chemistry==true)
    {
      if (config.compute_organic)	
	for (i=0;i<n;i++)
	  if (surrogate[i].is_organic)
	    {
	      double atot=surrogate[i].Ag;
              if (surrogate[i].hydrophobic)
                for (b=0;b<config.nbins;++b)
                  for (ilayer=0;ilayer<config.nlayer;++ilayer)
                    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                      atot+=surrogate[i].Ap_layer_init(b,ilayer,iphase);
	    
              if (surrogate[i].hydrophilic)
                for (b=0;b<config.nbins;++b)
                  atot+=surrogate[i].Aaq_bins_init(b);

              surrogate[i].Atot=atot;
            }

      if (config.compute_inorganic or config.chemistry)
	for (i=0;i<n;i++)
	  if (surrogate[i].is_inorganic_precursor and surrogate[i].is_solid==false)
	    {
	      surrogate[i].Atot=surrogate[i].Ag;
	      if (i==config.iNH3)
		surrogate[i].Atot+=sum(surrogate[config.iNH4p].Aaq_bins_init)/surrogate[config.iNH4p].MM*surrogate[i].MM; 
	      else if (i==config.iHNO3)
		surrogate[i].Atot+=sum(surrogate[config.iNO3m].Aaq_bins_init)/surrogate[config.iNO3m].MM*surrogate[i].MM;	      
	      else if (i==config.iHCl)
		surrogate[i].Atot+=sum(surrogate[config.iClm].Aaq_bins_init)/surrogate[config.iClm].MM*surrogate[i].MM; 	
	      else if (i==config.iH2SO4)
		surrogate[i].Atot+=sum(surrogate[config.iSO4mm].Aaq_bins_init)/surrogate[config.iSO4mm].MM*surrogate[i].MM+
		  sum(surrogate[config.iHSO4m].Aaq_bins_init)/surrogate[config.iHSO4m].MM*surrogate[i].MM;
	    }
    }
  
  
  if (config.chemistry==false)
    if (config.compute_organic)
      {
	Array <int, 1> ifound;
	ifound.resize(n);
	ifound=0;
  
	for (i=0;i<n;i++)
	  if (surrogate[i].hydrophilic)
	    {
	      int j=surrogate[i].iHyd;
	      if (j>-1 and ifound(i)==0)
		{
		  ifound(i)=1;
		  ifound(j)=1;
		  int k=surrogate[j].iHyd;		  
		  if (k>-1)
		    {
		      double atot=surrogate[i].Atot0/surrogate[i].MM+surrogate[j].Atot0/surrogate[j].MM+surrogate[k].Atot0/surrogate[k].MM;


		      
		      double atot2=surrogate[i].Ag/surrogate[i].MM+surrogate[j].Ag/surrogate[j].MM+surrogate[k].Ag/surrogate[k].MM;
		      if (surrogate[i].hydrophilic)
			for (b=0;b<config.nbins;++b)
			  atot2+=surrogate[i].Aaq_bins_init(b)/surrogate[i].MM+surrogate[j].Aaq_bins_init(b)/surrogate[j].MM+surrogate[k].Aaq_bins_init(b)/surrogate[k].MM;

		      if (surrogate[i].hydrophobic)
			for (b=0;b<config.nbins;++b)
			  for (ilayer=0;ilayer<config.nlayer;++ilayer)
			    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
			      atot2+=surrogate[i].Ap_layer_init(b,ilayer,iphase)/surrogate[i].MM+surrogate[j].Ap_layer_init(b,ilayer,iphase)/surrogate[j].MM+surrogate[k].Ap_layer_init(b,ilayer,iphase)/surrogate[k].MM;

		      //if (atot2<0.01*atot)
		      //cout << "weird " << atot2 << " " << atot << endl;
		      //atot2=atot;

		      atot2=max(min(atot2,2.*atot),0.5*atot);
		      atot2=atot;

		      if (atot2>0.)
			{			  
			  surrogate[i].Ag=surrogate[i].Ag/atot2*atot;
			  surrogate[j].Ag=surrogate[j].Ag/atot2*atot;
			  surrogate[k].Ag=surrogate[k].Ag/atot2*atot;

			  surrogate[i].Atot=surrogate[i].Ag;
			  surrogate[j].Atot=surrogate[j].Ag;
			  surrogate[k].Atot=surrogate[k].Ag;
		      
			  if (surrogate[i].hydrophilic)
			    for (b=0;b<config.nbins;++b)
			      {
				surrogate[i].Aaq_bins_init(b)=surrogate[i].Aaq_bins_init(b)/atot2*atot;
				surrogate[j].Aaq_bins_init(b)=surrogate[j].Aaq_bins_init(b)/atot2*atot;
				surrogate[k].Aaq_bins_init(b)=surrogate[k].Aaq_bins_init(b)/atot2*atot;

				surrogate[i].Atot+=surrogate[i].Aaq_bins_init(b);
				surrogate[j].Atot+=surrogate[j].Aaq_bins_init(b);
				surrogate[k].Atot+=surrogate[k].Aaq_bins_init(b);
			      }

			  if (surrogate[i].hydrophobic)
			    for (b=0;b<config.nbins;++b)
			      for (ilayer=0;ilayer<config.nlayer;++ilayer)
				for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
				  {				  
				    surrogate[i].Ap_layer_init(b,ilayer,iphase)=surrogate[i].Ap_layer_init(b,ilayer,iphase)/atot2*atot;
				    surrogate[j].Ap_layer_init(b,ilayer,iphase)=surrogate[j].Ap_layer_init(b,ilayer,iphase)/atot2*atot;
				    surrogate[k].Ap_layer_init(b,ilayer,iphase)=surrogate[k].Ap_layer_init(b,ilayer,iphase)/atot2*atot;

				    surrogate[i].Atot+=surrogate[i].Ap_layer_init(b,ilayer,iphase);
				    surrogate[j].Atot+=surrogate[j].Ap_layer_init(b,ilayer,iphase);
				    surrogate[k].Atot+=surrogate[k].Ap_layer_init(b,ilayer,iphase);
				  }
			}
		    }
		  else
		    {
		      double atot=surrogate[i].Atot0/surrogate[i].MM+surrogate[j].Atot0/surrogate[j].MM;
		      double atot2=surrogate[i].Ag/surrogate[i].MM+surrogate[j].Ag/surrogate[j].MM;
		      if (surrogate[i].hydrophilic)
			for (b=0;b<config.nbins;++b)
			  atot2+=surrogate[i].Aaq_bins_init(b)/surrogate[i].MM+surrogate[j].Aaq_bins_init(b)/surrogate[j].MM;

		      if (surrogate[i].hydrophobic)
			for (b=0;b<config.nbins;++b)
			  for (ilayer=0;ilayer<config.nlayer;++ilayer)
			    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
			      atot2+=surrogate[i].Ap_layer_init(b,ilayer,iphase)/surrogate[i].MM+surrogate[j].Ap_layer_init(b,ilayer,iphase)/surrogate[j].MM;

		      //if (atot2<0.01*atot)
		      //cout << "weird " << atot2 << " " << atot << endl;
		      //atot2=atot;

		      atot2=max(min(atot2,2.*atot),0.5*atot);
		      atot2=atot;

		      if (atot2>0.)
			{			  
			  surrogate[i].Ag=surrogate[i].Ag/atot2*atot;
			  surrogate[j].Ag=surrogate[j].Ag/atot2*atot;

			  surrogate[i].Atot=surrogate[i].Ag;
			  surrogate[j].Atot=surrogate[j].Ag;
			  
			  if (surrogate[i].hydrophilic)
			    for (b=0;b<config.nbins;++b)
			      {
				surrogate[i].Aaq_bins_init(b)=surrogate[i].Aaq_bins_init(b)/atot2*atot;
				surrogate[j].Aaq_bins_init(b)=surrogate[j].Aaq_bins_init(b)/atot2*atot;

				surrogate[i].Atot+=surrogate[i].Aaq_bins_init(b);
				surrogate[j].Atot+=surrogate[j].Aaq_bins_init(b);
			      }

			  if (surrogate[i].hydrophobic)
			    for (b=0;b<config.nbins;++b)
			      for (ilayer=0;ilayer<config.nlayer;++ilayer)
				for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
				  {				  
				    surrogate[i].Ap_layer_init(b,ilayer,iphase)=surrogate[i].Ap_layer_init(b,ilayer,iphase)/atot2*atot;
				    surrogate[j].Ap_layer_init(b,ilayer,iphase)=surrogate[j].Ap_layer_init(b,ilayer,iphase)/atot2*atot;

				    surrogate[i].Atot+=surrogate[i].Ap_layer_init(b,ilayer,iphase);
				    surrogate[j].Atot+=surrogate[j].Ap_layer_init(b,ilayer,iphase);
				  }
			}
		    }
		  
		  /*
		  else
		    {
		      double atot=surrogate[i].Atot1/surrogate[i].MM+surrogate[j].Atot1/surrogate[j].MM;
		      double atot2=surrogate[i].Ag/surrogate[i].MM+surrogate[j].Ag/surrogate[j].MM;
		      if (surrogate[i].hydrophilic)
			for (b=0;b<config.nbins;++b)
			  atot2+=surrogate[i].Aaq_bins_init(b)/surrogate[i].MM+surrogate[j].Aaq_bins_init(b)/surrogate[j].MM;

		      if (surrogate[i].hydrophobic)
			for (b=0;b<config.nbins;++b)
			  for (ilayer=0;ilayer<config.nlayer;++ilayer)
			    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
			      atot2+=surrogate[i].Ap_layer_init(b,ilayer,iphase)/surrogate[i].MM+surrogate[j].Ap_layer_init(b,ilayer,iphase)/surrogate[j].MM;
		      
		      atot2=max(min(atot2,2.*atot),0.5*atot);
		      atot2=atot;
		      
		      if (atot2>0.)
			{

			  surrogate[i].Ag=surrogate[i].Ag/atot2*atot;
			  surrogate[j].Ag=surrogate[j].Ag/atot2*atot;

			  surrogate[i].Atot=surrogate[i].Ag;
			  surrogate[j].Atot=surrogate[j].Ag;		   

			  if (surrogate[i].hydrophilic)
			    for (b=0;b<config.nbins;++b)
			      {
				surrogate[i].Aaq_bins_init(b)=surrogate[i].Aaq_bins_init(b)/atot2*atot;
				surrogate[j].Aaq_bins_init(b)=surrogate[j].Aaq_bins_init(b)/atot2*atot;			      

				surrogate[i].Atot+=surrogate[i].Aaq_bins_init(b);
				surrogate[j].Atot+=surrogate[j].Aaq_bins_init(b);			     
			      }

			  if (surrogate[i].hydrophobic)
			    for (b=0;b<config.nbins;++b)
			      for (ilayer=0;ilayer<config.nlayer;++ilayer)
				for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
				  {				  
				    surrogate[i].Ap_layer_init(b,ilayer,iphase)=surrogate[i].Ap_layer_init(b,ilayer,iphase)/atot2*atot;
				    surrogate[j].Ap_layer_init(b,ilayer,iphase)=surrogate[j].Ap_layer_init(b,ilayer,iphase)/atot2*atot;				 

				    surrogate[i].Atot+=surrogate[i].Ap_layer_init(b,ilayer,iphase);
				    surrogate[j].Atot+=surrogate[j].Ap_layer_init(b,ilayer,iphase);				 
				  }
			}
		    }*/
		}
	    }

	
	for (i=0;i<n;i++)
	  if (surrogate[i].is_organic)
	    {
	      double atot=0.; //surrogate[i].Ag;
	      if (surrogate[i].hydrophobic)
		for (b=0;b<config.nbins;++b)
		  for (ilayer=0;ilayer<config.nlayer;++ilayer)
		    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		      atot+=surrogate[i].Ap_layer_init(b,ilayer,iphase);

	      if (surrogate[i].hydrophilic)
		for (b=0;b<config.nbins;++b)
		  atot+=surrogate[i].Aaq_bins_init(b);

	
	      if (atot<-surrogate[i].Atot)
		{
		  surrogate[i].Ag=surrogate[i].Atot-atot;
		}
	      else if (surrogate[i].Atot>0)
		{
		  atot+=surrogate[i].Ag;
		  if (atot==0.)
		    {
		      cout << "error " << surrogate[i].Ag << " " << surrogate[i].Atot << " " << surrogate[i].name << endl;
		      exit(0);
		    }
		  //cout << surrogate[i].name << endl;

		  if (atot>0.)
		    {
		      surrogate[i].Ag*=surrogate[i].Atot/atot;
		      if (surrogate[i].hydrophobic)
			for (b=0;b<config.nbins;++b)
			  for (ilayer=0;ilayer<config.nlayer;++ilayer)
			    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
			      surrogate[i].Ap_layer_init(b,ilayer,iphase)*=surrogate[i].Atot/atot;
		  
		      if (surrogate[i].hydrophilic)
			for (b=0;b<config.nbins;++b)
			  surrogate[i].Aaq_bins_init(b)*=surrogate[i].Atot/atot;
		    }
		}
	    }
      }

  
  if (config.chemistry==false)
    if (config.compute_inorganic)
      for (i=0;i<n;++i)
	if (surrogate[i].is_inorganic_precursor)
	  {

	    double total=surrogate[i].Ag/surrogate[i].MM;
            double total0=surrogate[i].Ag0/surrogate[i].MM;
	    for (b=0;b<config.nbins;b++)	      
	      {
		if (i==config.iHCl)
		  {
		    total+=surrogate[config.iClm].Aaq_bins_init(b)/surrogate[config.iClm].MM;
		    total0+=surrogate[config.iClm].Aaq_bins_init0(b)/surrogate[config.iClm].MM;
		  }
		else if (i==config.iNH3)
		  {
		    total+=surrogate[config.iNH4p].Aaq_bins_init(b)/surrogate[config.iNH4p].MM;
		    total0+=surrogate[config.iNH4p].Aaq_bins_init0(b)/surrogate[config.iNH4p].MM;
		  }
		else if (i==config.iHNO3)
		  {
		    total+=surrogate[config.iNO3m].Aaq_bins_init(b)/surrogate[config.iNO3m].MM;
		    total0+=surrogate[config.iNO3m].Aaq_bins_init0(b)/surrogate[config.iNO3m].MM;
		  }		
		else if (i==config.iH2SO4)
		  {
		    total+=surrogate[config.iSO4mm].Aaq_bins_init(b)/surrogate[config.iSO4mm].MM
		      +surrogate[config.iHSO4m].Aaq_bins_init(b)/surrogate[config.iHSO4m].MM;
		    total0+=surrogate[config.iSO4mm].Aaq_bins_init0(b)/surrogate[config.iSO4mm].MM
		      +surrogate[config.iHSO4m].Aaq_bins_init0(b)/surrogate[config.iHSO4m].MM;
		  }
              }

            
            //if (total-total0>0.)
            //  surrogate[i].Ag=max(surrogate[i].Ag0-(total-total0)*surrogate[i].MM,0.);
            
            if (total>0)
              {
                surrogate[i].Ag=surrogate[i].Ag/total*total0;		
                for (b=0;b<config.nbins;b++)
                  {
                    if (i==config.iHCl)
                      {
                        surrogate[config.iClm].Aaq_bins_init(b)*=total0/total;
                      }
                    else if (i==config.iNH3)
                      {
                        surrogate[config.iNH4p].Aaq_bins_init(b)*=total0/total;
                      }
                    else if (i==config.iHNO3)
                      {
                        surrogate[config.iNO3m].Aaq_bins_init(b)*=total0/total;
                      }
                    else if (i==config.iH2SO4)
                      {
			surrogate[config.iSO4mm].Aaq_bins_init(b)*=total0/total;
			surrogate[config.iHSO4m].Aaq_bins_init(b)*=total0/total;
                      }
                        
                  }
              }
	      	  
	  }

  for (i=0;i<n;i++)
    if (i!=config.iH2O)
      {
	surrogate[i].Ag=factor*surrogate[i].Ag+(1.-factor)*surrogate[i].Ag1;
	for (b=0;b<config.nbins;b++)
	  if (surrogate[i].hydrophilic)
	    surrogate[i].Aaq_bins_init(b)=factor*surrogate[i].Aaq_bins_init(b)+(1.-factor)*surrogate[i].Aaq_bins(b);
	
	for (b=0;b<config.nbins;++b)
	  for (ilayer=0;ilayer<config.nlayer;++ilayer)
	    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
	      if (surrogate[i].hydrophobic)
		surrogate[i].Ap_layer_init(b,ilayer,iphase)=
		  factor*surrogate[i].Ap_layer_init(b,ilayer,iphase)+(1.-factor)*surrogate[i].Ap_layer(b,ilayer,iphase);

	for (b=0;b<config.nbins;b++)
	  if (surrogate[i].is_solid)
	    surrogate[i].Asol_bins_init(b)=factor*surrogate[i].Asol_bins_init(b)+(1.-factor)*surrogate[i].Asol_bins(b);
      }
  
  sum1=1.0;
  conceq=surrogate[config.iH2O].Ag;  
  for (b=0;b<config.nbins;++b) //mass of water absorbed by organics in the aqueous phase
    {	  	 	
      if (surrogate[config.iH2O].hydrophilic)
	if (LWCtot>config.LWClimit)  
	  {
	    conceq+=surrogate[config.iH2O].Aaq_bins_init(b);
	    sum1+=surrogate[config.iH2O].Kaq(b)*AQinit(b);
	  }
      
      if (surrogate[config.iH2O].hydrophobic) // and config.compute_organic)
	for (ilayer=0;ilayer<config.nlayer;++ilayer)
	  for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
	    {
	      conceq+=surrogate[config.iH2O].Ap_layer_init(b,ilayer,iphase);
	      sum1+=surrogate[config.iH2O].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase);                  
	    }
    }

  if (surrogate[config.iH2O].hydrophilic and config.hygroscopicity)
    if (LWCtot>config.LWClimit)
      for (b=0;b<config.nbins;++b) //mass of water absorbed by organics in the aqueous phase
	{	    	    	    
	  if (config.activity_model=="unifac" or config.compute_inorganic)
	    {
              
	      surrogate[config.iH2O].Aaq_bins_init(b)=factor*
		max(conceq*surrogate[config.iH2O].Kaq(b)*AQinit(b)/sum1 - LWC(b), 0.0) +
		(1.0-factor)*surrogate[config.iH2O].Aaq_bins_init(b);
	    }
	  else	      
	    {
	      Mads=AQinit(b)-LWC(b)-conc_inorganic(b);
	      surrogate[config.iH2O].Aaq_bins_init(b)=factor*
		surrogate[config.iH2O].Atot*surrogate[config.iH2O].Kaq(b)*Mads/sum1+
		(1.0-factor)*surrogate[config.iH2O].Aaq_bins_init(b);		
	    }
	}
  
  if (surrogate[config.iH2O].hydrophobic and config.compute_organic and config.hygroscopicity)
    for (b=0;b<config.nbins;++b)
      for (ilayer=0;ilayer<config.nlayer;++ilayer)
	for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
	  //if (surrogate[config.iH2O].time(b,ilayer,iphase)<tequilibrium)	      
	  surrogate[config.iH2O].Ap_layer_init(b,ilayer,iphase)=factor*
	    conceq*surrogate[config.iH2O].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)/sum1+
	    (1.0-factor)*surrogate[config.iH2O].Ap_layer_init(b,ilayer,iphase);


  for (b=0;b<config.nbins;++b)
    {
      conc_inorganic(b)=0.0;
      for (i=0;i<n;++i)
        if(surrogate[i].is_organic==false and surrogate[i].is_inorganic_precursor==false and i!=config.iH2O)
          conc_inorganic(b)+=surrogate[i].Aaq_bins_init(b);
    }
  
  //compute total concentration of the aqueous phase
  for (b=0;b<config.nbins;++b)
    {
      AQ(b)=LWC(b)+conc_inorganic(b);
      for (i=0;i<n;++i)
        if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophilic)
          AQ(b)+=surrogate[i].Aaq_bins_init(b);
    }

  //compute total concentrations of the organic phase
  if (config.compute_organic)
    {
      for (b=0;b<config.nbins;++b)
        for (ilayer=0;ilayer<config.nlayer;++ilayer)
          {
            for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
              {
                MO(b,ilayer,iphase)=0.0;
                for (i=0;i<n;++i)
                  if(surrogate[i].hydrophobic)
                    MO(b,ilayer,iphase)+=surrogate[i].Ap_layer_init(b,ilayer,iphase);
              }

            //make sure that concentrations in a bin are at least equal to MOmin
            sum1=0.0;
            for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
              sum1+=MO(b,ilayer,iphase);

            if (sum1>0.0)
              for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                MO(b,ilayer,iphase)=max(MO(b,ilayer,iphase),
                                        config.MOmin*config.Vlayer(ilayer)*MO(b,ilayer,iphase)/sum1);
            else
              for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                if(iphase==0)
                  MO(b,ilayer,iphase)=config.MOmin*config.Vlayer(ilayer);
                else
                  MO(b,ilayer,iphase)=0.0;
          }
    }
  else
    {
      for (b=0;b<config.nbins;++b)
        for (ilayer=0;ilayer<config.nlayer;++ilayer)
          {
            for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
              {
                MO(b,ilayer,iphase)=0.0;
                for (i=0;i<n;++i)
                  if(surrogate[i].hydrophobic and (surrogate[i].is_organic or i==config.iH2O))
                    MO(b,ilayer,iphase)+=surrogate[i].Ap_layer_init(b,ilayer,iphase);
              }
          }
            
    }
}






void twostep_aqorg_repart_ssh(model_config &config, vector<species>& surrogate, double &tequilibrium,
                              Array <double, 1> &AQinit,  Array <double, 1> &AQ,
                              Array <double, 1> &conc_inorganic, Array <double, 1> &ionic,
                              Array <double, 1> &ionic_organic, Array <double, 1> &organion,
                              Array <double, 1> &chp, Array <double, 1> &LWC,
                              Array <double, 3> &MOinit,  Array <double, 3> &MO, Array <double, 3> &MOW,
                              double &Temperature, double &RH, Array <double, 1> &MMaq,
                              bool compute_activity_coefficients, double factor, double t, double deltat,
                              int index, int b)
{
  //This routine computes the aqueous-phase and organic-phase concentrations being at equilibrium
  //FOR THE COUPLED SYSTEM
  //factor: weight factor to decrease variations in case of strong variations of composition
  //                between two iterations
  //compute_activity_coefficients: should the activity coefficients be recomputed?
  //tequilibrium: criterium to distinguish cases under equilibrium (characteristic time lower
  //                  than tequilibrium)
  int n=surrogate.size();
  int i,ilayer,iphase;
  //double Kaq;
  Array <double, 1> Kp,chp0,chp_save;
  Kp.resize(config.max_number_of_phases);
  chp0.resize(config.nbins);
  chp_save.resize(config.nbins);
  
  double conc_equilibrium;
  double sum1=0.0;
  bool is_equilibrium=false;
  double LWCtot=0.0;
  int index_b; //,index_layer;
  double Mads,conceq;
  int iter,index2;
  double error_chp,conc_org,total,total2;
  double tiny=1.0e-26;
  
  LWCtot+=LWC(b);

  chp0(b)=chp(b);
  chp_save(b)=chp(b);
      
  iter=0;
  error_chp=1000.;

  compute_conc_org_bins_ssh(config, surrogate, Temperature,  MMaq, AQinit, chp, ionic, LWC, b, conc_org);

  double inorganion=0.0;
  for (i=0;i<n;++i)
    if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)
      {
        surrogate[i].molality=surrogate[i].Aaq_bins_init(b)/surrogate[i].MM/conc_org*1000.0;           
        if (i!=config.iHp and i!=config.iOHm)
          inorganion-=surrogate[i].molality*surrogate[i].charge; 
      }

  organion(b)=0.0;      
  for (i=0;i<n;++i)
    if (surrogate[i].hydrophilic and surrogate[i].is_organic)
      if(surrogate[i].nonvolatile==false and (surrogate[i].aqt==1 or surrogate[i].aqt==2))
        {
          double fion1=0.;
          double fion2=0.;
          double Kploc=surrogate[i].Kp_eff_aqrealdyn_ssh(config, Temperature, ionic(b), chp(b),surrogate[config.iHp].LR(b),
                                                         surrogate[config.iHp].SRMR(b),MMaq(b), fion1, fion2, b);
            //surrogate[i].Kp_eff_aq_ssh(config, Temperature, ionic(b), chp(b),surrogate[config.iHp].LR(b),
            //                                     surrogate[config.iHp].SRMR(b),MMaq(b), fion1, fion2);
          //molality1: molality of ions HA- or A-  
          double molality1=surrogate[i].Aaq_bins_init(b)*fion1/surrogate[i].MM/conc_org*1000.0;
          //molality2: molality of ions A2-
          double molality2=surrogate[i].Aaq_bins_init(b)*fion2/surrogate[i].MM/conc_org*1000.0;
          //compute ionic_organic and organion
          organion(b)+=molality1+2*molality2;
        }

  if (config.compute_inorganic==false)
    {

      double chp_new=0.5*(organion(b)+inorganion+pow(pow(organion(b)+inorganion,2)+4*config.Ke/surrogate[config.iHp].gamma_aq_bins(b)/surrogate[config.iOHm].gamma_aq_bins(b),0.5));
      //cout << "new " << chp_new << endl;
      double change_ph=1.1;
      chp_new=max(max(min(chp_new,change_ph*chp(b)),1./change_ph*chp(b)),1.e-14);
      chp_new=min(chp_new,100.);
      //chp(b)=chp_new;
      chp(b)=factor*chp_new+(1.0-factor)*chp(b);
      //chp(b)=pow(chp_new,factor)*pow(chp(b),1-factor);
      if (chp(b)==0.0)
        chp(b)=pow(10.0,-5.6);
      //cout << b << " " << chp_new << " " << chp(b) << " " << AQinit(b) << " " << organion(b) << " " << inorganion << endl;
    }

  compute_ionic_strenght2_ssh(config, surrogate, Temperature, AQinit(b), conc_inorganic(b), ionic(b), chp(b),
                              organion(b), ionic_organic(b), conc_org, 1.0, config.Ke/surrogate[config.iHp].gamma_aq_bins(b)/surrogate[config.iOHm].gamma_aq_bins(b));

      
  if (compute_activity_coefficients)
    {
      if (config.compute_organic)
        {
          activity_coefficients_dyn_org_bins_ssh(config, surrogate, Temperature, MOW, b);
          compute_kp_org_bins_ssh(config, surrogate, MOinit, Temperature, MOW, b);
        }
     
      if (LWCtot>config.LWClimit)
        {
          activity_coefficients_dyn_aq_bins_ssh(config, surrogate, Temperature,AQinit,
                                                MOinit, conc_inorganic, ionic, ionic_organic,
                                                organion,chp,LWC,MMaq, factor, deltat, 1, b);
          compute_kp_aq_ssh(config, surrogate, Temperature, ionic, chp, AQinit, MMaq, b, b+1);
        }    
    }

  double apnew;
  for (i=0;i<n;++i)
    if (surrogate[i].is_organic or i==config.iH2O)
      if(surrogate[i].hydrophobic and surrogate[i].hydrophilic)   
        {
          double atot=0.;
          double ktot=0.;
          for (ilayer=0;ilayer<config.nlayer;++ilayer)
            for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
              {
                atot+=surrogate[i].Ap_layer_init(b,ilayer,iphase);
                ktot+=surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase);
              }
            
          if (LWCtot>config.LWClimit)
            {
              atot+=surrogate[i].Aaq_bins_init(b);
              ktot+=surrogate[i].Kaq(b)*AQinit(b);
            }

          if (ktot>0.)
            {
              for (ilayer=0;ilayer<config.nlayer;++ilayer)
                for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                  {
                    apnew=atot*surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)/ktot;
                    surrogate[i].Ap_layer_init(b,ilayer,iphase)=factor*apnew+(1.0-factor)*surrogate[i].Ap_layer_init(b,ilayer,iphase);
                  }
                
              if (LWCtot>config.LWClimit)
                {
                  apnew=atot*surrogate[i].Kaq(b)*AQinit(b)/ktot;
                  surrogate[i].Aaq_bins_init(b)=factor*apnew+(1.0-factor)*surrogate[i].Aaq_bins_init(b);
                }
               
            }
        }

  conc_inorganic(b)=0.0;
  for (i=0;i<n;++i)
    if(surrogate[i].is_organic==false and surrogate[i].is_inorganic_precursor==false and i!=config.iH2O)
      conc_inorganic(b)+=surrogate[i].Aaq_bins_init(b);
  
  //compute total concentration of the aqueous phase
  AQ(b)=LWC(b)+conc_inorganic(b);
  for (i=0;i<n;++i)
    if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophilic)
      AQ(b)+=surrogate[i].Aaq_bins_init(b);

  for (ilayer=0;ilayer<config.nlayer;++ilayer)
    {
      for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
        {
          MO(b,ilayer,iphase)=0.0;
          for (i=0;i<n;++i)
            if(surrogate[i].hydrophobic)
              MO(b,ilayer,iphase)+=surrogate[i].Ap_layer_init(b,ilayer,iphase);
        }
        
      //make sure that concentrations in a bin are at least equal to MOmin
      sum1=0.0;
      for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
        sum1+=MO(b,ilayer,iphase);
        
      if (sum1>0.0)
        for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
          MO(b,ilayer,iphase)=max(MO(b,ilayer,iphase),
                                  config.MOmin*config.Vlayer(ilayer)*MO(b,ilayer,iphase)/sum1);
      else
        for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
          if(iphase==0)
            MO(b,ilayer,iphase)=config.MOmin*config.Vlayer(ilayer);
          else
            MO(b,ilayer,iphase)=0.0;
    }
}

void equilibrium_inorg_ssh(model_config &config, vector<species>& surrogate, double &tequilibrium,
			   Array <double, 1> &AQinit, Array <double, 1> &conc_inorganic, 
			   Array <double, 1> &chp, Array <double, 1> &LWC, Array <double, 1> &ionic, 
			   double &Temperature, double &RH, Array <double, 1> &MMaq, double factor)
{
  int n=surrogate.size();
  int i,b;
  //double kelvin_effect=1.0;
  double conc_equilibrium;
  double sum=0.0;
  bool is_equilibrium=false;

  //This routine computes the aqueous-phase concentrations of inorganic compounds being at equilibrium
  //The inorganic aerosol is assumed metastable
  //compute_activity_coefficients: should the activity coefficients be recomputed?
  //tequilibrium: criterium to distinguish cases under equilibrium (characteristic time lower
  //                  than tequilibrium)

  Array <double,1> conc_org;
  conc_org.resize(config.nbins);
  for (b=0;b<config.nbins;b++)
    compute_conc_org_bins_ssh(config, surrogate, Temperature,  MMaq, AQinit, chp, ionic, LWC, b, conc_org(b));   

  for (i=0;i<n;++i)
    if(surrogate[i].is_inorganic_precursor and surrogate[i].is_solid==false)
      { 
        is_equilibrium=false;
        sum=1.0;
        conc_equilibrium=surrogate[i].Ag;
            
        for (b=0;b<config.nbins;++b)
          if (surrogate[i].time_aq(b) < tequilibrium)
            {
              is_equilibrium=true;
              if (i==config.iH2SO4)
                conc_equilibrium+=surrogate[config.iHSO4m].Aaq_bins_init(b)/surrogate[config.iHSO4m].MM*surrogate[i].MM
                  +surrogate[config.iSO4mm].Aaq_bins_init(b)/surrogate[config.iSO4mm].MM*surrogate[i].MM;
              else if (i==config.iNH3)
                conc_equilibrium+=surrogate[config.iNH4p].Aaq_bins_init(b)/surrogate[config.iNH4p].MM*surrogate[i].MM;
              else if (i==config.iHNO3)
                conc_equilibrium+=surrogate[config.iNO3m].Aaq_bins_init(b)/surrogate[config.iNO3m].MM*surrogate[i].MM;	      
              else if (i==config.iHCl)
                conc_equilibrium+=surrogate[config.iClm].Aaq_bins_init(b)/surrogate[config.iClm].MM*surrogate[i].MM;
                  
              //sum of Aaq/Ag ratios + 1
              sum+=surrogate[i].Kaq(b)*conc_org(b);
            }

        //compute gas phase concentrations
        if (is_equilibrium)
          {	
            
            if (i==config.iH2SO4)
              surrogate[i].Ag=0.0;
            else
              surrogate[i].Ag=factor*conc_equilibrium/sum+(1.0-factor)*surrogate[i].Ag;
            
            for (b=0;b<config.nbins;++b)
              {
                if (surrogate[i].time_aq(b)<tequilibrium)
                  {
                    if (i==config.iNH3)
                      {
                        if (surrogate[config.iNH4p].Aaq_bins_init(b)>0.0)
                          surrogate[config.iNH4p].Aaq_bins_init(b)=factor*surrogate[config.iNH4p].MM/surrogate[i].MM*conc_equilibrium*
                            surrogate[i].Kaq(b)*conc_org(b)/sum+(1.0-factor)*surrogate[config.iNH4p].Aaq_bins_init(b);
                        else
                          {
                            surrogate[config.iNH4p].Aaq_bins_init(b)=0.01*factor*surrogate[config.iNH4p].MM/surrogate[i].MM*conc_equilibrium*surrogate[i].Kaq(b)*conc_org(b)/sum;
                            surrogate[i].Ag+=99*surrogate[config.iNH4p].Aaq_bins_init(b);
                          }                     
                      }
                    else if (i==config.iHNO3)
                      {                        
                        if (surrogate[config.iNO3m].Aaq_bins_init(b)>0.0)
                          surrogate[config.iNO3m].Aaq_bins_init(b)=factor*surrogate[config.iNO3m].MM/surrogate[i].MM*conc_equilibrium*
                            surrogate[i].Kaq(b)*conc_org(b)/sum+(1.0-factor)*surrogate[config.iNO3m].Aaq_bins_init(b);
                        else
                          {
                            surrogate[config.iNO3m].Aaq_bins_init(b)=0.01*factor*surrogate[config.iNO3m].MM/surrogate[i].MM*conc_equilibrium*
                              surrogate[i].Kaq(b)*conc_org(b)/sum;
                            surrogate[i].Ag+=99*surrogate[config.iNO3m].Aaq_bins_init(b);
                          }
                      }
                    else if (i==config.iHCl)
                      {
                        if (surrogate[config.iClm].Aaq_bins_init(b)>0.0)
                          surrogate[config.iClm].Aaq_bins_init(b)=factor*surrogate[config.iClm].MM/surrogate[i].MM*conc_equilibrium*
                            surrogate[i].Kaq(b)*conc_org(b)/sum+(1.0-factor)*surrogate[config.iClm].Aaq_bins_init(b);
                        else
                          {
                            surrogate[config.iClm].Aaq_bins_init(b)=0.01*factor*surrogate[config.iClm].MM/surrogate[i].MM*conc_equilibrium*
                              surrogate[i].Kaq(b)*conc_org(b)/sum;
                            surrogate[i].Ag+=99*surrogate[config.iClm].Aaq_bins_init(b);
                          }
                      }
                  }
              }
          }

        if (i==config.iH2SO4)
          {
            double total;
            for (b=0;b<config.nbins;++b)
              {
                if (surrogate[i].time_aq(b)<tequilibrium)
                  total=conc_equilibrium/sum*surrogate[i].Kaq(b)*conc_org(b);
                else
                  total=(surrogate[config.iHSO4m].Aaq_bins_init(b)/surrogate[config.iHSO4m].MM+surrogate[config.iSO4mm].Aaq_bins_init(b)/surrogate[config.iSO4mm].MM)*surrogate[i].MM;
		double Keq=surrogate[i].Kequilibrium_ssh(Temperature)/chp(b)*surrogate[config.iHSO4m].gamma_aq_bins(b)/
		  (surrogate[config.iHp].gamma_aq_bins(b)*surrogate[config.iSO4mm].gamma_aq_bins(b));
                surrogate[config.iHSO4m].Aaq_bins_init(b)=(1.0-factor)*surrogate[config.iHSO4m].Aaq_bins_init(b)+factor*total*surrogate[config.iHSO4m].MM/surrogate[i].MM*1.0/(1.0+Keq);
                surrogate[config.iSO4mm].Aaq_bins_init(b)=(1.0-factor)*surrogate[config.iSO4mm].Aaq_bins_init(b)+factor*total*surrogate[config.iSO4mm].MM/surrogate[i].MM*Keq/(1.0+Keq);
              }
          }
      }

  for (b=0;b<config.nbins;++b)
    {
      conc_inorganic(b)=0.0;
      for (i=0;i<n;++i)
        if(surrogate[i].is_organic==false and surrogate[i].is_inorganic_precursor==false and i!=config.iH2O)
          conc_inorganic(b)+=surrogate[i].Aaq_bins_init(b);
    }
}

void equilibrium_aq_ssh(model_config &config, vector<species>& surrogate, double &tequilibrium,
			Array <double, 1> &AQinit,  Array <double, 1> &AQ,
			Array <double, 3> &MOinit,
			Array <double, 1> &conc_inorganic, Array <double, 1> &ionic,
			Array <double, 1> &ionic_organic, Array <double, 1> &organion,
			Array <double, 1> &chp, Array <double, 1> &LWC,
			double &Temperature, double &RH, Array <double, 1> &MMaq,
			bool compute_activity_coefficients, double factor)
{
  //This routine computes the aqueous-phase concentrations being at equilibrium
  //FOR THE UNCOUPLED SYSTEM
  //compute_activity_coefficients: should the activity coefficients be recomputed?
  //tequilibrium: criterium to distinguish cases under equilibrium (characteristic time lower
  //                  than tequilibrium)
  int n=surrogate.size();
  int i,b;
  double Mads;
  double conc_equilibrium;
  double sum=0.0;
  bool is_equilibrium=false; //indicates if there is concentrations at equilibrium

  //compute activity coefficients  
  if (compute_activity_coefficients)
    {
      activity_coefficients_dyn_aq_ssh(config, surrogate, Temperature,AQinit,MOinit,
				       conc_inorganic, ionic, ionic_organic,
				       organion,chp,LWC,MMaq, factor, 0., 0);
      
      compute_kp_aq_ssh(config, surrogate, Temperature, ionic, chp, AQinit, MMaq, 0, config.nbins);
    }
  
  if (config.compute_organic)
    for (i=0;i<n;++i)
      if(surrogate[i].is_organic and surrogate[i].hydrophilic)
        {
          is_equilibrium=false;
          sum=1.0;
          conc_equilibrium=surrogate[i].Atot;
		
          for (b=0;b<config.nbins;++b)
            if (surrogate[i].time_aq(b) < tequilibrium)
              is_equilibrium=true;
		
          for (b=0;b<config.nbins;++b)
            {
              //sum of concentrations at equilibrium
              if (surrogate[i].time_aq(b)>=tequilibrium)
                conc_equilibrium-=surrogate[i].Aaq_bins_init(b); 

              //sum of Aaq/Ag ratios +1
              if (surrogate[i].time_aq(b)<tequilibrium)
                sum+=surrogate[i].Kaq(b)*AQinit(b);
            }

          //compute gas phase concentrations
          if (is_equilibrium)
            {			
              surrogate[i].Ag=conc_equilibrium/sum;
              for (b=0;b<config.nbins;++b)
                if (surrogate[i].time_aq(b)<tequilibrium)
                  surrogate[i].Aaq_bins_init(b)=min(2.0*max(surrogate[i].Aaq_bins_init(b),1.0e-10),
                                                    factor*surrogate[i].Kaq(b)*AQinit(b)/sum*conc_equilibrium+(1.0-factor)*surrogate[i].Aaq_bins_init(b));
            }
          else
            {
              surrogate[i].Ag=surrogate[i].Atot;
              for (b=0;b<config.nbins;++b)
                surrogate[i].Ag-=surrogate[i].Aaq_bins_init(b);
            }
        }
  
  if (config.compute_inorganic)
    equilibrium_inorg_ssh(config, surrogate, tequilibrium, AQinit, conc_inorganic, chp, LWC, ionic,
			  Temperature, RH, MMaq, factor);

  sum=1.0;
  for (b=0;b<config.nbins;++b) //mass of water absorbed by organics in the aqueous phase	 	  
    if (surrogate[config.iH2O].hydrophilic)      
      if (surrogate[config.iH2O].time_aq(b)<tequilibrium)	    
	sum+=surrogate[config.iH2O].Kaq(b)*AQinit(b);
  
  if (surrogate[config.iH2O].hydrophilic and config.hygroscopicity)
    for (b=0;b<config.nbins;++b) //mass of water absorbed by organics in the aqueous phase
      if (surrogate[config.iH2O].time_aq(b)<tequilibrium)
	{	    
	  if (config.activity_model=="unifac" or config.compute_inorganic)
            {
              surrogate[config.iH2O].Aaq_bins_init(b)=factor*
                max(surrogate[config.iH2O].Atot*surrogate[config.iH2O].Kaq(b)*AQinit(b)/sum-LWC(b),0.0)+
                (1.0-factor)*surrogate[config.iH2O].Aaq_bins_init(b);
            }
	  else	      
	    {
              Mads=AQinit(b)-LWC(b)-conc_inorganic(b);
              surrogate[config.iH2O].Aaq_bins_init(b)=factor*
                surrogate[config.iH2O].Atot*surrogate[config.iH2O].Kaq(b)*Mads/sum+
                (1.0-factor)*surrogate[config.iH2O].Aaq_bins_init(b);
	    }
	}
  
  //compute total concentrations of the aqueous phase
  for (b=0;b<config.nbins;++b)
    {
      AQ(b)=LWC(b)+conc_inorganic(b);
      for (i=0;i<n;++i)
        if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophilic)
          AQ(b)+=surrogate[i].Aaq_bins_init(b);
    }
}

void equilibrium_tot_ssh(model_config &config, vector<species>& surrogate, double &tequilibrium,
			 Array <double, 1> &AQinit,  Array <double, 1> &AQ,
			 Array <double, 1> &conc_inorganic, Array <double, 1> &ionic,
			 Array <double, 1> &ionic_organic, Array <double, 1> &organion,
			 Array <double, 1> &chp, Array <double, 1> &LWC,
			 Array <double, 3> &MOinit,  Array <double, 3> &MO, Array <double, 3> &MOW,
			 double &Temperature, double &RH, Array <double, 1> &MMaq,
			 bool compute_activity_coefficients, double factor)
{
  //This routine computes the aqueous-phase and organic-phase concentrations being at equilibrium
  //FOR THE COUPLED SYSTEM
  //factor: weight factor to decrease variations in case of strong variations of composition
  //                between two iterations
  //compute_activity_coefficients: should the activity coefficients be recomputed?
  //tequilibrium: criterium to distinguish cases under equilibrium (characteristic time lower
  //                  than tequilibrium)
  int n=surrogate.size();
  int i,b,ilayer,iphase;
  //double Kaq;
  Array <double, 1> Kp;
  Kp.resize(config.max_number_of_phases);
  
  double conc_equilibrium;
  double sum=0.0;
  bool is_equilibrium=false;
  double LWCtot=0.0;
  int index_b; //,index_layer;
  double Mads,conceq;
  for (b=0;b<config.nbins;++b)
    LWCtot+=LWC(b);

  //cout << chp << endl;

  //compute activity coefficients
  if (compute_activity_coefficients)
    {
      if (config.compute_organic)
        {
          activity_coefficients_dyn_org_ssh(config, surrogate, Temperature, MOW);
          compute_kp_org_ssh(config, surrogate, MOinit, Temperature, MOW);
        }
     
      if (LWCtot>config.LWClimit)
        {
          activity_coefficients_dyn_aq_ssh(config, surrogate, Temperature,AQinit,
					   MOinit, conc_inorganic, ionic, ionic_organic,
					   organion,chp,LWC,MMaq, factor, 0., 0);
          compute_kp_aq_ssh(config, surrogate, Temperature, ionic, chp, AQinit, MMaq, 0, config.nbins);
        }    
    }
  
  if (config.compute_organic)
    for (i=0;i<n;++i)
      if(surrogate[i].is_organic)
        {	
          is_equilibrium=false; //indicates if there is concentrations at equilibrium
		
          sum=1.0;
          conc_equilibrium=surrogate[i].Atot;
          index_b=0;	 
          //compute conc_equilibrium: sum of concentrations at equilibrium
          // and compute sum: sum of Aaq/Ag + Ap/Ag ratios +1
          if(surrogate[i].hydrophilic and i!=config.iH2O and LWCtot>config.LWClimit)
            for (b=0;b<config.nbins;++b)
              {
                if (surrogate[i].time_aq(b)>=tequilibrium)
                  conc_equilibrium-=surrogate[i].Aaq_bins_init(b);
	     			  
                if (surrogate[i].time_aq(b)<tequilibrium)
                  {
                    is_equilibrium=true;
                    sum+=surrogate[i].Kaq(b)*AQinit(b);
                  }
              }

          if(surrogate[i].hydrophobic)
            for (b=0;b<config.nbins;++b)
              for (ilayer=0;ilayer<config.nlayer;++ilayer)
                for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                  { 
                    if (surrogate[i].time(b,ilayer,iphase)>=tequilibrium)
                      conc_equilibrium-=surrogate[i].Ap_layer_init(b,ilayer,iphase);
				  
                    if (surrogate[i].time(b,ilayer,iphase)<tequilibrium)
                      {
                        is_equilibrium=true;
                        sum+=surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase);
                      }
                  }

          //compute particulate-phase and gas-phase concentrations with a weight factor
          //new concentration= factor*(computed concentrations) + (1-factor)*(old concentrations)
          //
          if (is_equilibrium)
            {			
              surrogate[i].Ag=factor*conc_equilibrium/sum+(1.0-factor)*surrogate[i].Ag;
              if(LWCtot>config.LWClimit)
                {
                  if(surrogate[i].hydrophilic and i!=config.iH2O)
                    for (b=0;b<config.nbins;++b)
                      if (surrogate[i].time_aq(b)<tequilibrium)
                        surrogate[i].Aaq_bins_init(b)=min(2.0*max(surrogate[i].Aaq_bins_init(b),1.0e-10),factor*(surrogate[i].Kaq(b)*AQinit(b))/sum
							  *conc_equilibrium+(1.0-factor)*surrogate[i].Aaq_bins_init(b));
                }
              else
                for (b=0;b<config.nbins;++b)
                  surrogate[i].Aaq_bins_init(b)=0.0;

              if(surrogate[i].hydrophobic)
                for (b=0;b<config.nbins;++b)
                  for (ilayer=0;ilayer<config.nlayer;++ilayer)
                    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                      if (surrogate[i].time(b,ilayer,iphase)<tequilibrium)
                        surrogate[i].Ap_layer_init(b,ilayer,iphase)=min(2.0*max(surrogate[i].Ap_layer_init(b,ilayer,iphase),1.0e-10),
                                                                        factor*(surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase))/sum*
                                                                        conc_equilibrium+(1.0-factor)*surrogate[i].Ap_layer_init(b,ilayer,iphase));
            }
          else
            {
              surrogate[i].Ag=surrogate[i].Atot;
              if(surrogate[i].hydrophilic and i!=config.iH2O and LWCtot>config.LWClimit)
                for (b=0;b<config.nbins;++b)
                  surrogate[i].Ag-=surrogate[i].Aaq_bins_init(b);

              if(surrogate[i].hydrophobic)
                for (b=0;b<config.nbins;++b)
                  for (ilayer=0;ilayer<config.nlayer;++ilayer)
                    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                      surrogate[i].Ag-=surrogate[i].Ap_layer_init(b,ilayer,iphase);
            }
	}
  
  //compute absorption of water by organics in the aqueous phase
  //double kelvin_effect;
  if (config.compute_inorganic)
    equilibrium_inorg_ssh(config, surrogate, tequilibrium, AQinit, conc_inorganic, chp, LWC, ionic,
			  Temperature, RH, MMaq, factor);

  sum=1.0;
  conceq=surrogate[config.iH2O].Ag;  
  for (b=0;b<config.nbins;++b) //mass of water absorbed by organics in the aqueous phase
    {	  	 	
      if (surrogate[config.iH2O].hydrophilic)
	if (LWCtot>config.LWClimit)
	  if (surrogate[config.iH2O].time_aq(b)<tequilibrium)	   
	    { 
	      conceq+=surrogate[config.iH2O].Aaq_bins_init(b);
	      sum+=surrogate[config.iH2O].Kaq(b)*AQinit(b);
	    }
      
      if (surrogate[config.iH2O].hydrophobic and config.compute_organic)
	for (ilayer=0;ilayer<config.nlayer;++ilayer)
	  for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
	    if (surrogate[config.iH2O].time(b,ilayer,iphase)<tequilibrium)
	      {
		conceq+=surrogate[config.iH2O].Ap_layer_init(b,ilayer,iphase);
		sum+=surrogate[config.iH2O].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase);                  
	      }
    }

  if (surrogate[config.iH2O].hydrophilic and config.hygroscopicity)
    if (LWCtot>config.LWClimit)
      for (b=0;b<config.nbins;++b) //mass of water absorbed by organics in the aqueous phase
	if (surrogate[config.iH2O].time_aq(b)<tequilibrium)
	  {	    	    	    
	    if (config.activity_model=="unifac" or config.compute_inorganic)
	      {
		surrogate[config.iH2O].Aaq_bins_init(b)=factor*
		  max(conceq*surrogate[config.iH2O].Kaq(b)*AQinit(b)/sum - LWC(b), 0.0) +
		  (1.0-factor)*surrogate[config.iH2O].Aaq_bins_init(b);
	      }
	    else	      
	      {
		Mads=AQinit(b)-LWC(b)-conc_inorganic(b);
		surrogate[config.iH2O].Aaq_bins_init(b)=factor*
		  surrogate[config.iH2O].Atot*surrogate[config.iH2O].Kaq(b)*Mads/sum+
		  (1.0-factor)*surrogate[config.iH2O].Aaq_bins_init(b);		
	      }
	  }
  
  if (surrogate[config.iH2O].hydrophobic and config.compute_organic and config.hygroscopicity)
    for (b=0;b<config.nbins;++b)
      for (ilayer=0;ilayer<config.nlayer;++ilayer)
	for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
	  if (surrogate[config.iH2O].time(b,ilayer,iphase)<tequilibrium)	      
	    surrogate[config.iH2O].Ap_layer_init(b,ilayer,iphase)=factor*
              conceq*surrogate[config.iH2O].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)/sum+
	      (1.0-factor)*surrogate[config.iH2O].Ap_layer_init(b,ilayer,iphase);
  
  //compute total concentration of the aqueous phase
  for (b=0;b<config.nbins;++b)
    {
      AQ(b)=LWC(b)+conc_inorganic(b);
      for (i=0;i<n;++i)
        if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophilic)
          AQ(b)+=surrogate[i].Aaq_bins_init(b);
    }

  //compute total concentrations of the organic phase
  if (config.compute_organic)
    {
      for (b=0;b<config.nbins;++b)
        for (ilayer=0;ilayer<config.nlayer;++ilayer)
          {
            for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
              {
                MO(b,ilayer,iphase)=0.0;
                for (i=0;i<n;++i)
                  if(surrogate[i].hydrophobic)
                    MO(b,ilayer,iphase)+=surrogate[i].Ap_layer_init(b,ilayer,iphase);
              }

            //make sure that concentrations in a bin are at least equal to MOmin
            sum=0.0;
            for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
              sum+=MO(b,ilayer,iphase);

            if (sum>0.0)
              for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                MO(b,ilayer,iphase)=max(MO(b,ilayer,iphase),
                                        config.MOmin*config.Vlayer(ilayer)*MO(b,ilayer,iphase)/sum);
            else
              for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                if(iphase==0)
                  MO(b,ilayer,iphase)=config.MOmin*config.Vlayer(ilayer);
                else
                  MO(b,ilayer,iphase)=0.0;
          }
    }
  else
    for (b=0;b<config.nbins;++b)
      for (ilayer=0;ilayer<config.nlayer;++ilayer)
        for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
          MO(b,ilayer,iphase)=MOinit(b,ilayer,iphase);
}

void equilibrium_org_ssh(model_config &config, vector<species>& surrogate, double &tequilibrium,
			 Array <double, 3> &MOinit,  Array <double, 3> &MO,
			 double &Temperature, Array <double, 3> &MOW,
			 bool compute_activity_coefficients, double factor)
{
  //This routine computes the organic-phase concentrations being at equilibrium
  //FOR THE UNCOUPLED SYSTEM
  //compute_activity_coefficients: should the activity coefficients be recomputed?
  //tequilibrium: criterium to distinguish cases under equilibrium (characteristic time lower
  //                  than tequilibrium)
  int n=surrogate.size();
  int i,b,ilayer,iphase;
  //int index_b,index_layer;
  Array<double, 1> Kp;
  Kp.resize(config.max_number_of_phases);
  double conc_equilibrium;
  double sum=0.0;
  bool is_equilibrium=false;

  //compute activity coefficients
  if (config.compute_organic and compute_activity_coefficients)
    {
      activity_coefficients_dyn_org_ssh(config, surrogate, Temperature, MOW);
      compute_kp_org_ssh(config, surrogate, MOinit, Temperature, MOW);
    }
  
  if (config.compute_organic)
    {
      for (i=0;i<n;++i)
        if(surrogate[i].is_organic and surrogate[i].hydrophobic)
          {
            is_equilibrium=false;
            sum=1.0;
            conc_equilibrium=surrogate[i].Atot;
            
            for (b=0;b<config.nbins;++b)
              for (ilayer=0;ilayer<config.nlayer;++ilayer)
                for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                  {
                    //sum of concentrations at equilibrium
                    if (surrogate[i].time(b,ilayer,iphase)>=tequilibrium)
                      conc_equilibrium-=surrogate[i].Ap_layer_init(b,ilayer,iphase);
                    
                    //sum of Ap/Ag ratios + 1
                    if (surrogate[i].time(b,ilayer,iphase)<tequilibrium)
                      {
                        is_equilibrium=true;
                        sum+=surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase);
                      }
                  }	
            
            //compute Ag and Ap
            if (is_equilibrium)
              {			
                surrogate[i].Ag=conc_equilibrium/sum;
                for (b=0;b<config.nbins;++b)
                  for (ilayer=0;ilayer<config.nlayer;++ilayer)
                    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                      if (surrogate[i].time(b,ilayer,iphase)<tequilibrium)
                        surrogate[i].Ap_layer_init(b,ilayer,iphase)=
                          min(2.0*max(surrogate[i].Ap_layer_init(b,ilayer,iphase),1.0e-10),factor*surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)
                              /sum*conc_equilibrium
			      +(1.0-factor)*surrogate[i].Ap_layer_init(b,ilayer,iphase));
              }
            else
              {
                surrogate[i].Ag=surrogate[i].Atot;
                for (b=0;b<config.nbins;++b)
                  for (ilayer=0;ilayer<config.nlayer;++ilayer)
                    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                      surrogate[i].Ag-=surrogate[i].Ap_layer_init(b,ilayer,iphase);
              }
          }
      
      //compute organic mass concentrations
      for (b=0;b<config.nbins;++b)
        for (ilayer=0;ilayer<config.nlayer;++ilayer)
          {
            for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
              {
                MO(b,ilayer,iphase)=0.0;
                for (i=0;i<n;++i)
                  MO(b,ilayer,iphase)+=surrogate[i].Ap_layer_init(b,ilayer,iphase);
              }

            //To have organic masses at least equal to MOmin
            sum=0.0;
            for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
              sum+=MO(b,ilayer,iphase);
            
            if (sum>0.0)
              for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                MO(b,ilayer,iphase)=max(MO(b,ilayer,iphase),
                                        config.MOmin*config.Vlayer(ilayer)*MO(b,ilayer,iphase)/sum);
            else
              for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                if(iphase==0)
                  MO(b,ilayer,iphase)=config.MOmin*config.Vlayer(ilayer);
                else
                  MO(b,ilayer,iphase)=0.0;
          }

      sum=1.0;
      for (b=0;b<config.nbins;++b) //mass of water absorbed by organics in the aqueous phase
	{	  	 	  
	  if (surrogate[config.iH2O].hydrophobic and config.compute_organic)
	    for (ilayer=0;ilayer<config.nlayer;++ilayer)
	      for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		if (surrogate[config.iH2O].time(b,ilayer,iphase)<tequilibrium)
		  sum+=surrogate[config.iH2O].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase);                  
	}
  
      if (surrogate[config.iH2O].hydrophobic and config.compute_organic and config.hygroscopicity)
	for (b=0;b<config.nbins;++b)
	  for (ilayer=0;ilayer<config.nlayer;++ilayer)
	    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
	      if (surrogate[config.iH2O].time(b,ilayer,iphase)<tequilibrium)	      
		surrogate[config.iH2O].Ap_layer_init(b,ilayer,iphase)=factor*
		  surrogate[config.iH2O].Atot*surrogate[config.iH2O].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)/sum+
		  (1.0-factor)*surrogate[config.iH2O].Ap_layer_init(b,ilayer,iphase);

    }
  else
    {
      for (b=0;b<config.nbins;++b)
        for (ilayer=0;ilayer<config.nlayer;++ilayer)
          for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
            MO(b,ilayer,iphase)=MOinit(b,ilayer,iphase);
    }




}

void redistribution_ssh(model_config &config, vector<species>& surrogate,
			Array <double, 3> &MOinit,  Array <double, 3> &MO)
{
  //redistribute organic concentrations over the layer to ensure that:
  // MOlayer(ilayer)=alpha(ilayer)*MO
  // MOinit previous organic concentrations
  // MO new organic concentrations
  double MObins_init,MObins;
  int n=surrogate.size();
  int i,b,ilayer,ilayer2,iphase,jphase;
  double factor;
  double redistributed_mass;
  double sumphase,sumphase2;
  
  for (b=0;b<config.nbins;++b)
    {
      MObins=0.0;
      MObins_init=0.0;
      for (ilayer=0;ilayer<config.nlayer;++ilayer)
        for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
          {
            MObins+=MO(b,ilayer,iphase);
            MObins_init+=MOinit(b,ilayer,iphase);
          }
	  
      if (MObins > MObins_init) //condensation case
        for (ilayer=0;ilayer<config.nlayer-1;++ilayer)
          {	    
            sumphase=0.0;
            for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
              sumphase+=MO(b,ilayer,iphase);
	    			
            if (MObins*config.Vlayer(ilayer)-sumphase >0.0)
              {
                //MObins*config.Vlayer(ilayer)-sumphase: mass which has to be redistributed
                redistributed_mass=0.0; //mass whoch has been redistributed
                while (redistributed_mass < MObins*config.Vlayer(ilayer)-sumphase)
                  for (ilayer2=ilayer+1;ilayer2<config.nlayer;++ilayer2)
                    {
                      sumphase2=0.0;
                      for (iphase=0;iphase<config.nphase(b,ilayer2);++iphase)
                        sumphase2+=MO(b,ilayer2,iphase);
                      //fraction of the mass to be redistributed
                      if (sumphase2 == 0.0)
                        factor = 0.0;
                      else
                        factor=min((MObins*config.Vlayer(ilayer)-sumphase
                                    -redistributed_mass)/sumphase2,1.0);
                      redistributed_mass+=factor*sumphase2;

                      //redistritubed mass
                      for (iphase=0;iphase<config.nphase(b,ilayer2);++iphase)
                        {
                          jphase=min(iphase,config.nphase(b,ilayer)-1);
                          MO(b,ilayer,jphase)+=factor*MO(b,ilayer2,iphase);
                        }
					  
                      for (iphase=0;iphase<config.nphase(b,ilayer2);++iphase)
                        MO(b,ilayer2,iphase)-=factor*MO(b,ilayer2,iphase);

                      //redistribute concentrations of organic compounds
                      for (i=0;i<n;++i)
                        if (surrogate[i].hydrophobic)
                          for (iphase=0;iphase<config.nphase(b,ilayer2);++iphase)
                            {
                              jphase=min(iphase,config.nphase(b,ilayer)-1);
                              surrogate[i].Ap_layer_init(b,ilayer,jphase)+=factor*
                                surrogate[i].Ap_layer_init(b,ilayer2,iphase);
                              surrogate[i].Ap_layer_init(b,ilayer2,iphase)-=factor*
                                surrogate[i].Ap_layer_init(b,ilayer2,iphase);
                            }		     
                    }		
                sumphase+=redistributed_mass;
		
              }
            else if (MObins*config.Vlayer(ilayer)-sumphase <0.0)
              {
                if (sumphase == 0.0)
                  factor = 0.0;
                else
                  factor=(sumphase-MObins*config.Vlayer(ilayer))/sumphase;

                for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                  {		    
                    jphase=min(iphase,config.nphase(b,ilayer+1)-1);
                    MO(b,ilayer+1,jphase)+=factor*MO(b,ilayer,iphase);
                    MO(b,ilayer,iphase)-=factor*MO(b,ilayer,iphase);
                    for (i=0;i<n;++i)
                      if (surrogate[i].hydrophobic)
                        {
                          surrogate[i].Ap_layer_init(b,ilayer+1,jphase)+=factor*
                            surrogate[i].Ap_layer_init(b,ilayer,iphase);
                          surrogate[i].Ap_layer_init(b,ilayer,iphase)-=factor*
                            surrogate[i].Ap_layer_init(b,ilayer,iphase);
                        }
                  }
              }
          }
	  
      else if (MObins < MObins_init) //evaporation
        for (ilayer=config.nlayer-1;ilayer>0;ilayer--)
          {	    
            sumphase=0.0;
            for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
              sumphase+=MO(b,ilayer,iphase);
				    
            if (MObins*config.Vlayer(ilayer)-sumphase > 0.0)
              {
                //MObins*config.Vlayer(ilayer)-sumphase: mass which has to be redistributed
                redistributed_mass=0.0; //redistributed mass
                while (redistributed_mass < MObins*config.Vlayer(ilayer)-sumphase)
                  for (ilayer2=ilayer-1;ilayer2>=0;ilayer2--)
                    {
                      sumphase2=0.0;
                      for (iphase=0;iphase<config.nphase(b,ilayer2);++iphase)
                        sumphase2+=MO(b,ilayer2,iphase);
                      //fraction of the mass to be redistributed
                      if (sumphase2 == 0.0)
                        factor = 0.0;
                      else
                        factor=min((MObins*config.Vlayer(ilayer)-sumphase-redistributed_mass )
                                   /sumphase2,1.0);
                      redistributed_mass+=factor*sumphase2;

                      //redistrubed organic concentrations
                      for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                        {
                          jphase=min(iphase,config.nphase(b,ilayer)-1);
                          MO(b,ilayer,jphase)+=factor*MO(b,ilayer2,iphase);
                        }
					  
                      for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                        MO(b,ilayer2,iphase)-=factor*MO(b,ilayer2,iphase);

                      //redistribute concentrations of organic compounds
                      for (i=0;i<n;++i)
                        if (surrogate[i].hydrophobic)
                          for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                            {
                              jphase=min(iphase,config.nphase(b,ilayer)-1);
                              surrogate[i].Ap_layer_init(b,ilayer,jphase)+=factor*
                                surrogate[i].Ap_layer_init(b,ilayer2,iphase);
                              surrogate[i].Ap_layer_init(b,ilayer2,iphase)-=factor*
                                surrogate[i].Ap_layer_init(b,ilayer2,iphase);
                            }
                    }
                sumphase+=redistributed_mass;
              }
            else
              {
                if (sumphase == 0.0)
                  factor=0.0;
                else
                  factor=(sumphase-MObins*config.Vlayer(ilayer))/sumphase;
                for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                  {
                    jphase=min(iphase,config.nphase(b,ilayer-1)-1);
                    MO(b,ilayer-1,jphase)+=factor*MO(b,ilayer,iphase);
                    MO(b,ilayer,iphase)-=factor*MO(b,ilayer,iphase);
                    for (i=0;i<n;++i)
                      if (surrogate[i].hydrophobic)
                        {
                          surrogate[i].Ap_layer_init(b,ilayer-1,jphase)+=factor*
                            surrogate[i].Ap_layer_init(b,ilayer,iphase);
                          surrogate[i].Ap_layer_init(b,ilayer,iphase)-=factor*
                            surrogate[i].Ap_layer_init(b,ilayer,iphase);
                        }
                  }
              }
			
          }
    }
}



void dynamic_org_ssh(model_config &config, vector<species>& surrogate,
		     Array<double, 3> &MOinit,Array<double, 3> &MO,
		     Array<double, 1> &AQinit,
		     Array<double, 3> &MOW,
		     double &Temperature, 
		     double &DT2, double &tequilibrium,
		     bool compute_activity_coefficients)
{
  //compute the dynamic evolution of the organic phase concentrations with the roschem algorythm
  //FOR THE UNCOUPLED SYSTEM
  //compute_activity_coefficients: should the activity coefficients be recomputed?
  //tequilibrium: criterium to distinguish cases under equilibrium (characteristic time lower
  //                  than tequilibrium)
  int n=surrogate.size();
  int i,ilayer,iphase,b;
  double sumt=0.0;
  //int index_b,index_layer;
  double tiny=1.0e-26;
  double conc_available=0.0;
  double sum_rates;
  double gamma=1.7071;
  Array<double, 3> MOinit2;
  Array<double,1> Kp;
  MOinit2.resize(config.nbins,config.nlayer,config.max_number_of_phases);
  Kp.resize(config.max_number_of_phases);
  double sumconc,sumconc2;
  //double conc_equilibrium;
  
  //     gamma > .25     ->  A-stable  ( x(n+1)-x(n))*f(x(n)) > 0
  //     gamma < .25     ->  monotonic behaviour (x(n+1)-x*)(x(n)-x*) > 0
  //     gamma =  1+-sqrt(1/2) ->  L-stability
  
  //compute activity coefficients
  if (config.compute_organic)
    {
      if (compute_activity_coefficients)
        {
          activity_coefficients_dyn_org_ssh(config, surrogate, Temperature, MOW);
          compute_kp_org_ssh(config, surrogate, MOinit, Temperature, MOW);
        }

      //save initial concentrations in Ap0
      for (i=0;i<n;++i)
        if(surrogate[i].hydrophobic) 
          {
            surrogate[i].Ag0=surrogate[i].Ag;
            for (b=0;b<config.nbins;++b)		  
              for (ilayer=0;ilayer<config.nlayer;++ilayer)
                for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                  surrogate[i].Ap_layer_init0(b,ilayer,iphase)=
                    surrogate[i].Ap_layer_init(b,ilayer,iphase);
          }

      //compute first evaluation of kinetic rates
      flux_org_ssh(config, surrogate, MOinit, AQinit, tiny, 0);

      //compute the first evaluation of concentrations
      for (i=0;i<n;++i)
        if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophobic) 
          {
            conc_available=surrogate[i].Ag;
            sum_rates=0.0;
            for (b=0;b<config.nbins;++b)
              for (ilayer=0;ilayer<config.nlayer;++ilayer)
                for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                  if (surrogate[i].time(b,ilayer,iphase)<tequilibrium)
                    conc_available+=surrogate[i].Ap_layer_init(b,ilayer,iphase);
                  else
                    sum_rates+=surrogate[i].k1(b,ilayer,iphase,0);
		
            for (b=0;b<config.nbins;++b)
              for (ilayer=0;ilayer<config.nlayer;++ilayer)
                for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                  {
                    if (surrogate[i].time(b,ilayer,iphase)>=tequilibrium)
                      if (conc_available-surrogate[i].Ap_layer_init(b,ilayer,iphase)> tiny and
                          surrogate[i].k1(b,ilayer,iphase,0)>0.0)
                        surrogate[i].Jdn(b,ilayer,iphase,0)=-sum_rates/
                          (conc_available-surrogate[i].Ap_layer_init(b,ilayer,iphase));
                  }
		
            for (b=0;b<config.nbins;++b)		  
              for (ilayer=0;ilayer<config.nlayer;++ilayer)
                for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                  surrogate[i].k1(b,ilayer,iphase,0)=surrogate[i].k1(b,ilayer,iphase,0)*
                    DT2/(1.0-gamma*surrogate[i].Jdn(b,ilayer,iphase,0)*DT2);

            surrogate[i].Ag=surrogate[i].Atot;
            for (b=0;b<config.nbins;++b)		  
              for (ilayer=0;ilayer<config.nlayer;++ilayer)
                for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                  {
                    if (surrogate[i].time(b,ilayer,iphase)>=tequilibrium)
                      surrogate[i].Ap_layer_init(b,ilayer,iphase)=
                        max(0.01*surrogate[i].Ap_layer_init(b,ilayer,iphase),
                            surrogate[i].Ap_layer_init(b,ilayer,iphase)+
                            surrogate[i].k1(b,ilayer,iphase,0));
                    surrogate[i].Ag-=surrogate[i].Ap_layer_init(b,ilayer,iphase);
                  }
            surrogate[i].Ag=max(surrogate[i].Ag,0.01*surrogate[i].Ag0);
            surrogate[i].Ag1=surrogate[i].Ag;
          }

      //compute the first evaluation of the organic mass
      for (b=0;b<config.nbins;++b)		  
        for (ilayer=0;ilayer<config.nlayer;++ilayer)
          for (iphase=0;iphase<config.max_number_of_phases;++iphase)
            MOinit2(b,ilayer,iphase)=0.0;
  
      for (b=0;b<config.nbins;++b)		  
        for (ilayer=0;ilayer<config.nlayer;++ilayer)
          for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
            for (i=0;i<n;++i)
              MOinit2(b,ilayer,iphase)+=max(surrogate[i].Ap_layer_init(b,ilayer,iphase),
                                            config.MOmin*config.Vlayer(ilayer));

      //compute the new activity coefficients
      if (compute_activity_coefficients)
        {
          activity_coefficients_dyn_org_ssh(config, surrogate, Temperature, MOW);
          compute_kp_org_ssh(config, surrogate, MOinit2, Temperature, MOW);
        }
  
      //second estimation of kinetic rates
      flux_org_ssh(config, surrogate, MOinit2, AQinit, tiny, 1);

      //compute the new concentrations
      for (i=0;i<n;++i)
        if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophobic) 
          {
            conc_available=surrogate[i].Ag;
            sum_rates=0.0;
            for (b=0;b<config.nbins;++b)
              for (ilayer=0;ilayer<config.nlayer;++ilayer)
                for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                  if (surrogate[i].time(b,ilayer,iphase)<tequilibrium)
                    conc_available+=surrogate[i].Ap_layer_init(b,ilayer,iphase);
                  else
                    sum_rates+=surrogate[i].k1(b,ilayer,iphase,1);
		
            for (b=0;b<config.nbins;++b)		  
              for (ilayer=0;ilayer<config.nlayer;++ilayer)
                for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                  if (surrogate[i].time(b,ilayer,iphase)>=tequilibrium)
                    if (conc_available-surrogate[i].Ap_layer_init(b,ilayer,iphase)> tiny and
                        surrogate[i].k1(b,ilayer,iphase,1)>0.0)
                      surrogate[i].Jdn(b,ilayer,iphase,1)=-sum_rates/
                        (conc_available-surrogate[i].Ap_layer_init(b,ilayer,iphase));

            for (b=0;b<config.nbins;++b)		  
              for (ilayer=0;ilayer<config.nlayer;++ilayer)
                for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                  surrogate[i].k1(b,ilayer,iphase,1)=
                    (surrogate[i].k1(b,ilayer,iphase,1)*DT2
                     -gamma*DT2*(surrogate[i].Jdn(b,ilayer,iphase,0)+surrogate[i].Jdn(b,ilayer,iphase,1))
                     *surrogate[i].k1(b,ilayer,iphase,0))/
                    (1.0-gamma*surrogate[i].Jdn(b,ilayer,iphase,1)*DT2);

            surrogate[i].Ag=surrogate[i].Atot;
            for (b=0;b<config.nbins;++b)		  
              for (ilayer=0;ilayer<config.nlayer;++ilayer)
                for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                  {
                    if (surrogate[i].time(b,ilayer,iphase)>=tequilibrium)
                      surrogate[i].Ap_layer(b,ilayer,iphase)=
                        max(surrogate[i].Ap_layer_init0(b,ilayer,iphase)*0.01,
                            surrogate[i].Ap_layer_init0(b,ilayer,iphase)+0.5*
                            (surrogate[i].k1(b,ilayer,iphase,0)+surrogate[i].k1(b,ilayer,iphase,1)));
				
                    else
                      surrogate[i].Ap_layer(b,ilayer,iphase)=
                        surrogate[i].Ap_layer_init(b,ilayer,iphase);
                    surrogate[i].Ag-=surrogate[i].Ap_layer(b,ilayer,iphase);
                  }
          }
        else
          for (b=0;b<config.nbins;++b)		  
            for (ilayer=0;ilayer<config.nlayer;++ilayer)
              for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                surrogate[i].Ap_layer(b,ilayer,iphase)=0.0;

      //make sure that the sum of concentrations is not higher than the total concentration
      for (i=0;i<n;++i)
        if (surrogate[i].hydrophobic and surrogate[i].is_organic)
          {
            sumconc=0.0;
            sumconc2=0.0;
            for (b=0;b<config.nbins;++b)		  
              for (ilayer=0;ilayer<config.nlayer;++ilayer)
                for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                  if (surrogate[i].time(b,ilayer,iphase)>=tequilibrium)
                    sumconc+=surrogate[i].Ap_layer(b,ilayer,iphase);
                  else
                    sumconc2+=surrogate[i].Ap_layer(b,ilayer,iphase);
		
            if(sumconc>surrogate[i].Atot)
              {
                for (b=0;b<config.nbins;++b)		  
                  for (ilayer=0;ilayer<config.nlayer;++ilayer)
                    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                      if (surrogate[i].time(b,ilayer,iphase)>=tequilibrium)
                        surrogate[i].Ap_layer(b,ilayer,iphase)=
                          surrogate[i].Ap_layer(b,ilayer,iphase)/sumconc*
                          0.999*surrogate[i].Atot;
                      else
                        if (sumconc2 != 0.0) 
                          surrogate[i].Ap_layer(b,ilayer,iphase)=
                            surrogate[i].Ap_layer(b,ilayer,iphase)/sumconc2*
                            0.001*surrogate[i].Atot;
                surrogate[i].Ag=0.0;
              }
            else
              if(sumconc+sumconc2>surrogate[i].Atot)
                {
                  surrogate[i].Ag=0.0;
                  for (b=0;b<config.nbins;++b)		  
                    for (ilayer=0;ilayer<config.nlayer;++ilayer)
                      for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                        if (surrogate[i].time(b,ilayer,iphase)<tequilibrium)
                          surrogate[i].Ap_layer(b,ilayer,iphase)=
                            surrogate[i].Ap_layer(b,ilayer,iphase)/sumconc2*
                            (surrogate[i].Atot-sumconc);
                }
          }

      for (b=0;b<config.nbins;++b)
        for (ilayer=0;ilayer<config.nlayer;++ilayer)
          {
            for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
              {
                MO(b,ilayer,iphase)=0.0;
                for (i=0;i<n;++i)
                  MO(b,ilayer,iphase)+=surrogate[i].Ap_layer_init(b,ilayer,iphase);
              }
            sumt=0.0;
            for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
              sumt+=MO(b,ilayer,iphase);

            if (sumt>0.0)
              for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                MO(b,ilayer,iphase)=max(MO(b,ilayer,iphase),
                                        config.MOmin*config.Vlayer(ilayer)*MO(b,ilayer,iphase)/sumt);
            else
              for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                if(iphase==0)
                  MO(b,ilayer,iphase)=config.MOmin*config.Vlayer(ilayer);
                else
                  MO(b,ilayer,iphase)=0.0;
          }
    }
  else
    {
      for (b=0;b<config.nbins;++b)
        for (ilayer=0;ilayer<config.nlayer;++ilayer)
          for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
            {
              MO(b,ilayer,iphase)=MOinit(b,ilayer,iphase);
              for (i=0;i<n;++i)
                surrogate[i].Ap_layer(b,ilayer,iphase)=surrogate[i].Ap_layer_init(b,ilayer,iphase);
            }
    }
}
void dynamic_inorg_ssh(model_config &config, vector<species>& surrogate,
		       Array <double, 1> &conc_inorganic, 
		       Array <double, 1> &LWC, double &DT2, double &tequilibrium, int index)
{
  //compute the dynamic evolution of the inorganic compounds concentrations with the roschem algorythm
  //compute_activity_coefficients: should the activity coefficients be recomputed?
  //tequilibrium: criterium to distinguish cases under equilibrium (characteristic time lower
  //                  than tequilibrium)

  int n=surrogate.size();
  int i,b;
  //double sum=0.0;
  double tiny=1.0e-26;
  double conc_available=0.0;
  double sum_rates;
  double gamma=1.0+pow(2.,0.5)/2.;
  double redmax=0.01;

  for (b=0;b<config.nbins;++b)
    {
      LWC(b)=0.0;
      conc_inorganic(b)=surrogate[config.iHp].Aaq_bins_init(b)+surrogate[config.iOHm].Aaq_bins_init(b);
    }

  if (index==0) //first estimation of concentrations
    {      
      for (i=0;i<n;++i)
        if(surrogate[i].is_inorganic_precursor and surrogate[i].is_solid==false)
          {
            conc_available=surrogate[i].Ag;
            sum_rates=0.0;
            for (b=0;b<config.nbins;++b)
              {
                if (config.chemistry)
                  surrogate[i].k1_aq(b,0)+=surrogate[i].flux_chem_aq(b,0);
                
                if (surrogate[i].time_aq(b)<tequilibrium)
                  {
                    if (i==config.iH2SO4)
                      conc_available+=(surrogate[config.iHSO4m].Aaq_bins_init(b)/surrogate[config.iHSO4m].MM
                                       +surrogate[config.iSO4mm].Aaq_bins_init(b)/surrogate[config.iSO4mm].MM)*
                        surrogate[i].MM;
                    else if (i==config.iHNO3)
                      conc_available+=surrogate[config.iNO3m].Aaq_bins_init(b)/surrogate[config.iNO3m].MM*
                        surrogate[i].MM;
                    else if (i==config.iNH3)
                      conc_available+=surrogate[config.iNH4p].Aaq_bins_init(b)/surrogate[config.iNH4p].MM*
                        surrogate[i].MM;
                    else if (i==config.iHCl)
                      conc_available+=surrogate[config.iClm].Aaq_bins_init(b)/surrogate[config.iClm].MM*
                        surrogate[i].MM;
                  }
                else                
                  sum_rates+=surrogate[i].k1_aq(b,0); 
              }

            if (config.chemistry)
              surrogate[i].ktot1=surrogate[i].flux_chem_gas(0);
            else
              surrogate[i].ktot1=0.;
            
            for (b=0;b<config.nbins;++b)			  
              if (surrogate[i].time_aq(b)>=tequilibrium)
                {		 
                  double tmp=0.0;
                  if (i==config.iH2SO4)
                    {
                      tmp+=(surrogate[config.iHSO4m].Aaq_bins_init(b)/surrogate[config.iHSO4m].MM
                            +surrogate[config.iSO4mm].Aaq_bins_init(b)/surrogate[config.iSO4mm].MM)*
                        surrogate[i].MM;                    
                    }
                  else if (i==config.iHNO3)                    
                    tmp+=surrogate[config.iNO3m].Aaq_bins_init(b)/surrogate[config.iNO3m].MM*
                      surrogate[i].MM;                     
                  else if (i==config.iNH3)                    
                    tmp+=surrogate[config.iNH4p].Aaq_bins_init(b)/surrogate[config.iNH4p].MM*
                      surrogate[i].MM;                                       
                  else if (i==config.iHCl)                    
                    tmp+=surrogate[config.iClm].Aaq_bins_init(b)/surrogate[config.iClm].MM*
                      surrogate[i].MM;                                         

                  if (conc_available-tmp> tiny and surrogate[i].k1_aq(b,0)>0.0)
                    surrogate[i].Jdn_aq(b,0)=-sum_rates/(conc_available-tmp);
                }

            if (config.chemistry)
              for (b=0;b<config.nbins;++b)			             
                {		                  
                  if (i==config.iH2SO4)
                    {                     
                      surrogate[i].ktot1+=surrogate[config.iHSO4m].flux_chem_aq(b,0)/surrogate[config.iHSO4m].MM/surrogate[i].MM;
                      surrogate[i].ktot1+=surrogate[config.iSO4mm].flux_chem_aq(b,0)/surrogate[config.iSO4mm].MM/surrogate[i].MM;
                    }
                  else if (i==config.iHNO3)
                    surrogate[i].ktot1+=surrogate[config.iNO3m].flux_chem_aq(b,0)/surrogate[config.iNO3m].MM/surrogate[i].MM;                    
                  else if (i==config.iNH3)           
                    surrogate[i].ktot1+=surrogate[config.iNH4p].flux_chem_aq(b,0)/surrogate[config.iNH4p].MM/surrogate[i].MM;                    
                  else if (i==config.iHCl)            
                    surrogate[i].ktot1+=surrogate[config.iClm].flux_chem_aq(b,0)/surrogate[config.iClm].MM/surrogate[i].MM;                    
                }
                                 
            for (b=0;b<config.nbins;++b)
              if (surrogate[i].time_aq(b)>=tequilibrium)                     
                surrogate[i].k1_aq(b,0)=surrogate[i].k1_aq(b,0)*
                  DT2/(1.0-gamma*surrogate[i].Jdn_aq(b,0)*DT2);   

            if (config.chemistry)              
              surrogate[i].Atot=surrogate[i].Atot+surrogate[i].ktot1;                                          
            surrogate[i].Ag=surrogate[i].Atot;
            
            for (b=0;b<config.nbins;++b)
              {
                if (i==config.iH2SO4)
                  {
                    if (surrogate[i].time_aq(b)>=tequilibrium)
                      {			
                        double f=0.;
			if (surrogate[config.iHSO4m].Aaq_bins_init(b)/surrogate[config.iHSO4m].MM+surrogate[config.iSO4mm].Aaq_bins_init(b)/surrogate[config.iSO4mm].MM>0.0)
			  f=(surrogate[config.iHSO4m].Aaq_bins_init(b)/surrogate[config.iHSO4m].MM)/
			    (surrogate[config.iHSO4m].Aaq_bins_init(b)/surrogate[config.iHSO4m].MM+surrogate[config.iSO4mm].Aaq_bins_init(b)/surrogate[config.iSO4mm].MM);
                        surrogate[config.iHSO4m].Aaq_bins_init(b)=
                          max(redmax*surrogate[config.iHSO4m].Aaq_bins_init(b),
                              surrogate[config.iHSO4m].Aaq_bins_init(b)+f*surrogate[i].k1_aq(b,0)*surrogate[config.iHSO4m].MM/surrogate[i].MM+surrogate[config.iHSO4m].flux_chem_aq(b,0));
                        surrogate[config.iSO4mm].Aaq_bins_init(b)=
                          max(redmax*surrogate[config.iSO4mm].Aaq_bins_init(b),
                              surrogate[config.iSO4mm].Aaq_bins_init(b)+(1.0-f)*surrogate[i].k1_aq(b,0)*surrogate[config.iSO4mm].MM/surrogate[i].MM+surrogate[config.iSO4mm].flux_chem_aq(b,0));
                      }
                    else
                      {
                        surrogate[config.iHSO4m].Aaq_bins_init(b)=
                          max(redmax*surrogate[config.iHSO4m].Aaq_bins_init(b),
                              surrogate[config.iHSO4m].Aaq_bins_init(b)+surrogate[config.iHSO4m].flux_chem_aq(b,0));
                        surrogate[config.iSO4mm].Aaq_bins_init(b)=
                          max(redmax*surrogate[config.iSO4mm].Aaq_bins_init(b),
                              surrogate[config.iSO4mm].Aaq_bins_init(b)+surrogate[config.iSO4mm].flux_chem_aq(b,0));
                      }
                    
                    surrogate[i].Ag-=surrogate[config.iHSO4m].Aaq_bins_init(b)/surrogate[config.iHSO4m].MM*surrogate[i].MM
                      +surrogate[config.iSO4mm].Aaq_bins_init(b)/surrogate[config.iSO4mm].MM*surrogate[i].MM;
                    conc_inorganic(b)+=surrogate[config.iHSO4m].Aaq_bins_init(b)+surrogate[config.iSO4mm].Aaq_bins_init(b);

                  }

                if (i==config.iNH3)
                  {
                    if (surrogate[i].time_aq(b)>=tequilibrium)
                      surrogate[config.iNH4p].Aaq_bins_init(b)=
                        max(redmax*surrogate[config.iNH4p].Aaq_bins_init(b),
                            surrogate[config.iNH4p].Aaq_bins_init(b)+surrogate[i].k1_aq(b,0)*surrogate[config.iNH4p].MM/surrogate[i].MM+surrogate[config.iNH4p].flux_chem_aq(b,0));
                    else
                      surrogate[config.iNH4p].Aaq_bins_init(b)=
                        max(redmax*surrogate[config.iNH4p].Aaq_bins_init(b),
                            surrogate[config.iNH4p].Aaq_bins_init(b)+surrogate[config.iNH4p].flux_chem_aq(b,0));

                    surrogate[i].Ag-=surrogate[config.iNH4p].Aaq_bins_init(b)/surrogate[config.iNH4p].MM*surrogate[i].MM;
                    conc_inorganic(b)+=surrogate[config.iNH4p].Aaq_bins_init(b);
                  }
                
                if (i==config.iHNO3)
                  {                    
                    if (surrogate[i].time_aq(b)>=tequilibrium)
		      {
			surrogate[config.iNO3m].Aaq_bins_init(b)=
			  max(redmax*surrogate[config.iNO3m].Aaq_bins_init(b),
			      surrogate[config.iNO3m].Aaq_bins_init(b)+surrogate[i].k1_aq(b,0)*surrogate[config.iNO3m].MM/surrogate[i].MM+surrogate[config.iNO3m].flux_chem_aq(b,0));			
		      }
                    else
                      surrogate[config.iNO3m].Aaq_bins_init(b)=
                        max(redmax*surrogate[config.iNO3m].Aaq_bins_init(b),
                            surrogate[config.iNO3m].Aaq_bins_init(b)+surrogate[config.iNO3m].flux_chem_aq(b,0));

                    surrogate[i].Ag-=surrogate[config.iNO3m].Aaq_bins_init(b)/surrogate[config.iNO3m].MM*surrogate[i].MM;
                    conc_inorganic(b)+=surrogate[config.iNO3m].Aaq_bins_init(b);
                  }

                if (i==config.iHCl)
                  {
                    if (surrogate[i].time_aq(b)>=tequilibrium)
                      surrogate[config.iClm].Aaq_bins_init(b)=
                        max(redmax*surrogate[config.iClm].Aaq_bins_init(b),
                            surrogate[config.iClm].Aaq_bins_init(b)+surrogate[i].k1_aq(b,0)*surrogate[config.iClm].MM/surrogate[i].MM+surrogate[config.iClm].flux_chem_aq(b,0));
                    else
                      surrogate[config.iClm].Aaq_bins_init(b)=
                        max(redmax*surrogate[config.iClm].Aaq_bins_init(b),
                            surrogate[config.iClm].Aaq_bins_init(b)+surrogate[config.iClm].flux_chem_aq(b,0));

                    surrogate[i].Ag-=surrogate[config.iClm].Aaq_bins_init(b)/surrogate[config.iClm].MM*surrogate[i].MM;
                    conc_inorganic(b)+=surrogate[config.iClm].Aaq_bins_init(b);
                  } 
              }            
          }  
      for (i=0;i<n;i++)
	if (surrogate[i].is_inorganic_precursor and surrogate[i].is_solid==false)
	  { 
	    if (surrogate[i].Ag<0.0)	    
	      surrogate[i].Ag=0.0;
	    surrogate[i].Ag1=surrogate[i].Ag;
	  }
    
      for (b=0;b<config.nbins;++b)
        {
          conc_inorganic(b)=0.0;
          for (i=0;i<n;i++)
            if (surrogate[i].is_inorganic_precursor==false and surrogate[i].is_organic==false and i!=config.iH2O)
              conc_inorganic(b)+=surrogate[i].Aaq_bins_init(b);
        } 
    }
  else  //second estimation
    {      
      for (i=0;i<n;++i)
        if(surrogate[i].is_inorganic_precursor and surrogate[i].is_solid==false)
          {
            conc_available=surrogate[i].Ag;
            sum_rates=0.0;
            for (b=0;b<config.nbins;++b)			  
              {
                if (config.chemistry)
                  surrogate[i].k1_aq(b,1)+=surrogate[i].flux_chem_aq(b,1)/DT2;
                
                if (surrogate[i].time_aq(b)<tequilibrium)
                  {
                    if (i==config.iH2SO4)
                      conc_available+=0.0; 
                    else if (i==config.iHNO3)
                      conc_available+=surrogate[config.iNO3m].Aaq_bins_init(b)/surrogate[config.iNO3m].MM*
                        surrogate[i].MM;
                    else if (i==config.iNH3)
                      conc_available+=surrogate[config.iNH4p].Aaq_bins_init(b)/surrogate[config.iNH4p].MM*
                        surrogate[i].MM;
                    else if (i==config.iHCl)
                      conc_available+=surrogate[config.iClm].Aaq_bins_init(b)/surrogate[config.iClm].MM*
                        surrogate[i].MM;
                  }
                else
                  sum_rates+=surrogate[i].k1_aq(b,1);
              }

            if (config.chemistry)
              surrogate[i].ktot2=surrogate[i].flux_chem_gas(1);
            else
              surrogate[i].ktot2=0.0;
            
            for (b=0;b<config.nbins;++b)			  
              if (surrogate[i].time_aq(b)>=tequilibrium)
                {
                  double tmp=0.0;
                  if (i==config.iH2SO4)
                    tmp+=(surrogate[config.iHSO4m].Aaq_bins_init(b)/surrogate[config.iHSO4m].MM
                          +surrogate[config.iSO4mm].Aaq_bins_init(b)/surrogate[config.iSO4mm].MM)*
                      surrogate[i].MM;
                  else if (i==config.iHNO3)
                    tmp+=surrogate[config.iNO3m].Aaq_bins_init(b)/surrogate[config.iNO3m].MM*
                      surrogate[i].MM;
                  else if (i==config.iNH3)
                    tmp+=surrogate[config.iNH4p].Aaq_bins_init(b)/surrogate[config.iNH4p].MM*
                      surrogate[i].MM;
                  else if (i==config.iHCl)                    
                    tmp+=surrogate[config.iClm].Aaq_bins_init(b)/surrogate[config.iClm].MM*
                      surrogate[i].MM;                                        

                  if (conc_available-tmp> tiny and surrogate[i].k1_aq(b,1)>0.0)
                    surrogate[i].Jdn_aq(b,1)=-sum_rates/(conc_available-tmp);
                }

            if (config.chemistry)
              for (b=0;b<config.nbins;++b)			             
                {		                  
                  if (i==config.iH2SO4)
                    {                     
                      surrogate[i].ktot2+=surrogate[config.iHSO4m].flux_chem_aq(b,1)/surrogate[config.iHSO4m].MM/surrogate[i].MM;
                      surrogate[i].ktot2+=surrogate[config.iSO4mm].flux_chem_aq(b,1)/surrogate[config.iSO4mm].MM/surrogate[i].MM;
                    }
                  else if (i==config.iHNO3)
                    surrogate[i].ktot2+=surrogate[config.iNO3m].flux_chem_aq(b,1)/surrogate[config.iNO3m].MM/surrogate[i].MM;                    
                  else if (i==config.iNH3)           
                    surrogate[i].ktot2+=surrogate[config.iNH4p].flux_chem_aq(b,1)/surrogate[config.iNH4p].MM/surrogate[i].MM;                    
                  else if (i==config.iHCl)            
                    surrogate[i].ktot2+=surrogate[config.iClm].flux_chem_aq(b,1)/surrogate[config.iClm].MM/surrogate[i].MM;
                }
            
                            
            for (b=0;b<config.nbins;++b)              
              surrogate[i].k1_aq(b,1)=
                (surrogate[i].k1_aq(b,1)*DT2-gamma*DT2*(surrogate[i].Jdn_aq(b,0)+surrogate[i].Jdn_aq(b,1))
                 *surrogate[i].k1_aq(b,0))/(1.0-gamma*surrogate[i].Jdn_aq(b,1)*DT2);              

            if (config.chemistry)
              {
                surrogate[i].Atot=surrogate[i].Atot0+0.5*(surrogate[i].ktot1+surrogate[i].ktot2);                
              }
            surrogate[i].Ag=surrogate[i].Atot;
            /*
	      else
              for (b=0;b<config.nbins;++b)
              surrogate[i].k1_aq(b,1)=0.;                        */

            for (b=0;b<config.nbins;++b)
              if (i==config.iH2SO4)
                {
                  if (surrogate[i].time_aq(b)>=tequilibrium)
                    {
		      double f=0;
		      if(surrogate[config.iHSO4m].Aaq_bins_init(b)/surrogate[config.iHSO4m].MM+surrogate[config.iSO4mm].Aaq_bins_init(b)/surrogate[config.iSO4mm].MM>0.0)
			f=(surrogate[config.iHSO4m].Aaq_bins_init(b)/surrogate[config.iHSO4m].MM)/
			  (surrogate[config.iHSO4m].Aaq_bins_init(b)/surrogate[config.iHSO4m].MM+surrogate[config.iSO4mm].Aaq_bins_init(b)/surrogate[config.iSO4mm].MM);
                      surrogate[config.iHSO4m].Aaq_bins(b)=
                        max(redmax*surrogate[config.iHSO4m].Aaq_bins_init0(b),
                            surrogate[config.iHSO4m].Aaq_bins_init0(b)+f*0.5*(surrogate[i].k1_aq(b,0)+surrogate[i].k1_aq(b,1))*surrogate[config.iHSO4m].MM/surrogate[i].MM
                            +0.5*(surrogate[config.iHSO4m].flux_chem_aq(b,0)+surrogate[config.iHSO4m].flux_chem_aq(b,1)));
                      surrogate[config.iSO4mm].Aaq_bins(b)=
                        max(redmax*surrogate[config.iSO4mm].Aaq_bins_init0(b),
                            surrogate[config.iSO4mm].Aaq_bins_init0(b)+(1.0-f)*0.5*(surrogate[i].k1_aq(b,0)+surrogate[i].k1_aq(b,1))*surrogate[config.iSO4mm].MM/surrogate[i].MM
                            +0.5*(surrogate[config.iSO4mm].flux_chem_aq(b,0)+surrogate[config.iSO4mm].flux_chem_aq(b,1)));
                    }
                  else
                    {    
                      surrogate[config.iHSO4m].Aaq_bins(b)=
                        max(redmax*surrogate[config.iHSO4m].Aaq_bins_init0(b),
                            surrogate[config.iHSO4m].Aaq_bins_init0(b)+0.5*(surrogate[config.iHSO4m].flux_chem_aq(b,0)+surrogate[config.iHSO4m].flux_chem_aq(b,1)));
                      surrogate[config.iSO4mm].Aaq_bins(b)=
                        max(redmax*surrogate[config.iSO4mm].Aaq_bins_init0(b),
                            surrogate[config.iSO4mm].Aaq_bins_init0(b)+0.5*(surrogate[config.iSO4mm].flux_chem_aq(b,0)+surrogate[config.iSO4mm].flux_chem_aq(b,1)));
                    }
                    
                  surrogate[i].Ag-=surrogate[config.iHSO4m].Aaq_bins(b)/surrogate[config.iHSO4m].MM*surrogate[i].MM
                    +surrogate[config.iSO4mm].Aaq_bins(b)/surrogate[config.iSO4mm].MM*surrogate[i].MM;

                  conc_inorganic(b)+=surrogate[config.iHSO4m].Aaq_bins(b)+surrogate[config.iSO4mm].Aaq_bins(b);
                }
              else if (i==config.iNH3)
                {
                  if (surrogate[i].time_aq(b)>=tequilibrium)
                    surrogate[config.iNH4p].Aaq_bins(b)=
                      max(redmax*surrogate[config.iNH4p].Aaq_bins_init0(b),
                          surrogate[config.iNH4p].Aaq_bins_init0(b)+0.5*(surrogate[i].k1_aq(b,0)+surrogate[i].k1_aq(b,1))*surrogate[config.iNH4p].MM/surrogate[i].MM
                          +0.5*(surrogate[config.iNH4p].flux_chem_aq(b,0)+surrogate[config.iNH4p].flux_chem_aq(b,1)));
                  else
                    surrogate[config.iNH4p].Aaq_bins(b)=
                      max(redmax*surrogate[config.iNH4p].Aaq_bins_init0(b),
                          surrogate[config.iNH4p].Aaq_bins_init0(b)+0.5*(surrogate[config.iNH4p].flux_chem_aq(b,0)+surrogate[config.iNH4p].flux_chem_aq(b,1)));
                    
                  surrogate[i].Ag-=surrogate[config.iNH4p].Aaq_bins(b)/surrogate[config.iNH4p].MM*surrogate[i].MM;
                  conc_inorganic(b)+=surrogate[config.iNH4p].Aaq_bins(b);
                }
              else if (i==config.iHNO3)
                {
                  if (surrogate[i].time_aq(b)>=tequilibrium)
                    surrogate[config.iNO3m].Aaq_bins(b)=
                      max(redmax*surrogate[config.iNO3m].Aaq_bins_init0(b),
                          surrogate[config.iNO3m].Aaq_bins_init0(b)+0.5*(surrogate[i].k1_aq(b,0)+surrogate[i].k1_aq(b,1))*surrogate[config.iNO3m].MM/surrogate[i].MM
                          +0.5*(surrogate[config.iNO3m].flux_chem_aq(b,0)+surrogate[config.iNO3m].flux_chem_aq(b,1)));
                  else
                    surrogate[config.iNO3m].Aaq_bins(b)=
                      max(redmax*surrogate[config.iNO3m].Aaq_bins_init0(b),
                          surrogate[config.iNO3m].Aaq_bins_init0(b)+0.5*(surrogate[config.iNO3m].flux_chem_aq(b,0)+surrogate[config.iNO3m].flux_chem_aq(b,1)));

                  surrogate[i].Ag-=surrogate[config.iNO3m].Aaq_bins(b)/surrogate[config.iNO3m].MM*surrogate[i].MM;
                  conc_inorganic(b)+=surrogate[config.iNO3m].Aaq_bins(b);
                }

              else if (i==config.iHCl)
                {
                  if (surrogate[i].time_aq(b)>=tequilibrium)
                    surrogate[config.iClm].Aaq_bins(b)=
                      max(redmax*surrogate[config.iClm].Aaq_bins_init0(b),
                          surrogate[config.iClm].Aaq_bins_init0(b)+0.5*(surrogate[i].k1_aq(b,0)+surrogate[i].k1_aq(b,1))*surrogate[config.iClm].MM/surrogate[i].MM
                          +0.5*(surrogate[config.iClm].flux_chem_aq(b,0)+surrogate[config.iClm].flux_chem_aq(b,1)));
                  else
                    surrogate[config.iClm].Aaq_bins(b)=
                      max(redmax*surrogate[config.iClm].Aaq_bins_init0(b),
                          surrogate[config.iClm].Aaq_bins_init0(b)+0.5*(surrogate[config.iClm].flux_chem_aq(b,0)+surrogate[config.iClm].flux_chem_aq(b,1)));

                  surrogate[i].Ag-=surrogate[config.iClm].Aaq_bins(b)/surrogate[config.iClm].MM*surrogate[i].MM;
                  conc_inorganic(b)+=surrogate[config.iClm].Aaq_bins(b);
                } 
          }

      for (i=0;i<n;i++)
	if (surrogate[i].is_inorganic_precursor and surrogate[i].is_solid==false) 
	  if (surrogate[i].Ag<0.0)	    
	    surrogate[i].Ag=0.0;		 	   	    

      for (b=0;b<config.nbins;++b)
        {
          conc_inorganic(b)=0.0;
          for (i=0;i<n;i++)
            if (surrogate[i].is_inorganic_precursor==false and surrogate[i].is_organic==false and i!=config.iH2O)
              conc_inorganic(b)+=surrogate[i].Aaq_bins(b);
        } 
    

    }
}

void dynamic_tot_ssh(model_config &config, vector<species>& surrogate,
		     Array<double, 3> &MOinit,Array<double, 3> &MO,
		     Array<double, 3> &MOW,
		     Array<double, 1> &AQinit,Array<double, 1> &AQ,
		     Array <double, 1> &conc_inorganic,
		     Array <double, 1> &ionic,Array <double, 1> &ionic_organic,
		     Array <double, 1> &organion,Array <double, 1> &chp,
		     Array <double, 1> &chp1,Array <double, 1> &chp0,
		     Array <double, 1> &LWC,
		     Array<double, 1> &MMaq,
		     double &Temperature, double &RH, double &DT2, double &tequilibrium,
		     bool compute_activity_coefficients)
{
  //compute the dynamic evolution of the organic-phase and the aqueous-phase concentrations with
  //         the roschem algorythm
  //FOR THE COUPLED SYSTEM
  //compute_activity_coefficients: should the activity coefficients be recomputed?
  //tequilibrium: criterium to distinguish cases under equilibrium (characteristic time lower
  //                  than tequilibrium)
  int n=surrogate.size();
  int i,ilayer,iphase,b;
  //int index_b,index_layer;
  //double sumt=0.0;
  double LWCtot=0.0;
  double tiny=1.0e-26;
  double conc_available=0.0;
  double sum_rates;
  double gamma=1.0+pow(2.,0.5)/2.; //1.7071;
  //double Kaq;
  Array<double,1> Kp;
  Kp.resize(config.max_number_of_phases);
  Array<double, 3> MOinit2;
  Array<double, 1> AQinit2;
  Array <double, 1> chp2;
  //double conc_equilibrium;  
  double redmax=0.01;
  chp2.resize(config.nbins);
  AQinit2.resize(config.nbins);
  MOinit2.resize(config.nbins,config.nlayer,config.max_number_of_phases);
  
  //     gamma > .25     ->  A-stable  ( x(n+1)-x(n))*f(x(n)) > 0
  //     gamma < .25     ->  monotonic behaviour (x(n+1)-x*)(x(n)-x*) > 0
  //     gamma =  1+-sqrt(1/2) ->  L-stability

  LWCtot=sum(LWC);  
  chp0=chp;
  double error_chp=1;
  int iter=0;
  double conc_org;
  double sum_rates_gas;
  double sumt;
  double sumconc,sumconc2;
  double facph=1;
  
  //save initial concentrations 
  for (i=0;i<n;++i)
    {
      surrogate[i].Ag0=surrogate[i].Ag;
      surrogate[i].Ag1=surrogate[i].Ag;
      surrogate[i].Atot0=surrogate[i].Atot;
      surrogate[i].Ap_layer_init0=surrogate[i].Ap_layer_init;
      surrogate[i].Aaq_bins_init0=surrogate[i].Aaq_bins_init;      
    }

  while (iter<50 and error_chp>1.0e-3)
    {
      error_chp=0.0;
      iter++;
      if (config.compute_inorganic)	
	{
	  compute_ph_dyn_ssh(config, surrogate, Temperature, chp, AQinit, ionic, MMaq, LWC);
	  for (b=0;b<config.nbins;++b)
	    {
              compute_conc_org_bins_ssh(config, surrogate, Temperature,  MMaq, AQinit, chp, ionic, LWC, b, conc_org);              
              
	      //config.rho_aqueous=config.AQrho(b);
	      compute_ionic_strenght2_ssh(config, surrogate, Temperature, AQinit(b), conc_inorganic(b), ionic(b), chp(b),
					  organion(b), ionic_organic(b), conc_org, 1.0, config.Ke/surrogate[config.iHp].gamma_aq_bins(b)/surrogate[config.iOHm].gamma_aq_bins(b));
	      error_chp=max(error_chp,abs(chp0(b)-chp(b))/chp0(b));
	      chp0(b)=chp(b);	      
	    }      
	}
	  
      if (compute_activity_coefficients)
	{
	  if (config.compute_organic)
	    {
	      activity_coefficients_dyn_org_ssh(config, surrogate, Temperature, MOW);
	      compute_kp_org_ssh(config, surrogate, MOinit, Temperature, MOW);
	    }
	  
	  if (LWCtot>config.LWClimit)
	    {
	      activity_coefficients_dyn_aq_ssh(config, surrogate, Temperature,AQinit,
					       MOinit,conc_inorganic, ionic, ionic_organic,
					       organion,chp,LWC,MMaq, 1.0, 0., 0);
	      compute_kp_aq_ssh(config, surrogate, Temperature, ionic, chp, AQinit, MMaq, 0, config.nbins);   
	    }
	}
    }

  //compute kinetic rates
  flux_org_ssh(config, surrogate, MOinit, AQinit, tiny, 0);
  if (LWCtot>config.LWClimit)
    {
      flux_aq_ssh(config, surrogate, AQinit, LWC, MMaq, chp, ionic, MOinit, tiny, Temperature, 0);      
      //if (config.compute_inorganic and facph>0)
      //  correct_flux_ph_ssh(config, surrogate, Temperature, AQinit, MOinit, chp, chp2, MMaq, ionic, LWC, tiny, DT2/facph, 0);
    }
  
  if (config.chemistry)
    compute_flux_chem_ssh(config,surrogate,MOinit,MOW,AQinit,LWC,MMaq,chp,ionic,DT2,tiny,Temperature,0, RH);

  //compute the first evaluation of concentrations
  for (i=0;i<n;++i)
    if((config.compute_organic and surrogate[i].is_organic) or i==config.iH2O) 
      {
	conc_available=surrogate[i].Ag;
	sum_rates=0.0;
        sum_rates_gas=surrogate[i].flux_chem_gas(0)/DT2;
        if (config.chemistry==false)
          {        
            if (config.compute_organic)
              if (surrogate[i].hydrophobic)
                for (b=0;b<config.nbins;++b)
                  for (ilayer=0;ilayer<config.nlayer;++ilayer)
                    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)				
                      {
                        //surrogate[i].k1(b,ilayer,iphase,0)+=surrogate[i].flux_chem(b,ilayer,iphase,0)/DT2;
                        if (surrogate[i].time(b,ilayer,iphase)<tequilibrium)
                          conc_available+=surrogate[i].Ap_layer_init(b,ilayer,iphase);
                        else
                          {
                            sum_rates+=surrogate[i].k1(b,ilayer,iphase,0);
                            sum_rates_gas-=(surrogate[i].k1(b,ilayer,iphase,0)-surrogate[i].flux_chem(b,ilayer,iphase,0)/DT2);
                          }
                      }
	
            if (surrogate[i].hydrophilic and LWCtot>config.LWClimit and i!=config.iH2O)
              for (b=0;b<config.nbins;++b)
                {
                  //surrogate[i].k1_aq(b,0)+=surrogate[i].flux_chem_aq(b,0)/DT2;
                  if (surrogate[i].time_aq(b)<tequilibrium)
                    conc_available+=surrogate[i].Aaq_bins_init(b);
                  else
                    {
                      sum_rates+=surrogate[i].k1_aq(b,0);
                      sum_rates_gas-=surrogate[i].k1_aq(b,0)-surrogate[i].flux_chem_aq(b,0)/DT2;
                    }
                }
	
            if (config.compute_organic)
              if (surrogate[i].hydrophobic)
                for (b=0;b<config.nbins;++b)
                  for (ilayer=0;ilayer<config.nlayer;++ilayer)
                    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                      if (surrogate[i].time(b,ilayer,iphase)>=tequilibrium)
                        if (conc_available-surrogate[i].Ap_layer_init(b,ilayer,iphase)> tiny and
                            surrogate[i].k1(b,ilayer,iphase,0)>0.0)
                          surrogate[i].Jdn(b,ilayer,iphase,0)=-sum_rates/
                            (conc_available-surrogate[i].Ap_layer_init(b,ilayer,iphase));
	
            if (surrogate[i].hydrophilic and LWCtot>config.LWClimit) //and i!=config.iH2O)
              for (b=0;b<config.nbins;++b)			  
                if (surrogate[i].time_aq(b)>=tequilibrium)
                  if (conc_available-surrogate[i].Aaq_bins_init(b)> tiny and
                      surrogate[i].k1_aq(b,0)>0.0)
                    surrogate[i].Jdn_aq(b,0)=-sum_rates/
                      (conc_available-surrogate[i].Aaq_bins_init(b));    
          }

        surrogate[i].ktot1=surrogate[i].flux_chem_gas(0); // /(1.0-gamma*surrogate[i].Jdn_gas*DT2);
	if (config.compute_organic)
	  if (surrogate[i].hydrophobic)
	    for (b=0;b<config.nbins;++b)		  
	      for (ilayer=0;ilayer<config.nlayer;++ilayer)
		for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                  {
                    surrogate[i].k1(b,ilayer,iphase,0)=surrogate[i].k1(b,ilayer,iphase,0)*DT2/(1.0-gamma*surrogate[i].Jdn(b,ilayer,iphase,0)*DT2)+surrogate[i].flux_chem(b,ilayer,iphase,0);
                    //surrogate[i].flux_chem(b,ilayer,iphase,0)=surrogate[i].flux_chem(b,ilayer,iphase,0)/(1.0-gamma*surrogate[i].Jdn(b,ilayer,iphase,0);
                    surrogate[i].ktot1+=surrogate[i].flux_chem(b,ilayer,iphase,0); 
                  }
	
	if (surrogate[i].hydrophilic and LWCtot>config.LWClimit) //and i!=config.iH2O)
	  for (b=0;b<config.nbins;++b)
            {
              surrogate[i].k1_aq(b,0)=surrogate[i].k1_aq(b,0)*DT2/(1.0-gamma*surrogate[i].Jdn_aq(b,0)*DT2)+surrogate[i].flux_chem_aq(b,0);
              //surrogate[i].flux_chem_aq(b,0)=surrogate[i].flux_chem_aq(b,0)/(1.0-gamma*surrogate[i].Jdn_aq(b,0)*DT2);
              surrogate[i].ktot1+=surrogate[i].flux_chem_aq(b,0); 
            }                      

        surrogate[i].Atot=max(redmax*surrogate[i].Atot,surrogate[i].Atot0+surrogate[i].ktot1);	
	
	surrogate[i].Ag=surrogate[i].Atot;
	if (config.compute_organic)
	  if (surrogate[i].hydrophobic)
	    for (b=0;b<config.nbins;++b)		  
	      for (ilayer=0;ilayer<config.nlayer;++ilayer)
		for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		  {
		    if (surrogate[i].time(b,ilayer,iphase)>=tequilibrium)
		      surrogate[i].Ap_layer_init(b,ilayer,iphase)=
			max(redmax*surrogate[i].Ap_layer_init(b,ilayer,iphase),
			    surrogate[i].Ap_layer_init(b,ilayer,iphase)+
			    surrogate[i].k1(b,ilayer,iphase,0));
                    else
                      surrogate[i].Ap_layer_init(b,ilayer,iphase)=
			max(redmax*surrogate[i].Ap_layer_init(b,ilayer,iphase),
			    surrogate[i].Ap_layer_init(b,ilayer,iphase)+
			    surrogate[i].flux_chem(b,ilayer,iphase,0));
		    surrogate[i].Ag-=surrogate[i].Ap_layer_init(b,ilayer,iphase);
		  }
	
	if (surrogate[i].hydrophilic) // and i!=config.iH2O)
	  if(LWCtot>config.LWClimit)
	    for (b=0;b<config.nbins;++b)
	      {
		if (surrogate[i].time_aq(b)>=tequilibrium)
		  surrogate[i].Aaq_bins_init(b)=
		    max(redmax*surrogate[i].Aaq_bins_init(b),
			surrogate[i].Aaq_bins_init(b)+surrogate[i].k1_aq(b,0));
                else
                  surrogate[i].Aaq_bins_init(b)=
		    max(redmax*surrogate[i].Aaq_bins_init(b),
			surrogate[i].Aaq_bins_init(b)+surrogate[i].flux_chem_aq(b,0));
		
		surrogate[i].Ag-=surrogate[i].Aaq_bins_init(b);
	      }
	  else
	    for (b=0;b<config.nbins;++b)
	      surrogate[i].Aaq_bins_init(b)=0.0;
	else
	  for (b=0;b<config.nbins;++b)
	    surrogate[i].Aaq_bins_init(b)=surrogate[i].Aaq_bins_init0(b);
	
	surrogate[i].Ag=max(surrogate[i].Ag,redmax*surrogate[i].Ag0);
        
      }
    
  if (config.compute_inorganic)
    dynamic_inorg_ssh(config, surrogate, conc_inorganic, LWC, DT2, tequilibrium,0);

  for (i=0;i<n;++i)
    surrogate[i].Ag1=surrogate[i].Ag;

  //compute first evaluation of organic and aqueous masses 
  MOinit2=0.0;
  for (b=0;b<config.nbins;++b)		  
    for (ilayer=0;ilayer<config.nlayer;++ilayer)
      for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
        for (i=0;i<n;++i)
          MOinit2(b,ilayer,iphase)+=max(surrogate[i].Ap_layer_init(b,ilayer,iphase),
                                        config.MOmin*config.Vlayer(ilayer));
  
  for (b=0;b<config.nbins;++b)
    {
      AQinit2(b)=LWC(b)+conc_inorganic(b);
      for (i=0;i<n;++i)
        if(surrogate[i].is_organic or i==config.iH2O)
          AQinit2(b)+=surrogate[i].Aaq_bins_init(b);
    }
  
  chp1=chp;
  error_chp=1;
  iter=0;  
  while (iter<50 and error_chp>1.0e-3)
    {
      error_chp=0.0;
      iter++;
      if (config.compute_inorganic)
	{	  
	  compute_ph_dyn_ssh(config, surrogate, Temperature, chp, AQinit2, ionic, MMaq, LWC);	  
	  for (b=0;b<config.nbins;++b)
	    {
              compute_conc_org_bins_ssh(config, surrogate, Temperature,  MMaq, AQinit, chp, ionic, LWC, b, conc_org);
              /*
              conc_org=LWC(b);
              for (i=0;i<n;++i)
                if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophilic)
                  conc_org+=surrogate[i].Aaq_bins_init(b);
                  conc_org=max(conc_org,1.0e-5*AQinit(b)); //config.MOmin);*/

	      //config.rho_aqueous=config.AQrho(b);	      	      
	      compute_ionic_strenght2_ssh(config, surrogate, Temperature, AQinit2(b), conc_inorganic(b), ionic(b), chp(b),
					  organion(b), ionic_organic(b), conc_org, 1.0, config.Ke/surrogate[config.iHp].gamma_aq_bins(b)/surrogate[config.iOHm].gamma_aq_bins(b));	   
	      error_chp=max(error_chp,abs(chp1(b)-chp(b))/chp1(b));
	      chp1(b)=chp(b);	      
	    } 
	}
               
      if (compute_activity_coefficients)
	{
	  //compute activity coefficients
	  if (config.compute_organic)
	    {
	      activity_coefficients_dyn_org_ssh(config, surrogate, Temperature, MOW);
	      compute_kp_org_ssh(config, surrogate, MOinit2, Temperature, MOW);
	    }
	  
	  if (LWCtot>config.LWClimit)
	    {
	      activity_coefficients_dyn_aq_ssh(config, surrogate, Temperature,AQinit2,
					       MOinit2, conc_inorganic, ionic, ionic_organic,
					       organion,chp,LWC,MMaq,1.0, 0., 0);
	      compute_kp_aq_ssh(config, surrogate, Temperature, ionic, chp, AQinit2, MMaq, 0, config.nbins);
	    }
	}      
    }

  //compute the second evaluation of kinetic rates
  flux_org_ssh(config, surrogate, MOinit2, AQinit2, tiny, 1);
  if(LWCtot>config.LWClimit)
    {
      flux_aq_ssh(config, surrogate, AQinit2, LWC, MMaq, chp, ionic, MOinit2, tiny, Temperature, 1);           
      //if (config.compute_inorganic and facph>0)
      //correct_flux_ph_ssh(config, surrogate, Temperature, AQinit2, MOinit2, chp, chp2, MMaq, ionic, LWC, tiny, DT2/facph, 1);
    }

  if (config.chemistry)
    compute_flux_chem_ssh(config,surrogate,MOinit2,MOW,AQinit2,LWC,MMaq,chp,ionic,DT2,tiny,Temperature,1,RH);

  //compute the second evaluation of concentrations
  for (i=0;i<n;++i)
    if((config.compute_organic and surrogate[i].is_organic) or i==config.iH2O) 
      {
	conc_available=surrogate[i].Ag;
	sum_rates=0.0;
        if (config.chemistry==false)
          {
            if (config.compute_organic)
              if (surrogate[i].hydrophobic)
                for (b=0;b<config.nbins;++b)
                  for (ilayer=0;ilayer<config.nlayer;++ilayer)
                    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)		  
                      if (surrogate[i].time(b,ilayer,iphase)<tequilibrium)
                        conc_available+=surrogate[i].Ap_layer_init(b,ilayer,iphase);
                      else
                        sum_rates+=surrogate[i].k1(b,ilayer,iphase,1);
	
            if (surrogate[i].hydrophilic and LWCtot>config.LWClimit)// and i!=config.iH2O)
              for (b=0;b<config.nbins;++b)
                if (surrogate[i].time_aq(b)<tequilibrium)
                  conc_available+=surrogate[i].Aaq_bins_init(b);
                else
                  sum_rates+=surrogate[i].k1_aq(b,1);
	
        
            if (config.compute_organic)
              if (surrogate[i].hydrophobic)
                for (b=0;b<config.nbins;++b)		  
                  for (ilayer=0;ilayer<config.nlayer;++ilayer)
                    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                      if (surrogate[i].time(b,ilayer,iphase)>=tequilibrium)
                        if (conc_available-surrogate[i].Ap_layer_init(b,ilayer,iphase)> tiny and
                            surrogate[i].k1(b,ilayer,iphase,1)>0.0)
                          surrogate[i].Jdn(b,ilayer,iphase,1)=-sum_rates/
                            (conc_available-surrogate[i].Ap_layer_init(b,ilayer,iphase));
	
            if (surrogate[i].hydrophilic and LWCtot>config.LWClimit) // and i!=config.iH2O)
              for (b=0;b<config.nbins;++b)			  
                if (surrogate[i].time_aq(b)>=tequilibrium)
                  if (conc_available-surrogate[i].Aaq_bins_init(b)> tiny and
                      surrogate[i].k1_aq(b,1)>0.0)
                    surrogate[i].Jdn_aq(b,1)=-sum_rates/
                      (conc_available-surrogate[i].Aaq_bins_init(b));
          }
	
        surrogate[i].ktot2=surrogate[i].flux_chem_gas(1);
	if (config.compute_organic)
	  if (surrogate[i].hydrophobic)
	    for (b=0;b<config.nbins;++b)		  
	      for (ilayer=0;ilayer<config.nlayer;++ilayer)
		for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                  {
                    surrogate[i].k1(b,ilayer,iphase,1)=
                      (surrogate[i].k1(b,ilayer,iphase,1)*DT2
                       -gamma*DT2*(surrogate[i].Jdn(b,ilayer,iphase,0)+surrogate[i].Jdn(b,ilayer,iphase,1))
                       *surrogate[i].k1(b,ilayer,iphase,0))/
                      (1.0-gamma*surrogate[i].Jdn(b,ilayer,iphase,1)*DT2)+surrogate[i].flux_chem(b,ilayer,iphase,1);
                    /*surrogate[i].flux_chem(b,ilayer,iphase,1)=(surrogate[i].flux_chem(b,ilayer,iphase,1)
                      -gamma*DT2*(surrogate[i].Jdn(b,ilayer,iphase,0)+surr*/
                    surrogate[i].ktot2+=surrogate[i].flux_chem(b,ilayer,iphase,1);
                  }
	
	if (surrogate[i].hydrophilic and LWCtot>config.LWClimit) // and i!=config.iH2O)
	  for (b=0;b<config.nbins;++b)
            {
              surrogate[i].k1_aq(b,1)=surrogate[i].flux_chem_aq(b,1)+
                (surrogate[i].k1_aq(b,1)*DT2-gamma*DT2*(surrogate[i].Jdn_aq(b,0)+surrogate[i].Jdn_aq(b,1))
                 *surrogate[i].k1_aq(b,0))/(1.0-gamma*surrogate[i].Jdn_aq(b,1)*DT2);
              /*surrogate[i].flux_chem_aq(b,1)=(surrogate[i].flux_chem_aq(b,1)-gamma*DT2*(surrogate[i].Jdn_aq(b,0)+surrogate[i].Jdn_aq(b,1))
               *surrogate[i].flux_chem_aq(b,0))/(1.0-gamma*surrogate[i].Jdn_aq(b,1)*DT2);*/
              surrogate[i].ktot2+=surrogate[i].flux_chem_aq(b,1);
            }

        surrogate[i].Atot=max(redmax*surrogate[i].Atot0,surrogate[i].Atot0+0.5*(surrogate[i].ktot1+surrogate[i].ktot2));        
	surrogate[i].Ag=surrogate[i].Atot;
	if (config.compute_organic)
	  if (surrogate[i].hydrophobic)
	    for (b=0;b<config.nbins;++b)		  
	      for (ilayer=0;ilayer<config.nlayer;++ilayer)
		for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		  {
		    if (surrogate[i].time(b,ilayer,iphase)>=tequilibrium)
		      surrogate[i].Ap_layer(b,ilayer,iphase)=
                        max(surrogate[i].Ap_layer_init0(b,ilayer,iphase)*redmax,
                            surrogate[i].Ap_layer_init0(b,ilayer,iphase)+0.5*
                            (surrogate[i].k1(b,ilayer,iphase,0)+surrogate[i].k1(b,ilayer,iphase,1)));
		    
		    else
		      surrogate[i].Ap_layer(b,ilayer,iphase)=
                        max(surrogate[i].Ap_layer_init0(b,ilayer,iphase)*redmax,
                            surrogate[i].Ap_layer_init0(b,ilayer,iphase)+0.5*
                            (surrogate[i].flux_chem(b,ilayer,iphase,0)+surrogate[i].flux_chem(b,ilayer,iphase,1)));
		    surrogate[i].Ag-=surrogate[i].Ap_layer(b,ilayer,iphase);
		  }

	if (surrogate[i].hydrophilic and LWCtot>config.LWClimit) // and i!=config.iH2O)
	  for (b=0;b<config.nbins;++b)
	    {
	      if (surrogate[i].time_aq(b)>=tequilibrium)
		surrogate[i].Aaq_bins(b)=
		  max(surrogate[i].Aaq_bins_init0(b)*redmax,
		      surrogate[i].Aaq_bins_init0(b)+0.5*
		      (surrogate[i].k1_aq(b,0)+surrogate[i].k1_aq(b,1)));
	      
	      else
		surrogate[i].Aaq_bins(b)=
                  max(surrogate[i].Aaq_bins_init0(b)*redmax,
		      surrogate[i].Aaq_bins_init0(b)+0.5*
		      (surrogate[i].flux_chem_aq(b,0)+surrogate[i].flux_chem_aq(b,1)));
	      surrogate[i].Ag-=surrogate[i].Aaq_bins(b);
	    }
	else
	  for (b=0;b<config.nbins;++b)
	    surrogate[i].Aaq_bins(b)=surrogate[i].Aaq_bins_init(b);
      }
    else
      {		
	for (b=0;b<config.nbins;++b)		  
	  for (ilayer=0;ilayer<config.nlayer;++ilayer)
	    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
	      surrogate[i].Ap_layer(b,ilayer,iphase)=surrogate[i].Ap_layer_init(b,ilayer,iphase);
	
	for (b=0;b<config.nbins;++b)
	  surrogate[i].Aaq_bins(b)=surrogate[i].Aaq_bins_init(b);
      }

  if (config.compute_inorganic)
    dynamic_inorg_ssh(config, surrogate, conc_inorganic, LWC, DT2, tequilibrium,1);  

  //make sure that the sum of concentrations is not higher than the total concentration  
  if (config.chemistry==false)
    if (config.compute_organic)
      {
        for (i=0;i<n;++i)
          if (surrogate[i].is_organic)
            {
              sumconc=0.0;
              sumconc2=0.0;
              if (surrogate[i].hydrophobic)
                for (b=0;b<config.nbins;++b)		  
                  for (ilayer=0;ilayer<config.nlayer;++ilayer)
                    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                      if (surrogate[i].time(b,ilayer,iphase)>=tequilibrium)
                        sumconc+=surrogate[i].Ap_layer(b,ilayer,iphase);
                      else
                        sumconc2+=surrogate[i].Ap_layer(b,ilayer,iphase);
	  
              if (surrogate[i].hydrophilic and LWCtot>config.LWClimit)
                for (b=0;b<config.nbins;++b)
                  if (surrogate[i].time_aq(b)>=tequilibrium)
                    sumconc+=surrogate[i].Aaq_bins(b);
                  else
                    sumconc2+=surrogate[i].Aaq_bins(b);
		
              if(sumconc>surrogate[i].Atot)
                {
                  surrogate[i].Ag=0.0;
                  if (surrogate[i].hydrophobic)
                    for (b=0;b<config.nbins;++b)		  
                      for (ilayer=0;ilayer<config.nlayer;++ilayer)
                        for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                          if (surrogate[i].time(b,ilayer,iphase)>=tequilibrium)
                            surrogate[i].Ap_layer(b,ilayer,iphase)=
                              surrogate[i].Ap_layer(b,ilayer,iphase)/sumconc*
                              0.999*surrogate[i].Atot;
                          else
                            surrogate[i].Ap_layer(b,ilayer,iphase)=
                              surrogate[i].Ap_layer(b,ilayer,iphase)/sumconc2*
                              0.001*surrogate[i].Atot;
			
                  if (LWCtot>config.LWClimit)
                    if (surrogate[i].hydrophilic)
                      for (b=0;b<config.nbins;++b)
                        if (surrogate[i].time_aq(b)>=tequilibrium)
                          if (sumconc > 0.0)
                            surrogate[i].Aaq_bins(b)=surrogate[i].Aaq_bins(b)/sumconc*
                              0.999*surrogate[i].Atot;
			  else
			    surrogate[i].Aaq_bins(b)=surrogate[i].Aaq_bins(b)/sumconc2*
			      0.001*surrogate[i].Atot;
                }
              else
                if(sumconc+sumconc2>surrogate[i].Atot)
                  {
                    surrogate[i].Ag=0.0;
                    if (surrogate[i].hydrophobic)
                      for (b=0;b<config.nbins;++b)		  
                        for (ilayer=0;ilayer<config.nlayer;++ilayer)
                          for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                            if (surrogate[i].time(b,ilayer,iphase)<tequilibrium)
                              surrogate[i].Ap_layer(b,ilayer,iphase)=
                                surrogate[i].Ap_layer(b,ilayer,iphase)/sumconc2*
                                (surrogate[i].Atot-sumconc);
			  
                    if (LWCtot>config.LWClimit)
                      if (surrogate[i].hydrophilic)
                        for (b=0;b<config.nbins;++b)
                          surrogate[i].Aaq_bins(b)=surrogate[i].Aaq_bins(b)/sumconc2*
                            (surrogate[i].Atot-sumconc);
                  }
		
              surrogate[i].Ag=surrogate[i].Atot;
              if (surrogate[i].hydrophobic)
                for (b=0;b<config.nbins;++b)		  
                  for (ilayer=0;ilayer<config.nlayer;++ilayer)
                    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                      surrogate[i].Ag-=surrogate[i].Ap_layer(b,ilayer,iphase);
		
              if (LWCtot>config.LWClimit)
                if (surrogate[i].hydrophilic)
                  for (b=0;b<config.nbins;++b)
                    surrogate[i].Ag-=surrogate[i].Aaq_bins(b);
            }
        for (b=0;b<config.nbins;++b)
          for (ilayer=0;ilayer<config.nlayer;++ilayer)
            {
              for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                {
                  MO(b,ilayer,iphase)=0.0;
                  MOinit(b,ilayer,iphase)=0.0;
                  for (i=0;i<n;++i)
                    {
                      MO(b,ilayer,iphase)+=surrogate[i].Ap_layer(b,ilayer,iphase);
                      MOinit(b,ilayer,iphase)+=surrogate[i].Ap_layer_init(b,ilayer,iphase);
                    }
                }

              sumt=0.0;
              for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                sumt+=MO(b,ilayer,iphase);

              if (sumt>0.0)
                for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                  MO(b,ilayer,iphase)=max(MO(b,ilayer,iphase),
                                          config.MOmin*config.Vlayer(ilayer)*MO(b,ilayer,iphase)/sumt);
              else
                for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                  if(iphase==0)
                    MO(b,ilayer,iphase)=config.MOmin*config.Vlayer(ilayer);
                  else
                    MO(b,ilayer,iphase)=0.0;

              sumt=0.0;
              for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                sumt+=MOinit(b,ilayer,iphase);

              if (sumt>0.0)
                for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                  MOinit(b,ilayer,iphase)=max(MOinit(b,ilayer,iphase),
                                              config.MOmin*config.Vlayer(ilayer)*MOinit(b,ilayer,iphase)/sumt);
              else
                for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                  if(iphase==0)
                    MOinit(b,ilayer,iphase)=config.MOmin*config.Vlayer(ilayer);
                  else
                    MOinit(b,ilayer,iphase)=0.0;
            }
      }
    else
      for (b=0;b<config.nbins;++b)
        {
          for (ilayer=0;ilayer<config.nlayer;++ilayer)
            for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
              {
                MO(b,ilayer,iphase)=MOinit(b,ilayer,iphase);
                for (i=0;i<n;++i)
                  surrogate[i].Ap_layer(b,ilayer,iphase)=surrogate[i].Ap_layer_init(b,ilayer,iphase);
              }

          for (i=0;i<n;++i)
            if(surrogate[i].is_organic)
              surrogate[i].Aaq_bins(b)=surrogate[i].Aaq_bins_init(b);      
        }

 
  if (config.chemistry==false)
    if (config.compute_inorganic)
      for (i=0;i<n;++i)
	if (surrogate[i].is_inorganic_precursor)
	  {

	    double total=surrogate[i].Ag/surrogate[i].MM;
	    double total0=surrogate[i].Ag0/surrogate[i].MM;
	    for (b=0;b<config.nbins;b++)
	      if (surrogate[i].time_aq(b)>=config.tequilibrium)
		{
		  if (i==config.iHCl)
		    {
		      total+=surrogate[config.iClm].Aaq_bins(b)/surrogate[config.iClm].MM;
		      total0+=surrogate[config.iClm].Aaq_bins_init0(b)/surrogate[config.iClm].MM;
		    }
		  else if (i==config.iNH3)
		    {
		      total+=surrogate[config.iNH4p].Aaq_bins(b)/surrogate[config.iNH4p].MM;
		      total0+=surrogate[config.iNH4p].Aaq_bins_init0(b)/surrogate[config.iNH4p].MM;
		    }
		  else if (i==config.iHNO3)
		    {
		      total+=surrogate[config.iNO3m].Aaq_bins(b)/surrogate[config.iNO3m].MM;
		      total0+=surrogate[config.iNO3m].Aaq_bins_init0(b)/surrogate[config.iNO3m].MM;
		    }
		  else if (i==config.iH2SO4)
		    {
		      total+=surrogate[config.iSO4mm].Aaq_bins(b)/surrogate[config.iSO4mm].MM
			+surrogate[config.iHSO4m].Aaq_bins(b)/surrogate[config.iHSO4m].MM;
		      total0+=surrogate[config.iSO4mm].Aaq_bins_init0(b)/surrogate[config.iSO4mm].MM
			+surrogate[config.iHSO4m].Aaq_bins_init0(b)/surrogate[config.iHSO4m].MM;
		    }
		
		}

	    if (total>0)
	      {
		surrogate[i].Ag=surrogate[i].Ag/total*total0;		
		for (b=0;b<config.nbins;b++)
		  if (surrogate[i].time_aq(b)>=config.tequilibrium)
		    {
		      if (i==config.iHCl)
			{
			  surrogate[config.iClm].Aaq_bins(b)*=total0/total;
			}
		      else if (i==config.iNH3)
			{
			  surrogate[config.iNH4p].Aaq_bins(b)*=total0/total;
			}
		      else if (i==config.iHNO3)
			{
			  surrogate[config.iNO3m].Aaq_bins(b)*=total0/total;
			}
		      else if (i==config.iH2SO4)
			{
			  surrogate[config.iSO4mm].Aaq_bins(b)*=total0/total;
			  surrogate[config.iHSO4m].Aaq_bins(b)*=total0/total;
			}
				
		    }
	      }	  
	  }
  
  for (b=0;b<config.nbins;++b)
    {
      AQ(b)=LWC(b); //+conc_inorganic(b);
      AQinit(b)=LWC(b);
      for (i=0;i<n;++i)
        if(surrogate[i].hydrophilic) //surrogate[i].is_organic or i==config.iH2O)
	  {
	    AQ(b)+=surrogate[i].Aaq_bins(b);
	    AQinit(b)+=surrogate[i].Aaq_bins_init(b);
	  }

    }

  if (config.compute_inorganic)
    {
      compute_ph_dyn2_ssh(config, surrogate, Temperature, chp, AQ, ionic, MMaq, LWC);
      for (b=0;b<config.nbins;++b)
	{
          compute_conc_org_bins_ssh(config, surrogate, Temperature,  MMaq, AQinit, chp, ionic, LWC, b, conc_org);

	  //config.rho_aqueous=config.AQrho(b);
	  compute_ionic_strenght2_ssh(config, surrogate, Temperature, AQ(b), conc_inorganic(b), ionic(b), chp(b),
				      organion(b), ionic_organic(b), conc_org, 1.0, config.Ke/surrogate[config.iHp].gamma_aq_bins(b)/surrogate[config.iOHm].gamma_aq_bins(b));	 	  
	}    
    }
}

void dynamic_aq_ssh(model_config &config, vector<species>& surrogate,
		    Array<double, 1> &AQinit,Array<double, 1> &AQ,
		    Array<double, 3> &MOinit,
		    Array <double, 1> &conc_inorganic,
		    Array <double, 1> &ionic,Array <double, 1> &ionic_organic,
		    Array <double, 1> &organion,Array <double, 1> &chp,
		    Array <double, 1> &chp1, Array <double, 1> &chp0,
		    Array <double, 1> &LWC,
		    Array<double, 1> &MMaq,Array<double, 3> &MOW,
		    double &Temperature, double &RH, double &DT2, double &tequilibrium,
		    bool compute_activity_coefficients)
{
  //This routine computes the aqueous-phase concentrations being at equilibrium
  //FOR THE UNCOUPLED SYSTEM
  //compute_activity_coefficients: should the activity coefficients be recomputed?
  //tequilibrium: criterium to distinguish cases under equilibrium (characteristic time lower
  //                  than tequilibrium)
  int n=surrogate.size();
  int i,b;
  //int index_b;
  //double sum=0.0;
  double tiny=1.0e-26;
  double conc_available=0.0;
  double sum_rates;
  double gamma=1.7071;
  //double Kaq;
  Array<double, 1> AQinit2,chp2;
  AQinit2.resize(config.nbins);
  chp2.resize(config.nbins);
  //double conc_equilibrium;
  //     gamma > .25     ->  A-stable  ( x(n+1)-x(n))*f(x(n)) > 0
  //     gamma < .25     ->  monotonic behaviour (x(n+1)-x*)(x(n)-x*) > 0
  //     gamma =  1+-sqrt(1/2) ->  L-stability


  chp0=chp;
  double error_chp=1;
  int iter=0;
  double sumconc,sumconc2;
  double conc_org;
  double inorganion=0.0;

  while (iter<50 and error_chp>1.0e-3)
    {
      error_chp=0.0;
      iter++;
      if (config.compute_inorganic)	
	{
	  compute_ph_dyn_ssh(config, surrogate, Temperature, chp, AQinit, ionic, MMaq, LWC);
	  for (b=0;b<config.nbins;++b)
	    {
              compute_conc_org_bins_ssh(config, surrogate, Temperature,  MMaq, AQinit, chp, ionic, LWC, b, conc_org);

	      //config.rho_aqueous=config.AQrho(b);
	      compute_ionic_strenght2_ssh(config, surrogate, Temperature, AQinit(b), conc_inorganic(b), ionic(b), chp(b),
					  organion(b), ionic_organic(b), conc_org, 1.0, config.Ke/surrogate[config.iHp].gamma_aq_bins(b)/surrogate[config.iOHm].gamma_aq_bins(b));
	      error_chp=max(error_chp,abs(chp0(b)-chp(b))/chp0(b));
	      chp0(b)=chp(b);	      
	    }      
	}
	  
      if (compute_activity_coefficients)	         
	{
	  activity_coefficients_dyn_aq_ssh(config, surrogate, Temperature,AQinit,
					   MOinit,conc_inorganic, ionic, ionic_organic,
					   organion,chp,LWC,MMaq, 1.0, 0., 0);
	  compute_kp_aq_ssh(config, surrogate, Temperature, ionic, chp, AQinit, MMaq, 0, config.nbins);   	
	}
    }
  
  //save initial concentrations
  for (i=0;i<n;++i)
    if (surrogate[i].hydrophilic)
      {
        surrogate[i].Ag0=surrogate[i].Ag;
        for (b=0;b<config.nbins;++b)
          surrogate[i].Aaq_bins_init0(b)=surrogate[i].Aaq_bins_init(b);
      }

  //first evaluation of kinetic rates
  flux_aq_ssh(config, surrogate, AQinit, LWC, MMaq, chp, ionic, MOinit, tiny, Temperature, 0);
  //if (config.compute_inorganic)
  //  correct_flux_ph_ssh(config, surrogate, Temperature, AQinit, MOinit, chp, chp2, MMaq, ionic, LWC, tiny, DT2, 0);
  if (config.chemistry)
    compute_flux_chem_ssh(config,surrogate,MOinit,MOW,AQinit,LWC,MMaq,chp,ionic,DT2,tiny,Temperature,0,RH);

  //first evaluation of concentrations
  for (i=0;i<n;++i)
    if(((config.compute_organic and surrogate[i].is_organic) or i==config.iH2O) and surrogate[i].hydrophilic) 
      {
	conc_available=surrogate[i].Ag;
	sum_rates=0.0;
	for (b=0;b<config.nbins;++b)			  
	  if (surrogate[i].time_aq(b)<tequilibrium)
	    conc_available+=surrogate[i].Aaq_bins_init(b);
	  else
	    sum_rates+=surrogate[i].k1_aq(b,0);
	
	for (b=0;b<config.nbins;++b)			  
	  if (surrogate[i].time_aq(b)>=tequilibrium)
	    if (conc_available-surrogate[i].Aaq_bins_init(b)> tiny and
		surrogate[i].k1_aq(b,0)>0.0)
	      surrogate[i].Jdn_aq(b,0)=-sum_rates/
		(conc_available-surrogate[i].Aaq_bins_init(b));
	
	for (b=0;b<config.nbins;++b)
	  surrogate[i].k1_aq(b,0)=surrogate[i].k1_aq(b,0)*
	    DT2/(1.0-gamma*surrogate[i].Jdn_aq(b,0)*DT2);
	
	surrogate[i].Ag=surrogate[i].Atot;
	for (b=0;b<config.nbins;++b)
	  {
	    if (surrogate[i].time_aq(b)>=tequilibrium)
	      surrogate[i].Aaq_bins_init(b)=
		max(0.01*surrogate[i].Aaq_bins_init(b),
		    surrogate[i].Aaq_bins_init(b)+surrogate[i].k1_aq(b,0));
	    
	    surrogate[i].Ag-=surrogate[i].Aaq_bins_init(b);
	  }
	surrogate[i].Ag=max(surrogate[i].Ag,0.01*surrogate[i].Ag0);	
      }
  
  if (config.compute_inorganic)
    dynamic_inorg_ssh(config, surrogate, conc_inorganic, LWC, DT2, tequilibrium, 0); 

  for (i=0;i<n;++i)
    if (surrogate[i].hydrophilic)      
      surrogate[i].Ag1=surrogate[i].Ag;

  //first evalution of aqueous mass
  for (b=0;b<config.nbins;++b)
    {
      AQinit2(b)=LWC(b)+conc_inorganic(b);
      for (i=0;i<n;++i)
        if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophilic)
          AQinit2(b)+=surrogate[i].Aaq_bins_init(b);
    }

  if (config.compute_inorganic)
    {
      for (b=0;b<config.nbins;++b)
	{
          compute_conc_org_bins_ssh(config, surrogate, Temperature,  MMaq, AQinit, chp, ionic, LWC, b, conc_org);

	  inorganion=0.0;
	  for (i=0;i<n;++i)
	    if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)
	      {
		surrogate[i].molality=surrogate[i].Aaq_bins_init(b)/surrogate[i].MM/conc_org*1000.0;		
		if (i!=config.iHp and i!=config.iHp)
		  inorganion-=surrogate[i].molality*surrogate[i].charge;
	      }
	  
	  chp(b)=max(0.5*(organion(b)+inorganion+pow(pow(organion(b)+inorganion,2)+4*config.Ke/surrogate[config.iHp].gamma_aq_bins(b)/surrogate[config.iOHm].gamma_aq_bins(b),0.5)),1.e-14);
	  if (chp(b)==0.0)
	    chp(b)=pow(10.0,-5.6);	  
	  
	  //config.rho_aqueous=config.AQrho(b);
	  compute_ionic_strenght2_ssh(config, surrogate, Temperature, AQinit2(b), conc_inorganic(b), ionic(b), chp(b),
				      organion(b), ionic_organic(b), conc_org, 1.0, config.Ke/surrogate[config.iHp].gamma_aq_bins(b)/surrogate[config.iOHm].gamma_aq_bins(b));

	}      
    }

  //compute activity coefficients
  if (compute_activity_coefficients)
    {
      activity_coefficients_dyn_aq_ssh(config, surrogate, Temperature,AQinit2,
				       MOinit,conc_inorganic, ionic, ionic_organic,
				       organion,chp,LWC,MMaq, 1.0, 0., 0);
      compute_kp_aq_ssh(config, surrogate, Temperature, ionic, chp, AQinit2, MMaq, 0, config.nbins);
    }

  chp1=chp;
  error_chp=1;
  iter=0;  

  while (iter<50 and error_chp>1.0e-3)
    {
      error_chp=0.0;
      iter++;
      if (config.compute_inorganic)
	{	  
	  compute_ph_dyn_ssh(config, surrogate, Temperature, chp, AQinit2, ionic, MMaq, LWC);	  
	  for (b=0;b<config.nbins;++b)
	    {
              compute_conc_org_bins_ssh(config, surrogate, Temperature,  MMaq, AQinit, chp, ionic, LWC, b, conc_org);

	      //config.rho_aqueous=config.AQrho(b);	      	      
	      compute_ionic_strenght2_ssh(config, surrogate, Temperature, AQinit2(b), conc_inorganic(b), ionic(b), chp(b),
					  organion(b), ionic_organic(b), conc_org, 1.0, config.Ke/surrogate[config.iHp].gamma_aq_bins(b)/surrogate[config.iOHm].gamma_aq_bins(b));	   
	      error_chp=max(error_chp,abs(chp1(b)-chp(b))/chp1(b));
	      chp1(b)=chp(b);	      
	    } 
	}
               
      if (compute_activity_coefficients)
	{
	  //compute activity coefficients	  	  	
	  activity_coefficients_dyn_aq_ssh(config, surrogate, Temperature,AQinit2,
					   MOinit, conc_inorganic, ionic, ionic_organic,
					   organion,chp,LWC,MMaq,1.0,0., 0);
	  compute_kp_aq_ssh(config, surrogate, Temperature, ionic, chp, AQinit2, MMaq, 0, config.nbins);	    
	}      
    }  

  //second estimation of kinetic rates
  flux_aq_ssh(config, surrogate, AQinit2, LWC, MMaq, chp, ionic, MOinit, tiny, Temperature, 1);
  //if (config.compute_inorganic)
  //  correct_flux_ph(config, surrogate, Temperature, AQinit2, MOinit, chp, chp2, MMaq, ionic, LWC, tiny, DT2, 1);
  if (config.chemistry)
    compute_flux_chem_ssh(config,surrogate,MOinit,MOW,AQinit,LWC,MMaq,chp,ionic,DT2,tiny,Temperature,1,RH);

  //second estimation of concentrations
  for (i=0;i<n;++i)
    if(((config.compute_organic and surrogate[i].is_organic) or i==config.iH2O) and surrogate[i].hydrophilic)      
      {
	conc_available=surrogate[i].Ag;
	sum_rates=0.0;
	for (b=0;b<config.nbins;++b)			  
          if (surrogate[i].time_aq(b)<tequilibrium)
            conc_available+=surrogate[i].Aaq_bins_init(b);
          else
            sum_rates+=surrogate[i].k1_aq(b,1);
	
	for (b=0;b<config.nbins;++b)			  
	  if (surrogate[i].time_aq(b)>=tequilibrium)
	    if (conc_available-surrogate[i].Aaq_bins_init(b)> tiny and
                surrogate[i].k1_aq(b,1)>0.0)
	      surrogate[i].Jdn_aq(b,1)=-sum_rates/
		(conc_available-surrogate[i].Aaq_bins_init(b));
	
	for (b=0;b<config.nbins;++b)
	  surrogate[i].k1_aq(b,1)=
	    (surrogate[i].k1_aq(b,1)*DT2-gamma*DT2*(surrogate[i].Jdn_aq(b,0)+surrogate[i].Jdn_aq(b,1))
	     *surrogate[i].k1_aq(b,0))/(1.0-gamma*surrogate[i].Jdn_aq(b,1)*DT2);
	
	surrogate[i].Ag=surrogate[i].Atot;
	for (b=0;b<config.nbins;++b)
	  {
	    if (surrogate[i].time_aq(b)>=tequilibrium)
	      surrogate[i].Aaq_bins(b)=
		max(surrogate[i].Aaq_bins_init0(b)*0.01,
		    surrogate[i].Aaq_bins_init0(b)+0.5*
		    (surrogate[i].k1_aq(b,0)+surrogate[i].k1_aq(b,1)));
	    
	    else
	      surrogate[i].Aaq_bins(b)=surrogate[i].Aaq_bins_init(b);
	    surrogate[i].Ag-=surrogate[i].Aaq_bins(b);
	  }
      }
    else
      for (b=0;b<config.nbins;++b)
	surrogate[i].Aaq_bins(b)=surrogate[i].Aaq_bins_init(b);

  if (config.compute_inorganic)
    dynamic_inorg_ssh(config, surrogate, conc_inorganic, LWC, DT2, tequilibrium, 1); 
  
  //make sure that the sum of concentrations is not higher than the total concentration
  if (config.compute_organic)
    {
      for (i=0;i<n;++i)
        if(surrogate[i].is_organic and surrogate[i].hydrophilic)
          {
            sumconc=0.0;
            sumconc2=0.0;
            for (b=0;b<config.nbins;++b)
              if (surrogate[i].time_aq(b)>=tequilibrium)
                sumconc+=surrogate[i].Aaq_bins(b);
              else
                sumconc2+=surrogate[i].Aaq_bins(b);

            if(sumconc>surrogate[i].Atot)
              {
                surrogate[i].Ag=0.0;
                for (b=0;b<config.nbins;++b)
                  if (surrogate[i].time_aq(b)>=tequilibrium)
                    surrogate[i].Aaq_bins(b)=surrogate[i].Aaq_bins(b)/sumconc*
                      0.999*surrogate[i].Atot;
                  else
                    if (sumconc2 != 0.0)
                      surrogate[i].Aaq_bins(b)=surrogate[i].Aaq_bins(b)/sumconc2*
                        0.001*surrogate[i].Atot;
              }
            else
              if(sumconc+sumconc2>surrogate[i].Atot)
                {
                  surrogate[i].Ag=0.0;
                  for (b=0;b<config.nbins;++b)
                    surrogate[i].Aaq_bins(b)=surrogate[i].Aaq_bins(b)/sumconc2*
                      (surrogate[i].Atot-sumconc);
                }
          }  
    }
  else
    for (b=0;b<config.nbins;++b)
      for (i=0;i<n;++i)
        if(surrogate[i].is_organic)
          surrogate[i].Aaq_bins(b)=surrogate[i].Aaq_bins_init(b);

  if (config.chemistry==false)
    if (config.compute_inorganic)
      for (i=0;i<n;++i)
	if (surrogate[i].is_inorganic_precursor)
	  {

	    double total=surrogate[i].Ag/surrogate[i].MM;
	    double total0=surrogate[i].Ag0/surrogate[i].MM;
	    for (b=0;b<config.nbins;b++)
	      if (surrogate[i].time_aq(b)>=config.tequilibrium)
		{
		  if (i==config.iHCl)
		    {
		      total+=surrogate[config.iClm].Aaq_bins(b)/surrogate[config.iClm].MM;
		      total0+=surrogate[config.iClm].Aaq_bins_init0(b)/surrogate[config.iClm].MM;
		    }
		  else if (i==config.iNH3)
		    {
		      total+=surrogate[config.iNH4p].Aaq_bins(b)/surrogate[config.iNH4p].MM;
		      total0+=surrogate[config.iNH4p].Aaq_bins_init0(b)/surrogate[config.iNH4p].MM;
		    }
		  else if (i==config.iHNO3)
		    {
		      total+=surrogate[config.iNO3m].Aaq_bins(b)/surrogate[config.iNO3m].MM;
		      total0+=surrogate[config.iNO3m].Aaq_bins_init0(b)/surrogate[config.iNO3m].MM;
		    }
		  else if (i==config.iH2SO4)
		    {
		      total+=surrogate[config.iSO4mm].Aaq_bins(b)/surrogate[config.iSO4mm].MM
			+surrogate[config.iHSO4m].Aaq_bins(b)/surrogate[config.iHSO4m].MM;
		      total0+=surrogate[config.iSO4mm].Aaq_bins_init0(b)/surrogate[config.iSO4mm].MM
			+surrogate[config.iHSO4m].Aaq_bins_init0(b)/surrogate[config.iHSO4m].MM;
		    }
		
		}

	    if (total>0)
	      {
		surrogate[i].Ag=surrogate[i].Ag/total*total0;		
		for (b=0;b<config.nbins;b++)
		  if (surrogate[i].time_aq(b)>=config.tequilibrium)
		    {
		      if (i==config.iHCl)
			{
			  surrogate[config.iClm].Aaq_bins(b)*=total0/total;
			}
		      else if (i==config.iNH3)
			{
			  surrogate[config.iNH4p].Aaq_bins(b)*=total0/total;
			}
		      else if (i==config.iHNO3)
			{
			  surrogate[config.iNO3m].Aaq_bins(b)*=total0/total;
			}
		      else if (i==config.iH2SO4)
			{
			  surrogate[config.iSO4mm].Aaq_bins(b)*=total0/total;
			  surrogate[config.iHSO4m].Aaq_bins(b)*=total0/total;
			}
				
		    }
	      }
	  
	  }
  
  for (b=0;b<config.nbins;++b)
    {
      AQ(b)=LWC(b)+conc_inorganic(b);
      for (i=0;i<n;++i)
        if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophilic)
          AQ(b)+=surrogate[i].Aaq_bins(b);
    }

  if (config.compute_inorganic)
    {
      compute_ph_dyn2_ssh(config, surrogate, Temperature, chp, AQ, ionic, MMaq, LWC);
      for (b=0;b<config.nbins;++b)
	{	  
	  //config.rho_aqueous=config.AQrho(b);
          //compute_conc_org_bins_ssh(config, surrogate, Temperature,  MMaq, AQinit, chp, LWC, b, conc_org);
         
          double conc_org=LWC(b);
          for (i=0;i<n;i++)
            if (surrogate[i].is_organic or i==config.iH2O)
              if (surrogate[i].hydrophilic)
                if (surrogate[i].aqt==2 or surrogate[i].aqt==1)
                  {
                    double fion1,fion2;
                    double Kp=surrogate[i].Kp_eff_aqrealdyn_ssh(config, Temperature, ionic(b), chp(b),surrogate[config.iHp].LR(b),
                                                         surrogate[config.iHp].SRMR(b),MMaq(b), fion1, fion2,b);
                    conc_org+=surrogate[i].Aaq_bins(b)*(1-fion1-fion2);
                  }
                else
                  conc_org+=surrogate[i].Aaq_bins(b);
          conc_org=max(conc_org,1.0e-5*AQinit(b));
	  compute_ionic_strenght2_ssh(config, surrogate, Temperature, AQ(b), conc_inorganic(b), ionic(b), chp(b),
				      organion(b), ionic_organic(b), conc_org,1.0, config.Ke/surrogate[config.iHp].gamma_aq_bins(b)/surrogate[config.iOHm].gamma_aq_bins(b));	 	  
	}    
    }
}

void adapstep_ssh(model_config &config, vector<species>& surrogate, double &Temperature, double &tequilibrium,
		  double &deltat1, double &t, double &tend, double &deltatmin,
		  Array<double, 3> &MOinit, Array<double, 3> &MO, double LWCtot,
		  Array<double, 1> &AQinit, Array<double, 1> &AQ,
		  Array<double, 1> &LWC, Array<double, 1> &conc_inorganic,
		  Array<double, 1> &chp, Array<double, 1> &chp1, Array<double, 1> &chp0,
		  Array<double, 1> &Number)
{
  //estimate the new time step necessary 
  double n2err=0.0;
  double n2err2=0.0;
  double n2err_aq=0.0;
  double tinym=1.e-5;
  double tinym2=1.e-3;
  double R=1.0e1;
  int n=surrogate.size();
  int i,b,ilayer,iphase;
  int m=0;
  int maq=0;
  double sum1,sum2;
  
  for (i=0;i<n;++i)
    if(surrogate[i].is_organic or i==config.iH2O) 
      {
        if(LWCtot>config.LWClimit or surrogate[i].hydrophobic)
          if (surrogate[i].Ag1 > tinym or surrogate[i].Ag > tinym)
            {
              ++m;
              n2err+=pow((surrogate[i].Ag-surrogate[i].Ag1)/(surrogate[i].Ag1+tinym),2);
	      
	      if (config.explicit_representation or config.nlayer>10)
		n2err2=max(n2err2,abs(surrogate[i].Ag-surrogate[i].Ag1)/(surrogate[i].Ag1+tinym));	     
            }

       
        surrogate[i].Ap=0.0;
        if(surrogate[i].hydrophobic)
          for (b=0;b<config.nbins;++b)
            for (ilayer=0;ilayer<config.nlayer;++ilayer)
              for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                surrogate[i].Ap+=surrogate[i].Ap_layer(b,ilayer,iphase);


        if(surrogate[i].hydrophobic)
          for (ilayer=0;ilayer<config.nlayer;++ilayer)
            {
              surrogate[i].Ap=0.0;
              for (b=0;b<config.nbins;++b)
                for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                  surrogate[i].Ap+=surrogate[i].Ap_layer(b,ilayer,iphase);
			  
              for (b=0;b<config.nbins;++b)
		for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		  if (surrogate[i].time(b,ilayer,iphase)>=tequilibrium)
                    {
                      if(surrogate[i].Ap_layer_init(b,ilayer,iphase)>config.Vlayer(ilayer)*tinym or surrogate[i].Ap_layer(b,ilayer,iphase)>config.Vlayer(ilayer)*tinym)                       	       	
                        if (config.explicit_representation or config.nlayer>10)
                          n2err2=max(n2err2,abs(surrogate[i].Ap_layer_init(b,ilayer,iphase)-surrogate[i].Ap_layer(b,ilayer,iphase))/(config.Vlayer(ilayer)*tinym+surrogate[i].Ap_layer_init(b,ilayer,iphase)));			
			  		      
                      sum1=0.0;
                      sum2=0.0;
                      for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                        {
                          sum1+=surrogate[i].Ap_layer_init(b,ilayer,iphase);
                          sum2+=surrogate[i].Ap_layer(b,ilayer,iphase);
                        }
				  
                      if(sum1>config.Vlayer(ilayer)*tinym 
                         or sum2>config.Vlayer(ilayer)*tinym) 
                        {
                          ++m;
                          n2err+=pow(abs(sum1-sum2)/(config.Vlayer(ilayer)*tinym+sum1),2);
                        }		    
                    }
	    }

	
        if(surrogate[i].hydrophilic and LWCtot>config.LWClimit) 
          for (b=0;b<config.nbins;++b)
            if (surrogate[i].time_aq(b)>=tequilibrium and
                (surrogate[i].Aaq_bins(b) > tinym or surrogate[i].Aaq_bins_init(b)> tinym))
              {
                ++maq;
                n2err_aq+=pow((surrogate[i].Aaq_bins(b)-
                               surrogate[i].Aaq_bins_init(b))
                              /(surrogate[i].Aaq_bins_init(b)+tinym),2);
              }
      }

 
  for (ilayer=0;ilayer<config.nlayer;++ilayer)
    for (b=0;b<config.nbins;++b)
      {
        sum1=0.0;
        sum2=0.0;
        for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
          {
            sum1+=MOinit(b,ilayer,iphase);
	    sum2+=MO(b,ilayer,iphase);
          }
	sum1=max(sum1,config.Vlayer(ilayer)*config.MOmin);
	sum2=max(sum2,config.Vlayer(ilayer)*config.MOmin);

	if (sum1 > tinym)
	  if (Number(b) > tinym)
	    if(sum1>config.Vlayer(ilayer)*tinym2/Number(b) or sum2>config.Vlayer(ilayer)*tinym2/Number(b))
	      n2err2=max(n2err2,abs(sum1-sum2)/(sum1));
      }

  if (LWCtot>config.LWClimit) 
    for (b=0;b<config.nbins;++b)
      if(Number(b) > tinym)
	if(AQ(b)>tinym2/Number(b) or AQinit(b)>tinym2/Number(b))
	  n2err2=max(n2err2,abs(AQinit(b)-AQ(b))/(AQinit(b)));

  if (config.compute_inorganic and LWCtot>config.LWClimit)
    {
      for (i=0;i<n;++i)
	{
	  if(surrogate[i].is_inorganic_precursor and surrogate[i].is_solid==false) 
	    if (surrogate[i].Ag1 > tinym or surrogate[i].Ag > tinym)
	      {
		++m;
		n2err+=pow((surrogate[i].Ag-surrogate[i].Ag1)/(surrogate[i].Ag1+tinym),2);
	      }

	  if (surrogate[i].is_organic==false and surrogate[i].is_inorganic_precursor==false and config.iH2O!=i and config.iHp!=i and config.iOHm!=i)
	    {	    
	      for (b=0;b<config.nbins;++b)
		if (surrogate[i].Aaq_bins(b) > tinym or surrogate[i].Aaq_bins_init(b)> tinym)
		  {
		    ++maq;
		    n2err_aq+=pow((surrogate[i].Aaq_bins(b)-
				   surrogate[i].Aaq_bins_init(b))/(tinym+surrogate[i].Aaq_bins_init(b)),2.0);
		  }
	    }
	}

      for (b=0;b<config.nbins;++b)
	{
	  n2err2=max(n2err2,abs(chp(b)-chp1(b))/chp1(b));
	}
    }

  if (m>0)
    n2err=pow(n2err/(m*m),0.5);
  if (maq>0)
    n2err=max(pow(n2err_aq/(maq*maq),0.5),n2err);
  n2err=max(n2err,n2err2);
  n2err=min(n2err,config.EPSER*R*R);
  n2err=max(n2err,config.EPSER/(R*R));
  deltat1=deltat1*pow(config.EPSER/n2err,0.5);
  deltat1=max(deltatmin,deltat1);
  deltat1=min(tend-t,deltat1);
}

void compute_diameters_ssh(model_config &config, vector<species>& surrogate,
			   Array<double, 1> &Vsol, Array<double, 1> &number,  Array<double, 1> &LWC,
			   double LWCtot)
{
  //compute the new diameter of particles due to the growth of the organic phases
  //number: number of particles in a bin
  
  int i,b,ilayer,iphase;
  double volume;
  double pi=3.14159265358979323846;
  int n=surrogate.size();

  if (config.fixed_density)
    for (b=0;b<config.nbins;++b)
      {
	if (number(b) > 0.0)
	  {
	    volume=(Vsol(b)+LWC(b)*1.0e-9/config.density)/number(b);
	    for (ilayer=0;ilayer<config.nlayer;++ilayer)
	      for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		for (i=0;i<n;i++)
		  if (surrogate[i].hydrophobic)
		    volume+=surrogate[i].Ap_layer_init(b,ilayer,iphase)*1.0e-9
		      /(number(b)*config.density);
	  
	    for (i=0;i<n;i++)
	      if (surrogate[i].hydrophilic and LWCtot>config.LWClimit and surrogate[i].is_inorganic_precursor==false)
		volume+=surrogate[i].Aaq_bins_init(b)*1.0e-9/(number(b)*config.density);
	  
	    config.diameters(b)=pow(6.0/pi*volume,1.0/3.0)*1.0e6;
	    if(config.diameters(b) < 1e-4) config.diameters(b) = 1e-4;
	  }
      }
  else
    for (b=0;b<config.nbins;++b)
      {
	if (number(b) > 0.0)
	  {
	    volume=(Vsol(b)+LWC(b)*1.0e-9/config.AQrho(b))/number(b);
	    for (ilayer=0;ilayer<config.nlayer;++ilayer)
	      for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		for (i=0;i<n;i++)
		  if (surrogate[i].hydrophobic)
		    volume+=surrogate[i].Ap_layer_init(b,ilayer,iphase)*1.0e-9
		      /(number(b)*surrogate[i].rho);
	  
	    for (i=0;i<n;i++)
	      if (surrogate[i].hydrophilic and LWCtot>config.LWClimit and surrogate[i].is_inorganic_precursor==false)
		volume+=surrogate[i].Aaq_bins_init(b)*1.0e-9/(number(b)*config.AQrho(b));
	  
	    config.diameters(b)=pow(6.0/pi*volume,1.0/3.0)*1.0e6;
	    if(config.diameters(b) < 1e-4) config.diameters(b) = 1e-4;
	  }
      }
}

void phase_repartition_ssh(model_config &config,vector<species>& surrogate, double &Temperature,
			   Array<double, 3> &MOinit, Array<double, 3> &MO, Array<double, 3> &MOW,
			   double factor)
{
  //compute the separation of phases due to saturation by minimizing the gibbs energy of organic
  //phases
  int b,ilayer,iphase,i; //,jphase;
  int n=surrogate.size();
  double sumX,tot,temp;
  double sum2;
  
  for (b=0;b<config.nbins;++b)		  
    for (ilayer=0;ilayer<config.nlayer;++ilayer)
      {
        for (i=0;i<n;++i)
          if (surrogate[i].hydrophobic and surrogate[i].kp_from_experiment==false 
              and surrogate[i].time(b,ilayer,0)>=config.tequilibrium)
            {
              sumX=0.0;
              tot=0.0;
              sum2=0.0;
              for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                {
                  sumX+=MOinit(b,ilayer,iphase)/
                    (MOW(b,ilayer,iphase)*surrogate[i].gamma_org_layer(b,ilayer,iphase));
                  tot+=surrogate[i].Ap_layer_init(b,ilayer,iphase);
                  sum2+=surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase);
                }
			  
              for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                {
                  if (sumX>0.0)
                    {
                      temp=tot/sumX*MOinit(b,ilayer,iphase)/
                        (MOW(b,ilayer,iphase)*
                         surrogate[i].gamma_org_layer(b,ilayer,iphase));
                      surrogate[i].Ap_layer_init(b,ilayer,iphase)=
                        factor*temp+(1.0-factor)*
                        surrogate[i].Ap_layer_init(b,ilayer,iphase);
                    }
                  else
                    cout << "problem" << endl;
                }
            }
		
        for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
          {
            MO(b,ilayer,iphase)=0.0;
            for (i=0;i<n;++i)
              if (surrogate[i].hydrophobic)
                MO(b,ilayer,iphase)+=surrogate[i].Ap_layer_init(b,ilayer,iphase);
	    MO(b,ilayer,iphase)=max(MO(b,ilayer,iphase),config.MOmin*config.Vlayer(ilayer));
          }
		  
      }
}

void number_org_phases_ssh(model_config &config,vector<species>& surrogate, double &Temperature,
			   Array<double, 3> &MOinit, Array<double, 3> &MOW)
{
  //compute the separation of phases due to saturation by minimizing the gibbs energy of organic
  //phases
  int b,ilayer,iphase,i,jphase;
  int n=surrogate.size();
  int nh; //,h;
  int index;
  //bool stable_system=false;
  Array<double, 3> MOinit2,MOW2,mean_activity,stability;
  Array<bool, 2> is_stable;
  Array<double, 2> gibbs_energy,gibbs_energy_old;
  double sumX,error,tot,temp;
  MOinit2.resize(config.nbins,config.nlayer,config.max_number_of_phases);
  MOW2.resize(config.nbins,config.nlayer,config.max_number_of_phases);
  mean_activity.resize(config.nbins,config.nlayer,config.max_number_of_phases);
  stability.resize(config.nbins,config.nlayer,config.max_number_of_phases);
  is_stable.resize(config.nbins,config.nlayer);
  gibbs_energy.resize(config.nbins,config.nlayer);
  gibbs_energy_old.resize(config.nbins,config.nlayer);
  double factor, error_old,error_old2,sumMO;
  bool mixing,separation,low_concentrations;
  int nphase_old; //,maxi;
  
  for (b=0;b<config.nbins;++b)		  
    for (ilayer=0;ilayer<config.nlayer;++ilayer)
      {
        nphase_old=0;
        mixing=true;
        separation=true;
        while(config.nphase(b,ilayer)!=nphase_old)
          {
            //save initial conditions
            for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
              {
                MOinit2(b,ilayer,iphase)=MOinit(b,ilayer,iphase);
                MOW2(b,ilayer,iphase)=MOW(b,ilayer,iphase);
              }
			
            nphase_old=config.nphase(b,ilayer);
            for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
              {
                sumX=0.0;
                for(i=0;i<n;++i)
                  if (surrogate[i].hydrophobic)
                    {
                      surrogate[i].Ap_layer_init0(b,ilayer,iphase)=
                        surrogate[i].Ap_layer_init(b,ilayer,iphase);
                      surrogate[i].gamma_org_layer0(b,ilayer,iphase)=
                        surrogate[i].gamma_org_layer(b,ilayer,iphase);
                      surrogate[i].Xinit(b,ilayer,iphase)=
                        surrogate[i].Ap_layer_init(b,ilayer,iphase)
                        /surrogate[i].MM;
                      sumX+=surrogate[i].Xinit(b,ilayer,iphase);
                    }
				
                for(i=0;i<n;++i)
                  if (surrogate[i].hydrophobic)
                    surrogate[i].Xinit(b,ilayer,iphase)/=sumX;
				
                mean_activity(b,ilayer,iphase)=0.0;
				
                for(i=0;i<n;++i)
                  if (surrogate[i].hydrophobic and i!=config.iH2O)
                    mean_activity(b,ilayer,iphase)+=surrogate[i].Xinit(b,ilayer,iphase)
                      /(1.0-surrogate[config.iH2O].Xinit(b,ilayer,iphase))*
                      surrogate[i].gamma_org_layer(b,ilayer,iphase);
				
                stability(b,ilayer,iphase)=0.0;
                for (i=0;i<n;++i)
                  if (surrogate[i].hydrophobic)
                    if (surrogate[i].Xinit(b,ilayer,iphase) >0.0
                        and surrogate[i].Ap_layer_init(b,ilayer,iphase)>0.0)
                      stability(b,ilayer,iphase)+=surrogate[i].Ap_layer_init(b,ilayer,iphase)*
                        log(surrogate[i].Xinit(b,ilayer,iphase)*
                            surrogate[i].gamma_org_layer(b,ilayer,iphase));
				
              }
			
            //compute the gibbs energy before phase separation
            gibbs_energy_old(b,ilayer)=0.0;
            for (i=0;i<n;++i)
              if (surrogate[i].hydrophobic)
                for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                  if (surrogate[i].Xinit(b,ilayer,iphase)*
                      surrogate[i].gamma_org_layer(b,ilayer,iphase)>0.0)
                    gibbs_energy_old(b,ilayer)+=surrogate[i].Ap_layer_init(b,ilayer,iphase)
                      /surrogate[i].MM*
                      log(surrogate[i].Xinit(b,ilayer,iphase)*
                          surrogate[i].gamma_org_layer(b,ilayer,iphase));
			
            //mixing	    
            if(config.nphase(b,ilayer)>1 and mixing)
              {
                config.nphase(b,ilayer)--;
                for(i=0;i<n;++i)
                  if (surrogate[i].hydrophobic)
                    {
                      for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                        surrogate[i].Ap_layer_init(b,ilayer,iphase)=
                          surrogate[i].Ap_layer_init0(b,ilayer,iphase);

                      iphase=config.nphase(b,ilayer);
                      surrogate[i].Ap_layer_init(b,ilayer,iphase-1)+=
                        surrogate[i].Ap_layer_init0(b,ilayer,iphase);
					  
                      surrogate[i].Ap_layer_init(b,ilayer,iphase)=0.0;
                      for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                        surrogate[i].Xinit(b,ilayer,iphase)=surrogate[i].Ap_layer_init(b,ilayer,iphase)
                          /surrogate[i].MM;
                    }
				
                for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                  {
                    sumX=0.0;
                    for(i=0;i<n;++i)
                      if (surrogate[i].hydrophobic)
                        sumX+=surrogate[i].Xinit(b,ilayer,iphase);
					
                    for(i=0;i<n;++i)
                      if (surrogate[i].hydrophobic)
                        surrogate[i].Xinit(b,ilayer,iphase)/=sumX;
                  }
				
                //compute activity coefficients in bin=b, layer=ilayer
                activity_coefficients_dyn_sat_ssh(config, surrogate, Temperature, MOW, b, ilayer);
				
                //compute the new gubbs energy
                for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                  {
                    MOinit(b,ilayer,iphase)=0.0;
                    for (i=0;i<n;++i)
                      if (surrogate[i].hydrophobic)
                        MOinit(b,ilayer,iphase)+=surrogate[i].Ap_layer_init(b,ilayer,iphase);
                  }
				
                for (iphase=config.nphase(b,ilayer);iphase<config.max_number_of_phases;++iphase)
                  MOinit(b,ilayer,iphase)=0.0;
				
                gibbs_energy(b,ilayer)=0.0;
                for (i=0;i<n;++i)
                  if (surrogate[i].hydrophobic)
                    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                      if (surrogate[i].Xinit(b,ilayer,iphase)*
                          surrogate[i].gamma_org_layer(b,ilayer,iphase)>0.0)
                        gibbs_energy(b,ilayer)+=surrogate[i].Ap_layer_init(b,ilayer,iphase)
                          /surrogate[i].MM*
                          log(surrogate[i].Xinit(b,ilayer,iphase)*
                              surrogate[i].gamma_org_layer(b,ilayer,iphase));

                mixing=false;
                if (gibbs_energy(b,ilayer)<=gibbs_energy_old(b,ilayer)
                    or MOinit2(b,ilayer,config.nphase(b,ilayer))<1.0e-10)
                  {
                    mixing=true;
                    separation=false;
                    //compute activity coefficients in bin=b, layer=ilayer
                    //activity_coefficients_dyn_sat(config, surrogate, Temperature, MOW, b, ilayer);
                  }
				
                if (mixing==false)
                  config.nphase(b,ilayer)++;
              }
            else
              mixing=false;
			
			
            if (separation)
              {	
                for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                  for(i=0;i<n;++i)
                    if (surrogate[i].hydrophobic)
                      surrogate[i].Ap_layer_init(b,ilayer,iphase)=
                        surrogate[i].Ap_layer_init0(b,ilayer,iphase);
			
                if (config.nphase(b,ilayer) < config.max_number_of_phases)
                  {
                    config.nphase(b,ilayer)++;
                    //estimate in which phase each compound should more present
                    jphase=0;
                    for (iphase=0;iphase<config.nphase(b,ilayer)-1;++iphase)
                      if (iphase!=jphase and stability(b,ilayer,iphase)>stability(b,ilayer,jphase))
                        jphase=iphase;
					
                    for (i=0;i<n;++i)
                      if (surrogate[i].hydrophobic)
                        {
                          if(surrogate[i].gamma_org_layer(b,ilayer,jphase) <
                             mean_activity(b,ilayer,jphase))
                            surrogate[i].Ap_layer_init(b,ilayer,config.nphase(b,ilayer)-1)=0.0;
                          else
                            {
                              surrogate[i].Ap_layer_init(b,ilayer,config.nphase(b,ilayer)-1)=
                                surrogate[i].Ap_layer_init(b,ilayer,jphase);
                              surrogate[i].Ap_layer_init(b,ilayer,jphase)=0.0;
                            }
                          /*
                            for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                            surrogate[i].Xinit(b,ilayer,iphase)=
                            surrogate[i].Ap_layer_init(b,ilayer,iphase)/surrogate[i].MM;*/
                        }
					
                    error=10.0;
                    index = 0;
                    factor=1.0;
                    nh=1;
                    error_old=-1.0;
                    error_old2=-1.0;
                    //method to make sure that the system converge even if there is strong variations
                    //   of the composition
                    while (error>factor*config.relative_precision and index<config.max_iter and nh<100)
                      {
                        ++index;
                        if(index>3)
                          if ((abs(error_old2-error)/error<1.0e-3)
                              and (abs(error_old-error)/error>1.0e-3))
                            {
                              ++nh;
                              factor=1.0/nh;
                              error_old=-1.0;
                              error_old2=-1.0;
                            }
						
                        if(index>2)
                          error_old2=error_old;
                        if(index>1)
                          error_old=error;
                        error=0.0;

                        for (i=0;i<n;++i)
                          if (surrogate[i].hydrophobic)
                            for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                              surrogate[i].Xinit(b,ilayer,iphase)=
                                surrogate[i].Ap_layer_init(b,ilayer,iphase)/surrogate[i].MM;

                        for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                          {
                            sumX=0.0;
                            for(i=0;i<n;++i)
                              if (surrogate[i].hydrophobic)
                                sumX+=surrogate[i].Xinit(b,ilayer,iphase);
							
                            for(i=0;i<n;++i)
                              if (surrogate[i].hydrophobic)
                                surrogate[i].Xinit(b,ilayer,iphase)/=sumX;
                          }
						
                        for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                          {
                            MOinit(b,ilayer,iphase)=0.0;
                            for (i=0;i<n;++i)
                              if (surrogate[i].hydrophobic)
                                MOinit(b,ilayer,iphase)
                                  +=surrogate[i].Ap_layer_init(b,ilayer,iphase);
			    MOinit(b,ilayer,iphase)=max(MOinit(b,ilayer,iphase),config.MOmin*config.Vlayer(ilayer));
                          }
					
                        //compute activity coefficients in bin=b, layer=ilayer
                        activity_coefficients_dyn_sat_ssh(config, surrogate, Temperature, MOW, b, ilayer);

                        //estimate the repartition of organic compounds between phases 
                        for (i=0;i<n;++i)
                          if ((surrogate[i].is_organic or i==config.iH2O)
                              and surrogate[i].hydrophobic)
                            if(surrogate[i].kp_from_experiment==false)
                              {
                                sumX=0.0;
                                tot=0.0;
                                for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                                  {
                                    sumX+=MOinit(b,ilayer,iphase)/
                                      (MOW(b,ilayer,iphase)*
                                       surrogate[i].gamma_org_layer(b,ilayer,iphase));
                                    tot+=surrogate[i].Ap_layer_init(b,ilayer,iphase);
                                  }
								
                                for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                                  {
                                    if (sumX>0.0)
                                      {
                                        temp=tot/sumX*MOinit(b,ilayer,iphase)/
                                          (MOW(b,ilayer,iphase)*
                                           surrogate[i].gamma_org_layer(b,ilayer,iphase));
                                        if (temp>1.0e-10 or
                                            surrogate[i].Ap_layer_init(b,ilayer,iphase)>1.0e-10)
                                          error=max(error,(temp-
                                                           surrogate[i].Ap_layer_init(b,ilayer,iphase))/
                                                    surrogate[i].Ap_layer_init(b,ilayer,iphase));
                                        surrogate[i].Ap_layer_init(b,ilayer,iphase)=
                                          1.0/nh*temp+(1.0-1.0/nh)*
                                          surrogate[i].Ap_layer_init(b,ilayer,iphase);
                                      }
                                    else
                                      cout << "problem sep " << config.nphase(b,ilayer) << " " << MOinit << " " << MOW << " " << endl;
                                  }
                              }
                      }
					
                    //compute the new gubbs energy
                    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                      {
                        MOinit(b,ilayer,iphase)=0.0;
                        for (i=0;i<n;++i)
                          if (surrogate[i].hydrophobic)
                            MOinit(b,ilayer,iphase)+=surrogate[i].Ap_layer_init(b,ilayer,iphase);
                      }
					
                    for (iphase=config.nphase(b,ilayer);iphase<config.max_number_of_phases;++iphase)
                      MOinit(b,ilayer,iphase)=0.0;
				
                    gibbs_energy(b,ilayer)=0.0;
                    for (i=0;i<n;++i)
                      if (surrogate[i].hydrophobic)
                        for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                          if (surrogate[i].Xinit(b,ilayer,iphase)*
                              surrogate[i].gamma_org_layer(b,ilayer,iphase)>0.0)
                            gibbs_energy(b,ilayer)+=surrogate[i].Ap_layer_init(b,ilayer,iphase)
                              /surrogate[i].MM*
                              log(surrogate[i].Xinit(b,ilayer,iphase)*
                                  surrogate[i].gamma_org_layer(b,ilayer,iphase));

                    sumMO=0.0;
                    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                      sumMO+=MOinit(b,ilayer,iphase);

                    low_concentrations=false;
                    if (sumMO>0.0)
                      for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                        if (MOinit(b,ilayer,iphase)/sumMO<0.01)
                          low_concentrations=true;

                    separation=false;
                    if (gibbs_energy(b,ilayer)<gibbs_energy_old(b,ilayer)
                        and MOinit(b,ilayer,config.nphase(b,ilayer)-1)>1.0e-10
                        and low_concentrations==false)
                      {
                        separation=true;
                        mixing=false;
                      }
					
                    if (separation==false)
                      config.nphase(b,ilayer)--;
                  }
              }
            else
              separation=false;

			
            if (separation==false and mixing==false)
              {
                for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                  {
                    MOinit(b,ilayer,iphase)=MOinit2(b,ilayer,iphase);
                    MOW(b,ilayer,iphase)=MOW2(b,ilayer,iphase);
                  }
				
                for (iphase=config.nphase(b,ilayer);iphase<config.max_number_of_phases;++iphase)
                  {
                    MOinit(b,ilayer,iphase)=0.0;
                    MOW(b,ilayer,iphase)=200.0;
                  }
				
                for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                  for(i=0;i<n;++i)
                    if (surrogate[i].hydrophobic)
                      {
                        surrogate[i].Ap_layer_init(b,ilayer,iphase)=
                          surrogate[i].Ap_layer_init0(b,ilayer,iphase);
                        surrogate[i].gamma_org_layer(b,ilayer,iphase)=
                          surrogate[i].gamma_org_layer0(b,ilayer,iphase);
                      }
				
                for (iphase=config.nphase(b,ilayer);iphase<config.max_number_of_phases;++iphase)
                  for(i=0;i<n;++i)
                    if (surrogate[i].hydrophobic)
                      {
                        surrogate[i].Ap_layer_init(b,ilayer,iphase)=0.0;
                        surrogate[i].gamma_org_layer(b,ilayer,iphase)=1.0;
                      }
              }
          } 
      }
}


