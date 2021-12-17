//!!-----------------------------------------------------------------------
//!!     Copyright (C) 2019 CEREA (ENPC) - INERIS
//!!     SSH-aerosol is distributed under the GNU General Public License v3
//!!-----------------------------------------------------------------------

using namespace ssh_soap;

void compute_kp_org_ssh(model_config &config, vector<species>& surrogate,
			Array <double, 3> &MOinit, double &Temperature, Array<double, 3> &MOW)
{
  //compute partitioning constant of the organic phase by taking into account activity
  //coefficients and the kelvin effect
  double temp1,temp2,maxi,MOWsurf;
  double kelvin_effect=1.e0;
  int ilayer,i,b,iphase,jphase;
  int n=surrogate.size();
  for (b=0;b<config.nbins;++b)
    {      
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
			surrogate[i].Kp(b,ilayer,iphase)=pow(10.0,5)*200.0/
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
}

void compute_kp_aq_ssh(model_config &config, vector<species>& surrogate,
		       double &Temperature, Array <double, 1> &ionic,
		       Array <double, 1> &chp,Array<double, 1> &MMaq)
{
  //compute partitioning constant of the aqueous phase by taking into account activity
  //coefficients and the kelvin effect
  double kelvin_effect=1.;
  int b,i;
  int n=surrogate.size();
  double fion1,fion2;
  double R=8.314; //ideal gas constant (J/K/mol)
  double deltaH_over_RT0,deltaCp0_over_R;
  double T0=298.15;

  for (b=0;b<config.nbins;++b)      
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
	  if(surrogate[i].is_organic and surrogate[i].hydrophilic)
	    {	  
	      surrogate[i].gamma_LR=surrogate[i].LR(b);
	      if (surrogate[i].nonvolatile)
		surrogate[i].Kaq(b)=1000.0/surrogate[i].gamma_aq_bins(b)*18.0/MMaq(b);
	      else
		{
		  surrogate[i].Kaq(b)=
		    surrogate[i].Kp_eff_aqrealdyn_ssh(config, Temperature, ionic(b), chp(b),surrogate[config.iHp].LR(b),
						      surrogate[config.iHp].SRMR(b),MMaq(b),fion1,fion2,b)
		    /surrogate[i].gamma_aq_bins(b);
		
		  surrogate[i].fion1=fion1;
		  surrogate[i].fion2=fion2;
		}
	     
	      surrogate[i].Kaq(b)/=kelvin_effect;
	    }
	  else if (surrogate[i].is_inorganic_precursor and config.compute_inorganic and surrogate[i].is_solid==false)
	    {
	      if (surrogate[i].name=="H2SO4")
		surrogate[i].Kaq(b)=1.0e10;          
	      else if (surrogate[i].name=="NH3")
		{
		  deltaH_over_RT0=13.79;
		  deltaCp0_over_R=-5.39;
		  surrogate[i].Kaq(b)=
		    surrogate[i].Henry*exp(-deltaH_over_RT0*(T0/Temperature-1.0)-deltaCp0_over_R*(1+log(T0/Temperature)-T0/Temperature))
		    *R*Temperature/(1000.*1.0e6*1.013e5)
		    *(1.0+surrogate[i].Kequilibrium_ssh(Temperature)*chp(b)*surrogate[config.iHp].gamma_aq_bins(b)/surrogate[config.iNH4p].gamma_aq_bins(b));
		  surrogate[i].dKaq(b)=
		    surrogate[i].Henry*exp(-deltaH_over_RT0*(T0/Temperature-1.0)-deltaCp0_over_R*(1+log(T0/Temperature)-T0/Temperature))
		    *R*Temperature/(1000.*1.0e6*1.013e5)
		    *(surrogate[i].Kequilibrium_ssh(Temperature)*surrogate[config.iHp].gamma_aq_bins(b)/surrogate[config.iNH4p].gamma_aq_bins(b));
		  if (config.compute_kelvin_effect) //compute the kelvin effect
		    {
		      surrogate[i].Kaq(b)/=kelvin_effect;
		      surrogate[i].dKaq(b)/=kelvin_effect;
		    }
		}
	      else if (surrogate[i].name=="HNO3")
		{
		  deltaH_over_RT0=29.17;
		  deltaCp0_over_R=16.83;
		  surrogate[i].Kaq(b)=
		    surrogate[i].Henry*exp(-deltaH_over_RT0*(T0/Temperature-1.0)-deltaCp0_over_R*(1+log(T0/Temperature)-T0/Temperature))
		    *R*Temperature/(1000.*1.0e6*1.013e5)
		    *(1.0+surrogate[i].Kequilibrium_ssh(Temperature)/(chp(b)*surrogate[config.iHp].gamma_aq_bins(b)*surrogate[config.iNO3m].gamma_aq_bins(b)));
		  surrogate[i].dKaq(b)=-surrogate[i].Henry*exp(-deltaH_over_RT0*(T0/Temperature-1.0)-deltaCp0_over_R*(1+log(T0/Temperature)-T0/Temperature))
		    *R*Temperature/(1000.*1.0e6*1.013e5)
		    *(surrogate[i].Kequilibrium_ssh(Temperature)/(pow(chp(b),2)*surrogate[config.iHp].gamma_aq_bins(b)*surrogate[config.iClm].gamma_aq_bins(b)));
		  if (config.compute_kelvin_effect) //compute the kelvin effect
		    {
		      surrogate[i].Kaq(b)/=kelvin_effect;
		      surrogate[i].dKaq(b)/=kelvin_effect;
		    }
		}
	      else if (surrogate[i].name=="HCl")
		{
		  
		  deltaH_over_RT0=30.20;
		  deltaCp0_over_R=19.91;
		  surrogate[i].Kaq(b)=surrogate[i].Henry*exp(-deltaH_over_RT0*(T0/Temperature-1.0)-deltaCp0_over_R*(1+log(T0/Temperature)-T0/Temperature))
		    *R*Temperature/(1000.*1.0e6*1.013e5)
		    *(1.0+surrogate[i].Kequilibrium_ssh(Temperature)/(chp(b)*surrogate[config.iHp].gamma_aq_bins(b)*surrogate[config.iClm].gamma_aq_bins(b)));
		  surrogate[i].dKaq(b)=-surrogate[i].Henry*exp(-deltaH_over_RT0*(T0/Temperature-1.0)-deltaCp0_over_R*(1+log(T0/Temperature)-T0/Temperature))
		    *R*Temperature/(1000.*1.0e6*1.013e5)
		    *(surrogate[i].Kequilibrium_ssh(Temperature)/(pow(chp(b),2)*surrogate[config.iHp].gamma_aq_bins(b)*surrogate[config.iClm].gamma_aq_bins(b)));
		  if (config.compute_kelvin_effect) //compute the kelvin effect
		    {
		      surrogate[i].Kaq(b)/=kelvin_effect;
		      surrogate[i].dKaq(b)/=kelvin_effect;
		    }
		}
	    }
	  else if (i==config.iH2O)	    	    
	    surrogate[i].Kaq(b)=surrogate[i].kaqi/(MMaq(b)*surrogate[i].gamma_aq_bins(b))/kelvin_effect;			    	    
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
		    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		      surrogate[i].time(b,ilayer,iphase)=0.0;       		  
		}
	    }

	else
	  for (b=0;b<config.nbins;++b)
	    for (ilayer=0;ilayer<config.nlayer;++ilayer)
	      for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		surrogate[i].time(b,ilayer,iphase)=0.0;

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
		  for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		    surrogate[i].time(b,ilayer,iphase)=0.0;
		
		
	      }
    }  

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
				Array<double, 3> &MOinit)
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
            if (surrogate[i].name=="H2SO4")
              for (b=0;b<config.nbins;++b)
                surrogate[i].Aaq+=surrogate[config.iHSO4m].Aaq_bins_init(b)/surrogate[config.iHSO4m].MM*surrogate[i].MM
                  +surrogate[config.iSO4mm].Aaq_bins_init(b)/surrogate[config.iSO4mm].MM*surrogate[i].MM;
            else if (surrogate[i].name=="NH3")
              {
                for (b=0;b<config.nbins;++b)
                  surrogate[i].Aaq+=surrogate[config.iNH4p].Aaq_bins_init(b)/surrogate[config.iNH4p].MM*surrogate[i].MM;
              }
            else if (surrogate[i].name=="HNO3")
              {
                for (b=0;b<config.nbins;++b)
                  surrogate[i].Aaq+=surrogate[config.iNO3m].Aaq_bins_init(b)/surrogate[config.iNO3m].MM*surrogate[i].MM;
              }
            else if (surrogate[i].name=="HCl")
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
                if (surrogate[i].Aaq_bins_init(b)>0.0)
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
			
			AQ2=max(AQinit(b),(-b1+pow(delta,0.5))/(2.0*a1));
			f1=pow((AQ2+sumMO)/(AQinit(b)+sumMO),1.0/3.0);			
			surrogate[i].time_aq(b)=
			  (AQ2+sumMO)/AQ2*
			  (Kaq*AQ2*surrogate[i].tau_air(b)/f1)
			  /(1.0+Kaq*AQ2/surrogate[i].Aaq_bins_init(b)*surrogate[i].Aaq);
		      }

		  }
                else 
                  surrogate[i].time_aq(b)=0.0*config.tequilibrium;
              else
                surrogate[i].time_aq(b)=2.0*config.tequilibrium;
            }
        else if (surrogate[i].is_inorganic_precursor and config.compute_inorganic and surrogate[i].is_solid==false)
          for (b=0;b<config.nbins;++b)
            {
              conc_aq=0.0;
              if (surrogate[i].name=="H2SO4")
                conc_aq+=surrogate[config.iHSO4m].Aaq_bins_init(b)/surrogate[config.iHSO4m].MM*surrogate[i].MM
                  +surrogate[config.iSO4mm].Aaq_bins_init(b)/surrogate[config.iSO4mm].MM*surrogate[i].MM;
              else if (surrogate[i].name=="NH3")
                conc_aq+=surrogate[config.iNH4p].Aaq_bins_init(b)/surrogate[config.iNH4p].MM*surrogate[i].MM;
              else if (surrogate[i].name=="HNO3")
                conc_aq+=surrogate[config.iNO3m].Aaq_bins_init(b)/surrogate[config.iNO3m].MM*surrogate[i].MM;
              else if (surrogate[i].name=="HCl")
                conc_aq+=surrogate[config.iClm].Aaq_bins_init(b)/surrogate[config.iClm].MM*surrogate[i].MM;
              Kaq=surrogate[i].Kaq(b);
              if(Kaq<config.kp_low_volatility)
                if (conc_aq>0.0)
                  surrogate[i].time_aq(b)=(Kaq*AQinit(b)*surrogate[i].tau_air(b))
                    /(1.0+Kaq*AQinit(b)/conc_aq*surrogate[i].Aaq);
                else 
                  surrogate[i].time_aq(b)=0.0*config.tequilibrium;
              else
                surrogate[i].time_aq(b)=2.0*config.tequilibrium;

	      	      
              if (surrogate[i].name=="HCl" or surrogate[i].name=="HNO3" or surrogate[i].name=="NH3")		
		max_time_inorg(b)=max(max_time_inorg(b),surrogate[i].time_aq(b));	        		

              if (surrogate[i].name=="H2SO4")
                surrogate[i].time_aq(b)=2.0*config.tequilibrium;
            } 
    }
  else
    for (i=0;i<n;++i)
      if(surrogate[i].hydrophilic)
        for (b=0;b<config.nbins;++b)
          surrogate[i].time_aq(b)=0.0;

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
    if (surrogate[i].aq_type=="monoacid" or surrogate[i].aq_type=="diacid")
      for (b=0;b<config.nbins;b++)         	
	max_time_inorg(b)=max(max_time_inorg(b),(surrogate[i].fion1+surrogate[i].fion2)*surrogate[i].time_aq(b));       

  for (b=0;b<config.nbins;++b)
    surrogate[config.iH2O].time_aq(b)=-10.0; //H2O is forced to be at equilibrium
    
  for (i=0;i<n;i++)
    if (surrogate[i].aq_type=="monoacid" or surrogate[i].aq_type=="diacid")
      for (b=0;b<config.nbins;b++)         
	{	  
	  surrogate[i].time_aq(b)=(1.0-surrogate[i].fion1-surrogate[i].fion2)*surrogate[i].time_aq(b)
	    +(surrogate[i].fion1+surrogate[i].fion2)*max_time_inorg(b); //max(max_time_inorg(b),surrogate[i].time_aq(b));
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
			   Array<double, 1> &AQinit, Array<double, 1> &MMaq, 
			   Array<double, 1> &chp, double &deltat, double &tiny,int index)
{    
  int n=surrogate.size();
  int i,j,b,ilayer,iphase,jion;
  double gamma=1.+0.5*pow(2,0.5);  
  double sum_flux_gas=0.0;
  double Xmono=0.0;          
  double XH2O;
  double xmin=1.0e-10;

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
                  if (surrogate[i].is_monomer)
                    {
                      j=surrogate[i].ioligo;                       
                      double Keq2=surrogate[j].MM/surrogate[i].MM*pow(config.Keq_oligo,surrogate[i].moligo-1)*pow(max(Xmono,xmin),surrogate[i].moligo-2)/pow(max(XH2O,xmin),surrogate[i].moligo-1);                      
                      double flux=config.koligo*(surrogate[i].gamma_org_layer(b,ilayer,iphase)*surrogate[i].Ap_layer_init(b,ilayer,iphase)*Xmono-surrogate[j].gamma_org_layer(b,ilayer,iphase)*surrogate[j].Ap_layer_init(b,ilayer,iphase)/Keq2)*deltat;                                                       
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
            }
  
      for (b=0;b<config.nbins;++b)
        {        
          double Xmonoaq=0.0;
          for (i=0;i<n;i++)
            if (surrogate[i].is_monomer and surrogate[i].is_organic and surrogate[i].hydrophilic)  
              Xmonoaq+=surrogate[i].gamma_aq_bins(b)*surrogate[i].GAMMAinf*surrogate[i].Aaq_bins_init(b)*MMaq(b)/surrogate[i].MM/max(AQinit(b),1.0e-10);          

          double conc_org=0.0;//LWC(b);
          for (i=0;i<n;++i)
            if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophilic)
              conc_org+=surrogate[i].Aaq_bins_init(b);

          conc_org=max(conc_org,1.e-5*AQinit(b)); //config.MOmin);
          //conc_org=max(conc_org,config.MOmin);
      
          double XH2Oaq=surrogate[config.iH2O].gamma_aq_bins(b)*surrogate[config.iH2O].Aaq_bins_init(b)*MMaq(b)/surrogate[config.iH2O].MM/max(AQinit(b),1.0e-10);
          double xmin=1.0e-10;

          for (i=0;i<n;++i)
            {
              if (surrogate[i].is_organic)
                if (surrogate[i].is_organic and surrogate[i].hydrophilic) 
                  if (surrogate[i].is_monomer)
                    {
                      j=surrogate[i].ioligo;                                        
                      double Kaq2=surrogate[j].MM/surrogate[i].MM*pow(config.Keq_oligo,surrogate[i].moligo-1)*pow(max(Xmonoaq,xmin),surrogate[i].moligo-2)/pow(max(XH2Oaq,xmin),surrogate[i].moligo-1);                                            
                      double flux=config.koligo*(surrogate[i].gamma_aq_bins(b)*surrogate[i].GAMMAinf*surrogate[i].Aaq_bins_init(b)*Xmonoaq-surrogate[j].gamma_aq_bins(b)*surrogate[j].GAMMAinf*surrogate[j].Aaq_bins_init(b)/Kaq2)*deltat;                          
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
		      /*
			if (surrogate[i].name=="BiDER")
			{
			//cout << surrogate[i].name << " " << b << " " << flux << " " << config.koligo << " " << surrogate[i].gamma_aq_bins(b)*surrogate[i].GAMMAinf << " " << Xmonoaq <<  " " << config.koligo*surrogate[i].gamma_aq_bins(b)*surrogate[i].GAMMAinf*surrogate[i].Aaq_bins_init(b)*Xmonoaq << " " << deltat << " " << config.koligo*surrogate[j].gamma_aq_bins(b)*surrogate[j].GAMMAinf*surrogate[j].Aaq_bins_init(b)/Kaq2*deltat << endl;                          
			}*/
			
                    }

              if (surrogate[i].rion and surrogate[i].Aaq_bins_init(b)>0.0 and config.compute_inorganic)
                for (jion=0;jion<surrogate[i].nion;jion++)
                  {            
                    double molality=surrogate[surrogate[i].iion(jion)].gamma_aq_bins(b)*surrogate[surrogate[i].iion(jion)].Aaq_bins_init(b)/surrogate[i].MM/conc_org*1000.0;
                    double flux=surrogate[i].kion(jion)*molality*deltat*surrogate[i].GAMMAinf*surrogate[i].gamma_aq_bins(b)*surrogate[i].Aaq_bins_init(b);
                    double flux2=0.0;
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

                    if (surrogate[i].rion_catalyzed(jion)==false and molality>0.)
                      {
                        if (surrogate[i].iion(jion)==config.iHSO4m or surrogate[i].iion(jion)==config.iSO4mm)
                          {
                            double totsulf=surrogate[config.iHSO4m].Aaq_bins_init(b)/surrogate[config.iHSO4m].MM+surrogate[config.iSO4mm].Aaq_bins_init(b)/surrogate[config.iSO4mm].MM;
                            surrogate[config.iHSO4m].flux_chem_aq(b,index)-=surrogate[config.iHSO4m].Aaq_bins_init(b)/totsulf*flux/surrogate[i].MM;
                            surrogate[config.iSO4mm].flux_chem_aq(b,index)-=surrogate[config.iSO4mm].Aaq_bins_init(b)/totsulf*flux/surrogate[i].MM;
                            flux2-=surrogate[config.iHSO4m].flux_chem_aq(b,index)*surrogate[config.iSO4mm].MM/surrogate[config.iHSO4m].MM+surrogate[config.iSO4mm].flux_chem_aq(b,index);
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
                            flux2+=flux/surrogate[i].MM*surrogate[config.iNO3m].MM;
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
                            flux2+=flux/surrogate[i].MM*surrogate[config.iClm].MM;
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

                            flux2+=flux/surrogate[i].MM*surrogate[config.iNH4p].MM;
                          }
                      }                
                    surrogate[j].flux_chem_aq(b,index)+=fac2*(flux+flux2);		
                    surrogate[j].flux_chem_gas(index)+=(1.0-fac2)*(flux+flux2);
                  }
            }    
        }

  
      for (i=0;i<n;++i)
        {
      
          double sum_flux_gas=0.;
          if (surrogate[i].Ag>tiny)
            sum_flux_gas=surrogate[i].flux_chem_gas(index)/surrogate[i].Ag;

          for (b=0;b<config.nbins;++b)
            {
              for (ilayer=0;ilayer<config.nlayer;ilayer++)
                for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)       
                  if (surrogate[i].time(b,ilayer,iphase)>=config.tequilibrium)
                    {                        
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
                      if (surrogate[i].flux_chem(b,ilayer,iphase,index)/deltat<0.0 and surrogate[i].Ap_layer_init(b,ilayer,iphase)>tiny)
                        surrogate[i].Jdn(b,ilayer,iphase)=surrogate[i].flux_chem(b,ilayer,iphase,index)/deltat/
                          surrogate[i].Ap_layer_init(b,ilayer,iphase); 
                      else
                        surrogate[i].Jdn(b,ilayer,iphase)=0.0;                                    
                    }    
          
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
                  if (surrogate[i].is_monomer)
                    {
                      j=surrogate[i].ioligo;                       
                      double Keq2=surrogate[j].MM/surrogate[i].MM*pow(config.Keq_oligo,surrogate[i].moligo-1)*pow(max(Xmono,xmin),surrogate[i].moligo-2)/pow(max(XH2O,xmin),surrogate[i].moligo-1);                      
                      double flux=config.koligo*(surrogate[i].gamma_org_layer(b,ilayer,iphase)*surrogate[i].Ap_layer_init(b,ilayer,iphase)*Xmono-surrogate[j].gamma_org_layer(b,ilayer,iphase)*surrogate[j].Ap_layer_init(b,ilayer,iphase)/Keq2)*deltat;                                                       
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
            }
  
      for (b=0;b<config.nbins;++b)
        {        
          double Xmonoaq=0.0;
          for (i=0;i<n;i++)
            if (surrogate[i].is_monomer and surrogate[i].is_organic and surrogate[i].hydrophilic)  
              Xmonoaq+=surrogate[i].gamma_aq_bins(b)*surrogate[i].GAMMAinf*surrogate[i].Aaq_bins_init(b)*MMaq(b)/surrogate[i].MM/max(AQinit(b),1.0e-10);          

          double conc_org=0.0;//LWC(b);
          for (i=0;i<n;++i)
            if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophilic)
              conc_org+=surrogate[i].Aaq_bins_init(b);

          conc_org=max(conc_org,1.e-5*AQinit(b)); //config.MOmin);
          //conc_org=max(conc_org,config.MOmin);
      
          double XH2Oaq=surrogate[config.iH2O].gamma_aq_bins(b)*surrogate[config.iH2O].Aaq_bins_init(b)*MMaq(b)/surrogate[config.iH2O].MM/max(AQinit(b),1.0e-10);
          double xmin=1.0e-10;

          for (i=0;i<n;++i)
            if (surrogate[i].is_organic)
              {
                if (surrogate[i].is_organic and surrogate[i].hydrophilic) 
                  if (surrogate[i].is_monomer)
                    {
                      j=surrogate[i].ioligo;                                        
                      double Kaq2=surrogate[j].MM/surrogate[i].MM*pow(config.Keq_oligo,surrogate[i].moligo-1)*pow(max(Xmonoaq,xmin),surrogate[i].moligo-2)/pow(max(XH2Oaq,xmin),surrogate[i].moligo-1);                                            
                      double flux=config.koligo*(surrogate[i].gamma_aq_bins(b)*surrogate[i].GAMMAinf*surrogate[i].Aaq_bins_init(b)*Xmonoaq-surrogate[j].gamma_aq_bins(b)*surrogate[j].GAMMAinf*surrogate[j].Aaq_bins_init(b)/Kaq2)*deltat;                          
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
                  for (jion=0;jion<surrogate[i].nion;jion++)
                    {            
                      double molality=surrogate[surrogate[i].iion(jion)].gamma_aq_bins(b)*surrogate[surrogate[i].iion(jion)].Aaq_bins_init(b)/surrogate[i].MM/conc_org*1000.0;
                      double flux=surrogate[i].kion(jion)*molality*deltat*surrogate[i].GAMMAinf*surrogate[i].gamma_aq_bins(b)*surrogate[i].Aaq_bins_init(b);
                      double flux2=0.0;
                      double fac=1.0;
                      double fac2=1.0;
                      j=surrogate[i].iproduct(jion);
                      if (surrogate[i].time_aq(b)<config.tequilibrium)                                                  
                        if (surrogate[i].Ag+surrogate[i].Aaq_bins_init(b)>0.0)
                          fac=surrogate[i].Aaq_bins_init(b)/(surrogate[i].Ag+surrogate[i].Aaq_bins_init(b)); 

                      if (surrogate[j].time_aq(b)<config.tequilibrium)                        
                        if (surrogate[j].Ag+surrogate[j].Aaq_bins_init(b)>0.0)
                          fac2=surrogate[j].Aaq_bins_init(b)/(surrogate[j].Ag+surrogate[j].Aaq_bins_init(b));                                 

                      if (surrogate[i].rion_catalyzed(jion)==false and molality>0.)
                        {
                          if (surrogate[i].iion(jion)==config.iHSO4m or surrogate[i].iion(jion)==config.iSO4mm)
                            {
                              flux=flux/(1.0-gamma*surrogate[i].Jdn_aq(b,index)*deltat-gamma*(1.0-fac)*surrogate[i].Jdn_gas(index)*deltat-gamma*surrogate[i].Jdn_tot*deltat
                                         -gamma*surrogate[config.iHSO4m].Jdn_aq(b,index)/surrogate[config.iHSO4m].MM*surrogate[i].MM*deltat
                                         -gamma*surrogate[config.iSO4mm].Jdn_aq(b,index)/surrogate[config.iSO4mm].MM*surrogate[i].MM*deltat);
                              double totsulf=surrogate[config.iHSO4m].Aaq_bins_init(b)/surrogate[config.iHSO4m].MM+surrogate[config.iSO4mm].Aaq_bins_init(b)/surrogate[config.iSO4mm].MM;
                              surrogate[config.iHSO4m].flux_chem_aq(b,index)-=surrogate[config.iHSO4m].Aaq_bins_init(b)/totsulf*flux/surrogate[i].MM;
                              surrogate[config.iSO4mm].flux_chem_aq(b,index)-=surrogate[config.iSO4mm].Aaq_bins_init(b)/totsulf*flux/surrogate[i].MM;
                              flux2-=surrogate[config.iHSO4m].flux_chem_aq(b,index)*surrogate[config.iSO4mm].MM/surrogate[config.iHSO4m].MM+surrogate[config.iSO4mm].flux_chem_aq(b,index);                                                    
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
                              flux2+=flux/surrogate[i].MM*surrogate[config.iNO3m].MM;
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
                              flux2+=flux/surrogate[i].MM*surrogate[config.iClm].MM;
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

                              flux2+=flux/surrogate[i].MM*surrogate[config.iNH4p].MM;
                            }
                        }                

                      surrogate[i].flux_chem_aq(b,index)+=-fac*flux;
                      surrogate[i].flux_chem_gas(index)+=-(1.0-fac)*flux;
                      surrogate[j].flux_chem_aq(b,index)+=fac2*(flux+flux2);		
                      surrogate[j].flux_chem_gas(index)+=(1.0-fac2)*(flux+flux2);
                    }
    
              }
        } 
    }    
}



void prodloss_chem_ssh(model_config &config, vector<species>& surrogate,
			   Array<double, 3> &MOinit, Array<double, 3> &MOW, 
			   Array<double, 1> &AQinit, Array<double, 1> &MMaq, 
			   Array<double, 1> &chp, double &tiny,int index)
{    
  int n=surrogate.size();
  int i,j,b,ilayer,iphase,jion; 
  double Xmono=0.0;          
  double XH2O;
  double xmin=1.0e-10;

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
                  if (surrogate[i].is_monomer)
                    {
                      j=surrogate[i].ioligo;                       
                      double Keq2=surrogate[j].MM/surrogate[i].MM*pow(config.Keq_oligo,surrogate[i].moligo-1)*pow(max(Xmono,xmin),surrogate[i].moligo-2)/pow(max(XH2O,xmin),surrogate[i].moligo-1);                      
                      double flux1=config.koligo*surrogate[i].gamma_org_layer(b,ilayer,iphase)*surrogate[i].Ap_layer_init(b,ilayer,iphase)*Xmono;
		      double flux2=config.koligo*surrogate[j].gamma_org_layer(b,ilayer,iphase)*surrogate[j].Ap_layer_init(b,ilayer,iphase)/Keq2;
		      
                      surrogate[i].kprod(b,ilayer,iphase)+=flux2;
		      surrogate[j].kprod(b,ilayer,iphase)+=flux1;
		      surrogate[i].kloss(b,ilayer,iphase)+=config.koligo*surrogate[i].gamma_org_layer(b,ilayer,iphase)*Xmono;
		      surrogate[j].kloss(b,ilayer,iphase)+=config.koligo*surrogate[j].gamma_org_layer(b,ilayer,iphase)/Keq2;                      
                    }
            }
  
      for (b=0;b<config.nbins;++b)
        {        
          double Xmonoaq=0.0;
          for (i=0;i<n;i++)
            if (surrogate[i].is_monomer and surrogate[i].is_organic and surrogate[i].hydrophilic)  
              Xmonoaq+=surrogate[i].gamma_aq_bins(b)*surrogate[i].GAMMAinf*surrogate[i].Aaq_bins_init(b)*MMaq(b)/surrogate[i].MM/max(AQinit(b),1.0e-10);          

          double conc_org=0.0;//LWC(b);
          for (i=0;i<n;++i)
            if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophilic)
              conc_org+=surrogate[i].Aaq_bins_init(b);

          conc_org=max(conc_org,1.e-5*AQinit(b)); //config.MOmin);
          //conc_org=max(conc_org,config.MOmin);
      
          double XH2Oaq=surrogate[config.iH2O].gamma_aq_bins(b)*surrogate[config.iH2O].Aaq_bins_init(b)*MMaq(b)/surrogate[config.iH2O].MM/max(AQinit(b),1.0e-10);
          double xmin=1.0e-10;

          for (i=0;i<n;++i)
            {
              if (surrogate[i].is_organic)
                if (surrogate[i].is_organic and surrogate[i].hydrophilic) 
                  if (surrogate[i].is_monomer)
                    {
                      j=surrogate[i].ioligo;                                        
                      double Kaq2=surrogate[j].MM/surrogate[i].MM*pow(config.Keq_oligo,surrogate[i].moligo-1)*pow(max(Xmonoaq,xmin),surrogate[i].moligo-2)/pow(max(XH2Oaq,xmin),surrogate[i].moligo-1);                                            
                      double flux1=config.koligo*surrogate[i].gamma_aq_bins(b)*surrogate[i].GAMMAinf*surrogate[i].Aaq_bins_init(b)*Xmonoaq;
		      double flux2=config.koligo*surrogate[j].gamma_aq_bins(b)*surrogate[j].GAMMAinf*surrogate[j].Aaq_bins_init(b)/Kaq2;

		      surrogate[i].kprod_aq(b)+=flux2;
		      surrogate[j].kprod_aq(b)+=flux1;
		      surrogate[i].kloss_aq(b)+=config.koligo*surrogate[i].gamma_aq_bins(b)*surrogate[i].GAMMAinf*Xmonoaq;
		      surrogate[j].kloss_aq(b)+=config.koligo*surrogate[j].gamma_aq_bins(b)*surrogate[j].GAMMAinf/Kaq2;                                            
                    }

              if (surrogate[i].rion and surrogate[i].Aaq_bins_init(b)>0.0 and config.compute_inorganic)
                for (jion=0;jion<surrogate[i].nion;jion++)
                  {            
                    double molality=surrogate[surrogate[i].iion(jion)].gamma_aq_bins(b)*surrogate[surrogate[i].iion(jion)].Aaq_bins_init(b)/surrogate[i].MM/conc_org*1000.0;
                    double flux=surrogate[i].kion(jion)*molality*surrogate[i].GAMMAinf*surrogate[i].gamma_aq_bins(b)*surrogate[i].Aaq_bins_init(b);		    
                    double flux2=0.0;
		    j=surrogate[i].iproduct(jion);
		    surrogate[i].kloss_aq(b)+=surrogate[i].kion(jion)*molality*surrogate[i].GAMMAinf*surrogate[i].gamma_aq_bins(b);
                    if (surrogate[i].rion_catalyzed(jion)==false and molality>0.)
                      {                                                  
			surrogate[surrogate[i].iion(jion)].kloss_aq(b)+=flux/surrogate[i].MM*surrogate[surrogate[i].iion(jion)].MM/surrogate[surrogate[i].iion(jion)].Aaq_bins_init(b);
			flux2+=flux/surrogate[i].MM*surrogate[surrogate[i].iion(jion)].MM;
                      }                
                    surrogate[j].kprod_aq(b)+=flux+flux2;		                    
                  }
            }    
        }
    }
}

void prodloss_aq_ssh(model_config &config, vector<species>& surrogate, Array<double, 1> &AQinit,  Array<double, 1> &LWC, Array<double, 3> &MOinit, int index, double deltat)
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
      surrogate[i].deltat_exp=0.;
      surrogate[i].k1_gas(index)=0.0;
      surrogate[i].kloss_gas=0.0;
      surrogate[i].kprod_gas=0.0;
      for (b=0;b<config.nbins;b++)
	{
	  surrogate[i].k1_aq(b,index)=0.0;
	  surrogate[i].Jdn_aq(b,index)=0.0;
	  surrogate[i].kprod_aq(b)=0.0;
	  surrogate[i].kloss_aq(b)=0.0;
	}
      
      for (b=0;b<config.nbins;++b)		  
	for (ilayer=0;ilayer<config.nlayer;++ilayer)
	  for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
	    {
	      surrogate[i].kloss(b,ilayer,iphase)=0.0;
	      surrogate[i].kprod(b,ilayer,iphase)=0.0;
	      surrogate[i].k1(b,ilayer,iphase,index)=0.0;	      
	    }
	
    }
      

  if (sum(LWC)>config.LWClimit)
  for (b=0;b<config.nbins;++b)
    {
      conc_org=LWC(b);
      for (i=0;i<n;++i)
        if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophilic)
          conc_org+=surrogate[i].Aaq_bins_init(b);
      conc_org=max(conc_org,1.e-5*AQinit(b)); //config.MOmin);
      //conc_org=max(conc_org,config.MOmin);

      for (i=0;i<n;++i)
        if(surrogate[i].is_organic and surrogate[i].hydrophilic and config.compute_organic)
          {
            sum_mass=AQinit(b);
            for (ilayer=0;ilayer<config.nlayer;ilayer++)
              for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                sum_mass+=MOinit(b,ilayer,iphase);
          
            Kaq=surrogate[i].Kaq(b);          
	    if (Kaq > 1.e-20 and AQinit(b)> 1e-20 and sum_mass > 1e-20 and surrogate[i].tau_air(b) > 0.)
	      {
		if (index==0)
		  {
		    surrogate[i].k1_aq(b,index)=AQinit(b)/sum_mass*
		      (surrogate[i].Ag*Kaq*AQinit(b)-surrogate[i].Aaq_bins_init(b))/
		      (Kaq*AQinit(b)*surrogate[i].tau_air(b));
		    surrogate[i].k1_gas(index)-=AQinit(b)/sum_mass*
		      (surrogate[i].Ag*Kaq*AQinit(b)-surrogate[i].Aaq_bins_init(b))/
		      (Kaq*AQinit(b)*surrogate[i].tau_air(b));
		  }

		surrogate[i].kprod_aq(b)=surrogate[i].Ag*AQinit(b)/sum_mass/surrogate[i].tau_air(b);
		//cout << surrogate[i].tau_air << endl;
		//cout << surrogate[i].name << " " << AQinit(b) << " " << sum_mass << " " << Kaq << " " << AQinit(b) << " " << surrogate[i].tau_air(b) << endl;
		surrogate[i].kloss_aq(b)=AQinit(b)/sum_mass/(Kaq*AQinit(b)*surrogate[i].tau_air(b));
		surrogate[i].kprod_gas+=AQinit(b)/sum_mass*surrogate[i].Aaq_bins_init(b)/(Kaq*AQinit(b)*surrogate[i].tau_air(b));
		surrogate[i].kloss_gas+=AQinit(b)/sum_mass/surrogate[i].tau_air(b);
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
		if (index==0)
		  {
		    surrogate[i].k1_gas(index)-=surrogate[i].Ag/surrogate[i].tau_air(b);
		    surrogate[i].k1_aq(b,index)=surrogate[i].Ag/surrogate[i].tau_air(b)*surrogate[config.iSO4mm].MM/surrogate[i].MM;
		  }
		
		surrogate[config.iSO4mm].kprod_aq(b)=surrogate[i].Ag/surrogate[i].tau_air(b)*surrogate[config.iSO4mm].MM/surrogate[i].MM;
		surrogate[config.iSO4mm].kloss_aq(b)=0.0;
		surrogate[i].kprod_gas+=0.0;
		surrogate[i].kloss_gas+=1.0/surrogate[i].tau_air(b);		
              }
            else if (i==config.iNH3)
              {
                Kaq=surrogate[i].Kaq(b);
                //compute kinetic rate of absorption
                conc_aq=surrogate[config.iNH4p].Aaq_bins_init(b)/surrogate[config.iNH4p].MM*surrogate[i].MM;

		if (index==0)
		  {
		    surrogate[i].k1_gas(index)-=(surrogate[i].Ag*Kaq*conc_org-conc_aq)/(Kaq*conc_org*surrogate[i].tau_air(b))*AQinit(b)/sum_mass*surrogate[i].fac_corr_ph(b);
		    surrogate[config.iNH4p].k1_aq(b,index)=(surrogate[i].Ag*Kaq*conc_org-conc_aq)/(Kaq*conc_org*surrogate[i].tau_air(b))*AQinit(b)/sum_mass*surrogate[i].fac_corr_ph(b)
		      *surrogate[config.iNH4p].MM/surrogate[i].MM;
		  }

		surrogate[config.iNH4p].kprod_aq(b)=surrogate[i].Ag/surrogate[i].tau_air(b)*AQinit(b)/sum_mass*surrogate[i].fac_corr_ph(b)*surrogate[config.iNH4p].MM/surrogate[i].MM;
		surrogate[config.iNH4p].kloss_aq(b)=1.0/(Kaq*conc_org*surrogate[i].tau_air(b))*AQinit(b)/sum_mass*surrogate[i].fac_corr_ph(b);
		surrogate[i].kprod_gas+=conc_aq/(Kaq*conc_org*surrogate[i].tau_air(b))*AQinit(b)/sum_mass*surrogate[i].fac_corr_ph(b);
		surrogate[i].kloss_gas+=1.0/surrogate[i].tau_air(b)*AQinit(b)/sum_mass*surrogate[i].fac_corr_ph(b);

              }
            else if (i==config.iHNO3)
              {
                Kaq=surrogate[i].Kaq(b);
                //compute kinetic rate of absorption
                conc_aq=surrogate[config.iNO3m].Aaq_bins_init(b)/surrogate[config.iNO3m].MM*surrogate[i].MM;
		if (index==0)
		  {
		    surrogate[i].k1_gas(index)-=(surrogate[i].Ag*Kaq*conc_org-conc_aq)/(Kaq*conc_org*surrogate[i].tau_air(b))*AQinit(b)/sum_mass*surrogate[i].fac_corr_ph(b);
		    surrogate[config.iNO3m].k1_aq(b,index)=(surrogate[i].Ag*Kaq*conc_org-conc_aq)/(Kaq*conc_org*surrogate[i].tau_air(b))*AQinit(b)/sum_mass*surrogate[i].fac_corr_ph(b)
		      *surrogate[config.iNO3m].MM/surrogate[i].MM;
		  }

		surrogate[config.iNO3m].kprod_aq(b)=surrogate[i].Ag/surrogate[i].tau_air(b)*AQinit(b)/sum_mass*surrogate[i].fac_corr_ph(b)*surrogate[config.iNO3m].MM/surrogate[i].MM;
		surrogate[config.iNO3m].kloss_aq(b)=1.0/(Kaq*conc_org*surrogate[i].tau_air(b))*AQinit(b)/sum_mass*surrogate[i].fac_corr_ph(b);
		surrogate[i].kprod_gas+=conc_aq/(Kaq*conc_org*surrogate[i].tau_air(b))*AQinit(b)/sum_mass*surrogate[i].fac_corr_ph(b);
		surrogate[i].kloss_gas+=1.0/surrogate[i].tau_air(b)*AQinit(b)/sum_mass*surrogate[i].fac_corr_ph(b);


              }
            else if (i==config.iHCl)
              {
                Kaq=surrogate[i].Kaq(b);
                //compute kinetic rate of absorption
                conc_aq=surrogate[config.iClm].Aaq_bins_init(b)/surrogate[config.iClm].MM*surrogate[i].MM;
		if (index==0)
		  {
		    surrogate[i].k1_gas(index)-=(surrogate[i].Ag*Kaq*conc_org-conc_aq)/(Kaq*conc_org*surrogate[i].tau_air(b))*AQinit(b)/sum_mass*surrogate[i].fac_corr_ph(b);
		    surrogate[config.iClm].k1_aq(b,index)=(surrogate[i].Ag*Kaq*conc_org-conc_aq)/(Kaq*conc_org*surrogate[i].tau_air(b))*AQinit(b)/sum_mass*surrogate[i].fac_corr_ph(b)
		      *surrogate[config.iClm].MM/surrogate[i].MM;
		  }

		surrogate[config.iClm].kprod_aq(b)=surrogate[i].Ag/surrogate[i].tau_air(b)*AQinit(b)/sum_mass*surrogate[i].fac_corr_ph(b)*surrogate[config.iClm].MM/surrogate[i].MM;
		surrogate[config.iClm].kloss_aq(b)=1.0/(Kaq*conc_org*surrogate[i].tau_air(b))*AQinit(b)/sum_mass*surrogate[i].fac_corr_ph(b);
		surrogate[i].kprod_gas+=conc_aq/(Kaq*conc_org*surrogate[i].tau_air(b))*AQinit(b)/sum_mass*surrogate[i].fac_corr_ph(b);
		surrogate[i].kloss_gas+=1.0/surrogate[i].tau_air(b)*AQinit(b)/sum_mass*surrogate[i].fac_corr_ph(b);		
              }


          
          }
    }

  /*
  double redmax=0.0;
  double incmax=0.0;
  if (index==0)
    for (i=0;i<n;i++)
      if (surrogate[i].is_organic)
	{
	  surrogate[i].deltat_exp=0.5*deltat;
	  if (surrogate[i].k1_gas(index)<0.)
	    {
	      double time=redmax*surrogate[i].Ag/(-surrogate[i].k1_gas(index));
	      surrogate[i].deltat_exp=min(time,surrogate[i].deltat_exp);
	    }

	  if (surrogate[i].hydrophobic)
	    for (b=0;b<config.nbins;b++)
	      for (ilayer=0;ilayer<config.nlayer;++ilayer)
		for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		  if (surrogate[i].k1(b,ilayer,iphase,index)<0.)
		    {
		      double time=redmax*surrogate[i].Ap_layer_init(b,ilayer,iphase)/(-surrogate[i].k1(b,ilayer,iphase,index));
		      surrogate[i].deltat_exp=min(time,surrogate[i].deltat_exp);
		    }

	  if (surrogate[i].hydrophilic)
	    for (b=0;b<config.nbins;b++)
	      if (surrogate[i].k1_aq(b)<0.)
		{
		  double time=redmax*surrogate[i].Aaq_bins_init(b)/(-surrogate[i].k1_aq(b));
		  surrogate[i].deltat_exp=min(time,surrogate[i].deltat_exp);
		}  
	  
	}
      else if (surrogate[i].is_inorganic_precursor and config.compute_inorganic and surrogate[i].is_solid==false)
	{
	  surrogate[i].deltat_exp=0.5*deltat;
	  if (surrogate[i].k1_gas(index)<0.)
	    {
	      double time=redmax*surrogate[i].Ag/(-surrogate[i].k1_gas(index));
	      surrogate[i].deltat_exp=min(time,surrogate[i].deltat_exp);
	    }
	  else if (surrogate[i].k1_gas(index)>0.)
	    {
	      double time=incmax*surrogate[i].Ag/(surrogate[i].k1_gas(index));
	      surrogate[i].deltat_exp=min(time,surrogate[i].deltat_exp);
	    }
 
	  if (i==config.iHNO3)
	    {
	      for (b=0;b<config.nbins;b++)
		if (surrogate[config.iNO3m].k1_aq(b)<0.)
		  {
		    double time=redmax*surrogate[config.iNO3m].Aaq_bins_init(b)/(-surrogate[config.iNO3m].k1_aq(b));
		    surrogate[i].deltat_exp=min(time,surrogate[i].deltat_exp);
		  }
		else if (surrogate[config.iNO3m].k1_aq(b)>0.)
		  {
		    double time=incmax*surrogate[config.iNO3m].Aaq_bins_init(b)/surrogate[config.iNO3m].k1_aq(b);
		    surrogate[i].deltat_exp=min(time,surrogate[i].deltat_exp);
		  }
	      surrogate[config.iNO3m].deltat_exp=surrogate[i].deltat_exp;
	      
	    }
	  else if (i==config.iNH3)
	    {
	      for (b=0;b<config.nbins;b++)
		if (surrogate[config.iNH4p].k1_aq(b)<0.)
		  {
		    double time=redmax*surrogate[config.iNH4p].Aaq_bins_init(b)/(-surrogate[config.iNH4p].k1_aq(b));
		    surrogate[i].deltat_exp=min(time,surrogate[i].deltat_exp);
		  }
		else if (surrogate[config.iNH4p].k1_aq(b)>0.)
		  {
		    double time=incmax*surrogate[config.iNH4p].Aaq_bins_init(b)/surrogate[config.iNH4p].k1_aq(b);
		    surrogate[i].deltat_exp=min(time,surrogate[i].deltat_exp);
		  }
	      surrogate[config.iNH4p].deltat_exp=surrogate[i].deltat_exp;
	    }
	  else if (i==config.iHCl)
	    {
	      for (b=0;b<config.nbins;b++)
		if (surrogate[config.iClm].k1_aq(b)<0.)
		  {
		    double time=redmax*surrogate[config.iClm].Aaq_bins_init(b)/(-surrogate[config.iClm].k1_aq(b));
		    surrogate[i].deltat_exp=min(time,surrogate[i].deltat_exp);
		  }
		else if (surrogate[config.iClm].k1_aq(b)>0.)
		  {
		    double time=incmax*surrogate[config.iClm].Aaq_bins_init(b)/surrogate[config.iClm].k1_aq(b);
		    surrogate[i].deltat_exp=min(time,surrogate[i].deltat_exp);
		  }
	      surrogate[config.iClm].deltat_exp=surrogate[i].deltat_exp;
	    }
	  else if (i==config.iH2SO4)
	    {
	      //surrogate[i].deltat_exp=0.;
	      surrogate[config.iHSO4m].deltat_exp=surrogate[i].deltat_exp;
	      surrogate[config.iSO4mm].deltat_exp=surrogate[i].deltat_exp;
	    }
	  //cout << surrogate[i].name << " " << surrogate[i].deltat_exp << endl;
	  
	}
  */
  //exit(0);
 
}



void flux_aq_ssh(model_config &config, vector<species>& surrogate, Array<double, 1> &AQinit,  Array<double, 1> &LWC, Array<double, 3> &MOinit, double &tiny, int index)
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
      conc_org=LWC(b);
      for (i=0;i<n;++i)
        if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophilic)
          conc_org+=surrogate[i].Aaq_bins_init(b);
      conc_org=max(conc_org,1.e-5*AQinit(b)); //config.MOmin);
      //conc_org=max(conc_org,config.MOmin);

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
            if (surrogate[i].name=="H2SO4")
              {	      
                surrogate[i].k1_aq(b,index)=surrogate[i].Ag/surrogate[i].tau_air(b);      
                surrogate[i].Jdn_aq(b,index)=0.0;
              }
            else if (surrogate[i].name=="NH3")
              {
                Kaq=surrogate[i].Kaq(b);
                //compute kinetic rate of absorption
                conc_aq=surrogate[config.iNH4p].Aaq_bins_init(b)/surrogate[config.iNH4p].MM*surrogate[i].MM;
	      
                surrogate[i].k1_aq(b,index)=(surrogate[i].Ag*Kaq*conc_org-conc_aq)/(Kaq*conc_org*surrogate[i].tau_air(b))*AQinit(b)/sum_mass*surrogate[i].fac_corr_ph(b);	   
 
                surrogate[i].Jdn_aq(b,index)=0.0;
                if (surrogate[i].k1_aq(b,index)<0.0 and surrogate[config.iNH4p].Aaq_bins_init(b)>tiny)
                  surrogate[i].Jdn_aq(b,index)=surrogate[i].k1_aq(b,index)/surrogate[config.iNH4p].Aaq_bins_init(b)*surrogate[config.iNH4p].MM/surrogate[i].MM;
              }
            else if (surrogate[i].name=="HNO3")
              {
                Kaq=surrogate[i].Kaq(b);
                //compute kinetic rate of absorption
                conc_aq=surrogate[config.iNO3m].Aaq_bins_init(b)/surrogate[config.iNO3m].MM*surrogate[i].MM;
	      
                surrogate[i].k1_aq(b,index)=(surrogate[i].Ag*Kaq*conc_org-conc_aq)/(Kaq*conc_org*surrogate[i].tau_air(b))*AQinit(b)/sum_mass*surrogate[i].fac_corr_ph(b);     
              
                surrogate[i].Jdn_aq(b,index)=0.0;
                if (surrogate[i].k1_aq(b,index)<0.0 and surrogate[config.iNO3m].Aaq_bins_init(b)>tiny)
                  surrogate[i].Jdn_aq(b,index)=surrogate[i].k1_aq(b,index)/surrogate[config.iNO3m].Aaq_bins_init(b)*surrogate[config.iNO3m].MM/surrogate[i].MM;
              }
            else if (surrogate[i].name=="HCl")
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


void correct_flux_ph_ssh(model_config &config, vector<species>& surrogate, double &Temperature, Array<double, 1> &AQinit, Array <double, 3> &MOinit, 
			 Array<double, 1> &chp, 		     
			 Array<double, 1> &chpout, Array<double, 1> &MMaq, Array<double, 1> &ionic, Array<double, 1> &LWC, double tiny, double deltat, int index)
{
  int i,b;
  int n=surrogate.size();
  double a1,b1,c1,asulf; 
  double a2,b2,c2,d2;
  double Kaq,conc_aq;
  double Ke=1.010e-14*exp(-22.52*(298./Temperature-1.0)+26.92*(1+log(298./Temperature)-298./Temperature)); //Dissociation constant of water Ke=a(H+)*a(HO-)
  double Jhp,chp2,Jhp2;
  double epser=0.01;
  double tmini=min(config.deltatmin/100,1e-4);
  double t;
  double dt1;
  double delta,chpeq1,chpeq2;
  double organion,organion2,kelvin_effect;
  double sum_mass,faq;
  int ilayer,iphase;
  double conc_org,f;
  double fion1=0.0;
  double fion2=0.0;
  double ratio_gamma;
  double ratio_gamma1;
  double ratio_gamma2;
  double Kac,Kac1,Kac2;
  double Kh,dK,K,dJ,total;
  double gamma=1.0+pow(2.,0.5)/2.;
  
  tmini=deltat/10;
  //cout << " avant " << surrogate[config.iHNO3].k1_aq(0,index) << " " << index << endl;
  conc_org=LWC(0);
  for (i=0;i<n;++i)
    if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophilic)
      conc_org+=surrogate[i].Aaq_bins_init(0);
  //cout << index << " " << surrogate[config.iHNO3].Ag << " " << surrogate[config.iNO3m].Aaq_bins_init(0)/surrogate[config.iNO3m].MM*surrogate[config.iHNO3].MM/(surrogate[config.iHNO3].Kaq(0)*conc_org) << " " << surrogate[config.iHNO3].k1_aq(0,index) << endl;
  
  for (b=0;b<config.nbins;b++)
    //if (index==1)
    {
      t=0;
      a1=0.0;
      b1=0.0;
      c1=0.0;
      surrogate[config.iHp].k1_aq(b,index)=0.0;
      conc_org=LWC(b);
      for (i=0;i<n;++i)
        if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophilic)
          conc_org+=surrogate[i].Aaq_bins_init(b);
      
      sum_mass=AQinit(b);
      for (ilayer=0;ilayer<config.nlayer;ilayer++)
        for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
          sum_mass+=MOinit(b,ilayer,iphase);

      if (sum_mass>0.0)
        faq=AQinit(b)/sum_mass;
      else
        faq=0.0;

      kelvin_effect=1.0;
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
        if (surrogate[i].is_inorganic_precursor)
          {
            if (surrogate[i].name=="H2SO4")	      
              {
		
                f=1.0;
                if(surrogate[config.iSO4mm].Aaq_bins_init(b)/surrogate[config.iSO4mm].MM+
                   surrogate[config.iHSO4m].Aaq_bins_init(b)/surrogate[config.iHSO4m].MM>0.0)
                  f=surrogate[config.iSO4mm].Aaq_bins_init(b)/surrogate[config.iSO4mm].MM/
                    (surrogate[config.iSO4mm].Aaq_bins_init(b)/surrogate[config.iSO4mm].MM+
                     surrogate[config.iHSO4m].Aaq_bins_init(b)/surrogate[config.iHSO4m].MM);	       
                surrogate[config.iHp].k1_aq(b,index)+=(1.+f)*surrogate[i].k1_aq(b,index)/surrogate[i].MM;
		//Jhp2+=(1.+f)*surrogate[i].k1_aq(b,1)/surrogate[i].MM; 	
                surrogate[i].Agt=surrogate[i].Ag;
              }
            else if (surrogate[i].name=="NH3")
              {		
                Kaq=surrogate[i].Kaq(b)/chp(b);	       
                conc_aq=surrogate[config.iNH4p].Aaq_bins_init(b)/surrogate[config.iNH4p].MM;
                c1+=conc_aq/(Kaq*conc_org*surrogate[i].tau_air(b))*faq;
		//Jhp-=surrogate[i].k1_aq(b,0)/surrogate[i].MM;
		surrogate[config.iHp].k1_aq(b,index)-=surrogate[i].k1_aq(b,index)/surrogate[i].MM; 	
                surrogate[i].Agt=surrogate[i].Ag;	      
              }
            else if (surrogate[i].name=="HNO3")
              {
                Kaq=surrogate[i].Kaq(b)*chp(b);	       
                conc_aq=surrogate[config.iNO3m].Aaq_bins_init(b)/surrogate[config.iNO3m].MM;
		//Jhp+=surrogate[i].k1_aq(b,0)/surrogate[i].MM;
		surrogate[config.iHp].k1_aq(b,index)+=surrogate[i].k1_aq(b,index)/surrogate[i].MM; 
		
              }
            else if (surrogate[i].name=="HCl")
              {
                Kaq=surrogate[i].Kaq(b)*chp(b);	       
		//Jhp+=surrogate[i].k1_aq(b,0)/surrogate[i].MM;
		surrogate[config.iHp].k1_aq(b,index)+=surrogate[i].k1_aq(b,index)/surrogate[i].MM; 
              }
	    
          }
        else if (surrogate[i].is_organic and surrogate[i].time_aq(b)>=config.tequilibrium)
          if (surrogate[i].aq_type=="monoacid" or surrogate[i].aq_type=="diacid")
            {	      
              surrogate[i].Agt=surrogate[i].Ag;	     
              surrogate[i].Aaq=0.0;
            }
      
      double sensimax=0.;
      double sensicl=0;
      double sensino3=0.;
      double sensinh4=0.;
      if (config.iHCl>=0)
	sensicl=(surrogate[config.iClm].Aaq_bins_init(b)/surrogate[config.iClm].MM+surrogate[config.iHCl].Ag/surrogate[config.iHCl].MM)*
	  (surrogate[config.iHCl].dKaq(b)*AQinit(b)/(1+surrogate[config.iHCl].Kaq(b)*AQinit(b))-surrogate[config.iHCl].Kaq(b)*AQinit(b)/(pow(1+surrogate[config.iHCl].Kaq(b)*AQinit(b),2))*AQinit(b)*surrogate[config.iHCl].dKaq(b));
      if (config.iHNO3>=0)
	/*
	  sensino3=(surrogate[config.iNO3m].Aaq_bins_init(b)/surrogate[config.iNO3m].MM+surrogate[config.iHNO3].Ag/surrogate[config.iHNO3].MM)*
	  (surrogate[config.iHNO3].dKaq(b)*AQinit(b)/(1+surrogate[config.iHNO3].Kaq(b)*AQinit(b))-surrogate[config.iHNO3].Kaq(b)*AQinit(b)/(pow(1+surrogate[config.iHNO3].Kaq(b)*AQinit(b),2))*AQinit(b)*surrogate[config.iHNO3].dKaq(b));*/
	sensino3=-(surrogate[config.iHNO3].k1_aq(b,index)*deltat/(1.0-gamma*surrogate[config.iHNO3].Jdn_aq(b,0)*deltat) -surrogate[config.iHNO3].k1_aq(b,0))/surrogate[config.iHNO3].k1_aq(b,0);
      if (config.iNH3>=0)
	sensinh4=-(surrogate[config.iNH3].k1_aq(b,index)*deltat/(1.0-gamma*surrogate[config.iNH3].Jdn_aq(b,0)*deltat) -surrogate[config.iNH3].k1_aq(b,0))/surrogate[config.iNH3].k1_aq(b,0);
      /*(surrogate[config.iNH4p].Aaq_bins_init(b)/surrogate[config.iNH4p].MM+surrogate[config.iNH3].Ag/surrogate[config.iNH3].MM)*
	(surrogate[config.iNH3].dKaq(b)*AQinit(b)/(1+surrogate[config.iNH3].Kaq(b)*AQinit(b))-surrogate[config.iNH3].Kaq(b)*AQinit(b)/(pow(1+surrogate[config.iNH3].Kaq(b)*AQinit(b),2))*AQinit(b)*surrogate[config.iNH3].dKaq(b)); */
      
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

      
      //cout << surrogate[isensi].name << " " << sensimax << endl;
      /*
      //surrogate[isensi].k1_aq(b,index)*=0.1;
      surrogate[config.iHNO3].k1_aq(b,index)*=0.1;
      surrogate[config.iNH3].k1_aq(b,index)*=0.1;*/

      /*
	i=config.iHNO3;
	if (index==1 and surrogate[i].fac_corr_ph(b)>0.2)
	{
	surrogate[i].fac_corr_ph(b)=0.1;
	config.to_be_rejected=true;
	cout << "corr HNO3 " << surrogate[i].fac_corr_ph(b) << " " << surrogate[i].time_aq(b) << endl;
	}

	i=config.iNH3;
	if (index==1 and surrogate[i].fac_corr_ph(b)>0.2) // and surrogate[i].k1_aq(b,index)*surrogate[i].k1_aq(b,0)<0 and surrogate[i].fac_corr_ph(b)>0.2)
	{
	surrogate[i].fac_corr_ph(b)=0.1;
	config.to_be_rejected=true;
	cout << "corr NH3 " << surrogate[i].fac_corr_ph(b) << " " << surrogate[i].time_aq(b) << endl;
	}

	i=config.iHCl;
	if (index==10 and surrogate[i].k1_aq(b,index)*surrogate[i].k1_aq(b,0)<0 and surrogate[i].fac_corr_ph(b)>0.2)
	{
	surrogate[i].fac_corr_ph(b)=0.1;
	config.to_be_rejected=true;
	cout << "corr HCl " << surrogate[i].fac_corr_ph(b) << endl;
	}*/
      	
      if (index==1) // and surrogate[config.iHp].k1_aq(b,index)+ surrogate[config.iHp].k1_aq(b,0)<= 0.1*surrogate[config.iHp].k1_aq(b,0))
	//if (index==1 and surrogate[config.iHp].k1_aq(b,index)<0.2*surrogate[config.iHp].k1_aq(b,0))
	{
	  //cout << index << " " << surrogate[config.iHp].k1_aq(b,0) << " " << surrogate[config.iHp].k1_aq(b,index) << endl;
	  //surrogate[config.iHNO3].fac_corr_ph(b)*=0.9;
	  //surrogate[config.iNH3].fac_corr_ph(b)*=0.9; //-=surrogate[config.iHp].k1_aq(b,index)/2*surrogate[config.iHNO3].MM;

	  /*
	    if (surrogate[config.iHNO3].fac_corr_ph(b)>0.2 and surrogate[config.iNH3].fac_corr_ph(b)>0.2 and surrogate[config.iHCl].fac_corr_ph(b)>0.2)
	    {*/
	  //cout << surrogate[isensi].name << endl;

	  /*
	    surrogate[config.iHNO3].fac_corr_ph(b)=min(max(0.1,(surrogate[config.iHNO3].k1_aq(b,index)*deltat/(1.0-gamma*surrogate[config.iHNO3].Jdn_aq(b,0)*deltat)+surrogate[config.iHNO3].k1_aq(b,0))/2/surrogate[config.iHNO3].k1_aq(b,0))*surrogate[config.iHNO3].fac_corr_ph(b),1.);
	    surrogate[config.iNH3].fac_corr_ph(b)=min(max(0.1,(surrogate[config.iNH3].k1_aq(b,index)*deltat/(1.0-gamma*surrogate[config.iNH3].Jdn_aq(b,0)*deltat)+surrogate[config.iNH3].k1_aq(b,0))/2/surrogate[config.iNH3].k1_aq(b,0))*surrogate[config.iNH3].fac_corr_ph(b),1.);
	  */

	  //  }
	  /*i=config.iHNO3;
	  
	    if (surrogate[i].k1_aq(b,index)*surrogate[i].k1_aq(b,0)<0)
	    {
	    cout << "corr HNO3 " << surrogate[i].fac_corr_ph(b) << endl;
	      
	    if (surrogate[i].fac_corr_ph(b)>0.2)
	    {
	    surrogate[i].fac_corr_ph(b)=0.1;
	    config.to_be_rejected=true;
	    }
	      
	    }
	    //surrogate[i].fac_corr_ph(b)
	  
	    i=config.iHCl;
	    if (surrogate[i].k1_aq(b,index)*surrogate[i].k1_aq(b,0)<0)
	    surrogate[i].fac_corr_ph(b)=0.1;

	    i=config.iNH3;
	    if (surrogate[i].k1_aq(b,index)*surrogate[i].k1_aq(b,0)<0)
	    {
	    cout << "corr NH3 " << surrogate[i].fac_corr_ph(b) << endl;
	    surrogate[i].fac_corr_ph(b)=0.1;
	    }*/
	}
      
      //cout << index << " " << surrogate[config.iHp].k1_aq(b,index) << " " << surrogate[config.iHNO3].k1_aq(b,index) << " " << surrogate[config.iNH3].k1_aq(b,index) << " " << conc_org << " " << chp(0) << endl;
      
      //cout << surrogate[config.iHNO3].Ag << " " << surrogate[config.iNO3m].Aaq_bins_init(0)/surrogate[config.iNO3m].MM*surrogate[config.iHNO3].MM/(surrogate[config.iHNO3].Kaq(b)*conc_org) << endl;

      //cout << "apres " << Jhp << endl;
      /*
	chpout(b)=chp2;
	for (i=0;i<n;i++)
        if (surrogate[i].name=="NH3")
	{	   
	surrogate[i].k1_aq(b,index)=a2/deltat;
	    
	surrogate[i].Jdn_aq(b,index)=0.0;
	if (surrogate[i].k1_aq(b,index)<0.0 and surrogate[config.iNH4p].Aaq_bins_init(b)>tiny)
	surrogate[i].Jdn_aq(b,index)=surrogate[i].k1_aq(b,index)/surrogate[config.iNH4p].Aaq_bins_init(b)*surrogate[config.iNH4p].MM/surrogate[i].MM;
	}
        else if (surrogate[i].name=="HNO3")
	{	   
	surrogate[i].k1_aq(b,index)=b2/deltat;
	surrogate[i].Jdn_aq(b,index)=0.0;
	if (surrogate[i].k1_aq(b,index)<0.0 and surrogate[config.iNO3m].Aaq_bins_init(b)>tiny)
	surrogate[i].Jdn_aq(b,index)=surrogate[i].k1_aq(b,index)/surrogate[config.iNO3m].Aaq_bins_init(b)*surrogate[config.iNO3m].MM/surrogate[i].MM;
	}
        else if (surrogate[i].name=="HCl")
	{	   
	surrogate[i].k1_aq(b,index)=c2/deltat;
	surrogate[i].Jdn_aq(b,index)=0.0;
	if (surrogate[i].k1_aq(b,index)<0.0 and surrogate[config.iClm].Aaq_bins_init(b)>tiny)
	surrogate[i].Jdn_aq(b,index)=surrogate[i].k1_aq(b,index)/surrogate[config.iClm].Aaq_bins_init(b)*surrogate[config.iClm].MM/surrogate[i].MM;
	}    
        else if (surrogate[i].is_organic)
	if ((surrogate[i].aq_type=="monoacid" or surrogate[i].aq_type=="diacid") and surrogate[i].time_aq(b)>=config.tequilibrium)
	{
	if (surrogate[i].time_aq(b)>=config.tequilibrium)
	surrogate[i].k1_aq(b,index)=surrogate[i].Aaq/deltat;
	surrogate[i].Jdn_aq(b,index)=0.0;
	if (surrogate[i].k1_aq(b,index)<0.0 and surrogate[i].Aaq_bins_init(b)>tiny)
	surrogate[i].Jdn_aq(b,index)=surrogate[i].k1_aq(b,index)/surrogate[i].Aaq_bins_init(b);
	}*/
    }
 
  //cout << " apres " << surrogate[config.iHNO3].k1_aq(0,index) << " " << index << " " << b2 << " " << deltat << " " << tmini << " " << chp2 << " " << chp(0) << endl;
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
		      //if (surrogate[i].Ap_layer_init(b,ilayer,0)>0) cout << ilayer << " " << surrogate[i].tau_diffusion(b,ilayer,0) << " " << config.tequilibrium << endl;
		      for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
			if (surrogate[i].tau_diffusion(b,ilayer,iphase)<config.tequilibrium)
			  {
			    tau_interface=surrogate[i].tau_diffusion(b,ilayer,iphase);
			    ilayer_interface=ilayer;
			    //if (surrogate[i].Ap_layer_init(b,ilayer,iphase)>0) 
			    //  cout << surrogate[i].name << " " << ilayer_interface << endl;
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

		  /*
		  if (surrogate[i].name=="SOAlP")
		    cout << "Jinterface " << ktot1 << " " << surrogate[i].k1 << endl;*/
		     		  
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
		  
		  kcond=(1.0-AQinit(b)/sum_mass)*(kcond)/surrogate[i].tau_air(b);	  	      
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

		      surrogate[i].kprod(b,ilayer,iphase)+=surrogate[i].Ag*surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)/
			(sum/(1.0-AQinit(b)/sum_mass)*surrogate[i].tau_air(b));
		      surrogate[i].kloss(b,ilayer,iphase)+=1./(sum/(1.0-AQinit(b)/sum_mass)*surrogate[i].tau_air(b));
		      surrogate[i].kprod_gas+=surrogate[i].Ap_layer_init(b,ilayer,iphase)/
			(sum/(1.0-AQinit(b)/sum_mass)*surrogate[i].tau_air(b));
		      surrogate[i].kprod_gas+=surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)/
			    (sum/(1.0-AQinit(b)/sum_mass)*surrogate[i].tau_air(b));

		      if (index==0)
			surrogate[i].k1(b,ilayer,iphase,index)=
			  (surrogate[i].Ag*surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)
			   -surrogate[i].Ap_layer_init(b,ilayer,iphase))/
			  (sum/(1.0-AQinit(b)/sum_mass)*surrogate[i].tau_air(b));					      		     
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
	  /*
	  for (i=0;i<n;i++)
	    if (surrogate[i].name=="POAmP")
	      cout << surrogate[i].tau_diffusion << endl;
	    surrogate[i].tau_air=0.0001;*/
	  
	  for (i=0;i<n;i++)
	    for (b=0;b<config.nbins;b++)
	      for (ilayer=0;ilayer<config.nlayer-1;ilayer++)
		for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		  {
		    surrogate[i].kprod(b,ilayer,iphase)+=surrogate[i].Ap_layer_init(b,config.nlayer-1,iphase)/(surrogate[i].Kp(b,config.nlayer-1,iphase)*MOinit(b,config.nlayer-1,iphase))
		      *surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)/surrogate[i].tau_diffusion(b,ilayer,iphase);
		    surrogate[i].kloss(b,ilayer,iphase)+=1./surrogate[i].tau_diffusion(b,ilayer,iphase);
		    surrogate[i].kprod(b,config.nlayer-1,iphase,index)+=surrogate[i].Ap_layer_init(b,ilayer,iphase)/surrogate[i].tau_diffusion(b,ilayer,iphase);
		    surrogate[i].kloss(b,config.nlayer-1,iphase,index)+=1./(surrogate[i].Kp(b,config.nlayer-1,iphase)*MOinit(b,config.nlayer-1,iphase))
		      *surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)/surrogate[i].tau_diffusion(b,ilayer,iphase);
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
		      //if(surrogate[i].time(b,ilayer,iphase)>=config.tequilibrium)
		      {		      
			sum=0.0;
			for (jphase=0;jphase<config.nphase(b,ilayer);++jphase)
			  sum+=surrogate[i].Kp(b,ilayer,jphase)*MOinit(b,ilayer,jphase);
			
			surrogate[i].kprod(b,ilayer,iphase)+=surrogate[i].Ag*surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)
			  /(sum/(1.0-AQinit(b)/sum_mass)*surrogate[i].tau_air(b));
			surrogate[i].kloss(b,ilayer,iphase)+=1./(sum/(1.0-AQinit(b)/sum_mass)*surrogate[i].tau_air(b));
			surrogate[i].kloss_gas+=surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)
			  /(sum/(1.0-AQinit(b)/sum_mass)*surrogate[i].tau_air(b));
			surrogate[i].kprod_gas+=surrogate[i].Ap_layer_init(b,ilayer,iphase)/(sum/(1.0-AQinit(b)/sum_mass)*surrogate[i].tau_air(b));			
			if (sum > 0.0 and surrogate[i].tau_air(b) > 0.0 and index==0)
			  surrogate[i].k1(b,ilayer,iphase,index)=			
			    (surrogate[i].Ag*surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)
			     -surrogate[i].Ap_layer_init(b,ilayer,iphase))/
			    (sum/(1.0-AQinit(b)/sum_mass)*surrogate[i].tau_air(b));
			
			
		      }
	    }

	  /*
	  for (i=0;i<n;i++)
	    if (surrogate[i].name=="POAmP")
	      cout << surrogate[i].name << " " << surrogate[i].kloss_gas << " " << surrogate[i].kprod_gas << " " << surrogate[i].kprod << " " << surrogate[i].kloss << endl;
	  */
	  
	  /*
	  flux_org_ssh(config, surrogate, MOinit, AQinit, tiny, index);
	  for (i=0;i<n;i++)
	    for (b=0;b<config.nbins;b++)
	      for (ilayer=0;ilayer<config.nlayer;ilayer++)
		for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		  {
		    double diff=(surrogate[i].Ag*surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)
				 -surrogate[i].Ap_layer_init(b,ilayer,iphase));
		    double kloc=0.;
		    if (diff!=0.)
		      {
			kloc=surrogate[i].k1(b,ilayer,iphase,index)/diff;
		      }
		    if (kloc<0.)
		      {
			cout << "error kloc " << endl;
			exit(0);
		      }
		    double kloc2=0.*kloc+1.*surrogate[i].kloc(b,ilayer,iphase);
		    if (index==0)
		      surrogate[i].kloc(b,ilayer,iphase)=kloc;
		    kloc=kloc2; 
		    surrogate[i].kprod(b,ilayer,iphase)=kloc*surrogate[i].Ag*surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase);
		    surrogate[i].kloss(b,ilayer,iphase)=kloc;
		    surrogate[i].kloss_gas+=kloc*surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase);
		    surrogate[i].kprod_gas+=kloc*surrogate[i].Ap_layer_init(b,ilayer,iphase);
		  }
	  */
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
		      //if(surrogate[i].time(b,ilayer,iphase)>=config.tequilibrium)
		      {		      
			sum=0.0;
			for (jphase=0;jphase<config.nphase(b,ilayer);++jphase)
			  sum+=surrogate[i].Kp(b,ilayer,jphase)*MOinit(b,ilayer,jphase);

			surrogate[i].kprod(b,ilayer,iphase)=surrogate[i].Ag*surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)
			  /(sum/(1.0-AQinit(b)/sum_mass)*surrogate[i].tau_air(b));
			surrogate[i].kloss(b,ilayer,iphase)=1./(sum/(1.0-AQinit(b)/sum_mass)*surrogate[i].tau_air(b));
			surrogate[i].kloss_gas+=surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)
			  /(sum/(1.0-AQinit(b)/sum_mass)*surrogate[i].tau_air(b));
			surrogate[i].kprod_gas+=surrogate[i].Ap_layer_init(b,ilayer,iphase)/(sum/(1.0-AQinit(b)/sum_mass)*surrogate[i].tau_air(b));			
			if (sum > 0.0 and surrogate[i].tau_air(b) > 0.0 and index==0)
			  surrogate[i].k1(b,ilayer,iphase,index)=			
			    (surrogate[i].Ag*surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)
			     -surrogate[i].Ap_layer_init(b,ilayer,iphase))/
			    (sum/(1.0-AQinit(b)/sum_mass)*surrogate[i].tau_air(b));
			
			
		      }
	    }
	}
    }
}

void error_ph_bins_ssh(model_config &config, vector<species> &surrogate, int index_b, double Temperature, Array <double, 1> chp, 
		       double &error, double &derivative, Array <double, 1> AQinit, Array<double , 1> &ionic, Array<double, 1> &MMaq,
		       Array<double, 1> &LWC,
		       double &chp_new, double deltat)
{      
  int n=surrogate.size();
  int i;
  double inorganion=0.0;
  double organion=0.0;
  double total,dtotal;
  double sum_K_AQ;
  int b;
  double kelvin_effect=1.0;
  double Kaq;
  Array <double,1> conc_org;
  conc_org.resize(config.nbins);
  double ratio_gamma1=0.0;
  double ratio_gamma2=0.0;
  double Kac1=0.0;
  double Kac2=0.0;
  double Kh=0.0;
  double dK=0.0;
  double fion1=0.0;
  double fion2=0.0;
  double K, ratio_gamma, Kac;

  if (config.compute_kelvin_effect) //compute the kelvin effect    
    {
      kelvin_effect=2.0*config.surface_tension_aq*MMaq(index_b)/
	(8.314*Temperature*config.AQrho(index_b)*
	 0.5*config.diameters(index_b));
      if(kelvin_effect > 50.0)
	kelvin_effect = 50.0;
      kelvin_effect=exp(kelvin_effect);
    }
  for (b=0;b<config.nbins;b++)
    {
      conc_org(b)=LWC(b);  
      for (i=0;i<n;++i)
        if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophilic)
          conc_org(b)+=surrogate[i].Aaq_bins_init(b);
      conc_org(b)=max(conc_org(b),1.0e-5*AQinit(b)); //config.MOmin);
    }
  
  derivative=0.0;
  //cout << "ici " << deltat << endl;
  if (config.imethod==0)
    {
      //cout << "laaaaa" << endl;
      for (i=0;i<n;++i)
	if (surrogate[i].is_inorganic_precursor and surrogate[i].is_solid==false)
	  {
	    if (surrogate[i].name=="NH3")
	      {
		total=(surrogate[config.iNH4p].Aaq_bins_init0(index_b)
		       +0.5*deltat*(surrogate[config.iNH4p].k1_aq(index_b,0)+surrogate[config.iNH4p].kprod_aq(index_b)))/(1.0+0.5*deltat*surrogate[config.iNH4p].kloss_aq(index_b))
		  /surrogate[config.iNH4p].MM/conc_org(index_b)*1000.*surrogate[config.iNH4p].charge;

		dtotal=total/(1.0+0.5*deltat*surrogate[config.iNH4p].kloss_aq(index_b))*0.5*deltat*
		  surrogate[config.iNH4p].kloss_aq(index_b)/chp(index_b);
		derivative-=dtotal;
		inorganion-=total;
		//cout << "NH3: " << total << endl;
	      }
		
	    //cout << "NH3: " << inorganion << endl;
          
	
        else if (surrogate[i].name=="HNO3")
          {
	    total=(surrogate[config.iNO3m].Aaq_bins_init0(index_b)
		   +0.5*deltat*(surrogate[config.iNO3m].k1_aq(index_b,0)+surrogate[config.iNO3m].kprod_aq(index_b)))/(1.0+0.5*deltat*surrogate[config.iNO3m].kloss_aq(index_b))
	      /surrogate[config.iNO3m].MM/conc_org(index_b)*1000.*surrogate[config.iNO3m].charge;
	    
	    dtotal=-total/(1.0+0.5*deltat*surrogate[config.iNO3m].kloss_aq(index_b))*0.5*deltat*
	      surrogate[config.iNO3m].kloss_aq(index_b)/chp(index_b);
	    derivative-=dtotal;
	    inorganion-=total; 
	    //cout << "NO3: " << total << endl;  
          }
        else if (surrogate[i].name=="HCl")
          {
	    total=(surrogate[config.iClm].Aaq_bins_init0(index_b)
		   +0.5*deltat*(surrogate[config.iClm].k1_aq(index_b,0)+surrogate[config.iClm].kprod_aq(index_b)))/(1.0+0.5*deltat*surrogate[config.iClm].kloss_aq(index_b))
	      /surrogate[config.iClm].MM/conc_org(index_b)*1000.*surrogate[config.iClm].charge;
	    
	    dtotal=-total/(1.0+0.5*deltat*surrogate[config.iClm].kloss_aq(index_b))*0.5*deltat*
	      surrogate[config.iClm].kloss_aq(index_b)/chp(index_b);
	    derivative-=dtotal;
	    inorganion-=total; 
          }
        else if (surrogate[i].name=="H2SO4")
          {
            total=surrogate[config.iHSO4m].Aaq_bins_init(index_b)/surrogate[config.iHSO4m].MM
              +surrogate[config.iSO4mm].Aaq_bins_init(index_b)/surrogate[config.iSO4mm].MM;

	    K=surrogate[i].Kequilibrium_ssh(Temperature)*surrogate[config.iHSO4m].gamma_aq_bins(index_b)
              /(surrogate[config.iHp].gamma_aq_bins(index_b)*surrogate[config.iSO4mm].gamma_aq_bins(index_b));
	  
            derivative-=1000.*total/conc_org(index_b)*1.0/pow(1.0+K/chp(index_b),2.0)*K/(chp(index_b)*chp(index_b)); //HSO4m+SO4mm
            inorganion+=1000.*total/conc_org(index_b)*(2.0-1.0/(1.0+K/chp(index_b)));
	    //cout << "H2SO4: " << inorganion << " " << total << " " << conc_org << " " << endl;
          }
	else if (surrogate[i].is_organic and config.compute_organic)
	  {
	    /*
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
	      }*/
	  }

      
	  }



    }
  else
    for (i=0;i<n;++i)
      if (surrogate[i].is_inorganic_precursor and surrogate[i].is_solid==false)
	{
	  if (surrogate[i].name=="NH3")
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
	
	  else if (surrogate[i].name=="HNO3")
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
	  else if (surrogate[i].name=="HCl")
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
		  inorganion+=config.AQrho(index_b)*total*surrogate[i].Kaq(index_b)/sum_K_AQ;
		}            
	      else
		{
		  inorganion-=surrogate[config.iClm].Aaq_bins_init(index_b)/surrogate[config.iClm].MM/conc_org(index_b)*1000.*surrogate[config.iClm].charge;
		}
	      //cout << "HCL: " << inorganion << endl;
	    
	    }
	  else if (surrogate[i].name=="H2SO4")
	    {
	      total=surrogate[config.iHSO4m].Aaq_bins_init(index_b)/surrogate[config.iHSO4m].MM
		+surrogate[config.iSO4mm].Aaq_bins_init(index_b)/surrogate[config.iSO4mm].MM;

	      K=surrogate[i].Kequilibrium_ssh(Temperature)*surrogate[config.iHSO4m].gamma_aq_bins(index_b)
		/(surrogate[config.iHp].gamma_aq_bins(index_b)*surrogate[config.iSO4mm].gamma_aq_bins(index_b));
	  
	      derivative-=1000.*total/conc_org(index_b)*1.0/pow(1.0+K/chp(index_b),2.0)*K/(chp(index_b)*chp(index_b)); //HSO4m+SO4mm
	      inorganion+=1000.*total/conc_org(index_b)*(2.0-1.0/(1.0+K/chp(index_b)));
	      //cout << "H2SO4: " << inorganion << " " << total << " " << conc_org << " " << endl;
	    }
	}
      else if (surrogate[i].is_organic and config.compute_organic)
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

	}
  /*else 
    if (surrogate[i].name=="Na")*/
  if (config.iNa>=0) inorganion-=surrogate[config.iNa].Aaq_bins_init(index_b)/surrogate[config.iNa].MM/conc_org(index_b)*1000.*surrogate[config.iNa].charge;
  if (config.iMg>=0) inorganion-=surrogate[config.iMg].Aaq_bins_init(index_b)/surrogate[config.iMg].MM/conc_org(index_b)*1000.*surrogate[config.iMg].charge;
  if (config.iCa>=0) inorganion-=surrogate[config.iCa].Aaq_bins_init(index_b)/surrogate[config.iCa].MM/conc_org(index_b)*1000.*surrogate[config.iCa].charge;  
  if (config.iK>=0) inorganion-=surrogate[config.iK].Aaq_bins_init(index_b)/surrogate[config.iK].MM/conc_org(index_b)*1000.*surrogate[config.iK].charge;

  //cout << "final: " << inorganion << " " << organion << endl;

  derivative=derivative*0.5*(1.0+(organion+inorganion)/pow(pow(organion+inorganion,2)+4*config.Ke,0.5))-1.0;
  chp_new=0.5*(organion+inorganion+pow(pow(organion+inorganion,2)+4*config.Ke,0.5));
  error=chp_new-chp(index_b);
  //if (deltat>0)
  //  cout << "chp bins la " << chp_new << " " << inorganion << " " << chp(index_b) << endl;
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

  conc_org=LWC(index_b);
  for (i=0;i<n;++i)
    if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophilic)
      conc_org+=surrogate[i].Aaq_bins_init(index_b);
  conc_org=max(conc_org,1.e-5*AQinit(index_b)); //config.MOmin);

  derivative=0.0;
  for (i=0;i<n;++i)
    if (surrogate[i].is_inorganic_precursor and surrogate[i].is_solid==false)
      {
        if (surrogate[i].name=="NH3")
	  {
	    inorganion-=surrogate[config.iNH4p].Aaq_bins_init(index_b)/surrogate[config.iNH4p].MM/conc_org*1000.*surrogate[config.iNH4p].charge;	   
	  }
        else if (surrogate[i].name=="HNO3")
	  {
	    inorganion-=surrogate[config.iNO3m].Aaq_bins_init(index_b)/surrogate[config.iNO3m].MM/conc_org*1000.*surrogate[config.iNO3m].charge;	    
	  }
        else if (surrogate[i].name=="HCl")
	  {
	    inorganion-=surrogate[config.iClm].Aaq_bins_init(index_b)/surrogate[config.iClm].MM/conc_org*1000.*surrogate[config.iClm].charge;	  
	  }
        else if (surrogate[i].name=="H2SO4")
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
	if (surrogate[i].aq_type=="monoacid")
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
	else if (surrogate[i].aq_type=="diacid")
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
  /* else 
     if (surrogate[i].name=="Na")
     {
     inorganion-=surrogate[i].Aaq_bins_init(index_b)/surrogate[i].MM/conc_org*1000.*surrogate[i].charge;
     }*/

  if (config.iNa>=0) inorganion-=surrogate[config.iNa].Aaq_bins_init(index_b)/surrogate[config.iNa].MM/conc_org*1000.*surrogate[config.iNa].charge;
  if (config.iMg>=0) inorganion-=surrogate[config.iMg].Aaq_bins_init(index_b)/surrogate[config.iMg].MM/conc_org*1000.*surrogate[config.iMg].charge;
  if (config.iCa>=0) inorganion-=surrogate[config.iCa].Aaq_bins_init(index_b)/surrogate[config.iCa].MM/conc_org*1000.*surrogate[config.iCa].charge;  
  if (config.iK>=0) inorganion-=surrogate[config.iK].Aaq_bins_init(index_b)/surrogate[config.iK].MM/conc_org*1000.*surrogate[config.iK].charge;

  derivative=derivative*0.5*(1.0+(organion+inorganion)/pow(pow(organion+inorganion,2)+4*config.Ke,0.5))-1.0;
  chp_new=0.5*(organion+inorganion+pow(pow(organion+inorganion,2)+4*config.Ke,0.5));
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
      for (i=0;i<n;i++)
        if (surrogate[i].is_organic or i==config.iH2O)
          conc_org(b)+=surrogate[i].Aaq_bins_init(b);
      conc_org(b)=max(conc_org(b),1.e-5*AQinit(b)); //config.MOmin);
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
		  compute_kp_aq_ssh(config, surrogate, Temperature, ionic, chp2, MMaq);
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
        }
      nh++;
    }

  for (i=0;i<n;++i)
    if (surrogate[i].is_inorganic_precursor and surrogate[i].is_solid==false)
      {
        if (surrogate[i].name=="H2SO4")
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
  for (i=0;i<n;++i)
    if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophilic)
      conc_org+=surrogate[i].Aaq_bins(index_b);
  conc_org=max(conc_org,1.e-5*AQ(index_b)); //config.MOmin);

  derivative=0.0;
  for (i=0;i<n;++i)
    if (surrogate[i].is_inorganic_precursor and surrogate[i].is_solid==false)
      {
        if (surrogate[i].name=="NH3")
	  inorganion-=surrogate[config.iNH4p].Aaq_bins(index_b)/surrogate[config.iNH4p].MM/conc_org*1000.0*surrogate[config.iNH4p].charge;         
        else if (surrogate[i].name=="HNO3")
	  inorganion-=surrogate[config.iNO3m].Aaq_bins(index_b)/surrogate[config.iNO3m].MM/conc_org*1000.0*surrogate[config.iNO3m].charge;     
        else if (surrogate[i].name=="HCl")
	  inorganion-=surrogate[config.iClm].Aaq_bins(index_b)/surrogate[config.iClm].MM/conc_org*1000.0*surrogate[config.iClm].charge;
        else if (surrogate[i].name=="H2SO4")
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
	if (surrogate[i].aq_type=="monoacid")
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
	else if (surrogate[i].aq_type=="diacid")
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
  /*
    else 
    if (surrogate[i].name=="Na")
    inorganion-=surrogate[i].Aaq_bins(index_b)/surrogate[i].MM/conc_org*1000.*surrogate[i].charge;*/

  if (config.iNa>=0) inorganion-=surrogate[config.iNa].Aaq_bins(index_b)/surrogate[config.iNa].MM/conc_org*1000.*surrogate[config.iNa].charge;
  if (config.iMg>=0) inorganion-=surrogate[config.iMg].Aaq_bins(index_b)/surrogate[config.iMg].MM/conc_org*1000.*surrogate[config.iMg].charge;
  if (config.iCa>=0) inorganion-=surrogate[config.iCa].Aaq_bins(index_b)/surrogate[config.iCa].MM/conc_org*1000.*surrogate[config.iCa].charge;  
  if (config.iK>=0) inorganion-=surrogate[config.iK].Aaq_bins(index_b)/surrogate[config.iK].MM/conc_org*1000.*surrogate[config.iK].charge;


  derivative=derivative*0.5*(1.0+(organion+inorganion)/pow(pow(organion+inorganion,2)+4*config.Ke,0.5))-1.0;
  chp_new=0.5*(organion+inorganion+pow(pow(organion+inorganion,2)+4*config.Ke,0.5));
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
      for (i=0;i<n;i++)
        if (surrogate[i].is_organic or i==config.iH2O)
          conc_org(b)+=surrogate[i].Aaq_bins(b);
      conc_org(b)=max(conc_org(b),1.e-5*AQ(b)); //config.MOmin);
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
		  compute_kp_aq_ssh(config, surrogate, Temperature, ionic, chp2, MMaq);
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
        }
      nh++;
    }

  for (i=0;i<n;++i)
    if (surrogate[i].is_inorganic_precursor and surrogate[i].is_solid==false)
      {
        if (surrogate[i].name=="H2SO4")
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

void activity_coefficients_dyn_aq_ssh(model_config &config, vector<species>& surrogate,
				      double &Temperature, Array<double, 1> &AQinit,
				      Array<double, 3> &MOinit,
				      Array<double, 1> &conc_inorganic,
				      Array<double, 1> &ionic, Array<double, 1> &ionic_organic,
				      Array<double, 1> &organion, Array<double, 1> &chp,
				      Array<double, 1> &LWC,Array<double, 1> &MMaq, double factor, double deltat, int index2)
{
  //compute the activity coefficients with UNIFAC (short range interactions) for the organic phase
  //MOW: mean molar mass of the organic phase
  int i;
  int n=surrogate.size();
  double XH2O;
  Array<double, 1> X_unifac,gamma_unifac;
  int b;
  double AQ=0.0;
  double error_h=1000.0;
  double derivative_h;
  Array <double, 1> chp_save;
  chp_save.resize(config.nbins);
  Array <double, 1> chp2;
  chp2.resize(config.nbins);
  int b2;
  int index=0;
  double chp_new;
  double conc_org;	  
  
  if (config.use_global_dynamic_parameters)
    {
      double MMaqtemp=0.0;
      double AQ=0.0;
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
        {
          conc_org+=LWC(b);	  
          for (i=0;i<n;i++)
            if (surrogate[i].is_organic or i==config.iH2O)
              conc_org+=surrogate[i].Aaq_bins_init(b);
        }
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
          for (i=0;i<n;++i)
            if (surrogate[i].name=="H")
              if (AQ>0.0)
                surrogate[i].Aaq_bins_init(b)=AQinit(b)/AQ*surrogate[i].Aaq;
              else
                surrogate[i].Aaq_bins_init(b)=0.0;
        }	  
    }
  else
    {
      AQ=0.0;
      for (b=0;b<config.nbins;b++)
        chp_save(b)=chp(b);

      for (b=0;b<config.nbins;++b)
        {
          conc_inorganic(b)=0.0;
          for (i=0;i<n;i++)
            if (surrogate[i].is_inorganic_precursor==false and surrogate[i].is_organic==false and i!=config.iH2O)
              conc_inorganic(b)+=surrogate[i].Aaq_bins_init(b);
        } 
      
      for (b=0;b<config.nbins;++b)
        {
          AQ=conc_inorganic(b);
	  
          for (i=0;i<n;++i)
            if (surrogate[i].hydrophilic)	      
	      surrogate[i].Aaq=surrogate[i].Aaq_bins_init(b);	       	      
          surrogate[config.iH2O].Aaq+=LWC(b);

          for (i=0;i<n;++i)
            if (surrogate[i].hydrophilic)
              if(surrogate[i].is_organic or i==config.iH2O)
                AQ+=surrogate[i].Aaq;

          conc_org=LWC(b);	  
	  for (i=0;i<n;i++)
	    if (surrogate[i].is_organic or i==config.iH2O)
	      conc_org+=surrogate[i].Aaq_bins_init(b);
          conc_org=max(conc_org,1.e-5*AQinit(b)); //config.MOmin);
	  	
          //config.rho_aqueous=config.AQrho(b);          
          compute_ionic_strenght2_ssh(config, surrogate, Temperature, AQ, conc_inorganic(b), ionic(b), chp(b),
				      organion(b), ionic_organic(b), conc_org, factor);                    
          
          activity_coefficients_aq_ssh(config,surrogate,Temperature,0.0,MMaq(b),XH2O, conc_org);          
          if (config.compute_long_and_medium_range_interactions)
            activity_coefficients_LR_MR_ssh(config, surrogate, Temperature, 0.0, ionic(b));
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
        }

      if (factor>0.)
	for (b=0;b<config.nbins;++b)
	  {	  
	    chp2=chp(b);
	    if (config.compute_inorganic)
	      {
		error_h=1000.0;
		index=0;              
		while(abs(error_h/chp2(b))>1.0e-3 and index<1000)
		  {
		    index++;
		    
		    compute_kp_aq_ssh(config, surrogate, Temperature, ionic, chp2, MMaq);
		    if (config.imethod==0)
		      prodloss_aq_ssh(config, surrogate, AQinit, LWC, MOinit, index2, deltat);
		    error_ph_bins_ssh(config, surrogate, b, Temperature, chp2, error_h, derivative_h,AQinit,ionic,MMaq,LWC,chp_new,deltat);
		    if (chp2(b)==1.0e-14 and chp_new<=0.0)
		      error_h=0.0;
		    else if (chp_new <= 0.0 or chp_new!=chp_new)
		      {
			error_h=1.0e-14-chp2(b);
			chp2(b)=1.0e-14;
		      }
		    else
		      {
			if (chp2(b)-error_h/derivative_h>0.0 and derivative_h!=0.0)
			  chp2(b)=chp2(b)-error_h/derivative_h;
			else
			  {
			    if (chp2(b)+error_h<1.0e-14)
			      error_h=1.0e-14-chp2(b);
			    chp2(b)=max(chp2(b)+error_h,1.0e-14);
			  }
			  
		      }
		  }

		/*
		if (deltat>0.)
		  {
		    cout << "new ph: " << chp2(b) << " " << chp(b) << " " << index << " " << error_h << endl;
		    //exit(0);
		  }*/
		chp(b)=factor*min(chp2(b),100.)+(1.0-factor)*chp(b);	      

		conc_org=LWC(b);	  
		for (i=0;i<n;i++)
		  if (surrogate[i].is_organic or i==config.iH2O)
		    conc_org+=surrogate[i].Aaq_bins_init(b);
		conc_org=max(conc_org,1.e-5*AQinit(b)); //config.MOmin);
	      
		// Concentrations of H+ is not calculated to avoid numerical issues when pH given by ISORROPIA is too low
		surrogate[config.iHp].Aaq_bins_init(b)=chp(b)*conc_org/1000.0;
              
	      }	  
	  
	  }      
      compute_kp_aq_ssh(config, surrogate, Temperature, ionic, chp, MMaq);
      characteristic_time_aq_ssh(config, surrogate, Temperature, chp, LWC, AQinit, MOinit);
    } 

  compute_organion_ssh(config, surrogate, Temperature,  MMaq, AQinit, ionic,
		       chp, ionic_organic, organion, LWC);
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

  if (index==0)
    prodloss_aq_ssh(config, surrogate, AQinit, LWC, MOinit, 0, deltat);
    
  
  //cout << chp << endl;
  chp0=chp;
  chp_save=chp;

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
      conc_org=LWC(b);  
      for (i=0;i<n;++i)
        if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophilic)
          conc_org+=surrogate[i].Aaq_bins_init(b);
      conc_org=max(conc_org,1.e-5*AQinit(b)); //config.MOmin);

      double inorganion=0.0;
      for (i=0;i<n;++i)
        if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)
          {
            surrogate[i].molality=surrogate[i].Aaq_bins_init(b)/surrogate[i].MM/conc_org*1000.0;           
            if (surrogate[i].name!="H")
              inorganion-=surrogate[i].molality*surrogate[i].charge;
          }

      double chp_new=0.5*(organion(b)+inorganion+pow(pow(organion(b)+inorganion,2)+4*config.Ke,0.5));
      chp_new=max(min(chp_new,2.*chp(b)),0.5*chp(b));
      chp(b)=factor*chp_new+(1.0-factor)*chp(b);
      if (chp(b)==0.0)
        chp(b)=pow(10.0,-5.6);
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
					   MOinit, conc_inorganic, ionic, ionic_organic,
					   organion,chp,LWC,MMaq, factor, deltat, 1);
          compute_kp_aq_ssh(config, surrogate, Temperature, ionic, chp, MMaq);
        }    
    }

  //cout << chp << " " << index << endl;

  for (i=0;i<n;i++)
    {
      surrogate[i].kprod_gas=0.0;
      surrogate[i].kloss_gas=0.0;
    }

  //compute kinetic rates
  //prodloss_org_ssh(config, surrogate, MOinit, AQinit, tiny, index2);
  prodloss_aq_ssh(config, surrogate, AQinit, LWC, MOinit, 1, deltat);
  if (index==0)
    prodloss_org_ssh(config, surrogate, MOinit, AQinit, tiny, 0, deltat);
  else
    prodloss_org_ssh(config, surrogate, MOinit, AQinit, tiny, 1, deltat);

  if (config.chemistry)
    prodloss_chem_ssh(config, surrogate, MOinit, MOW, AQinit, MMaq, chp, tiny, 0);
  /*
    if (config.chemistry)
    compute_prodloss_chem_ssh(config,surrogate,MOinit,MOW,AQinit,MMaq,chp,DT2,tiny,index2);*/

  double apnew;
  for (i=0;i<n;++i)
      {
	if(surrogate[i].is_organic and config.compute_organic)
	  if(surrogate[i].hydrophobic)
            for (b=0;b<config.nbins;++b)
              for (ilayer=0;ilayer<config.nlayer;++ilayer)
                for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                  {		    	  		    
		    apnew=(surrogate[i].Ap_layer_init0(b,ilayer,iphase)
			   +surrogate[i].deltat_exp*surrogate[i].k1(b,ilayer,iphase,0)+(deltat-surrogate[i].deltat_exp)*
			   surrogate[i].kprod(b,ilayer,iphase))/(1.0+(deltat-surrogate[i].deltat_exp)*surrogate[i].kloss(b,ilayer,iphase));
		    if (surrogate[i].Ap_layer_init(b,ilayer,iphase)>tiny)
		      apnew=max(min(apnew,10.*surrogate[i].Ap_layer_init(b,ilayer,iphase)),0.1*surrogate[i].Ap_layer_init(b,ilayer,iphase));
		    surrogate[i].Ap_layer_init(b,ilayer,iphase)=apnew; //factor*apnew+(1.0-factor)*surrogate[i].Ap_layer_init(b,ilayer,iphase);
		  }

	if (LWCtot>config.LWClimit)
	  if((surrogate[i].is_organic and config.compute_organic) or ((i==config.iClm or i==config.iNO3m or i==config.iNH4p) and config.compute_inorganic))
	    if (surrogate[i].hydrophilic)
	      for (b=0;b<config.nbins;++b)
		{		    	  		    
		  apnew=(surrogate[i].Aaq_bins_init0(b)+surrogate[i].deltat_exp*surrogate[i].k1_aq(b,0)
			 +(deltat-surrogate[i].deltat_exp)*surrogate[i].kprod_aq(b))
		    /(1.0+(deltat-surrogate[i].deltat_exp)*surrogate[i].kloss_aq(b));
		  if (surrogate[i].Aaq_bins_init(b)>tiny)
		    apnew=max(min(apnew,10.*surrogate[i].Aaq_bins_init(b)),0.1*surrogate[i].Aaq_bins_init(b));
		  surrogate[i].Aaq_bins_init(b)=apnew; //factor*apnew+(1.0-factor)*surrogate[i].Aaq_bins_init(b);
		}

	if (LWCtot>config.LWClimit)
	  if (i==config.iSO4mm and config.compute_inorganic)
	    for (b=0;b<config.nbins;++b)
	      {
		double Keq=surrogate[config.iH2SO4].Kequilibrium_ssh(Temperature)/chp(b)*surrogate[config.iHSO4m].gamma_aq_bins(b)/
		  (surrogate[config.iHp].gamma_aq_bins(b)*surrogate[config.iSO4mm].gamma_aq_bins(b));	     
		total2=surrogate[config.iHSO4m].Aaq_bins_init(b)/surrogate[config.iHSO4m].MM+surrogate[config.iSO4mm].Aaq_bins_init(b)/surrogate[config.iSO4mm].MM;
		total=(surrogate[config.iHSO4m].Aaq_bins_init0(b)+surrogate[i].deltat_exp*surrogate[config.iHSO4m].k1_aq(b,0)
		       +(deltat-surrogate[i].deltat_exp)*surrogate[config.iHSO4m].kprod_aq(b))/
		  (1.0+(deltat-surrogate[i].deltat_exp)*surrogate[config.iHSO4m].kloss_aq(b))/surrogate[config.iHSO4m].MM;
		total+=(surrogate[config.iSO4mm].Aaq_bins_init0(b)+surrogate[i].deltat_exp*surrogate[config.iSO4mm].k1_aq(b,0)
			+(deltat-surrogate[i].deltat_exp)*surrogate[config.iSO4mm].kprod_aq(b))/
		  (1.0+(deltat-surrogate[i].deltat_exp)*surrogate[config.iSO4mm].kloss_aq(b))/surrogate[config.iSO4mm].MM;

		if (total2>tiny)
		  total=max(min(total,10.*total2),0.1*total2);	      
		surrogate[config.iHSO4m].Aaq_bins_init(b)=total*surrogate[config.iHSO4m].MM*1.0/(1.0+Keq); //(1.0-factor)*surrogate[config.iHSO4m].Aaq_bins_init(b)+factor*total*surrogate[config.iHSO4m].MM*1.0/(1.0+Keq);
		surrogate[config.iSO4mm].Aaq_bins_init(b)=total*surrogate[config.iSO4mm].MM*Keq/(1.0+Keq); //(1.0-factor)*surrogate[config.iSO4mm].Aaq_bins_init(b)+factor*total*surrogate[config.iSO4mm].MM*Keq/(1.0+Keq);
		
	      }
	
	if ((surrogate[i].is_organic and config.compute_organic) or (surrogate[i].is_inorganic_precursor and config.compute_inorganic))
	  {
	    apnew=(surrogate[i].Ag0+surrogate[i].deltat_exp*surrogate[i].k1_gas(0)+(deltat-surrogate[i].deltat_exp)*surrogate[i].kprod_gas)/
	      (1.0+(deltat-surrogate[i].deltat_exp)*surrogate[i].kloss_gas);
	    if (surrogate[i].Ag>tiny)
	      apnew=max(min(apnew,10.*surrogate[i].Ag),0.1*surrogate[i].Ag);
	    surrogate[i].Ag=apnew; //factor*apnew+(1.-factor)*surrogate[i].Ag;
	  }
      }

  /*
  for (i=0;i<n;i++)
    if (surrogate[i].Atot>0. and surrogate[i].name=="SOAlP")
      cout << surrogate[i].name << " " << surrogate[i].Ag << " " << surrogate[i].Ap_layer_init(0,config.nlayer-1,0) << " " << surrogate[i].kprod(0,config.nlayer-1,0) << " " << surrogate[i].kloss_gas << " " << config.diameters(0) << " " << surrogate[i].k1(0,0,0) << endl;*/

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

      if (config.compute_inorganic)
	for (i=0;i<n;i++)
	  if (surrogate[i].is_inorganic_precursor and surrogate[i].is_solid==false)
	    {
	      surrogate[i].Atot=surrogate[i].Ag;
	      if (surrogate[i].name=="NH3")
		surrogate[i].Atot+=sum(surrogate[config.iNH4p].Aaq_bins_init)/surrogate[config.iNH4p].MM*surrogate[i].MM; 
	      else if (surrogate[i].name=="HNO3")
		surrogate[i].Atot+=sum(surrogate[config.iNO3m].Aaq_bins_init)/surrogate[config.iNO3m].MM*surrogate[i].MM;
	      else if (surrogate[i].name=="HCl")
		surrogate[i].Atot+=sum(surrogate[config.iClm].Aaq_bins_init)/surrogate[config.iClm].MM*surrogate[i].MM; 	
	      else if (surrogate[i].name=="H2SO4")
		surrogate[i].Atot+=sum(surrogate[config.iSO4mm].Aaq_bins_init)/surrogate[config.iSO4mm].MM*surrogate[i].MM+
		  sum(surrogate[config.iHSO4m].Aaq_bins_init)/surrogate[config.iHSO4m].MM*surrogate[i].MM;
	    }
    } 

  if (config.chemistry==false)
    if (config.compute_organic)
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
      
      if (surrogate[config.iH2O].hydrophobic and config.compute_organic)
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

  /*
  for (i=0;i<n;i++)
    if (surrogate[i].name=="POAlP")
      {
	cout << surrogate[i].name << " " << surrogate[i].Ap_layer_init << " " << surrogate[i].Ag << endl;
	cout << surrogate[i].kprod << endl;
	cout << surrogate[i].kloss << endl;
	cout << surrogate[i].kprod_gas << " " << surrogate[i].kloss_gas << endl;
      }*/

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
    for (b=0;b<config.nbins;++b)
      for (ilayer=0;ilayer<config.nlayer;++ilayer)
        for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
          MO(b,ilayer,iphase)=MOinit(b,ilayer,iphase);
}


void equilibrium_inorg_ssh(model_config &config, vector<species>& surrogate, double &tequilibrium,
			   Array <double, 1> &AQinit, Array <double, 1> &conc_inorganic, 
			   Array <double, 1> &chp, Array <double, 1> &LWC,
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
    {
      conc_org(b)=LWC(b);  
      for (i=0;i<n;++i)
        if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophilic)
          conc_org(b)+=surrogate[i].Aaq_bins_init(b);
      conc_org(b)=max(conc_org(b),1.e-5*AQinit(b)); //config.MOmin);
    }
  

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
              if (surrogate[i].name=="H2SO4")
                conc_equilibrium+=surrogate[config.iHSO4m].Aaq_bins_init(b)/surrogate[config.iHSO4m].MM*surrogate[i].MM
                  +surrogate[config.iSO4mm].Aaq_bins_init(b)/surrogate[config.iSO4mm].MM*surrogate[i].MM;
              else if (surrogate[i].name=="NH3")
                conc_equilibrium+=surrogate[config.iNH4p].Aaq_bins_init(b)/surrogate[config.iNH4p].MM*surrogate[i].MM;
              else if (surrogate[i].name=="HNO3")
                conc_equilibrium+=surrogate[config.iNO3m].Aaq_bins_init(b)/surrogate[config.iNO3m].MM*surrogate[i].MM;
              else if (surrogate[i].name=="HCl")
                conc_equilibrium+=surrogate[config.iClm].Aaq_bins_init(b)/surrogate[config.iClm].MM*surrogate[i].MM;
                  
              //sum of Aaq/Ag ratios + 1
              sum+=surrogate[i].Kaq(b)*conc_org(b);
            }

        //compute gas phase concentrations
        if (is_equilibrium)
          {	
            
            if (surrogate[i].name=="H2SO4")
              surrogate[i].Ag=0.0;
            else
              surrogate[i].Ag=factor*conc_equilibrium/sum+(1.0-factor)*surrogate[i].Ag;
            
            for (b=0;b<config.nbins;++b)
              {
                if (surrogate[i].time_aq(b)<tequilibrium)
                  {
                    if (surrogate[i].name=="NH3")
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
                    else if (surrogate[i].name=="HNO3")
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
                    else if (surrogate[i].name=="HCl")
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

        if (surrogate[i].name=="H2SO4")
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
      
      compute_kp_aq_ssh(config, surrogate, Temperature, ionic, chp, MMaq);
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
    equilibrium_inorg_ssh(config, surrogate, tequilibrium, AQinit, conc_inorganic, chp, LWC,
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
          compute_kp_aq_ssh(config, surrogate, Temperature, ionic, chp, MMaq);
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
    equilibrium_inorg_ssh(config, surrogate, tequilibrium, AQinit, conc_inorganic, chp, LWC,
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
      conc_inorganic(b)=surrogate[config.iHp].Aaq_bins_init(b);
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
                    if (surrogate[i].name=="H2SO4")
                      conc_available+=(surrogate[config.iHSO4m].Aaq_bins_init(b)/surrogate[config.iHSO4m].MM
                                       +surrogate[config.iSO4mm].Aaq_bins_init(b)/surrogate[config.iSO4mm].MM)*
                        surrogate[i].MM;
                    else if (surrogate[i].name=="HNO3")
                      conc_available+=surrogate[config.iNO3m].Aaq_bins_init(b)/surrogate[config.iNO3m].MM*
                        surrogate[i].MM;
                    else if (surrogate[i].name=="NH3")
                      conc_available+=surrogate[config.iNH4p].Aaq_bins_init(b)/surrogate[config.iNH4p].MM*
                        surrogate[i].MM;
                    else if (surrogate[i].name=="HCl")
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
                  if (surrogate[i].name=="H2SO4")
                    {
                      tmp+=(surrogate[config.iHSO4m].Aaq_bins_init(b)/surrogate[config.iHSO4m].MM
                            +surrogate[config.iSO4mm].Aaq_bins_init(b)/surrogate[config.iSO4mm].MM)*
                        surrogate[i].MM;                    
                    }
                  else if (surrogate[i].name=="HNO3")                    
                    tmp+=surrogate[config.iNO3m].Aaq_bins_init(b)/surrogate[config.iNO3m].MM*
                      surrogate[i].MM;                     
                  else if (surrogate[i].name=="NH3")                    
                    tmp+=surrogate[config.iNH4p].Aaq_bins_init(b)/surrogate[config.iNH4p].MM*
                      surrogate[i].MM;                                       
                  else if (surrogate[i].name=="HCl")                    
                    tmp+=surrogate[config.iClm].Aaq_bins_init(b)/surrogate[config.iClm].MM*
                      surrogate[i].MM;                                         

                  if (conc_available-tmp> tiny and surrogate[i].k1_aq(b,0)>0.0)
                    surrogate[i].Jdn_aq(b,0)=-sum_rates/(conc_available-tmp);
                }

            if (config.chemistry)
              for (b=0;b<config.nbins;++b)			             
                {		                  
                  if (surrogate[i].name=="H2SO4")
                    {                     
                      surrogate[i].ktot1+=surrogate[config.iHSO4m].flux_chem_aq(b,0)/surrogate[config.iHSO4m].MM/surrogate[i].MM;
                      surrogate[i].ktot1+=surrogate[config.iSO4mm].flux_chem_aq(b,0)/surrogate[config.iSO4mm].MM/surrogate[i].MM;
                    }
                  else if (surrogate[i].name=="HNO3")
                    surrogate[i].ktot1+=surrogate[config.iNO3m].flux_chem_aq(b,0)/surrogate[config.iNO3m].MM/surrogate[i].MM;                    
                  else if (surrogate[i].name=="NH3")           
                    surrogate[i].ktot1+=surrogate[config.iNH4p].flux_chem_aq(b,0)/surrogate[config.iNH4p].MM/surrogate[i].MM;                    
                  else if (surrogate[i].name=="HCl")            
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
                if (surrogate[i].name=="H2SO4")
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

                if (surrogate[i].name=="NH3")
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
                
                if (surrogate[i].name=="HNO3")
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

                if (surrogate[i].name=="HCl")
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
                    if (surrogate[i].name=="H2SO4")
                      conc_available+=0.0; 
                    else if (surrogate[i].name=="HNO3")
                      conc_available+=surrogate[config.iNO3m].Aaq_bins_init(b)/surrogate[config.iNO3m].MM*
                        surrogate[i].MM;
                    else if (surrogate[i].name=="NH3")
                      conc_available+=surrogate[config.iNH4p].Aaq_bins_init(b)/surrogate[config.iNH4p].MM*
                        surrogate[i].MM;
                    else if (surrogate[i].name=="HCl")
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
                  if (surrogate[i].name=="H2SO4")
                    tmp+=(surrogate[config.iHSO4m].Aaq_bins_init(b)/surrogate[config.iHSO4m].MM
                          +surrogate[config.iSO4mm].Aaq_bins_init(b)/surrogate[config.iSO4mm].MM)*
                      surrogate[i].MM;
                  else if (surrogate[i].name=="HNO3")
                    tmp+=surrogate[config.iNO3m].Aaq_bins_init(b)/surrogate[config.iNO3m].MM*
                      surrogate[i].MM;
                  else if (surrogate[i].name=="NH3")
                    tmp+=surrogate[config.iNH4p].Aaq_bins_init(b)/surrogate[config.iNH4p].MM*
                      surrogate[i].MM;
                  else if (surrogate[i].name=="HCl")                    
                    tmp+=surrogate[config.iClm].Aaq_bins_init(b)/surrogate[config.iClm].MM*
                      surrogate[i].MM;                                        

                  if (conc_available-tmp> tiny and surrogate[i].k1_aq(b,1)>0.0)
                    surrogate[i].Jdn_aq(b,1)=-sum_rates/(conc_available-tmp);
                }

            if (config.chemistry)
              for (b=0;b<config.nbins;++b)			             
                {		                  
                  if (surrogate[i].name=="H2SO4")
                    {                     
                      surrogate[i].ktot2+=surrogate[config.iHSO4m].flux_chem_aq(b,1)/surrogate[config.iHSO4m].MM/surrogate[i].MM;
                      surrogate[i].ktot2+=surrogate[config.iSO4mm].flux_chem_aq(b,1)/surrogate[config.iSO4mm].MM/surrogate[i].MM;
                    }
                  else if (surrogate[i].name=="HNO3")
                    surrogate[i].ktot2+=surrogate[config.iNO3m].flux_chem_aq(b,1)/surrogate[config.iNO3m].MM/surrogate[i].MM;                    
                  else if (surrogate[i].name=="NH3")           
                    surrogate[i].ktot2+=surrogate[config.iNH4p].flux_chem_aq(b,1)/surrogate[config.iNH4p].MM/surrogate[i].MM;                    
                  else if (surrogate[i].name=="HCl")            
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
              if (surrogate[i].name=="H2SO4")
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
              else if (surrogate[i].name=="NH3")
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
              else if (surrogate[i].name=="HNO3")
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

              else if (surrogate[i].name=="HCl")
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
		     double &Temperature, double &DT2, double &tequilibrium,
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

  //for (i=config.iH2O+1;i<n;i++)
  //  cout << "Time :" << surrogate[i].name << " " << surrogate[i].time_aq(0) << endl;

  /*
    for (i=0;i<n;i++)
    if (surrogate[i].name=="BiDER")
    cout << surrogate[i].Atot+surrogate[surrogate[i].ioligo].Atot << " " << surrogate[surrogate[i].ioligo].Atot << " " << surrogate[surrogate[i].ioligo].Ag << " " << sum(surrogate[surrogate[i].ioligo].Aaq_bins_init) << endl;*/

  
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
              conc_org=LWC(b);
              for (i=0;i<n;++i)
                if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophilic)
                  conc_org+=surrogate[i].Aaq_bins_init(b);
              conc_org=max(conc_org,1.e-5*AQinit(b)); //config.MOmin);
              
              
	      //config.rho_aqueous=config.AQrho(b);
	      compute_ionic_strenght2_ssh(config, surrogate, Temperature, AQinit(b), conc_inorganic(b), ionic(b), chp(b),
					  organion(b), ionic_organic(b), conc_org, 1.0);
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
	      compute_kp_aq_ssh(config, surrogate, Temperature, ionic, chp, MMaq);   
	    }
	}
    }

  //compute kinetic rates
  flux_org_ssh(config, surrogate, MOinit, AQinit, tiny, 0);
  if (LWCtot>config.LWClimit)
    {
      flux_aq_ssh(config, surrogate, AQinit, LWC, MOinit, tiny, 0);      
      if (config.compute_inorganic and facph>0)
        correct_flux_ph_ssh(config, surrogate, Temperature, AQinit, MOinit, chp, chp2, MMaq, ionic, LWC, tiny, DT2/facph, 0);
    }

  if (config.chemistry)
    compute_flux_chem_ssh(config,surrogate,MOinit,MOW,AQinit,MMaq,chp,DT2,tiny,0);

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
              conc_org=LWC(b);
              for (i=0;i<n;++i)
                if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophilic)
                  conc_org+=surrogate[i].Aaq_bins_init(b);
              conc_org=max(conc_org,1.0e-5*AQinit(b)); //config.MOmin);

	      //config.rho_aqueous=config.AQrho(b);	      	      
	      compute_ionic_strenght2_ssh(config, surrogate, Temperature, AQinit2(b), conc_inorganic(b), ionic(b), chp(b),
					  organion(b), ionic_organic(b), conc_org, 1.0);	   
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
	      compute_kp_aq_ssh(config, surrogate, Temperature, ionic, chp, MMaq);
	    }
	}      
    }

  //compute the second evaluation of kinetic rates
  flux_org_ssh(config, surrogate, MOinit2, AQinit2, tiny, 1);
  if(LWCtot>config.LWClimit)
    {
      flux_aq_ssh(config, surrogate, AQinit2, LWC, MOinit2, tiny, 1);           
      if (config.compute_inorganic and facph>0)
	correct_flux_ph_ssh(config, surrogate, Temperature, AQinit2, MOinit2, chp, chp2, MMaq, ionic, LWC, tiny, DT2/facph, 1);
    }

  if (config.chemistry)
    compute_flux_chem_ssh(config,surrogate,MOinit2,MOW,AQinit2,MMaq,chp,DT2,tiny,1);

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
          conc_org=LWC(b);
          for (i=0;i<n;++i)
	    if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophilic)
	      conc_org+=surrogate[i].Aaq_bins(b);
          conc_org=max(conc_org,1.0e-5*AQinit(b)); //config.MOmin);


	  //config.rho_aqueous=config.AQrho(b);
	  compute_ionic_strenght2_ssh(config, surrogate, Temperature, AQ(b), conc_inorganic(b), ionic(b), chp(b),
				      organion(b), ionic_organic(b), conc_org, 1.0);	 	  
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
		    double &Temperature, double &DT2, double &tequilibrium,
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
              conc_org=LWC(b);
              for (i=0;i<n;++i)
                if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophilic)
                  conc_org+=surrogate[i].Aaq_bins_init(b);
              conc_org=max(conc_org,1.0e-5*AQinit(b)); //config.MOmin);

	      //config.rho_aqueous=config.AQrho(b);
	      compute_ionic_strenght2_ssh(config, surrogate, Temperature, AQinit(b), conc_inorganic(b), ionic(b), chp(b),
					  organion(b), ionic_organic(b), conc_org, 1.0);
	      error_chp=max(error_chp,abs(chp0(b)-chp(b))/chp0(b));
	      chp0(b)=chp(b);	      
	    }      
	}
	  
      if (compute_activity_coefficients)	         
	{
	  activity_coefficients_dyn_aq_ssh(config, surrogate, Temperature,AQinit,
					   MOinit,conc_inorganic, ionic, ionic_organic,
					   organion,chp,LWC,MMaq, 1.0, 0., 0);
	  compute_kp_aq_ssh(config, surrogate, Temperature, ionic, chp, MMaq);   	
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
  flux_aq_ssh(config, surrogate, AQinit, LWC, MOinit, tiny, 0);
  //if (config.compute_inorganic)
  //  correct_flux_ph_ssh(config, surrogate, Temperature, AQinit, MOinit, chp, chp2, MMaq, ionic, LWC, tiny, DT2, 0);
  if (config.chemistry)
    compute_flux_chem_ssh(config,surrogate,MOinit,MOW,AQinit,MMaq,chp,DT2,tiny,0);

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
	  conc_org=LWC(b);
	  for (i=0;i<n;++i)
	    if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophilic)
	      conc_org+=surrogate[i].Aaq_bins_init(b);
          conc_org=max(conc_org,1.0e-5*AQinit(b)); //config.MOmin);

	  inorganion=0.0;
	  for (i=0;i<n;++i)
	    if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)
	      {
		surrogate[i].molality=surrogate[i].Aaq_bins_init(b)/surrogate[i].MM/conc_org*1000.0;		
		if (surrogate[i].name!="H")
		  inorganion-=surrogate[i].molality*surrogate[i].charge;
	      }
	  
	  chp(b)=0.5*(organion(b)+inorganion+pow(pow(organion(b)+inorganion,2)+4*config.Ke,0.5));
	  if (chp(b)==0.0)
	    chp(b)=pow(10.0,-5.6);	  
	  
	  //config.rho_aqueous=config.AQrho(b);
	  compute_ionic_strenght2_ssh(config, surrogate, Temperature, AQinit2(b), conc_inorganic(b), ionic(b), chp(b),
				      organion(b), ionic_organic(b), conc_org, 1.0);

	}      
    }

  //compute activity coefficients
  if (compute_activity_coefficients)
    {
      activity_coefficients_dyn_aq_ssh(config, surrogate, Temperature,AQinit2,
				       MOinit,conc_inorganic, ionic, ionic_organic,
				       organion,chp,LWC,MMaq, 1.0, 0., 0);
      compute_kp_aq_ssh(config, surrogate, Temperature, ionic, chp, MMaq);
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
              conc_org=LWC(b);
              for (i=0;i<n;++i)
                if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophilic)
                  conc_org+=surrogate[i].Aaq_bins_init(b);
              conc_org=max(conc_org,1.0e-5*AQinit(b)); //config.MOmin);

	      //config.rho_aqueous=config.AQrho(b);	      	      
	      compute_ionic_strenght2_ssh(config, surrogate, Temperature, AQinit2(b), conc_inorganic(b), ionic(b), chp(b),
					  organion(b), ionic_organic(b), conc_org, 1.0);	   
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
	  compute_kp_aq_ssh(config, surrogate, Temperature, ionic, chp, MMaq);	    
	}      
    }  

  //second estimation of kinetic rates
  flux_aq_ssh(config, surrogate, AQinit2, LWC, MOinit, tiny, 1);
  //if (config.compute_inorganic)
  //  correct_flux_ph(config, surrogate, Temperature, AQinit2, MOinit, chp, chp2, MMaq, ionic, LWC, tiny, DT2, 1);
  if (config.chemistry)
    compute_flux_chem_ssh(config,surrogate,MOinit,MOW,AQinit,MMaq,chp,DT2,tiny,1);

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
          double conc_org=LWC(b);
	  for (i=0;i<n;++i)
	    if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophilic)
	      conc_org+=surrogate[i].Aaq_bins(b);
          conc_org=max(conc_org,1.0e-5*AQinit(b)); //config.MOmin);

	  compute_ionic_strenght2_ssh(config, surrogate, Temperature, AQ(b), conc_inorganic(b), ionic(b), chp(b),
				      organion(b), ionic_organic(b), conc_org,1.0);	 	  
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

	  if (surrogate[i].is_organic==false and surrogate[i].is_inorganic_precursor==false and config.iH2O!=i and config.iHp!=i)
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


