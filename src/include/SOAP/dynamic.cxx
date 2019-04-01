using namespace soap;

void compute_kp_org(model_config &config, vector<species>& surrogate,
                    Array <double, 3> &MOinit, double &Temperature, Array<double, 3> &MOW)
{
  //compute partitioning constant of the organic phase by taking into account activity
  //coefficients and the kelvin effect
  double temp1,temp2,maxi,MOWsurf;
  double kelvin_effect;
  int ilayer,i,b,iphase,jphase;
  int n=surrogate.size();
  for (i=0;i<n;++i)
    for (b=0;b<config.nbins;++b)
      if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophobic)
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
              temp2+=MOinit(b,config.nlayer-1,iphase)/MOW(b,config.nlayer-1,iphase);
            }
          if (temp1>0.0)
            MOWsurf=temp1/temp2;
          else
            MOWsurf=200.0;
		  
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
                        surrogate[i].Kp(b,ilayer,iphase)=surrogate[i].Kp_exp_org(Temperature);
                      else
                        surrogate[i].Kp(b,ilayer,iphase)=0.0;
                    else
                      surrogate[i].Kp(b,ilayer,iphase)=
                        surrogate[i].Kp_eff_org(Temperature,MOW(b,ilayer,iphase))/
                        surrogate[i].gamma_org_layer(b,ilayer,iphase);
					
                    if (config.compute_kelvin_effect and config.diameters(b) > 0.0) //compute the kelvin_effect
                      {
                        kelvin_effect=exp(2.0*config.surface_tension_org*1.0e-6*
                                          MOWsurf/(8.314*Temperature*config.rho_organic*
                                                   0.5*config.diameters(b)*1.0e-6));
                        surrogate[i].Kp(b,ilayer,iphase)/=kelvin_effect;
                      }
					
                  }
              else
                for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                  {
                    surrogate[i].Kp(b,ilayer,iphase)=
                      surrogate[i].Kp_eff_org(Temperature,MOW(b,ilayer,iphase))/
                      surrogate[i].gamma_org_layer(b,ilayer,iphase);
					
                    if (config.compute_kelvin_effect and config.diameters(b) > 0.0) //compute the kelvin_effect
                      {
                        kelvin_effect=exp(2.0*config.surface_tension_org*1.0e-6*
                                          MOWsurf/(8.314*Temperature*config.rho_organic*
                                                   0.5*config.diameters(b)*1.0e-6));
                        surrogate[i].Kp(b,ilayer,iphase)/=kelvin_effect;
                      }
					
                  }
            }
        }
}

void compute_kp_aq(model_config &config, vector<species>& surrogate,
                   double &Temperature, Array <double, 1> &ionic,
                   Array <double, 1> &chp,Array<double, 1> &MMaq)
{
  //compute partitioning constant of the aqueous phase by taking into account activity
  //coefficients and the kelvin effect
  double kelvin_effect;
  int b,i;
  int n=surrogate.size();
  double fion1,fion2;
  double R=8.314; //ideal gas constant (J/K/mol)
  double rho_h2o=1000.0; //volumic mass of H2O
  double MH2O=18.0;
  double deltaH_over_RT0,deltaCp0_over_R;
  double T0=298.15;

  for (i=0;i<n;++i)
    for (b=0;b<config.nbins;++b)
     if (MMaq(b) > 0.0)
      if(surrogate[i].is_organic and surrogate[i].hydrophilic)
        {	  
          surrogate[i].gamma_LR=surrogate[i].LR(b);
          if (surrogate[i].nonvolatile)
            surrogate[i].Kaq(b)=pow(10.0,3)/surrogate[i].gamma_aq_bins(b)*18.0/MMaq(b);
          else
            {
              surrogate[i].Kaq(b)=
                surrogate[i].Kp_eff_aq(config, Temperature, ionic(b), chp(b),surrogate[config.iHp].LR(b),
                                       surrogate[config.iHp].SRMR(b),MMaq(b),fion1,fion2)
                /surrogate[i].gamma_aq_bins(b);

              surrogate[i].fion1=fion1;
              surrogate[i].fion2=fion2;
            }
		  
          if (config.compute_kelvin_effect and config.diameters(b) > 0.0) //compute the kelvin effect
            {
              kelvin_effect=exp(2.0*config.surface_tension_aq*1.0e-6*MMaq(b)/
                                (8.314*Temperature*config.AQrho(b)*
                                 0.5*config.diameters(b)*1.0e-6));
              surrogate[i].Kaq(b)/=kelvin_effect;
            }
        }
      else if (surrogate[i].is_inorganic_precursor and config.compute_inorganic)
        {
          if (surrogate[i].name=="H2SO4")
            surrogate[i].Kaq(b)=1.0e10;          
          else if (surrogate[i].name=="NH3")
            {
              deltaH_over_RT0=13.79;
              deltaCp0_over_R=-5.39;
              surrogate[i].Kaq(b)=
                surrogate[i].Henry*exp(-deltaH_over_RT0*(T0/Temperature-1.0)-deltaCp0_over_R*(1+log(T0/Temperature)-T0/Temperature))
                *MH2O/MMaq(b)*R*Temperature/(rho_h2o*1.0e6*1.013e5)
                *(1.0+surrogate[i].Kequilibrium(Temperature)*chp(b)*surrogate[config.iHp].gamma_aq_bins(b)/surrogate[config.iNH4p].gamma_aq_bins(b));
              if (config.compute_kelvin_effect) //compute the kelvin effect
                {
                  kelvin_effect=exp(2.0*config.surface_tension_aq*1.0e-6*MMaq(b)/
                                    (8.314*Temperature*config.AQrho(b)*
                                     0.5*config.diameters(b)*1.0e-6));
                  surrogate[i].Kaq(b)/=kelvin_effect;
                }
            }
          else if (surrogate[i].name=="HNO3")
            {
              deltaH_over_RT0=29.17;
              deltaCp0_over_R=16.83;
              surrogate[i].Kaq(b)=
                surrogate[i].Henry*exp(-deltaH_over_RT0*(T0/Temperature-1.0)-deltaCp0_over_R*(1+log(T0/Temperature)-T0/Temperature))
                *R*Temperature/(rho_h2o*1.0e6*1.013e5)*MH2O/MMaq(b)
                *(1.0+surrogate[i].Kequilibrium(Temperature)/(chp(b)*surrogate[config.iHp].gamma_aq_bins(b)*surrogate[config.iNO3m].gamma_aq_bins(b)));
              if (config.compute_kelvin_effect) //compute the kelvin effect
                {
                  kelvin_effect=exp(2.0*config.surface_tension_aq*1.0e-6*MMaq(b)/
                                    (8.314*Temperature*config.AQrho(b)*
                                     0.5*config.diameters(b)*1.0e-6));
                  surrogate[i].Kaq(b)/=kelvin_effect;
                }
            }
          else if (surrogate[i].name=="HCl")
            {
              deltaH_over_RT0=30.20;
              deltaCp0_over_R=19.91;
              surrogate[i].Kaq(b)=surrogate[i].Henry*exp(-deltaH_over_RT0*(T0/Temperature-1.0)-deltaCp0_over_R*(1+log(T0/Temperature)-T0/Temperature))
                *R*Temperature/(rho_h2o*1.0e6*1.013e5)*MH2O/MMaq(b)
                *(1.0+surrogate[i].Kequilibrium(Temperature)/(chp(b)*surrogate[config.iHp].gamma_aq_bins(b)*surrogate[config.iClm].gamma_aq_bins(b)));
              if (config.compute_kelvin_effect) //compute the kelvin effect
                {
                  kelvin_effect=exp(2.0*config.surface_tension_aq*1.0e-6*MMaq(b)/
                                    (8.314*Temperature*config.AQrho(b)*
                                     0.5*config.diameters(b)*1.0e-6));
                  surrogate[i].Kaq(b)/=kelvin_effect;
                }
            }
	}
      else if (i==config.iH2O)
	{	  
	  surrogate[i].Kaq(b)=(760.0*8.202e-5*Temperature)/(MMaq(b)*1.0e6*surrogate[i].gamma_aq_bins(b)*surrogate[i].Psat(Temperature));
	  if (config.compute_kelvin_effect) //compute the kelvin effect
	    {
                  kelvin_effect=exp(2.0*config.surface_tension_aq*1.0e-6*MMaq(b)/
                                    (8.314*Temperature*config.AQrho(b)*
                                     0.5*config.diameters(b)*1.0e-6));
                  surrogate[i].Kaq(b)/=kelvin_effect;
	    }
	}
}

void characteristic_time(model_config &config, vector<species>& surrogate,
                         Array <double, 3> &MOinit, Array <double, 1> &AQinit,
                         double LWCtot)
{
  //compute characteristic time for each compound to reach an equilibrium between the gas phase and
  //the organic phase for each bin and each layer
  int i,b,ilayer,iphase,iphase2;
  int n=surrogate.size();
  double Kp,sum,sum2,sum3,sum4;
  double taumin;
  double time_dif;
  
  //compute total particulate concentrations
  if (config.explicit_representation)
    {      
      for (i=0;i<n;++i)
	if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophobic)
	  {
	    surrogate[i].Ap=0.0;
	    for (b=0;b<config.nbins;++b)
	      for (ilayer=0;ilayer<config.nlayer;++ilayer)
		for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		  surrogate[i].Ap+=surrogate[i].Ap_layer_init(b,ilayer,iphase);
	    for (b=0;b<config.nbins;++b)
	      {	      
	      time_dif=0.0;
	      double Vlayer_dif=0.0;
	
	      for (ilayer=config.nlayer-1;ilayer>=0;--ilayer)
		{
		  Vlayer_dif+=config.Vlayer(ilayer);
	  
		  sum=0.0;
		  for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		    sum+=surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase);

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
		    if(sum/sum3<config.kp_low_volatility)	     
		      for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)			
			{			    
			  surrogate[i].time(b,ilayer,iphase)=
			    (surrogate[i].tau_diffusion(b,ilayer,iphase)+sum/config.Vlayer(ilayer)*Vlayer_dif*surrogate[i].tau_air(b)*
			     (sum3/config.Vlayer(ilayer)+AQinit(b))/sum3*config.Vlayer(ilayer))
			    /(1.0+sum/sum2*surrogate[i].Ap);
			}
		    else 
		      for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
			surrogate[i].time(b,ilayer,iphase)=2.0*config.tequilibrium;
		  else
		    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		      surrogate[i].time(b,ilayer,iphase)=0.0;       		  
		}
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
	    surrogate[i].Ap=0.0;
	    for (b=0;b<config.nbins;++b)
	      for (ilayer=0;ilayer<config.nlayer;++ilayer)
		for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		  surrogate[i].Ap+=surrogate[i].Ap_layer_init(b,ilayer,iphase);

	    if (surrogate[i].Ap!=surrogate[i].Ap)
	      {
		cout << "nan" << "org " << surrogate[i].name << endl;
		throw string("NAN: org \"") + surrogate[i].name + "\".";
	      }
		
	    if (surrogate[i].hydrophilic and LWCtot>config.LWClimit)
	      for (b=0;b<config.nbins;++b)
		surrogate[i].Ap+=surrogate[i].Aaq_bins_init(b);

	    if (surrogate[i].Ap!=surrogate[i].Ap)
	      {
		cout << "nan" << "AQ " << surrogate[i].name << endl;
		throw string("NAN: aq \"") + surrogate[i].name + "\".";
	      }
	  }

      //compute the characteristic time
      for (i=0;i<n;++i)
	if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophobic)
	  for (b=0;b<config.nbins;++b)	  
	    for (ilayer=0;ilayer<config.nlayer;++ilayer)
	      {
		sum=0.0;
		for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		  sum+=surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase);

		sum2=0.0;
		for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		  sum2+=surrogate[i].Ap_layer_init(b,ilayer,iphase);

		sum3=0.0;
		for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		  sum3+=MOinit(b,ilayer,iphase);

		if (sum2>0.0)
		  if(sum/sum3<config.kp_low_volatility)
		    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		      {
			surrogate[i].time(b,ilayer,iphase)=
			  (surrogate[i].tau_diffusion(b,ilayer,iphase)+
			   sum/config.Vlayer(ilayer)*surrogate[i].tau_air(b)*
			   (sum3/config.Vlayer(ilayer)+AQinit(b))/sum3*config.Vlayer(ilayer))
			  /(1.0+sum/sum2*surrogate[i].Ap);

			if (i==config.iH2O)
			  {			    
			    double a1=sum/sum3;
			    double b1=1.0-sum-a1*surrogate[i].Atot;
			    double c1=-(sum3-sum2);
			    double delta=pow(b1,2.0)-4.0*a1*c1;			    
			    double MO2=max(sum3,(-b1+pow(delta,0.5))/(2.0*a1));			 
			    double f1=pow((AQinit(b)+MO2/config.Vlayer(ilayer))/(AQinit(b)+MOinit(b,ilayer,iphase)/config.Vlayer(ilayer)),1.0/3.0);			    
			    surrogate[i].time(b,ilayer,iphase)=
			      (surrogate[i].tau_diffusion(b,ilayer,iphase)+
			       MO2/sum3*sum/config.Vlayer(ilayer)*surrogate[i].tau_air(b)/f1*
			       (MO2/config.Vlayer(ilayer)+AQinit(b))/MO2*config.Vlayer(ilayer))
			      /(1.0+MO2/sum3*sum/sum2*surrogate[i].Ap);
			  }

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

void characteristic_time_aq(model_config &config, vector<species>& surrogate,
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
        else if(surrogate[i].is_inorganic_precursor and config.compute_inorganic)
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
      double conc_aq,sumMO;
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
		    double AQ2=AQinit(b); 
		    surrogate[i].time_aq(b)=
		      (AQ2+sumMO)/AQ2*
		      (Kaq*AQ2*surrogate[i].tau_air(b))
		      /(1.0+Kaq*AQ2/surrogate[i].Aaq_bins_init(b)*surrogate[i].Aaq);		    
		    if (i==config.iH2O)
		      {		       
			double a1=Kaq;
			double b1=1.0-Kaq*AQ2-Kaq*surrogate[i].Atot;
			double c1=-(AQinit(b)-surrogate[i].Aaq_bins_init(b));
			double delta=pow(b1,2.0)-4.0*a1*c1;
			
			AQ2=max(AQinit(b),(-b1+pow(delta,0.5))/(2.0*a1));
			double f1=pow((AQ2+sumMO)/(AQinit(b)+sumMO),1.0/3.0);			
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
        else if (surrogate[i].is_inorganic_precursor and config.compute_inorganic)
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
  
  for (i=0;i<n;i++)
    if (surrogate[i].aq_type=="monoacid" or surrogate[i].aq_type=="diacid")
      for (b=0;b<config.nbins;b++)         	
	max_time_inorg(b)=max(max_time_inorg(b),(surrogate[i].fion1+surrogate[i].fion2)*surrogate[i].time_aq(b));       

  for (b=0;b<config.nbins;++b)
    surrogate[config.iH2O].time_aq(b)=0.0; //H2O is forced to be at equilibrium
    
  for (i=0;i<n;i++)
    if (surrogate[i].aq_type=="monoacid" or surrogate[i].aq_type=="diacid")
      for (b=0;b<config.nbins;b++)         
	{	  
	  surrogate[i].time_aq(b)=(1.0-surrogate[i].fion1-surrogate[i].fion2)*surrogate[i].time_aq(b)
	    +(surrogate[i].fion1+surrogate[i].fion2)*max_time_inorg(b);
	}
}

void flux_aq(model_config &config, vector<species>& surrogate, Array<double, 1> &AQinit, Array<double, 3> &MOinit, double &tiny, int index)
{
  //compute kinetic rates for the absorption of a compound in the aqueous phase in a bin
  //index = 0 : first evaluation of rates
  //index = 1 : second evaluation of rates
  double Kaq;
  int n=surrogate.size();
  int i,b,ilayer,iphase;
  double sum_mass;
  
  for (i=0;i<n;++i)
    if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophilic)
      for (b=0;b<config.nbins;++b)
        {
          sum_mass=AQinit(b);
          for (ilayer=0;ilayer<config.nlayer;ilayer++)
            for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
              sum_mass+=MOinit(b,ilayer,iphase);
          
          Kaq=surrogate[i].Kaq(b);          
	  surrogate[i].k1_aq(b,index)=AQinit(b)/sum_mass*
	    (surrogate[i].Ag*Kaq*AQinit(b)-surrogate[i].Aaq_bins_init(b))/
	    (Kaq*AQinit(b)*surrogate[i].tau_air(b));
		  
          surrogate[i].Jdn_aq(b,index)=0.0;
          if (surrogate[i].k1_aq(b,index)<0.0 and surrogate[i].Aaq_bins_init(b)>tiny)
            surrogate[i].Jdn_aq(b,index)=surrogate[i].k1_aq(b,index)/
	    surrogate[i].Aaq_bins_init(b);
        }
    else if (surrogate[i].is_inorganic_precursor and config.compute_inorganic) //for inorganic compounds
      for (b=0;b<config.nbins;++b)
        {
	  sum_mass=AQinit(b);
          for (ilayer=0;ilayer<config.nlayer;ilayer++)
            for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
              sum_mass+=MOinit(b,ilayer,iphase);          

          double conc_aq;
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
	      
	      surrogate[i].k1_aq(b,index)=(surrogate[i].Ag*Kaq*AQinit(b)-conc_aq)/(Kaq*AQinit(b)*surrogate[i].tau_air(b))*AQinit(b)/sum_mass;	   
 
              surrogate[i].Jdn_aq(b,index)=0.0;
              if (surrogate[i].k1_aq(b,index)<0.0 and surrogate[config.iNH4p].Aaq_bins_init(b)>tiny)
                surrogate[i].Jdn_aq(b,index)=surrogate[i].k1_aq(b,index)/surrogate[config.iNH4p].Aaq_bins_init(b)*surrogate[config.iNH4p].MM/surrogate[i].MM;
            }
          else if (surrogate[i].name=="HNO3")
            {
              Kaq=surrogate[i].Kaq(b);
              //compute kinetic rate of absorption
              conc_aq=surrogate[config.iNO3m].Aaq_bins_init(b)/surrogate[config.iNO3m].MM*surrogate[i].MM;
	      
	      surrogate[i].k1_aq(b,index)=(surrogate[i].Ag*Kaq*AQinit(b)-conc_aq)/(Kaq*AQinit(b)*surrogate[i].tau_air(b))*AQinit(b)/sum_mass;     
              
              surrogate[i].Jdn_aq(b,index)=0.0;
              if (surrogate[i].k1_aq(b,index)<0.0 and surrogate[config.iNO3m].Aaq_bins_init(b)>tiny)
                surrogate[i].Jdn_aq(b,index)=surrogate[i].k1_aq(b,index)/surrogate[config.iNO3m].Aaq_bins_init(b)*surrogate[config.iNO3m].MM/surrogate[i].MM;
            }
          else if (surrogate[i].name=="HCl")
            {
              Kaq=surrogate[i].Kaq(b);
              //compute kinetic rate of absorption
              conc_aq=surrogate[config.iClm].Aaq_bins_init(b)/surrogate[config.iClm].MM*surrogate[i].MM;
	      
	      surrogate[i].k1_aq(b,index)=(surrogate[i].Ag*Kaq*AQinit(b)-conc_aq)/(Kaq*AQinit(b)*surrogate[i].tau_air(b))*AQinit(b)/sum_mass;       

              surrogate[i].Jdn_aq(b,index)=0.0;
              if (surrogate[i].k1_aq(b,index)<0.0 and surrogate[config.iClm].Aaq_bins_init(b)>tiny)
                surrogate[i].Jdn_aq(b,index)=surrogate[i].k1_aq(b,index)/surrogate[config.iClm].Aaq_bins_init(b)*surrogate[config.iClm].MM/surrogate[i].MM;
            }


          
        }
  
  for (i=0;i<n;++i)      
    if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophilic)
      {
        double sumkpositive=0.0;
        double sumknegative=0.0;
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


void correct_flux_ph(model_config &config, vector<species>& surrogate, double &Temperature, Array<double, 1> &AQinit, Array <double, 3> &MOinit, 
		     Array<double, 1> &chp, 		     
		     Array<double, 1> &chpout, Array<double, 1> &MMaq, Array<double, 1> &ionic, Array<double, 1> &LWC, double tiny, double deltat, int index)
{
  int i,b;
  int n=surrogate.size();
  double a1,b1,c1,asulf; 
  double a2,b2,c2,d2;
  double Kaq,conc_aq;
  double Ke=1.0e-14;
  double Jhp,chp2;
  double epser=0.001;
  double tmini=config.deltatmin/100.;
  double t;
  double dt1;
  double delta,chpeq1,chpeq2;
  double organion,organion2,kelvin_effect;
  double sum_mass,faq;
  int ilayer,iphase;

  for (b=0;b<config.nbins;b++)
    {
      t=0;
      a1=0.0;
      b1=0.0;
      c1=0.0;
      Jhp=0.0;

      double conc_org=LWC(b);
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
	  kelvin_effect=exp(2.0*config.surface_tension_aq*1.0e-6*MMaq(b)/
			    (8.314*Temperature*config.AQrho(b)*
			     0.5*config.diameters(b)*1.0e-6));
	}

      for (i=0;i<n;++i)
	if (surrogate[i].is_inorganic_precursor)
	  {
	    if (surrogate[i].name=="H2SO4")	      
	      {
		
		double f=1.0;
		if(surrogate[config.iSO4mm].Aaq_bins_init(b)/surrogate[config.iSO4mm].MM+
		   surrogate[config.iHSO4m].Aaq_bins_init(b)/surrogate[config.iHSO4m].MM>0.0)
		  f=surrogate[config.iSO4mm].Aaq_bins_init(b)/surrogate[config.iSO4mm].MM/
		    (surrogate[config.iSO4mm].Aaq_bins_init(b)/surrogate[config.iSO4mm].MM+
		     surrogate[config.iHSO4m].Aaq_bins_init(b)/surrogate[config.iHSO4m].MM);	       
		a1+=(1.+f)*surrogate[i].Ag/surrogate[i].tau_air(b)/surrogate[i].MM; 		
		surrogate[i].Agt=surrogate[i].Ag;
	      }
	    else if (surrogate[i].name=="NH3")
	      {		
		Kaq=surrogate[i].Kaq(b)/chp(b);	       
		conc_aq=surrogate[config.iNH4p].Aaq_bins_init(b)/surrogate[config.iNH4p].MM;
		c1+=conc_aq/(Kaq*AQinit(b)*surrogate[i].tau_air(b))*faq;
		a1-=surrogate[i].Ag/surrogate[i].tau_air(b)/surrogate[i].MM*faq;
		surrogate[i].Agt=surrogate[i].Ag;	      
	      }
	    else if (surrogate[i].name=="HNO3")
	      {
		Kaq=surrogate[i].Kaq(b)*chp(b);	       
		conc_aq=surrogate[config.iNO3m].Aaq_bins_init(b)/surrogate[config.iNO3m].MM;
		b1-=conc_aq/(Kaq*AQinit(b)*surrogate[i].tau_air(b))*faq;
		a1+=surrogate[i].Ag/surrogate[i].tau_air(b)/surrogate[i].MM*faq;	       
		surrogate[i].Agt=surrogate[i].Ag;
	      }
	    else if (surrogate[i].name=="HCl")
	      {
		Kaq=surrogate[i].Kaq(b)*chp(b);	       
		conc_aq=surrogate[config.iClm].Aaq_bins_init(b)/surrogate[config.iClm].MM;
		b1-=conc_aq/(Kaq*AQinit(b)*surrogate[i].tau_air(b))*faq;
		a1+=surrogate[i].Ag/surrogate[i].tau_air(b)/surrogate[i].MM*faq;	     
		surrogate[i].Agt=surrogate[i].Ag;
	      }
	    
	  }
	else if (surrogate[i].is_organic and surrogate[i].time_aq(b)>=config.tequilibrium)
	  if (surrogate[i].aq_type=="monoacid" or surrogate[i].aq_type=="diacid")
	    {	      
	      surrogate[i].Agt=surrogate[i].Ag;	     
	      surrogate[i].Aaq=0.0;
	    }	           

      chp2=chp(b);
      a2=0.0;
      b2=0.0;
      c2=0.0;
      d2=0.0;
      chpeq1=0.0;
      chpeq2=0.0;
     
      if (b1!=0.0 or c1!=0.0)
	while (t<deltat)
	  {
	    a1=0.0;
	    b1=0.0;
	    c1=0.0;
	    asulf=0.0;

	    organion=0.0;
	    organion2=0.0;
	    for (i=0;i<n;++i)
	      if (surrogate[i].is_organic and surrogate[i].time_aq(b)>=config.tequilibrium)
		{
		  if (surrogate[i].aq_type=="monoacid")
		    {
		      double fion1=0.0;
		      double fion2=0.0;
		      surrogate[i].gamma_LR=surrogate[i].LR(b);
		      double Kaq=surrogate[i].Kp_eff_aq(config, Temperature, ionic(b), chp2,surrogate[config.iHp].LR(b),
							surrogate[config.iHp].SRMR(b),MMaq(b), fion1, fion2)/kelvin_effect/surrogate[i].gamma_aq_bins(b);		      
		      organion=fion1*(surrogate[i].Agt-(surrogate[i].Aaq_bins_init(b)+surrogate[i].Ag-surrogate[i].Agt)/(Kaq*AQinit(b)))
			/surrogate[i].MM/conc_org*1000./surrogate[i].tau_air(b)*faq;			
		      organion2+=fion1*(1.0-fion1)*(surrogate[i].Aaq_bins_init(b)+surrogate[i].Ag-surrogate[i].Agt)/surrogate[i].MM/conc_org*1000./chp2;      
		    }	
		  else if (surrogate[i].aq_type=="diacid")
		    {
		      double fion1=0.0;
		      double fion2=0.0;
		      surrogate[i].gamma_LR=surrogate[i].LR(b);
		      double Kaq=surrogate[i].Kp_eff_aq(config, Temperature, ionic(b), chp2,surrogate[config.iHp].LR(b),
							surrogate[config.iHp].SRMR(b),MMaq(b), fion1, fion2)/kelvin_effect/surrogate[i].gamma_aq_bins(b);
		      
		      organion=(fion1+2.0*fion2)*(surrogate[i].Agt-(surrogate[i].Aaq_bins_init(b)+surrogate[i].Ag-surrogate[i].Agt)/(Kaq*AQinit(b)))
			/surrogate[i].MM/conc_org*1000./surrogate[i].tau_air(b)*faq;
		      organion2+=(fion1*(1.0-fion1)+4.0*fion2*(1.0-fion2)-4.0*fion2*fion1)*
			(surrogate[i].Aaq_bins_init(b)+surrogate[i].Ag-surrogate[i].Agt)/surrogate[i].MM/conc_org*1000./chp2;
		    }

		}
	      else if (surrogate[i].is_inorganic_precursor)
		{
		  if (surrogate[i].name=="H2SO4")	      
		    {
		      double K=surrogate[i].Kequilibrium(Temperature)*surrogate[config.iHSO4m].gamma_aq_bins(b)
			/(surrogate[config.iHp].gamma_aq_bins(b)*surrogate[config.iSO4mm].gamma_aq_bins(b));
		      a1+=(2.0-1.0/(1.0+K/chp2))*surrogate[i].Agt/surrogate[i].tau_air(b)/surrogate[i].MM; 		
		      
		      double total=surrogate[config.iSO4mm].Aaq_bins_init(b)/surrogate[config.iSO4mm].MM+
			surrogate[config.iHSO4m].Aaq_bins_init(b)/surrogate[config.iHSO4m].MM;	       		     
		      asulf=config.AQrho(b)*total/AQinit(b)*1.0/pow(1.0+K/chp2,2.0)*K/(chp2*chp2);
		    }
		  else if (surrogate[i].name=="NH3")
		    {		     
		      Kaq=surrogate[i].Kaq(b)/chp(b);	       
		      conc_aq=surrogate[config.iNH4p].Aaq_bins_init(b)/surrogate[config.iNH4p].MM+
			(surrogate[i].Ag-surrogate[i].Agt)/surrogate[i].MM;
		      c1+=conc_aq/(Kaq*AQinit(b)*surrogate[i].tau_air(b))*faq;
		      a1-=surrogate[i].Agt/surrogate[i].tau_air(b)/surrogate[i].MM*faq;
		    }
		  else if (surrogate[i].name=="HNO3")
		    {
		      Kaq=surrogate[i].Kaq(b)*chp(b);	       
		      conc_aq=surrogate[config.iNO3m].Aaq_bins_init(b)/surrogate[config.iNO3m].MM+
			(surrogate[i].Ag-surrogate[i].Agt)/surrogate[i].MM;
		      b1-=conc_aq/(Kaq*AQinit(b)*surrogate[i].tau_air(b))*faq;
		      a1+=surrogate[i].Agt/surrogate[i].tau_air(b)/surrogate[i].MM*faq;
		    }
		  else if (surrogate[i].name=="HCl")
		    {
		      Kaq=surrogate[i].Kaq(b)*chp(b);	       
		      conc_aq=surrogate[config.iClm].Aaq_bins_init(b)/surrogate[config.iClm].MM+
			(surrogate[i].Ag-surrogate[i].Agt)/surrogate[i].MM;
		      b1-=conc_aq/(Kaq*AQinit(b)*surrogate[i].tau_air(b))*faq;
		      a1+=surrogate[i].Agt/surrogate[i].tau_air(b)/surrogate[i].MM*faq;
		    }
		}

	    a1=a1/conc_org*1000.0;
	    b1=b1/conc_org*1000.0;
	    c1=c1/conc_org*1000.0;

	    Jhp=(a1*chp2+b1*pow(chp2,2)+c1+organion*chp2)/(Ke/chp2+chp2+organion2*chp2+chp2*asulf);
            Jhp=min(max(Jhp,-0.9*chp2/tmini),9.*chp2/tmini);

	    delta=pow(a1+organion,2)-4*b1*c1;
	    if (delta>0.0)
	      {
		chpeq1=max((-a1-organion-pow(delta,0.5))/(2.0*b1),0.0);
		chpeq2=max((-a1-organion+pow(delta,0.5))/(2.0*b1),0.0);
	      }
	    else if (delta==0.0)
	      {
		chpeq1=max(-(a1+organion)/(2.0*b1),0.0);
		chpeq2=0.0;
	      }
	    else
	      {
		chpeq1=0.0;
		chpeq2=0.0;
	      }
	    
	    if (Jhp!=0.0)
	      {
		dt1=abs(epser*chp2/Jhp);
		if ((chpeq1-chp2)/Jhp>0.0)
		  dt1=min(dt1,0.1*(chpeq1-chp2)/Jhp);
		
		if ((chpeq2-chp2)/Jhp>0.0)
		  dt1=min(dt1,0.1*(chpeq2-chp2)/Jhp);

		for (i=0;i<n;i++)
		  if (surrogate[i].is_inorganic_precursor)
		    {
		      if (surrogate[i].Agt>0.0 and surrogate[i].k1_aq(b,index)>0.0)						 
			dt1=min(dt1,0.1*abs(surrogate[i].Agt/surrogate[i].k1_aq(b,index)));			
			 
		      double conc_aq=0.0;
		      if (surrogate[i].name=="NH3")
			{			    
			  conc_aq=surrogate[config.iNH4p].Aaq_bins_init(b)/surrogate[config.iNH4p].MM*surrogate[i].MM+
			    surrogate[i].Ag-surrogate[i].Agt;
			}
		      else if (surrogate[i].name=="HNO3")
			{			  
			  conc_aq=surrogate[config.iNO3m].Aaq_bins_init(b)/surrogate[config.iNO3m].MM*surrogate[i].MM+
			    surrogate[i].Ag-surrogate[i].Agt;
			}
		      else if (surrogate[i].name=="HCl")
			{		
			  conc_aq=surrogate[config.iClm].Aaq_bins_init(b)/surrogate[config.iClm].MM*surrogate[i].MM+
			    surrogate[i].Ag-surrogate[i].Agt;		  
			}

		      if (conc_aq>0.0 and surrogate[i].k1_aq(b,index)>0.0)
			dt1=min(dt1,epser*abs(conc_aq/surrogate[i].k1_aq(b,index)));
		    }
		dt1=max(dt1,tmini);
		dt1=min(dt1,deltat-t);
	      }
	    else
	      dt1=deltat-t;
	    
	    for (i=0;i<n;i++)	      
	      if (surrogate[i].name=="H2SO4")
		{                  
		  surrogate[i].k1_aq(b,index)=surrogate[i].Agt/surrogate[i].tau_air(b);		  
                  d2+=surrogate[i].Agt-max(surrogate[i].Agt-surrogate[i].k1_aq(b,index)*dt1,0.01*surrogate[i].Agt);
		  surrogate[i].Agt=max(surrogate[i].Agt-surrogate[i].k1_aq(b,index)*dt1,0.01*surrogate[i].Agt);                  
		}
	      else if (surrogate[i].name=="NH3")
		{
		  Kaq=surrogate[i].Kaq(b)/chp(b);	       
		  conc_aq=surrogate[config.iNH4p].Aaq_bins_init(b)/surrogate[config.iNH4p].MM*surrogate[i].MM+
		    surrogate[i].Ag-surrogate[i].Agt;
		  
		  surrogate[i].k1_aq(b,index)=(surrogate[i].Agt*Kaq*chp2*AQinit(b)-conc_aq)/
		    (Kaq*chp2*AQinit(b)*surrogate[i].tau_air(b))*faq;
		  double dJ=conc_aq/(Kaq*AQinit(b)*surrogate[i].tau_air(b))*faq*Jhp/(chp2*chp2);
		  a2+=surrogate[i].k1_aq(b,index)*dt1+0.5*dJ*dt1*dt1;
		  
		  surrogate[i].Agt=max(surrogate[i].Agt-surrogate[i].k1_aq(b,index)*dt1-0.5*dJ*dt1*dt1,0.01*surrogate[i].Agt);		  
		}
	      else if (surrogate[i].name=="HNO3")
		{
		  Kaq=surrogate[i].Kaq(b)*chp(b);	       
		  conc_aq=surrogate[config.iNO3m].Aaq_bins_init(b)/surrogate[config.iNO3m].MM*surrogate[i].MM+
		    surrogate[i].Ag-surrogate[i].Agt;
		  
		  surrogate[i].k1_aq(b,index)=(surrogate[i].Agt*Kaq/chp2*AQinit(b)-conc_aq)/
		    (Kaq/chp2*AQinit(b)*surrogate[i].tau_air(b))*faq;
		  double dJ=-conc_aq/(Kaq*AQinit(b)*surrogate[i].tau_air(b))*faq*Jhp;
		  b2+=surrogate[i].k1_aq(b,index)*dt1+0.5*dJ*dt1*dt1;
		  
		  surrogate[i].Agt=max(surrogate[i].Agt-surrogate[i].k1_aq(b,index)*dt1-0.5*dJ*dt1*dt1,0.01*surrogate[i].Agt);		  
		}
	      else if (surrogate[i].name=="HCl")
		{
		  Kaq=surrogate[i].Kaq(b)*chp(b);	       
		  conc_aq=surrogate[config.iClm].Aaq_bins_init(b)/surrogate[config.iClm].MM*surrogate[i].MM+
		    surrogate[i].Ag-surrogate[i].Agt;
		  
		  surrogate[i].k1_aq(b,index)=(surrogate[i].Agt*Kaq/chp2*AQinit(b)-conc_aq)/
		    (Kaq/chp2*AQinit(b)*surrogate[i].tau_air(b))*faq;
		  double dJ=-conc_aq/(Kaq*AQinit(b)*surrogate[i].tau_air(b))*faq*Jhp;
		  c2+=surrogate[i].k1_aq(b,index)*dt1+0.5*dJ*dt1*dt1;
				  
		  surrogate[i].Agt=max(surrogate[i].Agt-surrogate[i].k1_aq(b,index)*dt1-0.5*dJ*dt1*dt1,0.01*surrogate[i].Agt);		 		  
		}
	    
	      else if (surrogate[i].is_organic and surrogate[i].time_aq(b)>=config.tequilibrium)
		{

		  if (surrogate[i].aq_type=="monoacid")
		    {	      
		      double fion1=0.0;
		      double fion2=0.0;
		      surrogate[i].gamma_LR=surrogate[i].LR(b);
		      double Kaq=surrogate[i].Kp_eff_aq(config, Temperature, ionic(b), chp2,surrogate[config.iHp].LR(b),
							surrogate[config.iHp].SRMR(b),MMaq(b), fion1, fion2)/kelvin_effect/surrogate[i].gamma_aq_bins(b);
		      
		      double ratio_gamma=pow(surrogate[config.iHp].LR(b),2.0)*surrogate[config.iHp].SRMR(b)/surrogate[i].gamma_LR;
		      double Kac=surrogate[i].Kacidity1/ratio_gamma;
		      double Kh=surrogate[i].Kpart_aq(Temperature, MMaq(b))/kelvin_effect/surrogate[i].gamma_aq_bins(b);

		      surrogate[i].k1_aq(b,index)=faq*
			(surrogate[i].Agt*Kaq*AQinit(b)-surrogate[i].Aaq_bins_init(b)-surrogate[i].Ag+surrogate[i].Agt)/
			(Kaq*AQinit(b)*surrogate[i].tau_air(b));
		      double dJ=-surrogate[i].Aaq_bins_init(b)/pow(Kaq*AQinit(b),2)/surrogate[i].tau_air(b)*faq*
				  Kh*AQinit(b)*Kac/(chp2*chp2)*Jhp;
		      surrogate[i].Aaq+=surrogate[i].k1_aq(b,index)*dt1+0.5*dJ*dt1*dt1;
		      surrogate[i].Agt=max(surrogate[i].Agt-surrogate[i].k1_aq(b,index)*dt1-0.5*dJ*dt1*dt1,0.01*surrogate[i].Agt);
		      
		    }
		  else if (surrogate[i].aq_type=="diacid")
		    {	      
		      double fion1=0.0;
		      double fion2=0.0;
		      surrogate[i].gamma_LR=surrogate[i].LR(b);
		      double Kaq=surrogate[i].Kp_eff_aq(config, Temperature, ionic(b), chp2,surrogate[config.iHp].LR(b),
							surrogate[config.iHp].SRMR(b),MMaq(b), fion1, fion2)/kelvin_effect/surrogate[i].gamma_aq_bins(b);
		      
		      double ratio_gamma1=pow(surrogate[config.iHp].LR(b),2.0)*surrogate[config.iHp].SRMR(b)/surrogate[i].gamma_LR;
		      double ratio_gamma2=pow(surrogate[config.iHp].LR(b),2.0)*surrogate[config.iHp].SRMR(b);
		      double Kac1=surrogate[i].Kacidity1/ratio_gamma1;
		      double Kac2=surrogate[i].Kacidity2/ratio_gamma2;
		      double Kh=surrogate[i].Kpart_aq(Temperature, MMaq(b))/kelvin_effect/surrogate[i].gamma_aq_bins(b);

		      surrogate[i].k1_aq(b,index)=faq*
			(surrogate[i].Agt*Kaq*AQinit(b)-surrogate[i].Aaq_bins_init(b)-surrogate[i].Ag+surrogate[i].Agt)/
			(Kaq*AQinit(b)*surrogate[i].tau_air(b));

		      double dK=-Kac1*Kh/pow(chp2,2.0)-2.0*Kac1*Kac2*Kh/pow(chp2,3.0);
		      double dJ=surrogate[i].Aaq_bins_init(b)/pow(Kaq*AQinit(b),2)/surrogate[i].tau_air(b)*faq*dK*AQinit(b)*Jhp;
		      surrogate[i].Aaq+=surrogate[i].k1_aq(b,index)*dt1+0.5*dJ*dt1*dt1;
		      surrogate[i].Agt=max(surrogate[i].Agt-surrogate[i].k1_aq(b,index)*dt1-0.5*dJ*dt1*dt1,0.01*surrogate[i].Agt);
		      
		    }
		}
	    
	    chp2=chp2+Jhp*dt1;
	    t+=dt1;
	  }
      
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
        else if (surrogate[i].name=="H2SO4")
          {	    
          }
	else if (surrogate[i].is_organic)
	  if ((surrogate[i].aq_type=="monoacid" or surrogate[i].aq_type=="diacid") and surrogate[i].time_aq(b)>=config.tequilibrium)
	    {
	      if (surrogate[i].time_aq(b)>=config.tequilibrium)
		surrogate[i].k1_aq(b,index)=surrogate[i].Aaq/deltat;
	      surrogate[i].Jdn_aq(b,index)=0.0;
	      if (surrogate[i].k1_aq(b,index)<0.0 and surrogate[i].Aaq_bins_init(b)>tiny)
		surrogate[i].Jdn_aq(b,index)=surrogate[i].k1_aq(b,index)/surrogate[i].Aaq_bins_init(b);
	    }
    }
}

void compute_morphology(model_config &config, Array<double,1> &vsol, Array<double, 1> & Number)
{
  int b,ilayer;
  //compute morpholgy
  for (b=0;b<config.nbins;++b)
   if (Number(b) > 0.0)  
    {
      config.dbound(b,0)=pow(3.0/(4.0*3.14159)*vsol(b),1.0/3.0)/Number(b);
      double vorg=(4.0/3.0*3.14159*pow(0.5e-6*config.diameters(b),3.0)-vsol(b)/Number(b));  
      for (ilayer=0;ilayer<config.nlayer;++ilayer)
	{
	  config.dbound(b,ilayer+1)=pow(pow(config.dbound(b,ilayer),3)+3.0/(4.0*3.14159)*vorg*config.Vlayer(ilayer),1.0/3.0);
	  config.Radius(b,ilayer)=0.5*(config.dbound(b,ilayer+1)+config.dbound(b,ilayer));	  
	}
    }
}

void flux_org(model_config &config, vector<species>& surrogate,
              Array<double, 3> &MOinit, Array<double, 1> &AQinit,
              double &tiny, int index)
{
  //compute kinetic rates for the absorption of a compound in the organic phase in a bin and a layer
  //index = 0 : first evaluation of rates)
  //index = 1 : second evaluation of rates
  double Kp;
  int n=surrogate.size();
  int i,b,ilayer,iphase,jphase;
  double sum,sum_mass;
  
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
			    (0.0
			     +sum/(1.0-AQinit(b)/sum_mass)*surrogate[i].tau_air(b));		
		      else
			  surrogate[i].k1(b,ilayer,iphase,index)=
			    (surrogate[i].Ag*surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase)
			     -surrogate[i].Ap_layer_init(b,ilayer,iphase))/
			    (0.0
			     +sum/(1.0-AQinit(b)/sum_mass)*surrogate[i].tau_air(b));					      		     
	
		      surrogate[i].Jdn(b,ilayer,iphase,index)=0.0;	  
		      if (surrogate[i].k1(b,ilayer,iphase,index)<0.0 and
			  surrogate[i].Ap_layer_init(b,ilayer,iphase)>tiny)
			surrogate[i].Jdn(b,ilayer,iphase,index)=surrogate[i].k1(b,ilayer,iphase,index)/
			  surrogate[i].Ap_layer_init(b,ilayer,iphase);
		    }

	    
	      for (ilayer=0;ilayer<config.nlayer;++ilayer)		    
		{
		  double F1=0.0;
		  double F2=0.0;
		  sum=0.0;
		  for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		    sum+=surrogate[i].Ap_layer_init(b,ilayer,iphase)/(surrogate[i].Kp(b,ilayer,iphase)*MOinit(b,ilayer,iphase));

		  if (ilayer>0)
		    {
		      double surf1=4.0*3.14159*pow(config.dbound(b,ilayer),2);
		      double sum1=0.0;
		      for (iphase=0;iphase<config.nphase(b,ilayer-1);++iphase)
			sum1+=surrogate[i].Ap_layer_init(b,ilayer-1,iphase)/(surrogate[i].Kp(b,ilayer-1,iphase)*MOinit(b,ilayer-1,iphase));

		      double dorg=(config.dbound(b,ilayer)-config.Radius(b,ilayer-1))/(config.Radius(b,ilayer)-config.Radius(b,ilayer-1))*surrogate[i].dif_org(b,ilayer-1)
			+(config.Radius(b,ilayer)-config.dbound(b,ilayer))/(config.Radius(b,ilayer)-config.Radius(b,ilayer-1))*surrogate[i].dif_org(b,ilayer);
			
		      F1=surf1*dorg*(sum1-sum)/(config.Radius(b,ilayer)-config.Radius(b,ilayer-1));
		    }

		  if (ilayer<config.nlayer-1)
		    {
		      double surf2=4.0*3.14159*pow(config.dbound(b,ilayer+1),2);
		      double sum2=0.0;
		      for (iphase=0;iphase<config.nphase(b,ilayer+1);++iphase)
			sum2+=surrogate[i].Ap_layer_init(b,ilayer+1,iphase)/(surrogate[i].Kp(b,ilayer+1,iphase)*MOinit(b,ilayer+1,iphase));

		      double dorg=(config.dbound(b,ilayer+1)-config.Radius(b,ilayer))/(config.Radius(b,ilayer+1)-config.Radius(b,ilayer))*surrogate[i].dif_org(b,ilayer)
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

		  bool eq_last_layer=true;
		  double sumkpositive=0.0;
		  double sumknegative=0.0;
		  double tau_interface=max(1.0,surrogate[i].tau_diffusion(b,0,0));
		  Array<double, 1> kpmo_interface;
		  double ap_interface=0.0;		  
		  double Vinterface=0.0;
		  double sumkpmo_interface=0.0;
		  int ilayer_interface=0;		  
		  kpmo_interface.resize(config.nlayer);
		  kpmo_interface=0.0;
		  Array <double,2> rJ_interface;
		  rJ_interface.resize(config.nlayer,config.max_number_of_phases);
		  rJ_interface=0.0;
		  for (ilayer=0;ilayer<config.nlayer-1;++ilayer)
		    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		      if (surrogate[i].tau_diffusion(b,ilayer,iphase)>config.tequilibrium and surrogate[i].tau_diffusion(b,ilayer+1,iphase)<config.tequilibrium)
			{
			  tau_interface=surrogate[i].tau_diffusion(b,ilayer,iphase);
			  ilayer_interface=ilayer+1;
			}

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

		  double Jinterface=(surrogate[i].Ag*sumkpmo_interface-ap_interface)/tau_interface;
	      
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
	      
		  double ktot1=0.0;
		  for (ilayer=0;ilayer<ilayer_interface;++ilayer)
		    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		      if (surrogate[i].time(b,ilayer,iphase)>=config.tequilibrium)		
			ktot1+=surrogate[i].k1(b,ilayer,iphase,index);
		     		    
		  for (ilayer=0;ilayer<config.nlayer;++ilayer)
		    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		      if (surrogate[i].time(b,ilayer,iphase)>=config.tequilibrium)			
			if (ilayer<ilayer_interface or
			    (Jinterface*ktot1<=0.0 and sumkpmo_interface*surrogate[i].tau_air(b)<surrogate[i].tau_diffusion(b,ilayer_interface-1,iphase)))
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
      	      
		  double kcond=0.0;	    	      
		  double ktot=sumkpositive+sumknegative;
	      
		  if (ktot>0.0)			
		    {
		      double sum1=0.0;
		      double sum2=0.0;
		      for (ilayer=0;ilayer<config.nlayer;++ilayer)		   
			if (surrogate[i].time(b,ilayer,0)>=config.tequilibrium)
			  {					  
			    double sumk=0.0;
			    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
			      sumk+=surrogate[i].k1(b,ilayer,iphase,index);
			
			    double sumap=0.0;
			    double kpmo=0.0;
	
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
		      double sum1=0.0;
		      double sum2=0.0;
		      for (ilayer=0;ilayer<config.nlayer;++ilayer)
			for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
			  if (surrogate[i].time(b,ilayer,0)>=config.tequilibrium)
			    {
			      double sumk=0.0;
			      for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
				sumk+=surrogate[i].k1(b,ilayer,iphase,index);
			  			 
			      double sumap=0.0;
			      double kpmo=0.0;
			  
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
		  double b2=0.0;
		  double c2=0.0;
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
		      int jlayer;
		      double sumk=0.0;
		      for (jlayer=ilayer_interface;jlayer<config.nlayer;jlayer++)
			for (jphase=0;jphase<config.nphase(b,jlayer);++jphase)
			  sumk+=surrogate[i].k1(b,jlayer,jphase,index);

		      for (ilayer=0;ilayer<config.nlayer;++ilayer)
			for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)		
			  if (surrogate[i].time(b,ilayer,iphase)>=config.tequilibrium)		  	  
			    if (surrogate[i].k1(b,ilayer,iphase,index)>0.0 and ktot>0.0 and kcond>0.0)
			      {			   			 
				double a;							       				
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
				      a=1.0/(1.0/(surrogate[i].k1(b,ilayer,iphase,index)/(sumkpositive-max(sumk,0.0))*(kcond-max(sumk,0.0)))+1.0/surrogate[i].k1(b,ilayer,iphase,index))-sumknegative*surrogate[i].k1(b,ilayer,iphase,index)/sumkpositive;
				      b2+=a;
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
		  sum_mass=max(sum_mass,config.MOmin);
	   
		  //compute kinetic rate of absorption	    
                  ilayer = 0;
                  for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                    if(surrogate[i].time(b,ilayer,iphase)>=config.tequilibrium)
                      {		      
                        sum=0.0;
                        for (jphase=0;jphase<config.nphase(b,ilayer);++jphase)
                          sum+=surrogate[i].Kp(b,ilayer,jphase)*MOinit(b,ilayer,jphase);
		      		      
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
	double sumkpositive=0.0;
	double sumknegative=0.0;
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
		
		//compute fluxes
		surrogate[i].Jdn(b,ilayer,iphase,index)=0.0; 
		if (surrogate[i].k1(b,ilayer,iphase,index)<0.0 and
		    surrogate[i].Ap_layer_init(b,ilayer,iphase)>tiny)
		  surrogate[i].Jdn(b,ilayer,iphase,index)=surrogate[i].k1(b,ilayer,iphase,index)/
		    surrogate[i].Ap_layer_init(b,ilayer,iphase);
	      }
      }
}

void error_ph_bins(model_config &config, vector<species> &surrogate, int index_b, double Temperature, Array <double, 1> chp, 
                   double &error, double &derivative, Array <double, 1> AQinit, Array<double , 1> &ionic, Array<double, 1> &MMaq,
		   Array<double, 1> &LWC,
		   double &chp_new)
{      
  int n=surrogate.size();
  int i;
  double inorganion=0.0;
  double organion=0.0;
  double total;
  double Ke=1.0e-14;
  double sum_K_AQ;
  int b;
  double kelvin_effect=1.0;
  if (config.compute_kelvin_effect) //compute the kelvin effect    
    kelvin_effect=exp(2.0*config.surface_tension_aq*1.0e-6*MMaq(index_b)/
		      (8.314*Temperature*config.AQrho(index_b)*
		       0.5*config.diameters(index_b)*1.0e-6));

  double conc_org=LWC(index_b);
  for (i=0;i<n;++i)
    if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophilic)
      conc_org+=surrogate[i].Aaq_bins_init(index_b);
  
  derivative=0.0;
  for (i=0;i<n;++i)
    if (surrogate[i].is_inorganic_precursor==true)
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
                      sum_K_AQ+=surrogate[i].Kaq(b)*AQinit(b);
                    }
                derivative-=1000.0*total*AQinit(index_b)/conc_org*
                  (surrogate[i].Kaq(index_b)/chp(index_b)/sum_K_AQ
                   -surrogate[i].Kaq(index_b)/pow(sum_K_AQ,2.0)*surrogate[i].Kaq(index_b)/chp(index_b)*AQinit(index_b));
                inorganion-=1000.0*total*surrogate[i].Kaq(index_b)*AQinit(index_b)/sum_K_AQ/conc_org; 
              }
            else
              inorganion-=surrogate[config.iNH4p].Aaq_bins_init(index_b)/surrogate[config.iNH4p].MM/conc_org*1000.0*surrogate[config.iNH4p].charge;
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
                      sum_K_AQ+=surrogate[i].Kaq(b)*AQinit(b);
                    }
                
                derivative+=1000.0/conc_org*AQinit(index_b)*total*
                  (-surrogate[i].Kaq(index_b)/chp(index_b)/sum_K_AQ
                   +surrogate[i].Kaq(index_b)/pow(sum_K_AQ,2.0)*surrogate[i].Kaq(index_b)/chp(index_b)*AQinit(index_b));
                inorganion+=1000.*total*surrogate[i].Kaq(index_b)*AQinit(index_b)/sum_K_AQ/conc_org;
              }
            else
              inorganion-=surrogate[config.iNO3m].Aaq_bins_init(index_b)/surrogate[config.iNO3m].MM/conc_org*1000.*surrogate[config.iNO3m].charge;
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
                      sum_K_AQ+=surrogate[i].Kaq(b)*AQinit(b);
                    }
                
                derivative+=1000.0*AQinit(index_b)/conc_org*total*
                  (-surrogate[i].Kaq(index_b)/chp(index_b)/(sum_K_AQ)
                   +surrogate[i].Kaq(index_b)/pow(sum_K_AQ,2.0)*surrogate[i].Kaq(index_b)/chp(index_b)*AQinit(index_b));
                inorganion+=config.AQrho(index_b)*total*surrogate[i].Kaq(index_b)/sum_K_AQ;
              }            
            else
              {
                inorganion-=surrogate[config.iClm].Aaq_bins_init(index_b)/surrogate[config.iClm].MM/conc_org*1000.*surrogate[config.iClm].charge;
              }
          }
        else if (surrogate[i].name=="H2SO4")
          {
            total=surrogate[config.iHSO4m].Aaq_bins_init(index_b)/surrogate[config.iHSO4m].MM
              +surrogate[config.iSO4mm].Aaq_bins_init(index_b)/surrogate[config.iSO4mm].MM;

	    double K=surrogate[i].Kequilibrium(Temperature)*surrogate[config.iHSO4m].gamma_aq_bins(index_b)
              /(surrogate[config.iHp].gamma_aq_bins(index_b)*surrogate[config.iSO4mm].gamma_aq_bins(index_b));
	  
            derivative-=1000.*total/conc_org*1.0/pow(1.0+K/chp(index_b),2.0)*K/(chp(index_b)*chp(index_b)); //HSO4m+SO4mm
            inorganion+=1000.*total/conc_org*(2.0-1.0/(1.0+K/chp(index_b)));
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

		double fion1=0.0;
		double fion2=0.0;
		surrogate[i].gamma_LR=surrogate[i].LR(index_b);
		double Kaq=surrogate[i].Kp_eff_aq(config, Temperature, ionic(index_b), chp(index_b),surrogate[config.iHp].LR(index_b),
						  surrogate[config.iHp].SRMR(index_b),MMaq(index_b), fion1, fion2)/kelvin_effect/surrogate[i].gamma_aq_bins(index_b);

		double ratio_gamma=pow(surrogate[config.iHp].LR(index_b),2.0)*surrogate[config.iHp].SRMR(index_b)/surrogate[i].gamma_LR;
		double Kac=surrogate[i].Kacidity1/ratio_gamma;
		double Kh=surrogate[i].Kpart_aq(Temperature, MMaq(index_b))/kelvin_effect/surrogate[i].gamma_aq_bins(index_b);
		double dK=-Kac*Kh/pow(chp(index_b),2.0);
		
                derivative+=1000.*AQinit(index_b)/conc_org*total*fion1*
		  (-(1.0-fion1)/chp(index_b)*surrogate[i].Kaq(index_b)/sum_K_AQ+dK/sum_K_AQ
		   -Kaq/pow(sum_K_AQ,2.0)*AQinit(index_b)*dK);                 
		organion+=1000.*fion1*total*surrogate[i].Kaq(index_b)*AQinit(index_b)/sum_K_AQ/conc_org;		
	      }
	    else
	      {
		total=surrogate[i].Aaq_bins_init(index_b)/surrogate[i].MM;
		double fion1=0.0;
		double fion2=0.0;
		surrogate[i].gamma_LR=surrogate[i].LR(index_b);
		double Kaq=surrogate[i].Kp_eff_aq(config, Temperature, ionic(index_b), chp(index_b),surrogate[config.iHp].LR(index_b),
						  surrogate[config.iHp].SRMR(index_b),MMaq(index_b), fion1, fion2)/kelvin_effect/surrogate[i].gamma_aq_bins(index_b);
		organion+=1000.*fion1*total/conc_org;
		derivative-=fion1*(1-fion1)/chp(index_b)*1000.*total/conc_org;
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

		double fion1=0.0;
		double fion2=0.0;
		surrogate[i].gamma_LR=surrogate[i].LR(index_b);
		double Kaq=surrogate[i].Kp_eff_aq(config, Temperature, ionic(index_b), chp(index_b),surrogate[config.iHp].LR(index_b),
						  surrogate[config.iHp].SRMR(index_b),MMaq(index_b), fion1, fion2)/kelvin_effect/surrogate[i].gamma_aq_bins(index_b);

		double ratio_gamma1=pow(surrogate[config.iHp].LR(index_b),2.0)*surrogate[config.iHp].SRMR(index_b)/surrogate[i].gamma_LR;
		double ratio_gamma2=pow(surrogate[config.iHp].LR(index_b),2.0)*surrogate[config.iHp].SRMR(index_b);
		double Kac1=surrogate[i].Kacidity1/ratio_gamma1;
		double Kac2=surrogate[i].Kacidity2/ratio_gamma2;
		double Kh=surrogate[i].Kpart_aq(Temperature, MMaq(index_b))/kelvin_effect/surrogate[i].gamma_aq_bins(index_b);
		double dK=-Kac1*Kh/pow(chp(index_b),2.0)-2.0*Kac1*Kac2*Kh/pow(chp(index_b),3.0);
		
		derivative+=1000.*total*AQinit(index_b)/conc_org*
		  ((fion1*(fion1-1.0)+4.0*fion1*fion2+4.0*fion2*(fion2-1))/chp(index_b)*surrogate[i].Kaq(index_b)/sum_K_AQ
		   +dK/sum_K_AQ-Kaq/pow(sum_K_AQ,2.0)*AQinit(index_b)*dK);                 
		organion+=1000.0*(fion1+2.0*fion2)*total*surrogate[i].Kaq(index_b)*AQinit(index_b)/sum_K_AQ/conc_org;		
	      }
	    else
	      {
		total=surrogate[i].Aaq_bins_init(index_b)/surrogate[i].MM;
		double fion1=0.0;
		double fion2=0.0;
		surrogate[i].gamma_LR=surrogate[i].LR(index_b);
		double Kaq=surrogate[i].Kp_eff_aq(config, Temperature, ionic(index_b), chp(index_b),surrogate[config.iHp].LR(index_b),
						  surrogate[config.iHp].SRMR(index_b),MMaq(index_b), fion1, fion2)/kelvin_effect/surrogate[i].gamma_aq_bins(index_b);
		organion+=1000.*(fion1+2.0*fion2)*total/conc_org;
		derivative+=(fion1*(fion1-1)+4.0*fion1*fion2+4.0*fion2*(fion2-1))/chp(index_b)*1000.*total/conc_org;
	      }
	  }

      }
    else 
      if (surrogate[i].name=="Na")
        inorganion-=surrogate[i].Aaq_bins_init(index_b)/surrogate[i].MM/conc_org*1000.*surrogate[i].charge;

  derivative=derivative*0.5*(1.0+(organion+inorganion)/pow(pow(organion+inorganion,2)+4*Ke,0.5))-1.0;
  chp_new=0.5*(organion+inorganion+pow(pow(organion+inorganion,2)+4*Ke,0.5));
  error=chp_new-chp(index_b);
}

void error_ph_dyn(model_config &config, vector<species> &surrogate, int index_b, double Temperature, Array <double, 1> chp, 
		  double &error, double &derivative, Array <double, 1> AQinit, Array<double , 1> &ionic, Array<double, 1> &MMaq,
		  Array<double, 1> &LWC,
		  double &chp_new)
{      
  int n=surrogate.size();
  int i;
  double inorganion=0.0;
  double organion=0.0;
  double total;
  double Ke=1.0e-14;
  double sum_K_AQ;
  int b;
  double kelvin_effect=1.0;
  if (config.compute_kelvin_effect) //compute the kelvin effect
    {
      kelvin_effect=exp(2.0*config.surface_tension_aq*1.0e-6*MMaq(index_b)/
			(8.314*Temperature*config.AQrho(index_b)*
			 0.5*config.diameters(index_b)*1.0e-6));
    }

  double conc_org=LWC(index_b);
  for (i=0;i<n;++i)
    if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophilic)
      conc_org+=surrogate[i].Aaq_bins_init(index_b);

  derivative=0.0;
  for (i=0;i<n;++i)
    if (surrogate[i].is_inorganic_precursor==true)
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

	    double K=surrogate[i].Kequilibrium(Temperature)*surrogate[config.iHSO4m].gamma_aq_bins(index_b)
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
	    double fion1=0.0;
	    double fion2=0.0;
	    surrogate[i].gamma_LR=surrogate[i].LR(index_b);
	    double Kaq=surrogate[i].Kp_eff_aq(config,Temperature, ionic(index_b), chp(index_b),surrogate[config.iHp].LR(index_b),
					      surrogate[config.iHp].SRMR(index_b),MMaq(index_b), fion1, fion2)/kelvin_effect/surrogate[i].gamma_aq_bins(index_b);
	    organion+=1000.*fion1*total/conc_org;
	    derivative-=fion1*(1-fion1)/chp(index_b)*1000.*total/conc_org;	   
	  }
	else if (surrogate[i].aq_type=="diacid")
	  {
	    total=surrogate[i].Aaq_bins_init(index_b)/surrogate[i].MM;
	    double fion1=0.0;
	    double fion2=0.0;
	    surrogate[i].gamma_LR=surrogate[i].LR(index_b);	    
	    double Kaq=surrogate[i].Kp_eff_aq(config, Temperature, ionic(index_b), chp(index_b),surrogate[config.iHp].LR(index_b),
					      surrogate[config.iHp].SRMR(index_b),MMaq(index_b), fion1, fion2)/kelvin_effect/surrogate[i].gamma_aq_bins(index_b);
	    organion+=1000.*(fion1+2.0*fion2)*total/conc_org;
	    derivative+=(fion1*(fion1-1.)+4.0*fion1*fion2+4.0*fion2*(fion2-1.))/chp(index_b)*1000.*total/conc_org;	   
	  }

      }
    else 
      if (surrogate[i].name=="Na")
	{
	  inorganion-=surrogate[i].Aaq_bins_init(index_b)/surrogate[i].MM/conc_org*1000.*surrogate[i].charge;
	}

  derivative=derivative*0.5*(1.0+(organion+inorganion)/pow(pow(organion+inorganion,2)+4*Ke,0.5))-1.0;
  chp_new=0.5*(organion+inorganion+pow(pow(organion+inorganion,2)+4*Ke,0.5));
  error=chp_new-chp(index_b);
}

void compute_ph_dyn(model_config &config, vector<species> &surrogate, double Temperature, Array <double, 1> &chp, 
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
  for (b=0;b<config.nbins;b++)
    for (i=0;i<n;i++)
      if (surrogate[i].is_organic or i==config.iH2O)
	conc_org(b)+=surrogate[i].Aaq_bins_init(b);

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
              Array <double, 1> chp2;
              chp2.resize(config.nbins);
              int b2;
              int index=0;
              for (b2=0;b2<config.nbins;b2++)
                chp2(b2)=chp(b2);
              double chp_new;
              
	      while(abs(error_h/chp2(b))>1.0e-3 and index<1000)
		{
		  index++;
		  compute_kp_aq(config, surrogate, Temperature, ionic, chp2, MMaq);
		  error_ph_dyn(config, surrogate, b, Temperature, chp2, error_h, derivative_h,AQinit,ionic,MMaq,LWC,chp_new);
		  if (chp2(b)==1.0e-14 and chp_new<=0.0)
		    error_h=0.0;
		  else if (chp_new <= 0.0 or chp_new!=chp_new)
		    chp2(b)=1.0e-14;
		  else
		    if (chp2(b)-error_h/derivative_h>0.0 and derivative_h!=0.0) 
		      chp2(b)=chp2(b)-error_h/derivative_h;
		    else
		      chp2(b)=chp2(b)+error_h;		  
		}	
	      error_tot=max(error_tot,abs(error_h)/chp2(b));
	      chp(b)=factor*chp2(b)+(1.0-factor)*chp(b);	      
            }	  
          surrogate[config.iHp].Aaq_bins_init(b)=chp(b)*conc_org(b)/1000.0;
        }
      nh++;
    }

  for (i=0;i<n;++i)
    if (surrogate[i].is_inorganic_precursor==true)
      {
        if (surrogate[i].name=="H2SO4")
	  for (b=0;b<config.nbins;++b)
	    {
	      total=surrogate[config.iHSO4m].Aaq_bins_init(b)/surrogate[config.iHSO4m].MM
		+surrogate[config.iSO4mm].Aaq_bins_init(b)/surrogate[config.iSO4mm].MM;
	      
	      double K=surrogate[i].Kequilibrium(Temperature)*surrogate[config.iHSO4m].gamma_aq_bins(b)
		/(surrogate[config.iHp].gamma_aq_bins(b)*surrogate[config.iSO4mm].gamma_aq_bins(b));
	      
	      surrogate[config.iHSO4m].Aaq_bins_init(b)=total*surrogate[config.iHSO4m].MM*1.0/(1.0+K/chp(b));
	      surrogate[config.iSO4mm].Aaq_bins_init(b)=total*surrogate[config.iSO4mm].MM*K/chp(b)/(1.0+K/chp(b));
	    }

      }
}	

void error_ph_dyn2(model_config &config, vector<species> &surrogate, int index_b, double Temperature, Array <double, 1> chp, 
		   double &error, double &derivative, Array <double, 1> AQ, Array<double , 1> &ionic, Array<double, 1> &MMaq,
		   Array<double, 1> &LWC, double &chp_new)
{      
  int n=surrogate.size();
  int i;
  double inorganion=0.0;
  double organion=0.0;
  double total;
  double Ke=1.0e-14;
  double sum_K_AQ;
  int b;
  double kelvin_effect=1.0;
  

  if (config.compute_kelvin_effect) //compute the kelvin effect
    {
      kelvin_effect=exp(2.0*config.surface_tension_aq*1.0e-6*MMaq(index_b)/
			(8.314*Temperature*config.AQrho(index_b)*
			 0.5*config.diameters(index_b)*1.0e-6));
    }

  double conc_org=LWC(index_b);
  for (i=0;i<n;++i)
    if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophilic)
      conc_org+=surrogate[i].Aaq_bins(index_b);

  derivative=0.0;
  for (i=0;i<n;++i)
    if (surrogate[i].is_inorganic_precursor==true)
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
	    double K=surrogate[i].Kequilibrium(Temperature)*surrogate[config.iHSO4m].gamma_aq_bins(index_b)
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
	    double Kaq=surrogate[i].Kp_eff_aq(config, Temperature, ionic(index_b), chp(index_b),surrogate[config.iHp].LR(index_b),
					      surrogate[config.iHp].SRMR(index_b),MMaq(index_b), fion1, fion2)/kelvin_effect/surrogate[i].gamma_aq_bins(index_b);
	    organion+=1000.*fion1*total/conc_org;
	    derivative-=fion1*(1-fion1)/chp(index_b)*1000.0*total/conc_org;	   
	  }
	else if (surrogate[i].aq_type=="diacid")
	  {
	    total=surrogate[i].Aaq_bins(index_b)/surrogate[i].MM;
	    double fion1=0.0;
	    double fion2=0.0;
	    surrogate[i].gamma_LR=surrogate[i].LR(index_b);
	    double Kaq=surrogate[i].Kp_eff_aq(config, Temperature, ionic(index_b), chp(index_b),surrogate[config.iHp].LR(index_b),
					      surrogate[config.iHp].SRMR(index_b),MMaq(index_b), fion1, fion2)/kelvin_effect/surrogate[i].gamma_aq_bins(index_b);
	    organion+=1000.0*(fion1+2.0*fion2)*total/conc_org;
	    derivative+=(fion1*(fion1-1)+4.0*fion1*fion2+4.0*fion2*(fion2-1))/chp(index_b)*1000.*total/conc_org;	   
	  }

      }
    else 
      if (surrogate[i].name=="Na")
        inorganion-=surrogate[i].Aaq_bins(index_b)/surrogate[i].MM/conc_org*1000.*surrogate[i].charge;


  derivative=derivative*0.5*(1.0+(organion+inorganion)/pow(pow(organion+inorganion,2)+4*Ke,0.5))-1.0;
  chp_new=0.5*(organion+inorganion+pow(pow(organion+inorganion,2)+4*Ke,0.5));
  error=chp_new-chp(index_b);
}
   
void compute_ph_dyn2(model_config &config, vector<species> &surrogate, double Temperature, Array <double, 1> &chp, 
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
  for (b=0;b<config.nbins;b++)
    for (i=0;i<n;i++)
      if (surrogate[i].is_organic or i==config.iH2O)
	conc_org(b)+=surrogate[i].Aaq_bins(b);

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
              Array <double, 1> chp2;
              chp2.resize(config.nbins);
              int b2;
              int index=0;
              for (b2=0;b2<config.nbins;b2++)
                chp2(b2)=chp(b2);
              double chp_new;
              
	      while(abs(error_h/chp2(b))>1.0e-3 and index<1000)
		{
		  index++;
		  compute_kp_aq(config, surrogate, Temperature, ionic, chp2, MMaq);
		  error_ph_dyn2(config, surrogate, b, Temperature, chp2, error_h, derivative_h,AQ,ionic,MMaq,LWC,chp_new);
		  if (chp2(b)==1.0e-14 and chp_new<=0.0)
		    error_h=0.0;
		  else if (chp_new <= 0.0 or chp_new!=chp_new)
		    chp2(b)=1.0e-14;
		  else
		    if (chp2(b)-error_h/derivative_h>0.0 and derivative_h!=0.0) 
		      chp2(b)=chp2(b)-error_h/derivative_h;
		    else
		      chp2(b)=chp2(b)+error_h;		  
		}	
	      error_tot=max(error_tot,abs(error_h)/chp2(b));
	      chp(b)=factor*chp2(b)+(1.0-factor)*chp(b);	      
            }	  
          surrogate[config.iHp].Aaq_bins(b)=chp(b)*conc_org(b)/1000.;
        }
      nh++;
    }

  for (i=0;i<n;++i)
    if (surrogate[i].is_inorganic_precursor==true)
      {
        if (surrogate[i].name=="H2SO4")
	  for (b=0;b<config.nbins;++b)
	    {
	      total=surrogate[config.iHSO4m].Aaq_bins(b)/surrogate[config.iHSO4m].MM
		+surrogate[config.iSO4mm].Aaq_bins(b)/surrogate[config.iSO4mm].MM;
	      
	      double K=surrogate[i].Kequilibrium(Temperature)*surrogate[config.iHSO4m].gamma_aq_bins(b)
		/(surrogate[config.iHp].gamma_aq_bins(b)*surrogate[config.iSO4mm].gamma_aq_bins(b));
	      
	      surrogate[config.iHSO4m].Aaq_bins(b)=total*surrogate[config.iHSO4m].MM*1.0/(1.0+K/chp(b));
	      surrogate[config.iSO4mm].Aaq_bins(b)=total*surrogate[config.iSO4mm].MM*K/chp(b)/(1.0+K/chp(b));
	    }
      }

}

void activity_coefficients_dyn_aq(model_config &config, vector<species>& surrogate,
                                  double &Temperature, Array<double, 1> &AQinit,
                                  Array<double, 3> &MOinit,
                                  Array<double, 1> &conc_inorganic,
                                  Array<double, 1> &ionic, Array<double, 1> &ionic_organic,
                                  Array<double, 1> &organion, Array<double, 1> &chp,
                                  Array<double, 1> &LWC,Array<double, 1> &MMaq, double factor)
{
  //compute the activity coefficients with UNIFAC (short range interactions) for the organic phase
  //MOW: mean molar mass of the organic phase
  int n_unifac,i;
  int n=surrogate.size();
  double sum,sumX_unifac,XH2O;
  Array<double, 1> X_unifac,gamma_unifac;
  int b;
  
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

      compute_ionic_strenght(config, surrogate,AQ, sum_inorganic, ionic_tmp, chp_tmp,
                             organion_tmp, ionic_organic_tmp, factor);

      activity_coefficients_aq(config,surrogate,Temperature,0.0,MMaqtemp,XH2O);
      if (config.compute_long_and_medium_range_interactions)
        activity_coefficients_LR_MR(config, surrogate, Temperature, 0.0, ionic_tmp);
						  
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
      double AQ=0.0;
      Array <double, 1> chp_save;
      chp_save.resize(config.nbins);
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
	  	 
          compute_ionic_strenght2(config, surrogate, AQ, conc_inorganic(b), ionic(b), chp(b),
				  organion(b), ionic_organic(b), factor);
	 	  
          activity_coefficients_aq(config,surrogate,Temperature,0.0,MMaq(b),XH2O);	 
          if (config.compute_long_and_medium_range_interactions)
            activity_coefficients_LR_MR(config, surrogate, Temperature, 0.0, ionic(b));
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
           
      for (b=0;b<config.nbins;++b)
        {
          if (config.compute_inorganic)
            {
              double error_h=1000.0;
              double derivative_h;
              Array <double, 1> chp2;
              chp2.resize(config.nbins);
              int b2;
              int index=0;
              for (b2=0;b2<config.nbins;b2++)
                chp2(b2)=chp_save(b2);
              double chp_new;
              
	      while(abs(error_h/chp2(b))>1.0e-3 and index<1000)
		{
		  index++;
		  compute_kp_aq(config, surrogate, Temperature, ionic, chp2, MMaq);
		  error_ph_bins(config, surrogate, b, Temperature, chp2, error_h, derivative_h,AQinit,ionic,MMaq,LWC,chp_new);
		  if (chp2(b)==1.0e-14 and chp_new<=0.0)
		    error_h=0.0;
		  else if (chp_new <= 0.0 or chp_new!=chp_new)
		    chp2(b)=1.0e-14;
		  else
		    if (chp2(b)-error_h/derivative_h>0.0 and derivative_h!=0.0) 
		      chp2(b)=chp2(b)-error_h/derivative_h;
		    else
		      chp2(b)=chp2(b)+error_h;
		  
		}
	      chp(b)=factor*chp2(b)+(1.0-factor)*chp(b);	      
            }	  
	  
	  double conc_org;	  
	  conc_org=LWC(b);	  
	  for (i=0;i<n;i++)
	    if (surrogate[i].is_organic or i==config.iH2O)
	      conc_org+=surrogate[i].Aaq_bins_init(b);

          surrogate[config.iHp].Aaq_bins_init(b)=chp(b)*conc_org/1000.0;
        }      
      compute_kp_aq(config, surrogate, Temperature, ionic, chp, MMaq);
      characteristic_time_aq(config, surrogate, LWC, AQinit, MOinit);
    } 

  compute_organion(config, surrogate, Temperature,  MMaq, AQinit, ionic,
                   chp, ionic_organic, organion, LWC);
}

void equilibrium_inorg(model_config &config, vector<species>& surrogate, double &tequilibrium,
                       Array <double, 1> &AQinit, Array <double, 1> &conc_inorganic, 
                       Array <double, 1> &chp, Array <double, 1> &LWC,
                       double &Temperature, double &RH, Array <double, 1> &MMaq, double factor)
{
  int n=surrogate.size();
  int i,b;
  double kelvin_effect=1.0;
  double conc_equilibrium;
  double sum=0.0;
  bool is_equilibrium=false;

  //This routine computes the aqueous-phase concentrations of inorganic compounds being at equilibrium
  //The inorganic aerosol is assumed metastable
  //compute_activity_coefficients: should the activity coefficients be recomputed?
  //tequilibrium: criterium to distinguish cases under equilibrium (characteristic time lower
  //                  than tequilibrium)

  for (i=0;i<n;++i)
    if(surrogate[i].is_inorganic_precursor and surrogate[i].hydrophilic)
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
              sum+=surrogate[i].Kaq(b)*AQinit(b);
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
                            surrogate[i].Kaq(b)*AQinit(b)/sum+(1.0-factor)*surrogate[config.iNH4p].Aaq_bins_init(b);
                        else
                          {
                            surrogate[config.iNH4p].Aaq_bins_init(b)=0.01*factor*surrogate[config.iNH4p].MM/surrogate[i].MM*conc_equilibrium*surrogate[i].Kaq(b)*AQinit(b)/sum;
                            surrogate[i].Ag+=99*surrogate[config.iNH4p].Aaq_bins_init(b);
                          }                     
                      }
                    else if (surrogate[i].name=="HNO3")
                      {
                        if (surrogate[config.iNO3m].Aaq_bins_init(b)>0.0)
                          surrogate[config.iNO3m].Aaq_bins_init(b)=factor*surrogate[config.iNO3m].MM/surrogate[i].MM*conc_equilibrium*
                            surrogate[i].Kaq(b)*AQinit(b)/sum+(1.0-factor)*surrogate[config.iNO3m].Aaq_bins_init(b);
                        else
                          {
                            surrogate[config.iNO3m].Aaq_bins_init(b)=0.01*factor*surrogate[config.iNO3m].MM/surrogate[i].MM*conc_equilibrium*
                            surrogate[i].Kaq(b)*AQinit(b)/sum;
                            surrogate[i].Ag+=99*surrogate[config.iNO3m].Aaq_bins_init(b);
                          }
                      }
                    else if (surrogate[i].name=="HCl")
                      {
                        if (surrogate[config.iClm].Aaq_bins_init(b)>0.0)
                          surrogate[config.iClm].Aaq_bins_init(b)=factor*surrogate[config.iClm].MM/surrogate[i].MM*conc_equilibrium*
                            surrogate[i].Kaq(b)*AQinit(b)/sum+(1.0-factor)*surrogate[config.iClm].Aaq_bins_init(b);
                        else
                          {
                            surrogate[config.iClm].Aaq_bins_init(b)=0.01*factor*surrogate[config.iClm].MM/surrogate[i].MM*conc_equilibrium*
                              surrogate[i].Kaq(b)*AQinit(b)/sum;
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
                  total=conc_equilibrium/sum*surrogate[i].Kaq(b)*AQinit(b);
                else
                  total=(surrogate[config.iHSO4m].Aaq_bins_init(b)/surrogate[config.iHSO4m].MM+surrogate[config.iSO4mm].Aaq_bins_init(b)/surrogate[config.iSO4mm].MM)*surrogate[i].MM;
		double Keq=surrogate[i].Kequilibrium(Temperature)/chp(b)*surrogate[config.iHSO4m].gamma_aq_bins(b)/
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

void equilibrium_aq(model_config &config, vector<species>& surrogate, double &tequilibrium,
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
  double Kaq;
  //compute activity coefficients
  if (compute_activity_coefficients)
    {
      activity_coefficients_dyn_aq(config, surrogate, Temperature,AQinit,MOinit,
                                   conc_inorganic, ionic, ionic_organic,
                                   organion,chp,LWC,MMaq, factor);
      compute_kp_aq(config, surrogate, Temperature, ionic, chp, MMaq);
    }
  
  double conc_equilibrium;
  double sum=0.0;
  bool is_equilibrium=false; //indicates if there is concentrations at equilibrium
  
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
  
  //compute absorption of water by organic compounds
  if (config.compute_inorganic)
    equilibrium_inorg(config, surrogate, tequilibrium, AQinit, conc_inorganic, chp, LWC,
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
	    surrogate[config.iH2O].Aaq_bins_init(b)=factor*
	      max(surrogate[config.iH2O].Atot*surrogate[config.iH2O].Kaq(b)*AQinit(b)/sum-LWC(b),0.0)+
	      (1.0-factor)*surrogate[config.iH2O].Aaq_bins_init(b);
	  else	      
	    {
		double Mads=AQinit(b)-LWC(b)-conc_inorganic(b);
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

void equilibrium_tot(model_config &config, vector<species>& surrogate, double &tequilibrium,
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
  double Kaq;
  Array <double, 1> Kp;
  Kp.resize(config.max_number_of_phases);
  
  double conc_equilibrium;
  double sum=0.0;
  bool is_equilibrium=false;
  double LWCtot=0.0;
  int index_b,index_layer;
  for (b=0;b<config.nbins;++b)
    LWCtot+=LWC(b);
  
  //compute activity coefficients
  if (compute_activity_coefficients)
    {
      if (config.compute_organic)
        {
          activity_coefficients_dyn_org(config, surrogate, Temperature, MOW);
          compute_kp_org(config, surrogate, MOinit, Temperature, MOW);
        }
	  
      if (LWCtot>config.LWClimit)
        {
          activity_coefficients_dyn_aq(config, surrogate, Temperature,AQinit,
                                       MOinit, conc_inorganic, ionic, ionic_organic,
                                       organion,chp,LWC,MMaq, factor);
          compute_kp_aq(config, surrogate, Temperature, ionic, chp, MMaq);
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
  double kelvin_effect;
  if (config.compute_inorganic)
    equilibrium_inorg(config, surrogate, tequilibrium, AQinit, conc_inorganic, chp, LWC,
                      Temperature, RH, MMaq, factor);

  sum=1.0;
  double conceq=surrogate[config.iH2O].Ag;
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
	      surrogate[config.iH2O].Aaq_bins_init(b)=factor*
		max(conceq*surrogate[config.iH2O].Kaq(b)*AQinit(b)/sum - LWC(b), 0.0) +
		(1.0-factor)*surrogate[config.iH2O].Aaq_bins_init(b);
	    else	      
	      {
		double Mads=AQinit(b)-LWC(b)-conc_inorganic(b);
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

void equilibrium_org(model_config &config, vector<species>& surrogate, double &tequilibrium,
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
  int index_b,index_layer;
  Array<double, 1> Kp;
  Kp.resize(config.max_number_of_phases);

  //compute activity coefficients
  if (config.compute_organic and compute_activity_coefficients)
    {
      activity_coefficients_dyn_org(config, surrogate, Temperature, MOW);
      compute_kp_org(config, surrogate, MOinit, Temperature, MOW);
    }
  
  double conc_equilibrium;
  double sum=0.0;
  bool is_equilibrium=false;
  
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

      double sum=1.0;
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

void redistribution(model_config &config, vector<species>& surrogate,
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



void dynamic_org(model_config &config, vector<species>& surrogate,
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
  int index_b,index_layer;
  double tiny=1.0e-26;
  double conc_available=0.0;
  double sum_rates;
  double gamma=1.7071;
  Array<double, 3> MOinit2;
  Array<double,1> Kp;
  MOinit2.resize(config.nbins,config.nlayer,config.max_number_of_phases);
  Kp.resize(config.max_number_of_phases);
  double conc_equilibrium;
  
  //     gamma > .25     ->  A-stable  ( x(n+1)-x(n))*f(x(n)) > 0
  //     gamma < .25     ->  monotonic behaviour (x(n+1)-x*)(x(n)-x*) > 0
  //     gamma =  1+-sqrt(1/2) ->  L-stability
  
  //compute activity coefficients
  if (config.compute_organic)
    {
      if (compute_activity_coefficients)
        {
          activity_coefficients_dyn_org(config, surrogate, Temperature, MOW);
          compute_kp_org(config, surrogate, MOinit, Temperature, MOW);
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
      flux_org(config, surrogate, MOinit, AQinit, tiny, 0);

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
          activity_coefficients_dyn_org(config, surrogate, Temperature, MOW);
          compute_kp_org(config, surrogate, MOinit2, Temperature, MOW);
        }
  
      //second estimation of kinetic rates
      flux_org(config, surrogate, MOinit2, AQinit, tiny, 1);

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
      double sumconc,sumconc2;
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

void dynamic_inorg(model_config &config, vector<species>& surrogate,
                   Array <double, 1> &conc_inorganic, 
                   Array <double, 1> &LWC, double &DT2, double &tequilibrium, int index)
{
  //compute the dynamic evolution of the inorganic compounds concentrations with the roschem algorythm
  //compute_activity_coefficients: should the activity coefficients be recomputed?
  //tequilibrium: criterium to distinguish cases under equilibrium (characteristic time lower
  //                  than tequilibrium)

  int n=surrogate.size();
  int i,b;
  double sum=0.0;
  double tiny=1.0e-26;
  double conc_available=0.0;
  double sum_rates;
  double gamma=1.7071;

  for (b=0;b<config.nbins;++b)
    {
      LWC(b)=0.0;
      conc_inorganic(b)=surrogate[config.iHp].Aaq_bins_init(b);
    }

  if (index==0) //first estimation of concentrations
    {      
      for (i=0;i<n;++i)
        if(surrogate[i].is_inorganic_precursor)
          {
            conc_available=surrogate[i].Ag;
            sum_rates=0.0;
            for (b=0;b<config.nbins;++b)			  
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

                  if (conc_available-tmp> tiny and surrogate[i].k1_aq(b,0)>0.0)
                    surrogate[i].Jdn_aq(b,0)=-sum_rates/(conc_available-tmp);
                }
            
            for (b=0;b<config.nbins;++b)
              if (surrogate[i].time_aq(b)>=tequilibrium)              
                surrogate[i].k1_aq(b,0)=surrogate[i].k1_aq(b,0)*
                  DT2/(1.0-gamma*surrogate[i].Jdn_aq(b,0)*DT2);
 
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
                          max(0.01*surrogate[config.iHSO4m].Aaq_bins_init(b),
                              surrogate[config.iHSO4m].Aaq_bins_init(b)+f*surrogate[i].k1_aq(b,0)*surrogate[config.iHSO4m].MM/surrogate[i].MM);
                        surrogate[config.iSO4mm].Aaq_bins_init(b)=
                          max(0.01*surrogate[config.iSO4mm].Aaq_bins_init(b),
                              surrogate[config.iSO4mm].Aaq_bins_init(b)+(1.0-f)*surrogate[i].k1_aq(b,0)*surrogate[config.iSO4mm].MM/surrogate[i].MM);
                      }
                    
                    surrogate[i].Ag-=surrogate[config.iHSO4m].Aaq_bins_init(b)/surrogate[config.iHSO4m].MM*surrogate[i].MM
                      +surrogate[config.iSO4mm].Aaq_bins_init(b)/surrogate[config.iSO4mm].MM*surrogate[i].MM;
                    conc_inorganic(b)+=surrogate[config.iHSO4m].Aaq_bins_init(b)+surrogate[config.iSO4mm].Aaq_bins_init(b);

                  }

                if (surrogate[i].name=="NH3")
                  {
                    if (surrogate[i].time_aq(b)>=tequilibrium)
                      surrogate[config.iNH4p].Aaq_bins_init(b)=
                        max(0.01*surrogate[config.iNH4p].Aaq_bins_init(b),
                            surrogate[config.iNH4p].Aaq_bins_init(b)+surrogate[i].k1_aq(b,0)*surrogate[config.iNH4p].MM/surrogate[i].MM);

                    surrogate[i].Ag-=surrogate[config.iNH4p].Aaq_bins_init(b)/surrogate[config.iNH4p].MM*surrogate[i].MM;
                    conc_inorganic(b)+=surrogate[config.iNH4p].Aaq_bins_init(b);
                  }
                
                if (surrogate[i].name=="HNO3")
                  {
                    if (surrogate[i].time_aq(b)>=tequilibrium)
                      surrogate[config.iNO3m].Aaq_bins_init(b)=
                        max(0.01*surrogate[config.iNO3m].Aaq_bins_init(b),
                            surrogate[config.iNO3m].Aaq_bins_init(b)+surrogate[i].k1_aq(b,0)*surrogate[config.iNO3m].MM/surrogate[i].MM);

                    surrogate[i].Ag-=surrogate[config.iNO3m].Aaq_bins_init(b)/surrogate[config.iNO3m].MM*surrogate[i].MM;
                    conc_inorganic(b)+=surrogate[config.iNO3m].Aaq_bins_init(b);
                  }

                if (surrogate[i].name=="HCl")
                  {
                    if (surrogate[i].time_aq(b)>=tequilibrium)
                      surrogate[config.iClm].Aaq_bins_init(b)=
                        max(0.01*surrogate[config.iClm].Aaq_bins_init(b),
                            surrogate[config.iClm].Aaq_bins_init(b)+surrogate[i].k1_aq(b,0)*surrogate[config.iClm].MM/surrogate[i].MM);

                    surrogate[i].Ag-=surrogate[config.iClm].Aaq_bins_init(b)/surrogate[config.iClm].MM*surrogate[i].MM;
                    conc_inorganic(b)+=surrogate[config.iClm].Aaq_bins_init(b);
                  } 
              }            
          }  
      for (i=0;i<n;i++)
	if (surrogate[i].is_inorganic_precursor)
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
        if(surrogate[i].is_inorganic_precursor)
          {
            conc_available=surrogate[i].Ag;
            sum_rates=0.0;
            for (b=0;b<config.nbins;++b)			  
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

            for (b=0;b<config.nbins;++b)
              surrogate[i].k1_aq(b,1)=
                (surrogate[i].k1_aq(b,1)*DT2-gamma*DT2*(surrogate[i].Jdn_aq(b,0)+surrogate[i].Jdn_aq(b,1))
                 *surrogate[i].k1_aq(b,0))/(1.0-gamma*surrogate[i].Jdn_aq(b,1)*DT2);

	    surrogate[i].Ag=surrogate[i].Atot;

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
                        max(0.01*surrogate[config.iHSO4m].Aaq_bins_init0(b),
                            surrogate[config.iHSO4m].Aaq_bins_init0(b)+f*surrogate[i].k1_aq(b,1)*surrogate[config.iHSO4m].MM/surrogate[i].MM);
                      surrogate[config.iSO4mm].Aaq_bins(b)=
                        max(0.01*surrogate[config.iSO4mm].Aaq_bins_init0(b),
                            surrogate[config.iSO4mm].Aaq_bins_init0(b)+(1.0-f)*surrogate[i].k1_aq(b,1)*surrogate[config.iSO4mm].MM/surrogate[i].MM);
                    }
                  else
                    {                        
                      surrogate[config.iHSO4m].Aaq_bins(b)=surrogate[config.iHSO4m].Aaq_bins_init(b);
                      surrogate[config.iSO4mm].Aaq_bins(b)=surrogate[config.iSO4mm].Aaq_bins_init(b); 
                    }
                    
                  surrogate[i].Ag-=surrogate[config.iHSO4m].Aaq_bins(b)/surrogate[config.iHSO4m].MM*surrogate[i].MM
                    +surrogate[config.iSO4mm].Aaq_bins(b)/surrogate[config.iSO4mm].MM*surrogate[i].MM;

                  conc_inorganic(b)+=surrogate[config.iHSO4m].Aaq_bins(b)+surrogate[config.iSO4mm].Aaq_bins(b);
                }
              else if (surrogate[i].name=="NH3")
                {
                  if (surrogate[i].time_aq(b)>=tequilibrium)
                    surrogate[config.iNH4p].Aaq_bins(b)=
                      max(0.01*surrogate[config.iNH4p].Aaq_bins_init0(b),
                          surrogate[config.iNH4p].Aaq_bins_init0(b)+surrogate[i].k1_aq(b,1)*surrogate[config.iNH4p].MM/surrogate[i].MM);
                  else
                    surrogate[config.iNH4p].Aaq_bins(b)=surrogate[config.iNH4p].Aaq_bins_init0(b);
                    
                  surrogate[i].Ag-=surrogate[config.iNH4p].Aaq_bins(b)/surrogate[config.iNH4p].MM*surrogate[i].MM;
                  conc_inorganic(b)+=surrogate[config.iNH4p].Aaq_bins(b);
                }
              else if (surrogate[i].name=="HNO3")
                {
                  if (surrogate[i].time_aq(b)>=tequilibrium)
                    surrogate[config.iNO3m].Aaq_bins(b)=
                      max(0.01*surrogate[config.iNO3m].Aaq_bins_init0(b),
                          surrogate[config.iNO3m].Aaq_bins_init0(b)+surrogate[i].k1_aq(b,1)*surrogate[config.iNO3m].MM/surrogate[i].MM);
                  else
                    surrogate[config.iNO3m].Aaq_bins(b)=surrogate[config.iNO3m].Aaq_bins_init(b);

                  surrogate[i].Ag-=surrogate[config.iNO3m].Aaq_bins(b)/surrogate[config.iNO3m].MM*surrogate[i].MM;
                  conc_inorganic(b)+=surrogate[config.iNO3m].Aaq_bins(b);
                }

              else if (surrogate[i].name=="HCl")
                {
                  if (surrogate[i].time_aq(b)>=tequilibrium)
                    surrogate[config.iClm].Aaq_bins(b)=
                      max(0.01*surrogate[config.iClm].Aaq_bins_init0(b),
                          surrogate[config.iClm].Aaq_bins_init0(b)+surrogate[i].k1_aq(b,0)*surrogate[config.iClm].MM/surrogate[i].MM);
                  else
                    surrogate[config.iClm].Aaq_bins(b)=surrogate[config.iClm].Aaq_bins_init(b);

                  surrogate[i].Ag-=surrogate[config.iClm].Aaq_bins(b)/surrogate[config.iClm].MM*surrogate[i].MM;
                  conc_inorganic(b)+=surrogate[config.iClm].Aaq_bins(b);
                } 
          }

      for (i=0;i<n;i++)
	if (surrogate[i].is_inorganic_precursor) 
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

void dynamic_tot(model_config &config, vector<species>& surrogate,
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
  int index_b,index_layer;
  //double sumt=0.0;
  double LWCtot=0.0;
  double tiny=1.0e-26;
  double conc_available=0.0;
  double sum_rates;
  double gamma=1.7071;
  double Kaq;
  Array<double,1> Kp;
  Kp.resize(config.max_number_of_phases);
  Array<double, 3> MOinit2;
  Array<double, 1> AQinit2;
  Array <double, 1> chp2;
  double conc_equilibrium;
  double Ke=1.0e-14;
  chp2.resize(config.nbins);
  AQinit2.resize(config.nbins);
  MOinit2.resize(config.nbins,config.nlayer,config.max_number_of_phases);
  
  //     gamma > .25     ->  A-stable  ( x(n+1)-x(n))*f(x(n)) > 0
  //     gamma < .25     ->  monotonic behaviour (x(n+1)-x*)(x(n)-x*) > 0
  //     gamma =  1+-sqrt(1/2) ->  L-stability
  
  for (b=0;b<config.nbins;++b)
    LWCtot+=LWC(b);
    
  chp0=chp;
  double error_chp=1;
  int iter=0;

  while (iter<50 and error_chp>1.0e-3)
    {
      error_chp=0.0;
      iter++;
      if (config.compute_inorganic)	
	{
	  compute_ph_dyn(config, surrogate, Temperature, chp, AQinit, ionic, MMaq, LWC);
	  for (b=0;b<config.nbins;++b)
	    {	  
	      compute_ionic_strenght2(config, surrogate, AQinit(b), conc_inorganic(b), ionic(b), chp(b),
				      organion(b), ionic_organic(b), 1.0);
	      error_chp=max(error_chp,abs(chp0(b)-chp(b))/chp0(b));
	      chp0(b)=chp(b);	      
	    }      
	}
	  
      if (compute_activity_coefficients)
	{
	  if (config.compute_organic)
	    {
	      activity_coefficients_dyn_org(config, surrogate, Temperature, MOW);
	      compute_kp_org(config, surrogate, MOinit, Temperature, MOW);
	    }
	  
	  if (LWCtot>config.LWClimit)
	    {
	      activity_coefficients_dyn_aq(config, surrogate, Temperature,AQinit,
					   MOinit,conc_inorganic, ionic, ionic_organic,
                                       organion,chp,LWC,MMaq, 1.0);
	      compute_kp_aq(config, surrogate, Temperature, ionic, chp, MMaq);   
	    }
	}
    }

  //save initial concentrations 
  for (i=0;i<n;++i)
    {
      surrogate[i].Ag0=surrogate[i].Ag;
      surrogate[i].Ag1=surrogate[i].Ag;
      for (b=0;b<config.nbins;++b)		  
        for (ilayer=0;ilayer<config.nlayer;++ilayer)
          for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
            surrogate[i].Ap_layer_init0(b,ilayer,iphase)=
              surrogate[i].Ap_layer_init(b,ilayer,iphase);
		
      for (b=0;b<config.nbins;++b)
        surrogate[i].Aaq_bins_init0(b)=surrogate[i].Aaq_bins_init(b);
    }

  //compute kinetic rates
  flux_org(config, surrogate, MOinit, AQinit, tiny, 0);
  if (LWCtot>config.LWClimit)
    {
      flux_aq(config, surrogate, AQinit, MOinit, tiny, 0);      
      if (config.compute_inorganic)
	correct_flux_ph(config, surrogate, Temperature, AQinit, MOinit, chp, chp2, MMaq, ionic, LWC, tiny, DT2, 0);
    }
 
  //compute the first evaluation of concentrations
  for (i=0;i<n;++i)
    if((config.compute_organic and surrogate[i].is_organic) or i==config.iH2O) 
      {
	conc_available=surrogate[i].Ag;
	sum_rates=0.0;
	if (config.compute_organic)
	  if (surrogate[i].hydrophobic)
	    for (b=0;b<config.nbins;++b)
	      for (ilayer=0;ilayer<config.nlayer;++ilayer)
		for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)				
		  if (surrogate[i].time(b,ilayer,iphase)<tequilibrium)
		    conc_available+=surrogate[i].Ap_layer_init(b,ilayer,iphase);
		  else
		    sum_rates+=surrogate[i].k1(b,ilayer,iphase,0);
	
	if (surrogate[i].hydrophilic and LWCtot>config.LWClimit and i!=config.iH2O)
	  for (b=0;b<config.nbins;++b)
	    if (surrogate[i].time_aq(b)<tequilibrium)
                  conc_available+=surrogate[i].Aaq_bins_init(b);
	    else
	      sum_rates+=surrogate[i].k1_aq(b,0);
	
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
	
	if (surrogate[i].hydrophilic and LWCtot>config.LWClimit)
	  for (b=0;b<config.nbins;++b)			  
	    if (surrogate[i].time(b)>=tequilibrium)
	      if (conc_available-surrogate[i].Aaq_bins_init(b)> tiny and
		  surrogate[i].k1_aq(b,0)>0.0)
		surrogate[i].Jdn_aq(b,0)=-sum_rates/
		  (conc_available-surrogate[i].Aaq_bins_init(b));
	
	if (config.compute_organic)
	  if (surrogate[i].hydrophobic)
	    for (b=0;b<config.nbins;++b)		  
	      for (ilayer=0;ilayer<config.nlayer;++ilayer)
		for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		  surrogate[i].k1(b,ilayer,iphase,0)=surrogate[i].k1(b,ilayer,iphase,0)*
		    DT2/(1.0-gamma*surrogate[i].Jdn(b,ilayer,iphase,0)*DT2);
	
	if (surrogate[i].hydrophilic and LWCtot>config.LWClimit)
	  for (b=0;b<config.nbins;++b)
	    surrogate[i].k1_aq(b,0)=surrogate[i].k1_aq(b,0)*
	      DT2/(1.0-gamma*surrogate[i].Jdn_aq(b,0)*DT2);
	
	surrogate[i].Ag=surrogate[i].Atot;
	if (config.compute_organic)
	  if (surrogate[i].hydrophobic)
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
	
	if (surrogate[i].hydrophilic)
	  if(LWCtot>config.LWClimit)
	    for (b=0;b<config.nbins;++b)
	      {
		if (surrogate[i].time_aq(b)>=tequilibrium)
		  surrogate[i].Aaq_bins_init(b)=
		    max(0.01*surrogate[i].Aaq_bins_init(b),
			surrogate[i].Aaq_bins_init(b)+surrogate[i].k1_aq(b,0));
		
		surrogate[i].Ag-=surrogate[i].Aaq_bins_init(b);
	      }
	  else
	    for (b=0;b<config.nbins;++b)
	      surrogate[i].Aaq_bins_init(b)=0.0;
	else
	  for (b=0;b<config.nbins;++b)
	    surrogate[i].Aaq_bins_init(b)=surrogate[i].Aaq_bins_init0(b);
	
	surrogate[i].Ag=max(surrogate[i].Ag,0.01*surrogate[i].Ag0);
        
      }
    
  if (config.compute_inorganic)
    dynamic_inorg(config, surrogate, conc_inorganic, LWC, DT2, tequilibrium,0);

  for (i=0;i<n;++i)
    surrogate[i].Ag1=surrogate[i].Ag;

  //compute first evaluation of organic and aqueous masses 
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
	  compute_ph_dyn(config, surrogate, Temperature, chp, AQinit2, ionic, MMaq, LWC);	  
	  for (b=0;b<config.nbins;++b)
	    {	  
	      compute_ionic_strenght2(config, surrogate, AQinit2(b), conc_inorganic(b), ionic(b), chp(b),
				      organion(b), ionic_organic(b), 1.0);	   
	      error_chp=max(error_chp,abs(chp1(b)-chp(b))/chp1(b));
	      chp1(b)=chp(b);	      
	    } 
	}
               
      if (compute_activity_coefficients)
	{
	  //compute activity coefficients
	  if (config.compute_organic)
	    {
	      activity_coefficients_dyn_org(config, surrogate, Temperature, MOW);
	      compute_kp_org(config, surrogate, MOinit2, Temperature, MOW);
	    }
	  
	  if (LWCtot>config.LWClimit)
	    {
	      activity_coefficients_dyn_aq(config, surrogate, Temperature,AQinit2,
					   MOinit2, conc_inorganic, ionic, ionic_organic,
					   organion,chp,LWC,MMaq,1.0);
	      compute_kp_aq(config, surrogate, Temperature, ionic, chp, MMaq);
	    }
	}      
    }

  //compute the second evaluation of kinetic rates
  flux_org(config, surrogate, MOinit2, AQinit2, tiny, 1);
  if(LWCtot>config.LWClimit)
    {
      flux_aq(config, surrogate, AQinit2, MOinit2, tiny, 1);           
      if (config.compute_inorganic)
	correct_flux_ph(config, surrogate, Temperature, AQinit2, MOinit2, chp, chp2, MMaq, ionic, LWC, tiny, DT2, 1);
    }

  //compute the second evaluation of concentrations
  for (i=0;i<n;++i)
    if((config.compute_organic and surrogate[i].is_organic) or i==config.iH2O) 
      {
	conc_available=surrogate[i].Ag;
	sum_rates=0.0;
	if (config.compute_organic)
	  if (surrogate[i].hydrophobic)
	    for (b=0;b<config.nbins;++b)
	      for (ilayer=0;ilayer<config.nlayer;++ilayer)
		for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)		  
		  if (surrogate[i].time(b,ilayer,iphase)<tequilibrium)
		    conc_available+=surrogate[i].Ap_layer_init(b,ilayer,iphase);
		  else
		    sum_rates+=surrogate[i].k1(b,ilayer,iphase,1);
	
	if (surrogate[i].hydrophilic and LWCtot>config.LWClimit)
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
	
	if (surrogate[i].hydrophilic and LWCtot>config.LWClimit)
	  for (b=0;b<config.nbins;++b)			  
	    if (surrogate[i].time_aq(b)>=tequilibrium)
	      if (conc_available-surrogate[i].Aaq_bins_init(b)> tiny and
		  surrogate[i].k1_aq(b,1)>0.0)
		surrogate[i].Jdn_aq(b,1)=-sum_rates/
		  (conc_available-surrogate[i].Aaq_bins_init(b));
	
	if (config.compute_organic)
	  if (surrogate[i].hydrophobic)
	    for (b=0;b<config.nbins;++b)		  
	      for (ilayer=0;ilayer<config.nlayer;++ilayer)
		for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		  surrogate[i].k1(b,ilayer,iphase,1)=
		    (surrogate[i].k1(b,ilayer,iphase,1)*DT2
		     -gamma*DT2*(surrogate[i].Jdn(b,ilayer,iphase,0)+surrogate[i].Jdn(b,ilayer,iphase,1))
		     *surrogate[i].k1(b,ilayer,iphase,0))/
		    (1.0-gamma*surrogate[i].Jdn(b,ilayer,iphase,1)*DT2);
	
	if (surrogate[i].hydrophilic and LWCtot>config.LWClimit)
	  for (b=0;b<config.nbins;++b)
	    surrogate[i].k1_aq(b,1)=
	      (surrogate[i].k1_aq(b,1)*DT2-gamma*DT2*(surrogate[i].Jdn_aq(b,0)+surrogate[i].Jdn_aq(b,1))
	       *surrogate[i].k1_aq(b,0))/(1.0-gamma*surrogate[i].Jdn_aq(b,1)*DT2);

	surrogate[i].Ag=surrogate[i].Atot;
	if (config.compute_organic)
	  if (surrogate[i].hydrophobic)
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

	if (surrogate[i].hydrophilic and LWCtot>config.LWClimit)
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
    dynamic_inorg(config, surrogate, conc_inorganic, LWC, DT2, tequilibrium,1);  

  //make sure that the sum of concentrations is not higher than the total concentration
  if (config.compute_organic)
    {
      double sumt;
      double sumconc,sumconc2;
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
  
  for (b=0;b<config.nbins;++b)
    {
      AQ(b)=LWC(b);
      AQinit(b)=LWC(b);
      for (i=0;i<n;++i)
        if(surrogate[i].hydrophilic)
	  {
	    AQ(b)+=surrogate[i].Aaq_bins(b);
	    AQinit(b)+=surrogate[i].Aaq_bins(b);
	  }

    }

  if (config.compute_inorganic)
    {
      compute_ph_dyn2(config, surrogate, Temperature, chp, AQ, ionic, MMaq, LWC);
      for (b=0;b<config.nbins;++b)
	{	  
	  compute_ionic_strenght2(config, surrogate, AQ(b), conc_inorganic(b), ionic(b), chp(b),
				  organion(b), ionic_organic(b), 1.0);	 	  
	}    
    }
}

void dynamic_aq(model_config &config, vector<species>& surrogate,
                Array<double, 1> &AQinit,Array<double, 1> &AQ,
                Array<double, 3> &MOinit,
                Array <double, 1> &conc_inorganic,
                Array <double, 1> &ionic,Array <double, 1> &ionic_organic,
                Array <double, 1> &organion,Array <double, 1> &chp,
		Array <double, 1> &chp1, Array <double, 1> &chp0,
                Array <double, 1> &LWC,
                Array<double, 1> &MMaq,
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
  int index_b;
  double sum=0.0;
  double tiny=1.0e-26;
  double conc_available=0.0;
  double sum_rates;
  double gamma=1.7071;
  double Kaq;
  Array<double, 1> AQinit2,chp2;
  AQinit2.resize(config.nbins);
  chp2.resize(config.nbins);
  double conc_equilibrium;
  double Ke=1.0e-14;
  
  //     gamma > .25     ->  A-stable  ( x(n+1)-x(n))*f(x(n)) > 0
  //     gamma < .25     ->  monotonic behaviour (x(n+1)-x*)(x(n)-x*) > 0
  //     gamma =  1+-sqrt(1/2) ->  L-stability


  chp0=chp;
  double error_chp=1;
  int iter=0;
  while (iter<50 and error_chp>1.0e-3)
    {
      error_chp=0.0;
      iter++;
      if (config.compute_inorganic)	
	{
	  compute_ph_dyn(config, surrogate, Temperature, chp, AQinit, ionic, MMaq, LWC);
	  for (b=0;b<config.nbins;++b)
	    {	  
	      compute_ionic_strenght2(config, surrogate, AQinit(b), conc_inorganic(b), ionic(b), chp(b),
				      organion(b), ionic_organic(b), 1.0);
	      error_chp=max(error_chp,abs(chp0(b)-chp(b))/chp0(b));
	      chp0(b)=chp(b);	      
	    }      
	}
	  
      if (compute_activity_coefficients)	         
	{
	  activity_coefficients_dyn_aq(config, surrogate, Temperature,AQinit,
				       MOinit,conc_inorganic, ionic, ionic_organic,
                                       organion,chp,LWC,MMaq, 1.0);
	  compute_kp_aq(config, surrogate, Temperature, ionic, chp, MMaq);   	
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
  flux_aq(config, surrogate, AQinit, MOinit, tiny, 0);
  if (config.compute_inorganic)
    correct_flux_ph(config, surrogate, Temperature, AQinit, MOinit, chp, chp2, MMaq, ionic, LWC, tiny, DT2, 0);

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
    dynamic_inorg(config, surrogate, conc_inorganic, LWC, DT2, tequilibrium, 0);

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
	  double conc_org=LWC(b);
	  for (i=0;i<n;++i)
	    if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophilic)
	      conc_org+=surrogate[i].Aaq_bins_init(b);

	  double inorganion=0.0;
	  for (i=0;i<n;++i)
	    if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)
	      {
		surrogate[i].molality=surrogate[i].Aaq_bins_init(b)/surrogate[i].MM/conc_org*1000.0;		
		if (surrogate[i].name!="H")
		  inorganion-=surrogate[i].molality*surrogate[i].charge;
	      }
	  
	  chp(b)=0.5*(organion(b)+inorganion+pow(pow(organion(b)+inorganion,2)+4*Ke,0.5));
	  if (chp(b)==0.0)
	    chp(b)=pow(10.0,-5.6);	  
	  
	  compute_ionic_strenght2(config, surrogate, AQinit2(b), conc_inorganic(b), ionic(b), chp(b),
				  organion(b), ionic_organic(b), 1.0);

	}      
    }

  //compute activity coefficients
  if (compute_activity_coefficients)
    {
      activity_coefficients_dyn_aq(config, surrogate, Temperature,AQinit2,
                                   MOinit,conc_inorganic, ionic, ionic_organic,
                                   organion,chp,LWC,MMaq, 1.0);
      compute_kp_aq(config, surrogate, Temperature, ionic, chp, MMaq);
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
	  compute_ph_dyn(config, surrogate, Temperature, chp, AQinit2, ionic, MMaq, LWC);	  
	  for (b=0;b<config.nbins;++b)
	    {	  
	      compute_ionic_strenght2(config, surrogate, AQinit2(b), conc_inorganic(b), ionic(b), chp(b),
				      organion(b), ionic_organic(b), 1.0);	   
	      error_chp=max(error_chp,abs(chp1(b)-chp(b))/chp1(b));
	      chp1(b)=chp(b);	      
	    } 
	}
               
      if (compute_activity_coefficients)
	{
	  //compute activity coefficients	  	  	
	  activity_coefficients_dyn_aq(config, surrogate, Temperature,AQinit2,
				       MOinit, conc_inorganic, ionic, ionic_organic,
				       organion,chp,LWC,MMaq,1.0);
	  compute_kp_aq(config, surrogate, Temperature, ionic, chp, MMaq);	    
	}      
    }  

  //second estimation of kinetic rates
  flux_aq(config, surrogate, AQinit2, MOinit, tiny, 1);
  if (config.compute_inorganic)
    correct_flux_ph(config, surrogate, Temperature, AQinit2, MOinit, chp, chp2, MMaq, ionic, LWC, tiny, DT2, 1);

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
    dynamic_inorg(config, surrogate, conc_inorganic, LWC, DT2, tequilibrium, 1);
  
  //make sure that the sum of concentrations is not higher than the total concentration
  if (config.compute_organic)
    {
      double sumconc,sumconc2;
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
  
  for (b=0;b<config.nbins;++b)
    {
      AQ(b)=LWC(b)+conc_inorganic(b);
      for (i=0;i<n;++i)
        if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophilic)
          AQ(b)+=surrogate[i].Aaq_bins(b);
    }

  if (config.compute_inorganic)
    {
      compute_ph_dyn2(config, surrogate, Temperature, chp, AQ, ionic, MMaq, LWC);
      for (b=0;b<config.nbins;++b)
	{	  
	  compute_ionic_strenght2(config, surrogate, AQ(b), conc_inorganic(b), ionic(b), chp(b),
				  organion(b), ionic_organic(b), 1.0);	 	  
	}    
    }
}

void adapstep(model_config &config, vector<species>& surrogate, double &tequilibrium,
              double &deltat1, double &t, double &tend, double &deltatmin,
              Array<double, 3> &MOinit, Array<double, 3> &MO, double LWCtot,
              Array<double, 1> &AQinit, Array<double, 1> &AQ,
              Array<double, 1> &LWC, Array<double, 1> &conc_inorganic,
	      Array<double, 1> &chp, Array<double, 1> &chp1, Array<double, 1> &chp0)
{
  //estimate the new time step necessary 
  double n2err=0.0;
  double n2err2=0.0;
  double n2err_aq=0.0;
  double tinym=1.0e-5;
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

        if(sum1>config.Vlayer(ilayer)*tinym or sum2>config.Vlayer(ilayer)*tinym)
          n2err2=max(n2err2,abs(sum1-sum2)/(sum1));
      }

  if (LWCtot>config.LWClimit) 
    for (b=0;b<config.nbins;++b)
      if(AQ(b)>tinym or AQinit(b)>tinym)
	n2err2=max(n2err2,abs(AQinit(b)-AQ(b))/(AQinit(b)));

  if (config.compute_inorganic and LWCtot>config.LWClimit)
    {
      for (i=0;i<n;++i)
	{
	  if(surrogate[i].is_inorganic_precursor) 
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

      double chp1a;
      double chp2; 
      double chp0a;
      double Ke=1.0e-14;
      for (b=0;b<config.nbins;++b)
	{
	  double inorg0=0.0;
	  double inorg1=0.0;
	  double inorg2=0.0;
	  for (i=0;i<n;++i)
	    if (surrogate[i].is_organic==false and surrogate[i].is_inorganic_precursor==false and config.iH2O!=i and config.iHp!=i)
	      {
		inorg0-=surrogate[i].Aaq_bins_init0(b)/surrogate[i].MM*surrogate[i].charge*config.AQrho(b)/AQinit(b);
		inorg1-=surrogate[i].Aaq_bins_init(b)/surrogate[i].MM*surrogate[i].charge*config.AQrho(b)/AQinit(b);
		inorg2-=surrogate[i].Aaq_bins(b)/surrogate[i].MM*surrogate[i].charge*config.AQrho(b)/AQ(b);	      
	      }
	  chp0a=0.5*(inorg0+pow(pow(inorg0,2)+4*Ke,0.5));
	  chp1a=0.5*(inorg1+pow(pow(inorg1,2)+4*Ke,0.5));
	  chp2=0.5*(inorg2+pow(pow(inorg2,2)+4*Ke,0.5));
	  n2err2=max(n2err2,abs(chp(b)-chp1(b))/chp1(b));
	}
    }

  if (m>0)
    n2err=max(pow(n2err/(m*m),0.5),n2err2);
  if (maq>0)
    n2err=max(pow(n2err_aq/(maq*maq),0.5),n2err);
  n2err=min(n2err,config.EPSER*R*R);
  n2err=max(n2err,config.EPSER/(R*R));
  deltat1=deltat1*pow(config.EPSER/n2err,0.5);
  deltat1=max(deltatmin,deltat1);
  deltat1=min(tend-t,deltat1);
}

void compute_diameters(model_config &config, vector<species>& surrogate,
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
      }
    }
}

void phase_repartition(model_config &config,vector<species>& surrogate, double &Temperature,
                       Array<double, 3> &MOinit, Array<double, 3> &MO, Array<double, 3> &MOW,
                       double factor)
{
  //compute the separation of phases due to saturation by minimizing the gibbs energy of organic
  //phases
  int b,ilayer,iphase,i,jphase;
  int n=surrogate.size();
  double sumX,error,tot,temp;
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

void number_org_phases(model_config &config,vector<species>& surrogate, double &Temperature,
                       Array<double, 3> &MOinit, Array<double, 3> &MOW)
{
  //compute the separation of phases due to saturation by minimizing the gibbs energy of organic
  //phases
  int b,ilayer,iphase,i,jphase;
  int n=surrogate.size();
  int h,nh;
  int index;
  bool stable_system=false;
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
  int nphase_old,maxi;
  
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
                activity_coefficients_dyn_sat(config, surrogate, Temperature, MOW, b, ilayer);
				
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
                        activity_coefficients_dyn_sat(config, surrogate, Temperature, MOW, b, ilayer);

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


