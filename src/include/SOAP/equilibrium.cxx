
void initialisation_eq(model_config &config, vector<species>& surrogate, double &Temperature, double &ionic, double &chp, bool all_hydrophobic)
{ 
  int n=surrogate.size();
  int i;
  double MOW=1.0;
  double gamma=pow(10,-0.511*pow(298.0/Temperature,1.5)*pow(ionic,0.5)/(1.0+pow(ionic,0.5)));

  for (i=0;i<n;i++)
    {
      if (surrogate[i].hydrophilic==false)
	surrogate[i].Aaq=0.0;
      if (surrogate[i].hydrophobic==false)
	surrogate[i].Ap=0.0;
    }

  if (config.compute_organic)
    for (i=0;i<n;i++)
      if (surrogate[i].is_organic)
        if (surrogate[i].nonvolatile)
          {            	    
            if(surrogate[i].hydrophobic and surrogate[i].hydrophilic)
              if (surrogate[i].Aaq==0.0)
                {
                  surrogate[i].Ap=0.99*surrogate[i].Atot;
                  surrogate[i].Aaq=0.01*surrogate[i].Atot;
                }
          }
        else
          {	    
            if(surrogate[i].hydrophobic or all_hydrophobic)
              {        
                if (surrogate[i].kp_from_experiment)
                  surrogate[i].kpi=surrogate[i].Kp_exp_org(Temperature);
                else if (surrogate[i].kp_from_experiment==false)
                  surrogate[i].kpi=surrogate[i].Kp_eff_org(Temperature, MOW);
              } 
            if (surrogate[i].hydrophilic)
              if (config.compute_aqueous_phase_properties or chp == 0.0)
                surrogate[i].kaqi=surrogate[i].Kpart_aq(Temperature,MOW);              
              else
                {
                  if (surrogate[i].aq_type=="diacid")
                    {
                      surrogate[i].kaqi=surrogate[i].Kpart_aq(Temperature,MOW)*
                        (1.0+surrogate[i].Kacidity1/(pow(gamma,2)*chp)*
                         (1.0+surrogate[i].Kacidity2/(pow(gamma,2)*chp)));
                      surrogate[i].fioni1=(surrogate[i].Kacidity1/(pow(gamma,2)*chp))/
                        (1.0+surrogate[i].Kacidity1/(pow(gamma,2)*chp)*(1.0+surrogate[i].Kacidity2/(pow(gamma,2)*chp)));
                      surrogate[i].fioni2=(surrogate[i].Kacidity1/(pow(gamma,2)*chp))*(surrogate[i].Kacidity2/(pow(gamma,2)*chp))/
                        (1.0+surrogate[i].Kacidity1/(pow(gamma,2)*chp)*(1.0+surrogate[i].Kacidity2/(pow(gamma,2)*chp)));
                    }
                  else if (surrogate[i].aq_type=="monoacid")
                    {
                      surrogate[i].kaqi=surrogate[i].Kpart_aq(Temperature,MOW)*(1.0+surrogate[i].Kacidity1/(pow(gamma,2)*chp));
                      surrogate[i].fioni1=(surrogate[i].Kacidity1/(pow(gamma,2)*chp))/(1.0+surrogate[i].Kacidity1/(pow(gamma,2)*chp));
                    }
                  else if (surrogate[i].aq_type=="aldehyde")
                    surrogate[i].kaqi=surrogate[i].Kpart_aq(Temperature,MOW)*(1.0+surrogate[i].Koligo_aq*pow(gamma*chp/pow(10,-surrogate[i].pHref),surrogate[i].beta));
                  else
                    surrogate[i].kaqi=surrogate[i].Kpart_aq(Temperature,MOW);
                }
          }
  
  if (config.compute_inorganic)
    {
      for (i=0;i<n;i++)
        if (surrogate[i].is_inorganic_precursor)          
          surrogate[i].keq=surrogate[i].Kequilibrium(Temperature);

      Kpideal_inorganic(config, surrogate, Temperature);      
    }
}

void error_org(model_config &config, vector<species>& surrogate,double &MOinit,double &MOW,
               double &Temperature, double &error, double &derivative, double &RH,
               bool all_hydrophobic, bool compute_activity_coefficients)
{
  //Compute the value of the function error in the organic phase
  //error=MOinit-MO
  //MOinit: initial concentrations of the organic phase (µg/m3)
  //MO=sum of Ap=sum(Atot*Kpi*MOinit/(1+Kpi*MOinit)
  //Kp: partitioning constant of a species
  // Ap/Ag=Kp*MOinit
  //error=0.0 when the solution is found
  //derivative=d(error)/d(MOinit)1
  //all_hydrophobic:  do all the compounds condense on the organic?
  // (If low mass of particulate water)
  //MOW: mean molar mass of the organic phase (g/mol)
  int n=surrogate.size();
  int i;
  double MO=0.0;
  double Kp;

  MOinit=max(MOinit,config.MOmin);

  //compute activity coefficients
  if(config.compute_organic)
    {
      if (compute_activity_coefficients)
        activity_coefficients_org(config, surrogate, all_hydrophobic, Temperature, MOW);

      derivative=1.0;
      for (i=0;i<n;++i)
        if ((surrogate[i].hydrophobic and surrogate[i].is_organic) or
            (all_hydrophobic and surrogate[i].is_organic))
          if (surrogate[i].nonvolatile)
            {
              //If nonvolatile, the compounds is entirely in the organic phase
              surrogate[i].Ap=surrogate[i].Atot;
              MO+=surrogate[i].Ap;
            }
          else
            if (surrogate[i].kp_from_experiment)
              {
                Kp=surrogate[i].kpi;
                surrogate[i].Ap=surrogate[i].Atot*Kp*MOinit/(1+Kp*MOinit);
                MO+=surrogate[i].Ap;
                derivative+=surrogate[i].Atot*pow(Kp,2)*MOinit/(pow(1+Kp*MOinit,2))
                  -surrogate[i].Atot*Kp/(1+Kp*MOinit);
              }
            else
              {
                Kp=surrogate[i].kpi/MOW/surrogate[i].gamma_org;
                surrogate[i].Ap=surrogate[i].Atot*Kp*MOinit/(1+Kp*MOinit);
                MO+=surrogate[i].Ap;
                derivative+=surrogate[i].Atot*pow(Kp,2)*MOinit/(pow(1+Kp*MOinit,2))
                  -surrogate[i].Atot*Kp/(1+Kp*MOinit);
              }
  
      if (config.iH2O>=0 and config.hygroscopicity)
        if (surrogate[config.iH2O].hydrophobic) //Can H2O condense on the organic phase
          hygroscopicity_org(config, surrogate, Temperature, MOW, RH, MOinit, MO, derivative, 1.0);

      MO=max(MO,config.MOmin);
      
      error=MOinit-MO;
    }
  else
    {
      MO=MOinit;
      error=0.0;
    }
}

void error_ph(model_config &config, vector<species> &surrogate, double Temperature, double &chp, 
              double organion, double &error, double &derivative, double AQinit, double LWC)
{ 
  //This routine is used to compute the electroneutrality conditions with a 
  //method of newton raphson. The routine computes the error between two 
  //iterations and derivative of the error.
  //In this routine, organic ions are taken into account but are assumed not to
  //strongly impact the pH (the derivative of organic ions concentrations) do
  //not have to be taken into account. 
  int n=surrogate.size();
  int i;
  double inorganion=0.0;
  double total;
  double Ke=1.0e-14;
    
  double conc_org=LWC;
  for (i=0;i<n;++i)
    if (surrogate[i].is_organic or i==config.iH2O)      
      conc_org+=surrogate[i].Aaq;
      
  conc_org=max(conc_org,config.MOmin);

  derivative=-1.0;
  for (i=0;i<n;++i)
    if (surrogate[i].is_inorganic_precursor==true) 
      //adding concentrations of inorganic ions depending on pH
      {
        if (surrogate[i].name=="NH3")
          {
            total=surrogate[i].Ag/surrogate[i].MM+surrogate[config.iNH4p].Aaq/surrogate[config.iNH4p].MM; //gas + particle concentration
            derivative-=1000.*total*AQinit/conc_org*
              (surrogate[i].Kaq_inorg/chp/(1.0+surrogate[i].Kaq_inorg*AQinit)
               -surrogate[i].Kaq_inorg/pow(1.0+surrogate[i].Kaq_inorg*AQinit,2.0)*surrogate[i].Kaq_inorg/chp*AQinit);
            inorganion-=1000.*total*surrogate[i].Kaq_inorg/(1.0+surrogate[i].Kaq_inorg*AQinit); //concentration of NH3+
          }
        else if (surrogate[i].name=="HNO3")
          {
            total=surrogate[i].Ag/surrogate[i].MM+surrogate[config.iNO3m].Aaq/surrogate[config.iNO3m].MM; //gas + particle concentration
            derivative+=1000.*total*AQinit/conc_org*
              (-surrogate[i].Kaq_inorg/chp/(1.0+surrogate[i].Kaq_inorg*AQinit)
               +surrogate[i].Kaq_inorg/pow(1.0+surrogate[i].Kaq_inorg*AQinit,2.0)*surrogate[i].Kaq_inorg/chp*AQinit);
            inorganion+=1000.*total*surrogate[i].Kaq_inorg/(1.0+surrogate[i].Kaq_inorg*AQinit); //concentration of NO3-
          }
        else if (surrogate[i].name=="HCl")
          {
            total=surrogate[i].Ag/surrogate[i].MM+surrogate[config.iClm].Aaq/surrogate[config.iClm].MM; //gas + particle concentration
            derivative+=1000.*total*AQinit/conc_org*
              (-surrogate[i].Kaq_inorg/chp/(1.0+surrogate[i].Kaq_inorg*AQinit)
               +surrogate[i].Kaq_inorg/pow(1.0+surrogate[i].Kaq_inorg*AQinit,2.0)*surrogate[i].Kaq_inorg/chp*AQinit);
            inorganion+=1000.*total*surrogate[i].Kaq_inorg/(1.0+surrogate[i].Kaq_inorg*AQinit); //concentration of Cl-
          }
        else if (surrogate[i].name=="H2SO4")
          {
            total=surrogate[i].Ag/surrogate[i].MM+surrogate[config.iHSO4m].Aaq/surrogate[config.iHSO4m].MM
                  +surrogate[config.iSO4mm].Aaq/surrogate[config.iSO4mm].MM; //gas + particle concentration
            double K=surrogate[i].keq*surrogate[config.iHSO4m].gamma_aq
              /(surrogate[config.iHp].gamma_aq*surrogate[config.iSO4mm].gamma_aq);
	    
            derivative-=1000.*total/conc_org*1.0/pow(1.0+K/chp,2.0)*K/(chp*chp); //HSO4-+SO4--
            inorganion+=1000.*total/AQinit*(2.0-1.0/(1.0+K/chp));
            //concentration of HSO4- + 2*SO4--
          }
        else
          inorganion-=surrogate[i].Aaq/surrogate[i].MM/AQinit*1000.*surrogate[i].charge;
      }
    else 
      //adding concentrations of nonvolatile inorganic ions (not depending on pH)
      if (surrogate[i].name=="Na") 
        inorganion-=surrogate[i].Aaq/surrogate[i].MM/AQinit*1000.*surrogate[i].charge;

  inorganion*=AQinit/conc_org;

  derivative=derivative*0.5*(1.0+(organion+inorganion)/pow(pow(organion+inorganion,2)+4*Ke,0.5))-1.0;
  error=0.5*(organion+inorganion+pow(pow(organion+inorganion,2)+4*Ke,0.5))-chp;
}

void error_aq(model_config &config, vector<species>& surrogate,
              double &AQinit,double &LWC, double &conc_inorganic,
              double &ionic, double &chp, double &MMaq,
              double &Temperature, double &error, double &derivative, double &RH,
              double &organion, double &ionic_organic, double factor,
              bool compute_activity_coefficients)
{
  //Compute the value of the function error in the aqueous phase
  //error=AQinit-AQ
  //AQinit: initial concentrations of the aqueous phase (µg/m3)
  //AQ=sum of Aaq=sum(Atot*Kaqi*AQinit/(1+Kaqi*AQinit)+LWC+conc_inorganic
  //Kaq: partitioning constant of a species
  // Aaq/Ag=Kaq*AQinit
  //error=0.0 when the solution is found
  //derivative=d(error)/d(AQinit)
  
  //conc_inorganic: concentrations of inorganic ions
  //LWC liquid water content due to inorganic ions
  //ionic: ionic strength (mol/kg)
  //ionic_organic: ionic_strength due organic ions
  //organion: sum( -charge of organic ion molality of the organic ion)
  int n=surrogate.size();
  int i;
  double AQ;
  double Kp;
  double XH2O;
  //MMaq: Mean molar mass of the aqueous phase

  AQinit=max(AQinit,config.MOmin);

  compute_ionic_strenght2(config, surrogate, AQinit, conc_inorganic, ionic, chp, organion, ionic_organic, factor);

  //initialize AQ
  AQ=0.0;

  //compute acitivity coefficients and MMaq
  if (compute_activity_coefficients)
    {
      activity_coefficients_aq(config, surrogate, Temperature, LWC, MMaq, XH2O);
      if (config.compute_long_and_medium_range_interactions)
        activity_coefficients_LR_MR(config, surrogate, Temperature, LWC, ionic);
    }

  double conc_org=LWC;
  for (i=0;i<n;++i)
    if (surrogate[i].is_organic or i==config.iH2O)      
      conc_org+=surrogate[i].Aaq;      
  conc_org=max(conc_org,config.MOmin);

  //pH computation
  if (config.compute_aqueous_phase_properties)
    {
      //If inorganic ion concentrations are computed by SOAP, used a method of
      //newton raphson
      if (config.compute_inorganic)
        {
          double error_h=1000.0;
          double derivative_h;
          double chp2=chp;
          while(abs(error_h/chp2)>1.0e-3)
            {
              Kpreal_inorganic(config, surrogate, chp2, MMaq);
              error_ph(config, surrogate, Temperature, chp2, organion, error_h, derivative_h,AQinit,LWC);
              if (chp2-error_h/derivative_h>0.0 and derivative_h!=0.0)
		chp2=chp2-error_h/derivative_h;
	      else
		{		 
		  if (chp2+error_h<=0.0)
		    chp2=0.99*chp2;
		  else
		    chp2=chp2+error_h;
		}	    	      
            }
          chp=factor*chp2+(1.0-factor)*chp;
        }
      else
        {
          double Ke=1.0e-14;
          double inorganion=0.0;
          for (i=0;i<n;++i)
            if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)
              if (surrogate[i].name!="H")
                inorganion-=surrogate[i].molality*surrogate[i].charge;
      
          chp=factor*0.5*(organion+inorganion+pow(pow(organion+inorganion,2)+4*Ke,0.5))+(1.0-factor)*chp;
        }
    }

  derivative=1.0;
  double organion_tmp=0.0;
  double fion1,fion2,molality1,molality2;
  double ionic_organic_tmp=0.0;
  if (config.compute_organic)
    {
      for (i=0;i<n;++i)
        if (surrogate[i].hydrophilic and surrogate[i].is_organic)
          if (surrogate[i].nonvolatile)
            {
              surrogate[i].Aaq=surrogate[i].Atot;
              AQ+=surrogate[i].Aaq;
            }
          else
            {
              fion1=0.0;
              fion2=0.0;
	      Kp=surrogate[i].Kp_eff_aqreal(config, Temperature, ionic, chp,surrogate[config.iHp].gamma_LR,
                                            surrogate[config.iHp].gamma_SRMR,MMaq,fion1,fion2)
                /surrogate[i].gamma_aq;
              surrogate[i].Aaq=factor*surrogate[i].Atot*Kp*AQinit/(1+Kp*AQinit)+
                (1.0-factor)*surrogate[i].Aaq;
              AQ+=surrogate[i].Aaq;

              derivative+=surrogate[i].Atot*pow(Kp,2)*AQinit/(pow(1+Kp*AQinit,2))
                -surrogate[i].Atot*Kp/(1+Kp*AQinit);
              //molality1: molality of ions HA- or A-  
              molality1=surrogate[i].Aaq*fion1/surrogate[i].MM/conc_org*1000.0;
              //molality2: molality of ions A2-
              molality2=surrogate[i].Aaq*fion2/surrogate[i].MM/conc_org*1000.0;
              //compute ionic_organic and organion
              ionic_organic_tmp+=0.5*molality1+0.5*molality2*4;
              organion_tmp+=molality1+2*molality2;
            }
    }
  else
    for (i=0;i<n;++i)
      if (surrogate[i].hydrophilic and surrogate[i].is_organic)
        AQ+=surrogate[i].Aaq;

  if (config.compute_inorganic)
    {
      //compute the concentrations of inorganics
      double total;
      Kpreal_inorganic(config, surrogate, chp, MMaq);
      for (i=0;i<n;i++)
        if (surrogate[i].is_inorganic_precursor)
	  {
            if (surrogate[i].name=="H2SO4") //sulfate
	      {
                total=surrogate[i].Ag+surrogate[config.iHSO4m].Aaq/surrogate[config.iHSO4m].MM*surrogate[i].MM
                  +surrogate[config.iSO4mm].Aaq/surrogate[config.iSO4mm].MM*surrogate[i].MM; //total concentrations 
                double Keq=surrogate[i].keq/chp*surrogate[config.iHSO4m].gamma_aq/surrogate[config.iHp].gamma_aq/surrogate[config.iSO4mm].gamma_aq;
                surrogate[i].Ag=(1.0-factor)*surrogate[i].Ag;
                surrogate[config.iHSO4m].Aaq=factor*total*surrogate[config.iHSO4m].MM/surrogate[i].MM/(1.0+Keq)
                  +(1.0-factor)*surrogate[config.iHSO4m].Aaq; //HSO4-
                surrogate[config.iSO4mm].Aaq=factor*total*surrogate[config.iSO4mm].MM/surrogate[i].MM*Keq/
                  (1.0+Keq)+(1.0-factor)*surrogate[config.iSO4mm].Aaq; //SO4--
	      }

            if (surrogate[i].name=="NH3") //ammoniac
	      {
                total=surrogate[i].Ag+surrogate[config.iNH4p].Aaq/surrogate[config.iNH4p].MM*surrogate[i].MM; //total concentrations
                surrogate[i].Ag=factor*total/(1.0+surrogate[i].Kaq_inorg*AQinit)+(1.0-factor)*surrogate[i].Ag;
                surrogate[config.iNH4p].Aaq=factor*total*surrogate[i].Kaq_inorg*AQinit/(1.0+surrogate[i].Kaq_inorg*AQinit)*surrogate[config.iNH4p].MM/surrogate[i].MM
                  +(1.0-factor)*surrogate[config.iNH4p].Aaq; //NH4+
                derivative+=total*surrogate[config.iNH4p].MM/surrogate[i].MM*pow(surrogate[i].Kaq_inorg,2)*AQinit/(pow(1.+surrogate[i].Kaq_inorg*AQinit,2))
                  -total*surrogate[config.iNH4p].MM/surrogate[i].MM*surrogate[i].Kaq_inorg/(1.+surrogate[i].Kaq_inorg*AQinit);
	      }

            if (surrogate[i].name=="HNO3") //nitrate
	      {
                total=surrogate[i].Ag+surrogate[config.iNO3m].Aaq/surrogate[config.iNO3m].MM*surrogate[i].MM; //total concentrations
                surrogate[i].Ag=factor*total/(1.0+surrogate[i].Kaq_inorg*AQinit)+(1.0-factor)*surrogate[i].Ag;
                surrogate[config.iNO3m].Aaq=factor*total*surrogate[i].Kaq_inorg*AQinit/(1.0+surrogate[i].Kaq_inorg*AQinit)*surrogate[config.iNO3m].MM/surrogate[i].MM
                  +(1.0-factor)*surrogate[config.iNO3m].Aaq; //NO3-
                derivative+=total*surrogate[config.iNO3m].MM/surrogate[i].MM*pow(surrogate[i].Kaq_inorg,2)*AQinit/(pow(1.+surrogate[i].Kaq_inorg*AQinit,2))
                -total*surrogate[config.iNO3m].MM/surrogate[i].MM*surrogate[i].Kaq_inorg/(1.+surrogate[i].Kaq_inorg*AQinit);
	      }

            if (surrogate[i].name=="HCl") //chloride
	      {
                total=surrogate[i].Ag+surrogate[config.iClm].Aaq/surrogate[config.iClm].MM*surrogate[i].MM; //total concentrations
                surrogate[i].Ag=factor*total/(1.0+surrogate[i].Kaq_inorg*AQinit)+(1.0-factor)*surrogate[i].Ag;
                surrogate[config.iClm].Aaq=factor*total*surrogate[i].Kaq_inorg*AQinit/(1.0+surrogate[i].Kaq_inorg*AQinit)*surrogate[config.iClm].MM/surrogate[i].MM
                  +(1.0-factor)*surrogate[config.iClm].Aaq; //Cl-
                derivative+=total*surrogate[config.iClm].MM/surrogate[i].MM*pow(surrogate[i].Kaq_inorg,2)*AQinit/(pow(1.+surrogate[i].Kaq_inorg*AQinit,2))
                  -total*surrogate[config.iClm].MM/surrogate[i].MM*surrogate[i].Kaq_inorg/(1.+surrogate[i].Kaq_inorg*AQinit);
	      }
	  }

      conc_inorganic=0.0;
      for (i=0;i<n;i++)
        if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)
          conc_inorganic+=surrogate[i].Aaq; //compute total inorganic concentration

      AQ+=conc_inorganic;
      if (config.iH2O>=0)
        //compute hygroscopicity
        {
	  LWC=0.0; //In the case where SOAP compute inorganics, LWC is only used for the initialization
	  hygroscopicity_tot(config, surrogate, Temperature, RH, AQinit, conc_inorganic,MMaq,AQ,derivative,factor);
        }

    }
  else
    {
      AQ+=LWC+conc_inorganic;  
      if (config.iH2O>=0) //compute hygroscopicity
        if (surrogate[config.iH2O].hydrophilic and config.hygroscopicity) //Can H2O condense on the organic phase
          hygroscopicity_aq(config, surrogate, Temperature, RH, AQinit, LWC,
                            conc_inorganic,MMaq,AQ,derivative,factor);
    }

  AQ=max(AQ,config.MOmin);
  error=AQinit-AQ;
  ionic_organic=ionic_organic_tmp;
  organion=organion_tmp;
}

void error_coupled(model_config &config, vector<species>& surrogate,
                   double &MOinit, double &MOW, double &MMaq,
                   double &AQinit, double &LWC, double &conc_inorganic,
                   double &ionic, double &ionic_organic, double &chp, double &organion, 
                   double &Temperature, double &RH,
                   double &error1, double &deriv_error1_MO, double &deriv_error1_AQ,
                   double &error2, double &deriv_error2_MO, double &deriv_error2_AQ,
                   double factor, bool compute_activity_coefficients)
{
  //compute the errors when the system is coupled
  //system coupled: if at least one compound condense on both the organic and the aqueous phases
  // the concentrations of both phases has to be solved simultaneously
  //error1=MOinit-MO
  //error2=AQinit-AQ
  //AQinit: initial concentrations of the aqueous phase (µg/m3)
  //MOinit: initial concentrations of the aqueous phase (µg/m3)
  //MO=sum of Ap=sum(Atot*Kpi*MOinit/(1+Kaqi*AQinit+Kpi*MOinit)
  //AQ=sum of Aaq=sum(Atot*Kaqi*AQinit/(1+Kaqi*AQinit+Kpi*MOinit)+LWC+conc_inorganic
  //Kp and Kaq: partitioning constants of a species in the two phases
  // Ap/Ag=Kp*MOinit and Aaq/Ag=Kaq*AQinit 
  //error1=0.0 and error2=0.0 when the solution is found
  //derivative_error1_MO=d(error1)/d(MOinit)
  //derivative_error1_AQ=d(error1)/d(AQinit)
  //derivative_error2_MO=d(error2)/d(MOinit)
  //derivative_error2_AQ=d(error2)/d(AQinit)

  //conc_inorganic: concentrations of inorganic ions (µg/m3)
  //LWC liquid water content due to inorganic ions (µg/m3)
  //ionic: ionic strength (mol/kg)
  //ionic_organic: ionic_strength due organic ions (mol/kg)
  //organion: sum( -charge of organic ion *molality of the organic ion)
  //MMaq: mean molar mass of the aqueous phase (g/mol)
  //MOW: mean molar mass of the organic phase (g/mol)
  int n=surrogate.size();
  int i;
  double AQ,MO;
  double Kp;
  double XH2O;
  double chp2=chp;

  compute_ionic_strenght2(config, surrogate, AQinit, conc_inorganic, ionic, chp2,
                          organion, ionic_organic, factor);
  //initialize AQ and MO
  AQ=0.0;
  MO=0.0;
  deriv_error1_MO=1.0;
  deriv_error1_AQ=0.0;
  deriv_error2_MO=0.0;
  deriv_error2_AQ=1.0;

  AQinit=max(AQinit,config.MOmin);
  MOinit=max(MOinit,config.MOmin);

  double conc_org=LWC;
  for (i=0;i<n;++i)
    if (surrogate[i].is_organic or i==config.iH2O)      
      conc_org+=surrogate[i].Aaq;
      
  conc_org=max(conc_org,config.MOmin);

  //compute acitivity coefficients, MOW and MMaq
  if (compute_activity_coefficients)
    {
      if (config.compute_organic)
        activity_coefficients_org(config, surrogate, false, Temperature, MOW);
      activity_coefficients_aq(config, surrogate, Temperature, LWC, MMaq, XH2O);
      if (config.compute_long_and_medium_range_interactions)
        activity_coefficients_LR_MR(config, surrogate, Temperature, LWC, ionic);
    }
  
  if (config.compute_aqueous_phase_properties)
    {
      //If inorganic ion concentrations are computed by SOAP, used a method of
      //newton raphson
      if (config.compute_inorganic)
	{
	  double error_h=1000.0;
	  double derivative_h;
	  chp2=chp;
	  int index=0;
	  while(abs(error_h/chp2)>1.0e-3 and index<10000)
	    {
	      Kpreal_inorganic(config, surrogate, chp2, MMaq);
	      error_ph(config, surrogate, Temperature, chp2, organion, error_h, derivative_h,AQinit,LWC);
	      if (chp2-error_h/derivative_h>0.0 and derivative_h!=0.0)
		chp2=chp2-error_h/derivative_h;
	      else
		chp2=chp2+error_h;

	      chp2=max(chp2,1.0e-20);	      
	      index++;
	    }
	  chp=factor*chp2+(1.0-factor)*chp;
	  Kpreal_inorganic(config, surrogate, chp, MMaq);
	}
      else
	{
	  double Ke=1.0e-14;
	  double inorganion=0.0;
	  for (i=0;i<n;++i)
	    if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)
	      if (surrogate[i].name!="H")
		inorganion-=surrogate[i].molality*surrogate[i].charge;
      
	  chp=factor*0.5*(organion+inorganion+pow(pow(organion+inorganion,2)+4*Ke,0.5))+(1.0-factor)*chp;
	}
    }
  
  double organion_tmp=0.0;
  double fion1,fion2,molality1,molality2;
  double ionic_organic_tmp=0.0;
  double Kp_aq,Kp_org;
  if (config.compute_organic)
    {
      for (i=0;i<n;++i)
        if (surrogate[i].hydrophobic and surrogate[i].is_organic
            and surrogate[i].hydrophilic==false)
          //compute absorption for compounds which are only hydrophobic
          if (surrogate[i].nonvolatile)
            {
              surrogate[i].Ap=surrogate[i].Atot;
              surrogate[i].Aaq=0.0;
              MO+=surrogate[i].Ap;
            }
          else
            if (surrogate[i].kp_from_experiment)
              {
                Kp=surrogate[i].kpi;
                surrogate[i].Ap=factor*surrogate[i].Atot*Kp*MOinit/(1+Kp*MOinit)+(1.0-factor)*surrogate[i].Ap;
                surrogate[i].Aaq=0.0;
                MO+=surrogate[i].Ap;
                deriv_error1_MO+=surrogate[i].Atot*pow(Kp,2)*MOinit/(pow(1+Kp*MOinit,2))
                  -surrogate[i].Atot*Kp/(1+Kp*MOinit);
              }
            else
              {
                Kp=surrogate[i].kpi/MOW/surrogate[i].gamma_org;
                surrogate[i].Aaq=0.0;
                surrogate[i].Ap=factor*surrogate[i].Atot*Kp*MOinit/(1+Kp*MOinit)+(1.0-factor)*surrogate[i].Ap;
                MO+=surrogate[i].Ap;
                deriv_error1_MO+=surrogate[i].Atot*pow(Kp,2)*MOinit/(pow(1+Kp*MOinit,2))
                  -surrogate[i].Atot*Kp/(1+Kp*MOinit);
              }
        else if (surrogate[i].hydrophobic==false and surrogate[i].is_organic
                 and surrogate[i].hydrophilic)
          //compute absorption for compounds which are only hydrophilic
          if (surrogate[i].nonvolatile)
            {
              surrogate[i].Aaq=surrogate[i].Atot;
              AQ+=surrogate[i].Aaq;
              surrogate[i].Ap=0.0;
            }
  
          else
            {
              fion1=0.0;
              fion2=0.0;
              Kp=surrogate[i].Kp_eff_aqreal(config, Temperature, ionic, chp,surrogate[config.iHp].gamma_LR,
                                            surrogate[config.iHp].gamma_SRMR,MMaq, fion1, fion2)
                /surrogate[i].gamma_aq;
              surrogate[i].Aaq=factor*surrogate[i].Atot*Kp*AQinit/(1+Kp*AQinit)+(1.0-factor)*surrogate[i].Aaq;
              surrogate[i].Ap=0.0;
              AQ+=surrogate[i].Aaq;
              deriv_error2_AQ+=surrogate[i].Atot*pow(Kp,2)*AQinit/(pow(1+Kp*AQinit,2))
                -surrogate[i].Atot*Kp/(1+Kp*AQinit);
              //molality1: molality of ions HA- or A-  
              molality1=surrogate[i].Aaq*fion1/surrogate[i].MM/conc_org*1000.0;
              //molality2: molality of ions A2-
              molality2=surrogate[i].Aaq*fion2/surrogate[i].MM/conc_org*1000.0;
              //compute ionic_organic and organion
              ionic_organic_tmp+=0.5*molality1+0.5*molality2*4;
              organion_tmp+=molality1+2*molality2;
            }
        else if (surrogate[i].hydrophobic and surrogate[i].is_organic
                 and surrogate[i].hydrophilic)
          //compute absorption for compounds which are both hydrophobic and hydrophilic
          if (surrogate[i].nonvolatile)
            if (AQinit==0.0)
              {
                surrogate[i].Ap=surrogate[i].Atot;
                surrogate[i].Aaq=0.0;
                MO+=surrogate[i].Ap;
              }
            else
              {
                //If the compounds is the nonvolatile the concentrations has to respect the
                // following conditions:
                //   activity_organic_phase = activity_aqueous_phase
                //   Aaq/Ap = (AQinit/MOinit)*gamma_org/gamma_aq)
                surrogate[i].Aaq=factor*surrogate[i].Atot/
                  (1.0+MOinit/AQinit*surrogate[i].gamma_aq/surrogate[i].gamma_org)+(1.0-factor)*surrogate[i].Aaq;
                surrogate[i].Ap=surrogate[i].Atot-surrogate[i].Aaq;
                MO+=surrogate[i].Ap;
                AQ+=surrogate[i].Aaq;
                deriv_error1_MO+=-surrogate[i].Atot/
                  (pow(1+MOinit/AQinit*surrogate[i].gamma_aq/surrogate[i].gamma_org,2))
                  *(surrogate[i].gamma_aq/surrogate[i].gamma_org)/AQinit;
                deriv_error1_AQ+=surrogate[i].Atot/
                  (pow(1+MOinit/AQinit*surrogate[i].gamma_aq/surrogate[i].gamma_org,2))
                  *(surrogate[i].gamma_aq/surrogate[i].gamma_org)*MOinit/(pow(AQinit,2));
                deriv_error2_AQ+=-surrogate[i].Atot/
                  (pow(1+MOinit/AQinit*surrogate[i].gamma_aq/surrogate[i].gamma_org,2))
                  *(surrogate[i].gamma_aq/surrogate[i].gamma_org)*MOinit/(pow(AQinit,2));
                deriv_error2_MO+=surrogate[i].Atot/
                  (pow(1+MOinit/AQinit*surrogate[i].gamma_aq/surrogate[i].gamma_org,2))
                  *(surrogate[i].gamma_aq/surrogate[i].gamma_org)/AQinit;
              }
          else
            {
              fion1=0.0;
              fion2=0.0;
              Kp_aq=surrogate[i].Kp_eff_aqreal(config, Temperature, ionic, chp,
                                               surrogate[config.iHp].gamma_LR,
                                               surrogate[config.iHp].gamma_SRMR,
                                               MMaq, fion1, fion2)
                /surrogate[i].gamma_aq;
              if (surrogate[i].kp_from_experiment)
                Kp_org=surrogate[i].kpi;
              else
                Kp_org=surrogate[i].kpi/MOW/surrogate[i].gamma_org;
          
              surrogate[i].Aaq=factor*surrogate[i].Atot*Kp_aq*AQinit/(1+Kp_aq*AQinit+Kp_org*MOinit)+(1.0-factor)*surrogate[i].Aaq;
              surrogate[i].Ap=factor*surrogate[i].Atot*Kp_org*MOinit/(1+Kp_aq*AQinit+Kp_org*MOinit)+(1.0-factor)*surrogate[i].Ap;
              AQ+=surrogate[i].Aaq;
              MO+=surrogate[i].Ap;
              deriv_error1_MO+=-surrogate[i].Atot*Kp_org/(1+Kp_aq*AQinit+Kp_org*MOinit)
                +surrogate[i].Atot*pow(Kp_org,2)*MOinit/(pow(1+Kp_aq*AQinit+Kp_org*MOinit,2));
              deriv_error1_AQ+=surrogate[i].Atot*Kp_org*Kp_aq*MOinit/
                (pow(1+Kp_aq*AQinit+Kp_org*MOinit,2));
              deriv_error2_MO+=surrogate[i].Atot*Kp_org*Kp_aq*AQinit/
                (pow(1+Kp_aq*AQinit+Kp_org*MOinit,2));
              deriv_error2_AQ+=-surrogate[i].Atot*Kp_aq/(1+Kp_aq*AQinit+Kp_org*MOinit)
                +surrogate[i].Atot*pow(Kp_aq,2)*AQinit/(pow(1+Kp_aq*AQinit+Kp_org*MOinit,2));
              //molality1: molality of ions HA- or A-
              molality1=surrogate[i].Aaq*fion1/surrogate[i].MM/conc_org*1000.0;
              //molality2: molality of ions A2-
              molality2=surrogate[i].Aaq*fion2/surrogate[i].MM/conc_org*1000.0;
              //compute ionic_organic and organion
              ionic_organic_tmp+=0.5*molality1+0.5*molality2*4;
              organion_tmp+=molality1+2*molality2;
            }
    }
  else
    {
      MO=MOinit;
      for (i=0;i<n;++i)
        if (surrogate[i].hydrophilic and surrogate[i].is_organic)
          AQ+=surrogate[i].Aaq;
    } 

  //compute inorganic ion concentrations
  if (config.compute_inorganic)
    {
      double total;
      conc_inorganic=0.0;

      for (i=0;i<n;i++)
        if (surrogate[i].is_inorganic_precursor)
	  {
            if (surrogate[i].name=="H2SO4") //sulfate
	      {
                total=surrogate[i].Ag+surrogate[config.iHSO4m].Aaq/surrogate[config.iHSO4m].MM*surrogate[i].MM
                  +surrogate[config.iSO4mm].Aaq/surrogate[config.iSO4mm].MM*surrogate[i].MM; //gas+particle concentrations
                double Keq=surrogate[i].keq/chp*surrogate[config.iHSO4m].gamma_aq/surrogate[config.iHp].gamma_aq/surrogate[config.iSO4mm].gamma_aq;
                surrogate[i].Ag=(1.0-factor)*surrogate[i].Ag;
                surrogate[config.iHSO4m].Aaq=factor*total*surrogate[config.iHSO4m].MM/surrogate[i].MM/(1.0+Keq)
                  +(1.0-factor)*surrogate[config.iHSO4m].Aaq; //HSO4-
                surrogate[config.iSO4mm].Aaq=factor*total*surrogate[config.iSO4mm].MM/surrogate[i].MM*Keq/
                  (1.0+Keq)+(1.0-factor)*surrogate[config.iSO4mm].Aaq; //SO4--
	      }
            else if (surrogate[i].name=="NH3") //ammoniac
	      {
                total=surrogate[i].Ag+surrogate[config.iNH4p].Aaq/surrogate[config.iNH4p].MM*surrogate[i].MM; //gas+particle concentrations
                surrogate[i].Ag=factor*total/(1.0+surrogate[i].Kaq_inorg*AQinit)+(1.0-factor)*surrogate[i].Ag; 
                surrogate[config.iNH4p].Aaq=factor*total*surrogate[i].Kaq_inorg*AQinit/(1.0+surrogate[i].Kaq_inorg*AQinit)*surrogate[config.iNH4p].MM/surrogate[i].MM
                  +(1.0-factor)*surrogate[config.iNH4p].Aaq; //NH4+
                deriv_error2_AQ+=total*surrogate[config.iNH4p].MM/surrogate[i].MM*pow(surrogate[i].Kaq_inorg,2)*AQinit/(pow(1.+surrogate[i].Kaq_inorg*AQinit,2))
                  -total*surrogate[config.iNH4p].MM/surrogate[i].MM*surrogate[i].Kaq_inorg/(1.+surrogate[i].Kaq_inorg*AQinit);
	      }
            else if (surrogate[i].name=="HNO3") //nitrate
	      {
                total=surrogate[i].Ag+surrogate[config.iNO3m].Aaq/surrogate[config.iNO3m].MM*surrogate[i].MM; //gas+particle concentrations
                surrogate[i].Ag=factor*total/(1.0+surrogate[i].Kaq_inorg*AQinit)+(1.0-factor)*surrogate[i].Ag;
                surrogate[config.iNO3m].Aaq=factor*total*surrogate[i].Kaq_inorg*AQinit/(1.0+surrogate[i].Kaq_inorg*AQinit)*surrogate[config.iNO3m].MM/surrogate[i].MM
                  +(1.0-factor)*surrogate[config.iNO3m].Aaq; //NO3-
                deriv_error2_AQ+=total*surrogate[config.iNO3m].MM/surrogate[i].MM*pow(surrogate[i].Kaq_inorg,2)*AQinit/(pow(1.+surrogate[i].Kaq_inorg*AQinit,2))
                  -total*surrogate[config.iNO3m].MM/surrogate[i].MM*surrogate[i].Kaq_inorg/(1.+surrogate[i].Kaq_inorg*AQinit);
	      }
            else if (surrogate[i].name=="HCl") //chloride
	      {
                total=surrogate[i].Ag+surrogate[config.iClm].Aaq/surrogate[config.iClm].MM*surrogate[i].MM; //gas+particle concentrations
                surrogate[i].Ag=factor*total/(1.0+surrogate[i].Kaq_inorg*AQinit)+(1.0-factor)*surrogate[i].Ag;
                surrogate[config.iClm].Aaq=factor*total*surrogate[i].Kaq_inorg*AQinit/(1.0+surrogate[i].Kaq_inorg*AQinit)*surrogate[config.iClm].MM/surrogate[i].MM
                  +(1.0-factor)*surrogate[config.iClm].Aaq; //Cl-
                deriv_error2_AQ+=total*surrogate[config.iClm].MM/surrogate[i].MM*pow(surrogate[i].Kaq_inorg,2)*AQinit/(pow(1.+surrogate[i].Kaq_inorg*AQinit,2))
                  -total*surrogate[config.iClm].MM/surrogate[i].MM*surrogate[i].Kaq_inorg/(1.+surrogate[i].Kaq_inorg*AQinit);
	      }
	  }

      for (i=0;i<n;i++)
        if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)
          conc_inorganic+=surrogate[i].Aaq; //compute total inorganic concentrations

      AQ+=conc_inorganic;
      if (config.iH2O>=0)
        {
          //compute hygroscopicity
	  LWC=0.0; //In the case where SOAP compute inorganics, LWC is only used for the initialization
	  if (config.hygroscopicity)
	    hygroscopicity_coupled(config, surrogate, Temperature, MOW, MMaq, RH, LWC, conc_inorganic,
				   MOinit, MO, AQinit, AQ, deriv_error1_MO, 
				   deriv_error1_AQ, deriv_error2_MO, deriv_error2_AQ, factor);
        }
    }
  else
    {
      AQ+=LWC+conc_inorganic;  
      if (config.iH2O>=0 and config.hygroscopicity)        
	if (config.hygroscopicity)
	  hygroscopicity_coupled(config, surrogate, Temperature, MOW, MMaq, RH, LWC, conc_inorganic,
				 MOinit, MO, AQinit, AQ, deriv_error1_MO, 
				 deriv_error1_AQ, deriv_error2_MO, deriv_error2_AQ, factor);        
    }

  AQ=max(AQ,config.MOmin);
  MO=max(MO,config.MOmin);
  
  error1=MOinit-MO;
  error2=AQinit-AQ;
  ionic_organic=ionic_organic_tmp;
  organion=organion_tmp;

  

}

void newton_raphson_coupled(double &MO, double &AQ,
                            double &error1, double &deriv_error1_MO,double &deriv_error1_AQ,
                            double &error2, double &deriv_error2_MO,double &deriv_error2_AQ)
{
  //method of newton raphson used to solved the coupled system
  //    Jacobian = | d(error1)/d(MO) d(error1)/d(AQ) |
  //               | d(error2)/d(MO) d(error2)/d(AQ) |
  //    Xi = | MO | at iteration i
  //         | AQ |
  //
  //    Y = |error1|
  //        |error2|
  //
  //    Xi+1 = Xi + inverse(Jacobian) Y
  //
  //    inverse(Jacobian) = 1/determinant * | d(error2)/d(AQ)  -d(error1)/d(AQ) |
  //                                        | -d(error2)/d(MO) d(error1)/d(MO)  |  
  
  double MO_save=MO;
  double AQ_save=AQ;
  double determinant = deriv_error1_MO*deriv_error2_AQ- deriv_error1_AQ*deriv_error2_MO;
  if (determinant > 0.0)
    {
      double Jacobian_inverse[2][2];
      Jacobian_inverse[0][0]=deriv_error2_AQ/determinant;
      Jacobian_inverse[0][1]=-deriv_error1_AQ/determinant;
      Jacobian_inverse[1][0]=-deriv_error2_MO/determinant;
      Jacobian_inverse[1][1]=deriv_error1_MO/determinant;
      MO-=Jacobian_inverse[0][0]*error1+Jacobian_inverse[0][1]*error2;
      AQ-=Jacobian_inverse[1][0]*error1+Jacobian_inverse[1][1]*error2;
    }
  else
    {
      MO-=error1;
      AQ-=error2;
    }
  if (MO <=0.0 or AQ<=0.0)
    {
      MO=MO_save-error1;
      AQ=AQ_save-error2;
    }
}

void init_saturation2(model_config &config, vector<species>& surrogate,
                      Array<double, 1> &MOinit, 
                      double &AQinit, double LWC, double &conc_inorganic)
{
  //a method to initialize saturation:
  // add a new organic phase with MO = 1 µg/m3 (but keep the compounds in old phases)
  // do not guess in which phase the compound should mainly be present
  AQinit=LWC+conc_inorganic;
  int i,j;
  int n=surrogate.size();
  int nphase=MOinit.size();
  for (j=0;j<nphase-1;++j)
    MOinit(j)=0.0;

  MOinit(nphase-1)=1.0;
  
  for (i=0;i<n;++i)
    {
      AQinit+=surrogate[i].Aaq;
      surrogate[i].Aaq_old=surrogate[i].Aaq;
      surrogate[i].gamma_aq_old=surrogate[i].gamma_aq;
      for (j=0;j<nphase-1;++j)
        {
          surrogate[i].Ap_sat(j)=surrogate[i].Ap_sat_old(j);
          MOinit(j)+=surrogate[i].Ap_sat(j);
        }
      surrogate[i].Ap_sat(nphase-1)=0.0;
    } 
}


void init_saturation(model_config &config, vector<species>& surrogate,
                     Array<double, 1> &MOinit, 
                     double &AQinit, double &LWC, double &conc_inorganic,
                     bool all_hydrophobic)
{
  //a method to initialize saturation
  // guess in which phase the compound should mainly be present
  AQinit=LWC;
  int i,j;
  int n=surrogate.size();
  int nphase=MOinit.size();
  Array<double, 1> mean_activity,stability;
  mean_activity.resize(nphase-1);
  stability.resize(nphase-1);
  
  for (i=0;i<n;++i)
    if (surrogate[i].hydrophilic)
      {
        AQinit+=surrogate[i].Aaq;
        surrogate[i].Aaq_old=surrogate[i].Aaq;
        surrogate[i].gamma_aq_old=surrogate[i].gamma_aq;
      }
  
  int jphase=0; //index of the less stable phase
  
  for (j=0;j<nphase-1;++j)
    {
      double sumX=0.0;
      mean_activity(j)=0.0;
      for (i=0;i<n;++i)
        if ((surrogate[i].is_organic or i==config.iH2O)
            and (surrogate[i].hydrophobic or all_hydrophobic))
          {
            surrogate[i].Xorg_sat_old(j)=surrogate[i].Ap_sat_old(j)/surrogate[i].MM;
            sumX+=surrogate[i].Xorg_sat_old(j);
          }
	  
      for (i=0;i<n;++i)
        if ((surrogate[i].is_organic or i==config.iH2O)
            and (surrogate[i].hydrophobic or all_hydrophobic))
          surrogate[i].Xorg_sat_old(j)/=sumX;

      for (i=0;i<n;++i)
        if ((surrogate[i].is_organic or i==config.iH2O)
            and (surrogate[i].hydrophobic or all_hydrophobic))
          mean_activity(j)+=surrogate[i].Xorg_sat_old(j)*surrogate[i].gamma_org_sat_old(j);

	  
      stability(j)=0.0;
      for (i=0;i<n;++i)
        if ((surrogate[i].is_organic or i==config.iH2O)
            and (surrogate[i].hydrophobic or all_hydrophobic))
          if (surrogate[i].Xorg_sat_old(j)>0.0 and surrogate[i].Ap_sat_old(j)>0.0)
            {
              stability(j)+=surrogate[i].Ap_sat_old(j)*
                log(surrogate[i].Xorg_sat_old(j)*surrogate[i].gamma_org_sat_old(j));
            }
	  
      if (j!=jphase and stability(j)>stability(jphase))
        jphase=j;
    }
  
  double stability_aq=0.0;
  double mean_activity_aq=0.0;
  if (AQinit>0.0)
    {
      double sumX=0.0;
      for (i=0;i<n;++i)
        if (surrogate[i].hydrophilic)
          {
            surrogate[i].Xaq_old=surrogate[i].Aaq_old/surrogate[i].MM;
            sumX+=surrogate[i].Xaq_old;
          }
	  
      for (i=0;i<n;++i)
        if (surrogate[i].hydrophilic)
          surrogate[i].Xaq_old/=sumX;
	  
      for (i=0;i<n;++i)
        if (surrogate[i].hydrophilic)
          mean_activity_aq+=surrogate[i].Xaq_old*surrogate[i].gamma_aq_old*surrogate[i].GAMMAinf;
  
	  
      for (i=0;i<n;++i)
        if (surrogate[i].hydrophilic)
          if (surrogate[i].Xaq_old>0.0 and surrogate[i].Aaq_old>0.0)
            {
              stability_aq+=surrogate[i].Aaq_old*
                log(surrogate[i].Xaq_old*surrogate[i].gamma_aq_old*surrogate[i].GAMMAinf);
            }
    }
  
  if (stability_aq<=stability(jphase) or (AQinit==0.0 or config.coupled_phases==false))
    {      
      for (i=0;i<n;++i)
        if ((surrogate[i].is_organic or i==config.iH2O)
            and (surrogate[i].hydrophobic or all_hydrophobic))
          for (j=0;j<nphase-1;++j)
            if (j==jphase)
              {
                if (surrogate[i].gamma_org_sat_old(jphase) < mean_activity(jphase))
                  {
                    surrogate[i].Ap_sat(jphase)=surrogate[i].Ap_sat_old(jphase);
                    surrogate[i].Ap_sat(nphase-1)=0.0;
                  }
                else
                  {
                    surrogate[i].Ap_sat(jphase)=0.0;
                    surrogate[i].Ap_sat(nphase-1)=surrogate[i].Ap_sat_old(jphase);
                  }
              }
            else
              surrogate[i].Ap_sat(j)=surrogate[i].Ap_sat_old(j);
        else
          for (j=0;j<nphase-1;++j)
            surrogate[i].Ap_sat(j)=0.0;
	  
      for (j=0;j<nphase;++j)
        {
          MOinit(j)=0.0;
          for (i=0;i<n;++i)
            if ((surrogate[i].is_organic or i==config.iH2O)
                and (surrogate[i].hydrophobic or all_hydrophobic))
              MOinit(j)+=surrogate[i].Ap_sat(j);
        }
	  
      if (MOinit(jphase)-surrogate[config.iH2O].Ap_sat(jphase)<1.0e-6) //only water in this phase
        {
          double max_activity=-1.e6;
          int imax=0;
		  
          for (i=0;i<n;++i)
            if ((surrogate[i].is_organic or i==config.iH2O)
                and (surrogate[i].hydrophobic or all_hydrophobic))
              if (surrogate[i].Xorg_sat_old(jphase)>0.0 and (i != config.iH2O) and 
                  surrogate[i].Xorg_sat_old(jphase)*surrogate[i].gamma_org_sat_old(jphase)
                  >max_activity)
                {
                  max_activity=surrogate[i].Xorg_sat_old(jphase)*
                    surrogate[i].gamma_org_sat_old(jphase);
                  imax=i;
                }
		  
          surrogate[imax].Ap_sat(jphase)=surrogate[imax].Ap_sat_old(jphase);
          surrogate[imax].Ap_sat(nphase-1)=0.0;
          for (j=0;j<nphase;++j)
            {
              MOinit(j)=0.0;
              for (i=0;i<n;++i)
                if ((surrogate[i].is_organic or i==config.iH2O)
                    and (surrogate[i].hydrophobic or all_hydrophobic))
                  MOinit(j)+=surrogate[i].Ap_sat(j);
            }
        }
	  
      if (MOinit(nphase-1)-surrogate[config.iH2O].Ap_sat(nphase-1)<1.0e-6) //only water in this phase
        {
          double min_activity=1.e6;
          int imin=0;
          for (i=0;i<n;++i)
            if (surrogate[i].is_organic and (surrogate[i].hydrophobic or all_hydrophobic))
              if (surrogate[i].Xorg_sat_old(jphase)>0.0 and 
                  surrogate[i].Xorg_sat_old(jphase)*surrogate[i].gamma_org_sat_old(jphase)
                  <min_activity)
                {
                  min_activity=surrogate[i].Xorg_sat_old(jphase)*
                    surrogate[i].gamma_org_sat_old(jphase);
                  imin=i;
                }
	  
          surrogate[imin].Ap_sat(nphase-1)=surrogate[imin].Ap_sat_old(jphase);
          surrogate[imin].Ap_sat(jphase)=0.0;
          for (j=0;j<nphase;++j)
            {
              MOinit(j)=0.0;
              for (i=0;i<n;++i)
                if ((surrogate[i].is_organic or i==config.iH2O)
                    and (surrogate[i].hydrophobic or all_hydrophobic))
                  MOinit(j)+=surrogate[i].Ap_sat(j);
            }
        }
    }
  else
    {      
      AQinit=LWC;      
      for (i=0;i<n;++i)
        {
          if (surrogate[i].hydrophobic and surrogate[i].hydrophilic)
            {
              if (surrogate[i].is_organic)
                {
                  if (surrogate[i].gamma_aq_old*surrogate[i].GAMMAinf < mean_activity_aq)
                    {
                      surrogate[i].Aaq=surrogate[i].Aaq_old;
                      surrogate[i].Ap_sat(nphase-1)=0.0;
                    }
                  else
                    {
                      surrogate[i].Aaq=0.0;
                      surrogate[i].Ap_sat(nphase-1)=surrogate[i].Aaq_old;
                    }
                }
              else
                {
                  surrogate[i].Aaq=surrogate[i].Aaq_old;
                  surrogate[i].Ap_sat(nphase-1)=0.0;
                }
		  
              for (j=0;j<nphase-1;++j)
                surrogate[i].Ap_sat(j)=surrogate[i].Ap_sat_old(j);
            }
          else
            {
              for (j=0;j<nphase-1;++j)
                surrogate[i].Ap_sat(j)=surrogate[i].Ap_sat_old(j);
              surrogate[i].Ap_sat(nphase-1)=0.0;
            }
		  
          if (surrogate[i].hydrophilic)
            AQinit+=surrogate[i].Aaq;
        }
      for (j=0;j<nphase;++j)
        {
          MOinit(j)=0.0;
          for (i=0;i<n;++i)
            if ((surrogate[i].is_organic or i==config.iH2O)
                and (surrogate[i].hydrophobic or all_hydrophobic))
              MOinit(j)+=surrogate[i].Ap_sat(j);
        }
    }

  for (i=0;i<n;++i)
    {
      double maximum=0.0;
      if(nphase>2)
        surrogate[i].jmain_phase_old=surrogate[i].jmain_phase;
      else
        surrogate[i].jmain_phase_old=0;
	  
      surrogate[i].jmain_phase=0;
      if (surrogate[i].hydrophobic and surrogate[i].is_organic)
        for (j=0;j<nphase;++j)
          {
            if (surrogate[i].Ap_sat(j)>maximum)
              {
                maximum=surrogate[i].Ap_sat(j);
                surrogate[i].jmain_phase=j;
              }
          }
    }
  
  for(i=0;i<n;++i)
    if((surrogate[i].is_organic==false and i!=config.iH2O)or
       (surrogate[i].hydrophobic==false and all_hydrophobic==false))
      for (j=0;j<nphase;++j)
        surrogate[i].Ap_sat(j)=0.0;

  for(i=0;i<n;++i)
    if (surrogate[i].hydrophilic==false or all_hydrophobic and surrogate[i].is_organic)
      surrogate[i].Aaq=0.0;
  
}

void error_saturation(model_config &config, vector<species>& surrogate,
                      Array <double, 1> &MOinit, Array <double, 1> &MOW, double &MMaq,
                      double &AQinit, double &LWC, double &conc_inorganic,
                      double &ionic, double &ionic_organic, double &chp, double &organion, 
                      double &Temperature, double &RH,
                      Array <double, 1> &error, Array <double, 2> &Jacobian, double factor,
                      bool compute_activity_coefficients)
{
  int n=surrogate.size();
  int nphase=MOinit.size();
  int i,j,k;
  double AQ;
  Array <double, 1> MO;
  MO.resize(nphase);
  double Kp;
  double XH2O;

  //initialize AQ and MO
  AQ=0.0;
  for (i=0;i<nphase;++i)
    {
      MO(i)=0.0;  
      MOinit(i)=max(MOinit(i),config.MOmin);
    }

  AQinit=max(AQinit,config.MOmin);

  double conc_org=LWC;
  for (i=0;i<n;++i)
    if (surrogate[i].is_organic or i==config.iH2O)     
      conc_org+=surrogate[i].Aaq;
      
  conc_org=max(conc_org,config.MOmin);

  compute_ionic_strenght2(config, surrogate, AQinit, conc_inorganic, ionic, chp, organion, ionic_organic, factor);  
  
  for (i=0;i<nphase+1;++i)
    {
      for (j=0;j<nphase+1;++j)
        Jacobian(i,j)=0.0;
      Jacobian(i,i)=1.0;
    }

  //compute acitivity coefficients, MOW and MMaq
  if (compute_activity_coefficients)
    {
      activity_coefficients_saturation(config, surrogate, false, Temperature, MOW);
      activity_coefficients_aq(config, surrogate, Temperature, LWC, MMaq, XH2O);
      if (config.compute_long_and_medium_range_interactions)
        activity_coefficients_LR_MR(config, surrogate, Temperature, LWC, ionic);
    }

  //compute pH
  if (config.compute_aqueous_phase_properties)
    {
      //If inorganic ion concentrations are computed by SOAP, used a method of
      //newton raphson
      if (config.compute_inorganic)
        {
          double error_h=1000.0;
          double derivative_h;
          double chp2=chp;
          while(abs(error_h/chp2)>1.0e-3)
            {
              Kpreal_inorganic(config, surrogate, chp2, MMaq);
              error_ph(config, surrogate, Temperature, chp2, organion, error_h, derivative_h,AQinit,LWC);
              if (chp2-error_h/derivative_h>0.0 and derivative_h!=0.0)
		chp2=chp2-error_h/derivative_h;
	      else
                chp2=chp2+error_h;    
            }
          chp=factor*chp2+(1.0-factor)*chp;
        }
      else
        {
          double Ke=1.0e-14;
          double inorganion=0.0;
          for (i=0;i<n;++i)
            if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)
              if (surrogate[i].name!="H")
                inorganion-=surrogate[i].molality*surrogate[i].charge;
      
          chp=factor*0.5*(organion+inorganion+pow(pow(organion+inorganion,2)+4*Ke,0.5))+(1.0-factor)*chp;
        }
    }

  double organion_tmp=0.0;
  double fion1,fion2,molality1,molality2;
  double ionic_organic_tmp=0.0;
  double Kp_aq;
  Array<double, 1> Kp_org;
  Kp_org.resize(nphase);
  double sum;
  
  for (i=0;i<n;++i)
    if (surrogate[i].hydrophobic and surrogate[i].is_organic
        and surrogate[i].hydrophilic==false)
      //compute absorption for compounds which are only hydrophobic
      if (surrogate[i].nonvolatile)
        {
          sum=0.0;
          for (j=0;j<nphase;++j)
            sum+=MOinit(j)/surrogate[i].gamma_org_sat(j);
		  
          surrogate[i].Aaq=0.0;
          for (j=0;j<nphase;++j)
            {
              surrogate[i].Ap_sat(j)=factor*surrogate[i].Atot*
                (MOinit(j)/surrogate[i].gamma_org_sat(j))/sum+(1.0-factor)*surrogate[i].Ap_sat(j);
              MO(j)+=surrogate[i].Ap_sat(j);
              for (k=0;k<nphase;++k)
                if (j==k)
                  Jacobian(j,k)-=surrogate[i].Atot/(sum*surrogate[i].gamma_org_sat(j))
                    -surrogate[i].Atot*MOinit(j)/(pow(sum*surrogate[i].gamma_org_sat(j),2));
                else
                  Jacobian(j,k)+=surrogate[i].Atot*MOinit(j)/
                    (pow(sum,2)*surrogate[i].gamma_org_sat(j)*surrogate[i].gamma_org_sat(k));	  
            }
        }
      else
        if (surrogate[i].kp_from_experiment)
          {
            Kp=surrogate[i].kpi;
            for (j=0;j<nphase;++j)
              if (j!=surrogate[i].jmain_phase)
                surrogate[i].Ap_sat(j)=0.0;
			
            surrogate[i].Ap_sat(surrogate[i].jmain_phase)=factor*surrogate[i].Atot*Kp*
              MOinit(surrogate[i].jmain_phase)/(1+Kp*MOinit(surrogate[i].jmain_phase))+
              (1.0-factor)*surrogate[i].Ap_sat(surrogate[i].jmain_phase);
			
            surrogate[i].Aaq=0.0;
            MO(surrogate[i].jmain_phase)+=surrogate[i].Ap_sat(surrogate[i].jmain_phase);
            Jacobian(surrogate[i].jmain_phase,surrogate[i].jmain_phase)+=
              surrogate[i].Atot*pow(Kp,2)*MOinit(surrogate[i].jmain_phase)
              /(pow(1+Kp*MOinit(surrogate[i].jmain_phase),2))
              -surrogate[i].Atot*Kp/(1+Kp*MOinit(surrogate[i].jmain_phase));
          }
        else
          {
            for (j=0;j<nphase;++j)
              Kp_org(j)=surrogate[i].kpi/MOW(j)/surrogate[i].gamma_org_sat(j);
            surrogate[i].Aaq=0.0;
            sum=1.0;
            for (j=0;j<nphase;++j)
              sum+=Kp_org(j)*MOinit(j);
			
            for (j=0;j<nphase;++j)
              {
                surrogate[i].Ap_sat(j)=factor*surrogate[i].Atot*Kp_org(j)*MOinit(j)/sum+
                  (1.0-factor)*surrogate[i].Ap_sat(j);
                MO(j)+=surrogate[i].Ap_sat(j);
                for (k=0;k<nphase;++k)
                  if (j==k)
                    Jacobian(j,k)-=surrogate[i].Atot*Kp_org(j)/sum
                      -surrogate[i].Atot*pow(Kp_org(j),2)*MOinit(j)/pow(sum,2);
                  else
                    Jacobian(j,k)+=surrogate[i].Atot*
                      Kp_org(j)*MOinit(j)*Kp_org(k)/pow(sum,2);
              }
          }
    else if (surrogate[i].hydrophobic==false and surrogate[i].is_organic
             and surrogate[i].hydrophilic)
      //compute absorption for compounds which are only hydrophilic
      if (surrogate[i].nonvolatile)
        {
          surrogate[i].Aaq=surrogate[i].Atot;
          AQ+=surrogate[i].Aaq;
          for (j=0;j<nphase;++j)
            surrogate[i].Ap_sat(j)=0.0;
        }

      else
        {
          fion1=0.0;
          fion2=0.0;
          Kp=surrogate[i].Kp_eff_aqreal(config, Temperature, ionic, chp,surrogate[config.iHp].gamma_LR,
                                        surrogate[config.iHp].gamma_SRMR,MMaq, fion1, fion2)
            /surrogate[i].gamma_aq;
          surrogate[i].Aaq=factor*surrogate[i].Atot*Kp*AQinit/(1+Kp*AQinit)
            +(1.0-factor)*surrogate[i].Aaq;
          for (j=0;j<nphase;++j)
            surrogate[i].Ap_sat(j)=0.0;
          AQ+=surrogate[i].Aaq;
          Jacobian(nphase,nphase)+=
            surrogate[i].Atot*pow(Kp,2)*AQinit/(pow(1+Kp*AQinit,2))
            -surrogate[i].Atot*Kp/(1+Kp*AQinit);
          //molality1: molality of ions HA- or A-  
          molality1=surrogate[i].Aaq*fion1/surrogate[i].MM/conc_org*1000.0;
          //molality2: molality of ions A2-
          molality2=surrogate[i].Aaq*fion2/surrogate[i].MM/conc_org*1000.0;
          //compute ionic_organic and organion
          ionic_organic_tmp+=0.5*molality1+0.5*molality2*4;
          organion_tmp+=molality1+2*molality2;
        }
    else if (surrogate[i].hydrophobic and surrogate[i].is_organic
             and surrogate[i].hydrophilic)
      //compute absorption for compounds which are both hydrophobic and hydrophilic
      if (surrogate[i].nonvolatile)
        {
          sum=AQinit/surrogate[i].gamma_aq;
          for (j=0;j<nphase;++j)
            sum+=MOinit(j)/surrogate[i].gamma_org_sat(j);
		  
          surrogate[i].Aaq=factor*surrogate[i].Atot*AQinit/(surrogate[i].gamma_aq*sum)
            +(1.0-factor)*surrogate[i].Aaq;
          AQ+=surrogate[i].Aaq;
		  
          for (j=0;j<nphase;++j)
            {
              surrogate[i].Ap_sat(j)=(1.0-factor)*surrogate[i].Ap_sat(j)+factor*surrogate[i].Atot*
                (MOinit(j)/surrogate[i].gamma_org_sat(j))/sum;
              MO(j)+=surrogate[i].Ap_sat(j);
              for (k=0;k<nphase;++k)
                if (j==k)
                  Jacobian(j,k)-=surrogate[i].Atot/(sum*surrogate[i].gamma_org_sat(j))
                    -surrogate[i].Atot*MOinit(j)/(pow(sum*surrogate[i].gamma_org_sat(j),2));
                else
                  Jacobian(j,k)+=surrogate[i].Atot*MOinit(j)/
                    (pow(sum,2)*surrogate[i].gamma_org_sat(j)*surrogate[i].gamma_org_sat(k));
              Jacobian(j,nphase)+=surrogate[i].Atot*MOinit(j)/
                (pow(sum,2)*surrogate[i].gamma_org_sat(j)*surrogate[i].gamma_aq);
            }
		  
          for (k=0;k<nphase;++k)
            Jacobian(nphase,k)+=surrogate[i].Atot*AQinit/
              (surrogate[i].gamma_aq*pow(sum,2)*surrogate[i].gamma_org_sat(k));

          Jacobian(nphase,nphase)-=surrogate[i].Atot/(sum*surrogate[i].gamma_aq)
            -surrogate[i].Atot*AQinit/(pow(sum*surrogate[i].gamma_aq,2));
        }
      else
        {
          fion1=0.0;
          fion2=0.0;
          Kp_aq=surrogate[i].Kp_eff_aqreal(config, Temperature, ionic, chp,
                                           surrogate[config.iHp].gamma_LR,
                                           surrogate[config.iHp].gamma_SRMR,
                                           MMaq, fion1, fion2)
            /surrogate[i].gamma_aq;

          for (j=0;j<nphase;++j)
            {
              if (surrogate[i].kp_from_experiment)
                if (j==surrogate[i].jmain_phase)
                  Kp_org(j)=surrogate[i].kpi;
                else
                  Kp_org(j)=0.0;
              else
                Kp_org(j)=surrogate[i].kpi/MOW(j)/surrogate[i].gamma_org_sat(j);
            }

          sum=1.0+Kp_aq*AQinit;
          for (j=0;j<nphase;++j)
            sum+=Kp_org(j)*MOinit(j);
		  
          surrogate[i].Aaq=factor*surrogate[i].Atot*Kp_aq*AQinit/sum
            +(1.0-factor)*surrogate[i].Aaq;
          AQ+=surrogate[i].Aaq;
		  
          for (j=0;j<nphase;++j)
            {
              surrogate[i].Ap_sat(j)=factor*surrogate[i].Atot*Kp_org(j)*MOinit(j)/sum
                +(1.0-factor)*surrogate[i].Ap_sat(j);
              MO(j)+=surrogate[i].Ap_sat(j);
              for (k=0;k<nphase;++k)
                if (j==k)
                  Jacobian(j,k)-=surrogate[i].Atot*Kp_org(j)/sum
                    -surrogate[i].Atot*pow(Kp_org(j),2)*MOinit(j)/pow(sum,2);
                else
                  Jacobian(j,k)+=surrogate[i].Atot*
                    Kp_org(j)*MOinit(j)*Kp_org(k)/pow(sum,2);
              Jacobian(j,nphase)+=surrogate[i].Atot*
                Kp_org(j)*MOinit(j)*Kp_aq/pow(sum,2);
            }
		  
          for (k=0;k<nphase;++k)
            Jacobian(nphase,k)+=surrogate[i].Atot*AQinit*Kp_aq*Kp_org(k)/pow(sum,2);

          Jacobian(nphase,nphase)-=surrogate[i].Atot*Kp_aq/sum
            -surrogate[i].Atot*pow(Kp_aq,2)*AQinit/(pow(sum,2));
		  
          //molality1: molality of ions HA- or A-
          molality1=surrogate[i].Aaq*fion1/surrogate[i].MM/conc_org*1000.0;
          //molality2: molality of ions A2-
          molality2=surrogate[i].Aaq*fion2/surrogate[i].MM/conc_org*1000.0;
          //compute ionic_organic and organion
          ionic_organic_tmp+=0.5*molality1+0.5*molality2*4;
          organion_tmp+=molality1+2*molality2;
        }
  
  //compute inorganic ion concentrations
  if (config.compute_inorganic)
    {
      double total;
      conc_inorganic=0.0;
      Kpreal_inorganic(config, surrogate, chp, MMaq);
      for (i=0;i<n;i++)
        if (surrogate[i].is_inorganic_precursor)
          if (surrogate[i].name=="H2SO4") //sulfate 
            {
              total=surrogate[i].Ag+surrogate[config.iHSO4m].Aaq/surrogate[config.iHSO4m].MM*surrogate[i].MM
                +surrogate[config.iSO4mm].Aaq/surrogate[config.iSO4mm].MM*surrogate[i].MM;
              double Keq=surrogate[i].keq/chp*surrogate[config.iHSO4m].gamma_aq/surrogate[config.iHp].gamma_aq/surrogate[config.iSO4mm].gamma_aq;
              surrogate[i].Ag=(1.0-factor)*surrogate[i].Ag;
              surrogate[config.iHSO4m].Aaq=factor*total*surrogate[config.iHSO4m].MM/surrogate[i].MM/(1.0+Keq)
                +(1.0-factor)*surrogate[config.iHSO4m].Aaq;
              surrogate[config.iSO4mm].Aaq=factor*total*surrogate[config.iSO4mm].MM/surrogate[i].MM*Keq/
                (1.0+Keq)+(1.0-factor)*surrogate[config.iSO4mm].Aaq;
            }
          else if (surrogate[i].name=="NH3") //ammonium
            {
              total=surrogate[i].Ag+surrogate[config.iNH4p].Aaq/surrogate[config.iNH4p].MM*surrogate[i].MM;
              surrogate[i].Ag=factor*total/(1.0+surrogate[i].Kaq_inorg*AQinit)+(1.0-factor)*surrogate[i].Ag;
              surrogate[config.iNH4p].Aaq=factor*total*surrogate[i].Kaq_inorg*AQinit/(1.0+surrogate[i].Kaq_inorg*AQinit)*surrogate[config.iNH4p].MM/surrogate[i].MM
                +(1.0-factor)*surrogate[config.iNH4p].Aaq;
              Jacobian(nphase,nphase)+=total*surrogate[config.iNH4p].MM/surrogate[i].MM*pow(surrogate[i].Kaq_inorg,2)*AQinit/(pow(1.+surrogate[i].Kaq_inorg*AQinit,2))
                -total*surrogate[config.iNH4p].MM/surrogate[i].MM*surrogate[i].Kaq_inorg/(1.+surrogate[i].Kaq_inorg*AQinit);
            }
          else if (surrogate[i].name=="HNO3") //nitrate
            {
              total=surrogate[i].Ag+surrogate[config.iNO3m].Aaq/surrogate[config.iNO3m].MM*surrogate[i].MM;
              surrogate[i].Ag=factor*total/(1.0+surrogate[i].Kaq_inorg*AQinit)+(1.0-factor)*surrogate[i].Ag;
              surrogate[config.iNO3m].Aaq=factor*total*surrogate[i].Kaq_inorg*AQinit/(1.0+surrogate[i].Kaq_inorg*AQinit)*surrogate[config.iNO3m].MM/surrogate[i].MM
                +(1.0-factor)*surrogate[config.iNO3m].Aaq;
              Jacobian(nphase,nphase)+=total*surrogate[config.iNO3m].MM/surrogate[i].MM*pow(surrogate[i].Kaq_inorg,2)*AQinit/(pow(1.+surrogate[i].Kaq_inorg*AQinit,2))
                -total*surrogate[config.iNO3m].MM/surrogate[i].MM*surrogate[i].Kaq_inorg/(1.+surrogate[i].Kaq_inorg*AQinit);
            }
      
          else if (surrogate[i].name=="HCl") //chloride
            {
              total=surrogate[i].Ag+surrogate[config.iClm].Aaq/surrogate[config.iClm].MM*surrogate[i].MM;
              surrogate[i].Ag=factor*total/(1.0+surrogate[i].Kaq_inorg*AQinit)+(1.0-factor)*surrogate[i].Ag;
              surrogate[config.iClm].Aaq=factor*total*surrogate[i].Kaq_inorg*AQinit/(1.0+surrogate[i].Kaq_inorg*AQinit)*surrogate[config.iClm].MM/surrogate[i].MM
                +(1.0-factor)*surrogate[config.iClm].Aaq;
              Jacobian(nphase,nphase)+=total*surrogate[config.iClm].MM/surrogate[i].MM*pow(surrogate[i].Kaq_inorg,2)*AQinit/(pow(1.+surrogate[i].Kaq_inorg*AQinit,2))
                -total*surrogate[config.iClm].MM/surrogate[i].MM*surrogate[i].Kaq_inorg/(1.+surrogate[i].Kaq_inorg*AQinit);
            }

      for (i=0;i<n;i++)
        if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)
          conc_inorganic+=surrogate[i].Aaq; //compute total inorganic concentrations
  
      AQ+=conc_inorganic;
      if (config.iH2O>=0)
        {
	  LWC=0.0;
	  
	  if (surrogate[config.iH2O].hydrophilic)
	    hygroscopicity_org_sat(config, surrogate, Temperature, MOW, MMaq, RH, LWC, conc_inorganic, MOinit, MO, AQinit, AQ, Jacobian,
				   true, factor);
	  else
	    hygroscopicity_org_sat(config, surrogate, Temperature, MOW, MMaq, RH, LWC, conc_inorganic, MOinit, MO, AQinit, AQ, Jacobian,
				   false, factor);	      
        }
      
    }
  else
    {	
      AQ+=LWC+conc_inorganic;  
      if (config.iH2O>=0 and config.hygroscopicity)
	if (surrogate[config.iH2O].hydrophilic)
	  hygroscopicity_org_sat(config, surrogate, Temperature, MOW, MMaq, RH, LWC, conc_inorganic, MOinit, MO, AQinit, AQ, Jacobian,
				 true, factor);
	else
	  hygroscopicity_org_sat(config, surrogate, Temperature, MOW, MMaq, RH, LWC, conc_inorganic, MOinit, MO, AQinit, AQ, Jacobian,
				 false, factor);	       
    }
  
  for (i=0;i<n;++i)
    if (surrogate[i].is_organic==false and i!=config.iH2O)
      for (j=0;j<nphase;++j)
        surrogate[i].Ap_sat(j)=0.0; //non-organic ions are not present in the organic phase
  
  for (j=0;j<nphase;++j)
    {
      MO(j)=max(MO(j),config.MOmin);
      error(j)=MOinit(j)-MO(j);
    }

  AQ=max(AQ,config.MOmin);
  error(nphase)=AQinit-AQ;
  ionic_organic=ionic_organic_tmp;
  organion=organion_tmp;
}

void error_saturation_hydrophobic(model_config &config, vector<species>& surrogate,
                                  Array <double, 1> &MOinit, Array <double, 1> &MOW, 
                                  double &Temperature, double &RH,
                                  Array <double, 1> &error, Array <double, 2> &Jacobian,
                                  double factor, bool compute_activity_coefficients)
{
  int n=surrogate.size();
  int nphase=MOinit.size();
  int i,j,k;
  Array <double, 1> MO;
  MO.resize(nphase);
  double Kp;

  //initialize AQ and MO
  for (i=0;i<nphase;++i)
    {
      MO(i)=0.0;
      MOinit(i)=max(MOinit(i),config.MOmin);
    }
  
  for (i=0;i<nphase+1;++i)
    {
      for (j=0;j<nphase+1;++j)
        Jacobian(i,j)=0.0;
      Jacobian(i,i)=1.0;
    }
  
  //compute acitivity coefficients and MOW
  if (compute_activity_coefficients)
    activity_coefficients_saturation(config, surrogate, true, Temperature, MOW);
  
  Array<double, 1> Kp_org;
  Kp_org.resize(nphase);
  double sum;
  
  for (i=0;i<n;++i)
    if (surrogate[i].is_organic)
      //compute absorption for compounds which are only hydrophobic
      if (surrogate[i].nonvolatile)
        {
          sum=0.0;
          for (j=0;j<nphase;++j)
            sum+=MOinit(j)/surrogate[i].gamma_org_sat(j);
		  
          surrogate[i].Aaq=0.0;
          for (j=0;j<nphase;++j)
            {
              surrogate[i].Ap_sat(j)=factor*surrogate[i].Atot*
                (MOinit(j)/surrogate[i].gamma_org_sat(j))/sum+(1.0-factor)*surrogate[i].Ap_sat(j);
              MO(j)+=surrogate[i].Ap_sat(j);
              for (k=0;k<nphase;++k)
                if (j==k)
                  Jacobian(j,k)-=surrogate[i].Atot/(sum*surrogate[i].gamma_org_sat(j))
                    -surrogate[i].Atot*MOinit(j)/(pow(sum*surrogate[i].gamma_org_sat(j),2));
                else
                  Jacobian(j,k)+=surrogate[i].Atot*MOinit(j)/
                    (pow(sum,2)*surrogate[i].gamma_org_sat(j)*surrogate[i].gamma_org_sat(k));	  
            }
        }
      else
        if (surrogate[i].kp_from_experiment)
          {
            Kp=surrogate[i].kpi;
            for (j=0;j<nphase;++j)
              if(j!=surrogate[i].jmain_phase)
                surrogate[i].Ap_sat(j)=0.0;
            surrogate[i].Ap_sat(surrogate[i].jmain_phase)=factor*surrogate[i].Atot*Kp*
              MOinit(surrogate[i].jmain_phase)/(1+Kp*MOinit(surrogate[i].jmain_phase))
              +(1.0-factor)*surrogate[i].Ap_sat(surrogate[i].jmain_phase);
            surrogate[i].Aaq=0.0;
            MO(surrogate[i].jmain_phase)+=surrogate[i].Ap_sat(surrogate[i].jmain_phase);
            Jacobian(surrogate[i].jmain_phase,surrogate[i].jmain_phase)+=
              surrogate[i].Atot*pow(Kp,2)*MOinit(surrogate[i].jmain_phase)
              /(pow(1+Kp*MOinit(surrogate[i].jmain_phase),2))
              -surrogate[i].Atot*Kp/(1+Kp*MOinit(surrogate[i].jmain_phase));
          }
        else
          {
            for (j=0;j<nphase;++j)
              Kp_org(j)=surrogate[i].kpi/surrogate[i].gamma_org_sat(j);

            surrogate[i].Aaq=0.0;
            sum=1.0;
            for (j=0;j<nphase;++j)
              sum+=Kp_org(j)*MOinit(j);
			
            for (j=0;j<nphase;++j)
              {
                surrogate[i].Ap_sat(j)=factor*surrogate[i].Atot*Kp_org(j)*MOinit(j)/sum
                  +(1.0-factor)*surrogate[i].Ap_sat(j);
                MO(j)+=surrogate[i].Ap_sat(j);
                for (k=0;k<nphase;++k)
                  if (j==k)
                    Jacobian(j,k)+=surrogate[i].Atot*pow(Kp_org(j),2)*MOinit(j)/pow(sum,2)
                      -surrogate[i].Atot*Kp_org(j)/sum;
                  else
                    Jacobian(j,k)+=surrogate[i].Atot*
                      Kp_org(j)*MOinit(j)*Kp_org(k)/pow(sum,2);
              }
          }
  
  for (i=0;i<n;++i)
    if (surrogate[i].is_organic==false and i!=config.iH2O)
      for (j=0;j<nphase;++j)
        surrogate[i].Ap_sat(j)=0.0;  //non-organic ions are not present in the organic phase
  
  if (config.iH2O>=0 and config.hygroscopicity)
    if (surrogate[config.iH2O].hydrophobic) //Can H2O condense on the organic phase
      if (config.iH2O>=0 and config.hygroscopicity)
	{
	  double AQinit=0.0;
	  double AQ=0.0;
	  double MMaq=18.0;
	  hygroscopicity_org_sat(config, surrogate, Temperature, MOW, MMaq, RH, 0.0, 0.0, MOinit, MO, AQinit, AQ, Jacobian,
				 false, factor);	  
	}

  for (j=0;j<nphase;++j)
    {
      MO(j)=max(MO(j),config.MOmin);
      error(j)=MOinit(j)-MO(j);      
    }
  error(nphase)=0.0;
  
}

void stability(model_config &config, vector<species>& surrogate, double &Temperature,
               bool &stable_system, bool all_hydrophobic, Array <double, 1> MO)
{
  //This routine compute the gibbs energy before and after phase separation
  //If gibbs energy is lower before phase separation, the system is stable and phase
  // separation does not occur 
  double stability_param=0.0;
  double stability_param_old=0.0;
  int nphase=surrogate[0].Ap_sat.size();
  int n=surrogate.size();
  int i,j;
  double sum,sumMo;
  bool low_concentrations;

  sumMo=0.0;
  for (j=0;j<nphase;++j)
    sumMo+=MO(j);

  low_concentrations=false;
  if (sumMo>0.0)
    for (j=0;j<nphase;++j)
      if(MO(j)/sumMo<0.01)
        low_concentrations=true;
  
  for (i=0;i<n;++i)
    if (surrogate[i].is_organic)
      if(surrogate[i].nonvolatile==false)
        {
          surrogate[i].Ag=surrogate[i].Atot;
          surrogate[i].Ag_old=surrogate[i].Atot;
          if (surrogate[i].hydrophilic and all_hydrophobic==false)
            {
              surrogate[i].Ag-=surrogate[i].Aaq;
              surrogate[i].Ag_old-=surrogate[i].Aaq_old;
            }

          if (surrogate[i].hydrophobic or all_hydrophobic)
            for (j=0;j<nphase;++j)
              surrogate[i].Ag-=surrogate[i].Ap_sat(j);
		  
          if (surrogate[i].hydrophobic or all_hydrophobic)
            for (j=0;j<nphase-1;++j)
              surrogate[i].Ag_old-=surrogate[i].Ap_sat_old(j);
		  
          surrogate[i].partial_pressure=surrogate[i].Ag/surrogate[i].MM*1.e-6*8.314*Temperature;
          surrogate[i].partial_pressure_old=surrogate[i].Ag_old/
            surrogate[i].MM*1.e-6*8.314*Temperature;
        }
      else
        {
          surrogate[i].Ag=0.0;
          surrogate[i].partial_pressure=0.0;
          surrogate[i].Ag_old=0.0;
          surrogate[i].partial_pressure_old=0.0;
        }
    else if (surrogate[i].is_inorganic_precursor)
      if(surrogate[i].nonvolatile==false)
        {	  
          surrogate[i].partial_pressure=surrogate[i].Ag/surrogate[i].MM*1.e-6*8.314*Temperature;
          surrogate[i].partial_pressure_old=surrogate[i].Ag_old/
            surrogate[i].MM*1.e-6*8.314*Temperature;
        }
      else
        {
          surrogate[i].Ag=0.0;
          surrogate[i].partial_pressure=0.0;
          surrogate[i].Ag_old=0.0;
          surrogate[i].partial_pressure_old=0.0;
        }
    else
      {
        surrogate[i].Ag=0.0;
        surrogate[i].partial_pressure=0.0;
        surrogate[i].Ag_old=0.0;
        surrogate[i].partial_pressure_old=0.0;
      }

  for (j=0;j<nphase;++j)
    {
      sum=0.0;
      for (i=0;i<n;++i)
        if (surrogate[i].hydrophobic or all_hydrophobic)
          {
            surrogate[i].Xorg_sat(j)=surrogate[i].Ap_sat(j)/surrogate[i].MM;
            sum+=surrogate[i].Xorg_sat(j);
          }
      for (i=0;i<n;++i)
        if (surrogate[i].hydrophobic or all_hydrophobic)
          surrogate[i].Xorg_sat(j)/=sum;
    }

  for (j=0;j<nphase-1;++j)
    {
      sum=0.0;
      for (i=0;i<n;++i)
        if (surrogate[i].hydrophobic or all_hydrophobic)
          {
            surrogate[i].Xorg_sat_old(j)=surrogate[i].Ap_sat_old(j)/surrogate[i].MM;
            sum+=surrogate[i].Xorg_sat_old(j);
          }
      for (i=0;i<n;++i)
        if (surrogate[i].hydrophobic or all_hydrophobic)
          surrogate[i].Xorg_sat_old(j)/=sum;
    }
  
  sum=0.0;
  for (i=0;i<n;++i)
    if (surrogate[i].hydrophilic and all_hydrophobic==false)
      {
        surrogate[i].Xaq=surrogate[i].Aaq/surrogate[i].MM;
        sum+=surrogate[i].Xaq;
      }

  if (sum>0.0)
    for (i=0;i<n;++i)
      if (surrogate[i].hydrophilic and all_hydrophobic==false)
        surrogate[i].Xaq/=sum;

  sum=0.0;
  for (i=0;i<n;++i)
    if (surrogate[i].hydrophilic and all_hydrophobic==false)
      {
        surrogate[i].Xaq_old=surrogate[i].Aaq_old/surrogate[i].MM;
        sum+=surrogate[i].Xaq_old;
      }

  if (sum>0.0)
    for (i=0;i<n;++i)
      if (surrogate[i].hydrophilic and all_hydrophobic==false)
        surrogate[i].Xaq_old/=sum;
  
  for (i=0;i<n;++i)
    {
      if (surrogate[i].is_organic)
        if (surrogate[i].partial_pressure>0.0)
          if (surrogate[i].kp_from_experiment)
            stability_param+=surrogate[i].Ag/surrogate[i].MM*
              log(surrogate[i].Xorg_sat(surrogate[i].jmain_phase)*
                  surrogate[i].gamma_org_sat(surrogate[i].jmain_phase));
          else
            stability_param+=surrogate[i].Ag/surrogate[i].MM*
              log(surrogate[i].partial_pressure/surrogate[i].Psat(Temperature));

      if (surrogate[i].hydrophilic and all_hydrophobic==false)
        if (surrogate[i].Xaq*surrogate[i].gamma_aq*surrogate[i].GAMMAinf>0.0)
          stability_param+=surrogate[i].Aaq*
            log(surrogate[i].Xaq*surrogate[i].gamma_aq*surrogate[i].GAMMAinf)/
            surrogate[i].MM;

      if ((surrogate[i].hydrophobic or all_hydrophobic) and
          (surrogate[i].is_organic or i==config.iH2O))
        for (j=0;j<nphase;++j)
          if (surrogate[i].Xorg_sat(j)*surrogate[i].gamma_org_sat(j)>0.0)
            stability_param+=surrogate[i].Ap_sat(j)/surrogate[i].MM*
              log(surrogate[i].Xorg_sat(j)*surrogate[i].gamma_org_sat(j));
	  
      if (surrogate[i].is_organic)
        if (surrogate[i].partial_pressure_old>0.0)
          if (surrogate[i].kp_from_experiment)
            stability_param_old+=surrogate[i].Ag_old/surrogate[i].MM*
              log(surrogate[i].Xorg_sat_old(surrogate[i].jmain_phase_old)*
                  surrogate[i].gamma_org_sat_old(surrogate[i].jmain_phase_old));
          else
            stability_param_old+=surrogate[i].Ag_old/surrogate[i].MM*
              log(surrogate[i].partial_pressure_old/surrogate[i].Psat(Temperature));

      if (surrogate[i].hydrophilic and all_hydrophobic==false)
        if (surrogate[i].Xaq_old*surrogate[i].gamma_aq_old*surrogate[i].GAMMAinf>0.0)
          stability_param_old+=surrogate[i].Aaq_old/surrogate[i].MM
            *log(surrogate[i].Xaq_old*surrogate[i].gamma_aq_old*surrogate[i].GAMMAinf);
	  
      if ((surrogate[i].hydrophobic or all_hydrophobic) and
          (surrogate[i].is_organic or i==config.iH2O))
        for (j=0;j<nphase-1;++j)
          if (surrogate[i].Xorg_sat_old(j)*surrogate[i].gamma_org_sat_old(j)>0.0)
            stability_param_old+=surrogate[i].Ap_sat_old(j)/surrogate[i].MM*
              log(surrogate[i].Xorg_sat_old(j)*surrogate[i].gamma_org_sat_old(j));
    }

  if ((stability_param_old<=stability_param
       or abs(stability_param_old-stability_param)<0.0001*abs(stability_param_old))
      or low_concentrations)
    {
      stable_system=true;
      nphase=nphase-1;
      for (i=0;i<n;++i)
        {
          surrogate[i].Ag=surrogate[i].Ag_old;
          surrogate[i].Aaq=surrogate[i].Aaq_old;
          surrogate[i].Ap_sat.resize(nphase);
          surrogate[i].gamma_org_sat_old.resize(nphase);
          surrogate[i].gamma_aq=surrogate[i].gamma_aq_old;
          for (j=0;j<nphase;++j)
            if ((surrogate[i].hydrophobic or all_hydrophobic) and
                (surrogate[i].is_organic or i==config.iH2O))
              {
                surrogate[i].gamma_org_sat(j)=surrogate[i].gamma_org_sat_old(j);
                surrogate[i].Ap_sat(j)=surrogate[i].Ap_sat_old(j);
              }
            else
              {
                surrogate[i].gamma_org_sat(j)=1.0;
                surrogate[i].Ap_sat(j)=0.0;
              }
        }
    }
}


void mineur(int &n, int &i, int &j, Array<double, 2> &Jacobian, Array<double, 2> &com_Jacobian)
{
  com_Jacobian.resize(n-1,n-1);
  int index1,index2;
  int ic,jc;
  for (index1=0;index1<n;++index1)
    for (index2=0;index2<n;++index2)
      if (index1!=i and index2 !=j)
        {
          if (index1 < i)
            ic=index1;
          else
            ic=index1-1;
          if (index2 < j)
            jc=index2;
          else
            jc=index2-1;
          com_Jacobian(ic,jc)=Jacobian(index1,index2);
		  
        }
}


double determinant(int n, Array<double, 2> &Jacobian)
{
  double value=0.0;
  int i,j,l,i2,j2;
  int nl;
  int m=n;
  Array<double, 3> Matrix,Matrix_old;
  Array<double, 2> com_Jacobian,com_Jacobian_old;
  Array<double, 1> Factor,Factor_old;
  
  if (n>2)
    {
      nl=1;
      l=0;
      Matrix_old.resize(1,m,m);
      for (i=0;i<m;++i)
        for (j=0;j<m;++j)
          Matrix_old(l,i,j)=Jacobian(i,j);

      Factor_old.resize(nl);
      Factor_old(0)=1.0;
	  
      while (m>2)
        {
          int index=1;
          for (i=n;i>=m;i--)
            index*=i;

          Factor.resize(index);
          Matrix.resize(index,m-1,m-1);
          int iter=0;
          com_Jacobian_old.resize(m,m);
			
          for (l=0;l<nl;++l)
            {
              for (i=0;i<m;++i)
                for (j=0;j<m;++j)
                  com_Jacobian_old(i,j)=Matrix_old(l,i,j);
			  
              i=0;
              for (j=0;j<m;++j)
                {
                  mineur(m, i,j, com_Jacobian_old, com_Jacobian);
                  Factor(iter)=Factor_old(l)*pow(-1.0,i+j)*com_Jacobian_old(i,j);
                  for (i2=0;i2<m-1;++i2)
                    for (j2=0;j2<m-1;++j2)
                      Matrix(iter,i2,j2)=com_Jacobian(i2,j2);
                  ++iter;	
                }
            }
          Matrix_old.resize(index,m-1,m-1);
          for (l=0;l<index;++l)
            for (i=0;i<m-1;++i)
              for (j=0;j<m-1;++j)
                Matrix_old(l,i,j)=Matrix(l,i,j);
          m--;
          nl=index;
          Factor_old.resize(nl);
          for (l=0;l<nl;++l)
            Factor_old(l)=Factor(l);
        }

      value=0.0;
      for (l=0;l<nl;++l)
        value+=Factor(l)*(Matrix_old(l,0,0)*Matrix_old(l,1,1)-Matrix_old(l,1,0)*Matrix_old(l,0,1));
    }
  else if (n==2)
    value=Jacobian(0,0)*Jacobian(1,1)-Jacobian(0,1)*Jacobian(1,0);
  else
    value=Jacobian(0,0);
  
  return value;
  
}

void compute_comatrice(int &n, Array<double, 2> &Jacobian, Array<double, 2> &comatrice)
{
  int i,j;
  Array<double, 2> com_Jacobian;
  comatrice.resize(n,n);
  for (i=0;i<n;++i)
    for (j=0 ; j<n; ++j)
      {
        mineur(n,i,j,Jacobian,com_Jacobian);
        comatrice(i,j)=pow(-1.0,i+j)*determinant(n-1,com_Jacobian);
      }
}

void inversion_matrice(int &n, Array<double, 2> &Jacobian, Array<double, 2> &inverse)
{
  int i,j;
  Array<double, 2> comatrice;
  compute_comatrice(n,Jacobian,comatrice);
  double det=determinant(n,Jacobian);
  inverse.resize(n,n);
  
  for (i=0;i<n;++i)
    for (j=0;j<n;++j)
      inverse(i,j)=comatrice(j,i)/det;
}

void newton_raphson_sat(model_config &config, Array<double, 1> &MO_sat, double &AQ_sat,
                        Array<double, 1> &error, Array<double, 2> &Jacobian,
                        bool point_fixe)
{
  //method of newton raphson used to solved the coupled system
  //    Jacobian = | d(error1)/d(MO) d(error1)/d(AQ) |
  //               | d(error2)/d(MO) d(error2)/d(AQ) |
  //    Xi = | MO | at iteration i
  //         | AQ |
  //
  //    Y = |error1|
  //        |error2|
  //
  //    Xi+1 = Xi + inverse(Jacobian) Y
  //
  //    inverse(Jacobian) = 1/determinant * | d(error2)/d(AQ)  -d(error1)/d(AQ) |
  //                                        | -d(error2)/d(MO) d(error1)/d(MO)  |  
  int n=error.size();
  int i,j;
  Array<double, 2> inverse;
  double AQ_sat_old=AQ_sat;
  Array<double, 1> MO_sat_old;
  MO_sat_old.resize(n-1);
  for (j=0;j<n-1;++j)
    MO_sat_old(j)=MO_sat(j);
  
  if (determinant(n,Jacobian)!=0.0 and point_fixe==false)
    {
      inversion_matrice(n, Jacobian, inverse);
      for (i=0;i<n-1;++i)
        {
          for (j=0;j<n;++j)
            MO_sat(i)-=inverse(i,j)*error(j);
        }
      for (j=0;j<n;++j)
        AQ_sat-=inverse(n-1,j)*error(j);
      bool negative=false;
      for (i=0;i<n-1;++i)
        if (MO_sat(i)<=0.0)
          negative=true;
      if (AQ_sat<=0.0)
        negative=true;
      if (negative)
        {
          for (j=0;j<n-1;++j)
            MO_sat(j)=MO_sat_old(j)-error(j);
          AQ_sat=AQ_sat_old-error(n-1);

        }
	  
    }
  else
    {
      for (j=0;j<n-1;++j)
        MO_sat(j)-=error(j);
      AQ_sat-=error(n-1);
    }
  
}


void saturation(model_config &config, vector<species>& surrogate,
                bool all_hydrophobic,
                double &LWC, double &ionic, double &conc_inorganic,
                double &ionic_organic, double &chp, double &organion,
                double &Temperature, double &RH)
{
  //Determine the number of organic phase for which is system is stable and
  // compute the partitioning of compounds for this number of phases
  int nphase=2;
  int i,j;
  int n=surrogate.size();
  bool stable_system=false;
  int index_iter;
  double MMaq;

  if (config.initialized_saturation)
    if (surrogate[0].Ap_sat.size()>1)
      for (i=0;i<n;++i)
	{
	  surrogate[i].Ap_sat_save.resize(surrogate[0].Ap_sat.size());
	  for (j=0;j<surrogate[0].Ap_sat.size();++j)
	    surrogate[i].Ap_sat_save(j)=surrogate[i].Ap_sat(j);
	  surrogate[i].Aaq_save=surrogate[i].Aaq;
	}
    else
      for (i=0;i<n;++i)
	surrogate[i].Ap_sat_save.resize(1);
			  
  while (stable_system==false and nphase <= config.max_number_of_phases)
    {
      Array<double, 1> MO_sat, MOW,error_sat;
      Array<double, 2> Jacobian;
      double AQ_sat_old=0.0;
      double AQ_sat_old2=0.0;
      MO_sat.resize(nphase);
      MOW.resize(nphase);
      error_sat.resize(nphase+1);
      Jacobian.resize(nphase+1,nphase+1);
				  
      for (i=0;i<n;++i)
        {
          surrogate[i].gamma_org_sat_old.resize(nphase-1);
          surrogate[i].Ap_sat_old.resize(nphase-1);
          surrogate[i].Xorg_sat_old.resize(nphase-1);
        }
			  
      for (i=0;i<n;++i)
        if (surrogate[i].is_organic or i==config.iH2O)
          {
            if (nphase==2)
              {
                surrogate[i].gamma_org_sat_old(0)=surrogate[i].gamma_org;
                surrogate[i].Ap_sat_old(0)=surrogate[i].Ap;
              }
            else
              for (j=0;j<nphase-1;++j)
                {
                  surrogate[i].gamma_org_sat_old(j)=surrogate[i].gamma_org_sat(j);
                  surrogate[i].Ap_sat_old(j)=surrogate[i].Ap_sat(j);
                }
          }
        else
          for (j=0;j<nphase-1;++j)
            {
              surrogate[i].gamma_org_sat_old(j)=1.0;
              surrogate[i].Ap_sat_old(j)=0.0;
            }
	  
      for (i=0;i<n;++i)
        {
          surrogate[i].gamma_org_sat.resize(nphase);
          surrogate[i].Ap_sat.resize(nphase);
          surrogate[i].Xorg_sat.resize(nphase);
          surrogate[i].Ag_old=surrogate[i].Ag;
        }
				  
      double AQ_sat;
	  
      index_iter=0;
      bool reached_precision=false;
      bool point_fixe=false;
      if (config.initialized_saturation and surrogate[0].Ap_sat_save.size()>1)
	{	  
	  for (i=0;i<n;++i)
	    {
	      for (j=0;j<nphase;++j)
		surrogate[i].Ap_sat(j)=surrogate[i].Ap_sat_save(j);
	      if (nphase<surrogate[0].Ap_sat_save.size())
		for (j=nphase;j<surrogate[0].Ap_sat_save.size();++j)
		  surrogate[i].Ap_sat(nphase-1)+=surrogate[i].Ap_sat_save(j);

	      surrogate[i].Aaq=surrogate[i].Aaq_save;
	    }

	  AQ_sat=0.0;
	  for (i=0;i<n;++i)
	    AQ_sat+=surrogate[i].Aaq;

	  for (j=0;j<nphase;++j)
	    {
	      MO_sat(j)=0.0;
	      for (i=0;i<n;++i)
		MO_sat(j)+=surrogate[i].Ap_sat(j);	      
	    }
	}
      else
	init_saturation(config,surrogate,MO_sat,AQ_sat,LWC,conc_inorganic,all_hydrophobic);      
      
      double error_tot=0.0;
      int nh;
      if (config.compute_inorganic)
        nh=max(config.nh_inorg_init,max(config.nh_org_init,config.nh_aq_init));
      else
        nh=max(config.nh_org_init,config.nh_aq_init);

      Array<double, 2> vec_error;
      Array<double, 1> vec_error_chp;
      vec_error.resize(config.max_iter,nphase+1);
      vec_error_chp.resize(config.max_iter);
      bool non_convergence;
      int iiter=0;
      bool compute_activity_coefficients=true;
      double error3=1000.0;
      double chp_old;
      while ((index_iter < config.max_iter) and reached_precision==false)
        {
          if (config.first_evaluation_activity_coefficients==false or index_iter==0)
            compute_activity_coefficients=true;
          else
            compute_activity_coefficients=false;

          if (iiter>20)
            {
              non_convergence=false;
              for (i=max(index_iter-20,0);i<index_iter-1;i++)
		{
		  for (j=0;j<nphase+1;++j)
		    if ((abs((vec_error(i,j)-abs(error_sat(j)))/error_sat(j))<1.0e-3 and abs(error_sat(j))*nh>config.precision) or abs(error_sat(j)>1.0))
		      non_convergence=true;
		  
		  if ((abs((vec_error_chp(i)-abs(error3))/error3)<1.0e-3 and abs(error3)*nh>1.0e-3) or abs(error3)>1.0)
		    non_convergence=true;
		}
              
              if (non_convergence and nh<config.nh_max)
                {
                  ++nh;
                  for (i=max(index_iter-10,0);i<index_iter;i++)
                    for (j=0;j<nphase+1;++j)
                      vec_error(i,j)=-1.0;
                  iiter=0;
                }
            }
          iiter++;
	  
	  chp_old=chp; 
          if (all_hydrophobic)
            error_saturation_hydrophobic(config,surrogate,MO_sat,MOW,
                                         Temperature,RH,error_sat,Jacobian,1.0/nh, compute_activity_coefficients);
          else
	    error_saturation(config,surrogate,MO_sat,MOW,MMaq,AQ_sat,LWC,conc_inorganic,
			     ionic,ionic_organic,chp,organion,Temperature,RH,
			     error_sat,Jacobian,1.0/nh, compute_activity_coefficients);

	 
	  if (config.compute_inorganic)
	    error3=(chp_old-chp)/chp_old;
	  else
	    error3=0.0;

	  vec_error_chp(index_iter)=abs(error3);
          error_tot=0.0;
          for (j=0;j<nphase+1;++j)
            {
              error_tot=max(error_tot,abs(error_sat(j)));
              vec_error(index_iter,j)=abs(error_sat(j));
            }
		  
          for (j=0;j<nphase;++j)
            MO_sat(j)-=error_sat(j);
          AQ_sat-=error_sat(nphase); 

          reached_precision=true;
	  if (error3*nh>1.0e-3 or error_tot*nh > config.precision)
            reached_precision=false;
          ++index_iter;	  
        }        
      stability(config,surrogate,Temperature,stable_system,all_hydrophobic,MO_sat);
      ++nphase;
    }
  
  nphase=surrogate[0].Ap_sat.size();
  for (i=0;i<n;++i)
    {
      if (surrogate[i].hydrophobic==false and all_hydrophobic==false)
        for (j=0;j<nphase;++j)
          surrogate[i].Ap_sat(j)=0.0;
      if (surrogate[i].hydrophobic==false and all_hydrophobic)
        if (surrogate[i].is_organic==false and i!=config.iH2O)
          for (j=0;j<nphase;++j)
            surrogate[i].Ap_sat(j)=0.0;
	  
      if (surrogate[i].hydrophilic and all_hydrophobic)
	surrogate[i].Aaq=0.0;
      if (surrogate[i].hydrophilic==false)
	surrogate[i].Aaq=0.0; 
    }

  for (i=0;i<n;++i)
    {
      surrogate[i].Ap=0.0;
      for (j=0;j<nphase;++j)
        surrogate[i].Ap+=surrogate[i].Ap_sat(j);
    }
}
