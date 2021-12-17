//!!-----------------------------------------------------------------------
//!!     Copyright (C) 2019 CEREA (ENPC) - INERIS
//!!     SSH-aerosol is distributed under the GNU General Public License v3
//!!-----------------------------------------------------------------------

//#include "properties.cxx"
//using namespace ssh_soap;

void initialisation_eq_ssh(model_config &config, vector<species>& surrogate, double &Temperature, double &RH, double &ionic, double &chp, double &AQinit, bool all_hydrophobic)
{ 
  int n=surrogate.size();
  int i;
  double MOW=1.0;
  if (Temperature == 0.0)
    throw string("Error: Temperature is zero.");

  double gamma=pow(10,-0.511*pow(298.0/Temperature,1.5)*pow(ionic,0.5)/(1.0+pow(ionic,0.5)));
  config.initAQ=AQinit;
  config.ionicinit=ionic;
  config.chpinit=chp;
  for (i=0;i<n;i++)
    {
      surrogate[i].Aaqinit=surrogate[i].Aaq;
      surrogate[i].Aginit=surrogate[i].Ag;
    }

  for (i=0;i<n;i++)
    {
      surrogate[i].gamma_aq=1.0;
      if (surrogate[i].hydrophilic==false)
	surrogate[i].Aaq=0.0;
      if (surrogate[i].hydrophobic==false)
	surrogate[i].Ap=0.0;

    }

  config.Ke=1.010e-14*exp(-22.52*(298./Temperature-1.0)+26.92*(1+log(298./Temperature)-298./Temperature)); //Dissociation constant of water Ke=a(H+)*a(HO-)
  if (config.compute_inorganic)
    {
      double inorganion=0.0;
      double conc_org=0.0;      
      AQinit=0.0;

      surrogate[config.iSO4mm].Aaq+=surrogate[config.iH2SO4].Ag/surrogate[config.iH2SO4].MM*surrogate[config.iSO4mm].MM;
      surrogate[config.iH2SO4].Ag=0.;

      /* double conc_inorg=0.0;
         for (i=0;i<n;i++)        
         if (surrogate[i].is_inorganic_precursor==false and surrogate[i].is_organic==false and i!=config.iH2O and i!=config.iHp)           
         conc_inorg+=surrogate[i].Aaq/surrogate[i].MM;

         aH2O=RH/(1-RH)*surrogate[i].Aaq/surrogate[i].MM;*/


      
      //if (surrogate[config.iH2O].Aaq<1.0e-5*config.MOmin)
      //  {
      /*
        double mol=0.0;
        for (i=0;i<n;++i)
        if (i!=config.iH2O)
        {
        mol+=surrogate[i].Aaq/surrogate[i].MM;
        }
         
        mol=max(mol,config.MOmin);
        double RH2=min(RH,0.97);
        double ah2o=mol*RH2/(1.0-RH2)*18.0;
        //        }
      
        //if (ah2o>10*surrogate[config.iH2O].Aaq)
        surrogate[config.iH2O].Aaq=ah2o;

        double anh3=max(surrogate[config.iNH4p].Aaq/surrogate[config.iNH4p].MM-surrogate[config.iNO3m].Aaq/surrogate[config.iNO3m].MM,0.0);
        double sulf=surrogate[config.iSO4mm].Aaq/surrogate[config.iSO4mm].MM+surrogate[config.iHSO4m].Aaq/surrogate[config.iHSO4m].MM;
        surrogate[config.iSO4mm].Aaq=min(sulf,anh3*0.5)*surrogate[config.iSO4mm].MM;
        surrogate[config.iHSO4m].Aaq=(sulf-surrogate[config.iSO4mm].Aaq/surrogate[config.iSO4mm].MM)*surrogate[config.iHSO4m].MM;*/
      for (i=0;i<n;++i)
        {
          AQinit+=surrogate[i].Aaq;
          conc_org=0.;
          if (surrogate[i].is_organic or i==config.iH2O)      
            conc_org+=surrogate[i].Aaq;
        }
      
      for (i=0;i<n;++i)
        if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)          
	  conc_org=max(conc_org,surrogate[i].Aaq/surrogate[i].MM/config.molalmax*1000.0);      

      if (AQinit<=config.MOmin and config.solids==false)
        {
          //cout << "AQinit is zero" << endl;
          //exit(0);
          
          //If possible put some NO3 in the particle associate with ammonia
          double in_part=0.0;          
          in_part=min(0.1*surrogate[config.iHNO3].Ag/surrogate[config.iHNO3].MM,0.1*surrogate[config.iNH3].Ag/surrogate[config.iNH3].MM);
          surrogate[config.iNO3m].Aaq+=in_part*surrogate[config.iNO3m].MM;
          surrogate[config.iHNO3].Ag-=in_part*surrogate[config.iHNO3].MM;
          surrogate[config.iNH4p].Aaq+=0.999*in_part*surrogate[config.iNH4p].MM;
          surrogate[config.iNH3].Ag-=0.999*in_part*surrogate[config.iNH3].MM;

          AQinit=0.0;
          for (i=0;i<n;++i)
            {
              AQinit+=surrogate[i].Aaq;
              conc_org=0.0;
              if (surrogate[i].is_organic or i==config.iH2O)      
                conc_org+=surrogate[i].Aaq;
            }
          
          for (i=0;i<n;++i)
            if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)          
              conc_org=max(conc_org,surrogate[i].Aaq/surrogate[i].MM/config.molalmax*1000.0);
          
          if (AQinit<=config.MOmin)
            {
              //cout << "AQinit is zero" << endl;
              //exit(0);
          
              //If possible put some NO3 in the particle associate with ammonia
              double in_part=0.0;
              in_part=min(0.1*surrogate[config.iHCl].Ag/surrogate[config.iHCl].MM,0.1*surrogate[config.iNH3].Ag/surrogate[config.iNH3].MM);
              surrogate[config.iClm].Aaq+=in_part*surrogate[config.iClm].MM;
              surrogate[config.iHCl].Ag-=in_part*surrogate[config.iHCl].MM;
              surrogate[config.iNH4p].Aaq+=0.9*in_part*surrogate[config.iNH4p].MM;
              surrogate[config.iNH3].Ag-=0.9*in_part*surrogate[config.iNH3].MM;

              AQinit=0.0;
              conc_org=0.;
              for (i=0;i<n;++i)
                {
                  AQinit+=surrogate[i].Aaq;
                  if (surrogate[i].is_organic or i==config.iH2O)      
                    conc_org+=surrogate[i].Aaq;
                }

              for (i=0;i<n;++i)
                if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)          
                  conc_org=max(conc_org,surrogate[i].Aaq/surrogate[i].MM/config.molalmax*1000.0);
            }
        }

      AQinit=max(AQinit,config.MOmin);
      conc_org=max(conc_org,0.01*AQinit); //config.MOmin);

      //Ensure that the initial aerosol is at least neutral to prevent numerical issues
      for (i=0;i<n;++i)
        if (surrogate[i].is_inorganic_precursor==false and surrogate[i].is_organic==false and i!=config.iH2O and i!=config.iHp) 
          {
            inorganion-=surrogate[i].Aaq/surrogate[i].MM/conc_org*1000.*surrogate[i].charge;          
          }

      if (inorganion<0.) 
        {
          //Put enough NO3 in the particle to neutralize the aerosol (if possible)
          surrogate[config.iNO3m].Aaq+=min(-inorganion*conc_org/1000,surrogate[config.iHNO3].Ag/surrogate[config.iHNO3].MM)*surrogate[config.iNO3m].MM;
          surrogate[config.iHNO3].Ag-=min(-inorganion*conc_org/1000*surrogate[config.iHNO3].MM,surrogate[config.iHNO3].Ag);

          //Recompute inorganion
          inorganion=0.0;
          for (i=0;i<n;++i)
            if (surrogate[i].is_inorganic_precursor==false and surrogate[i].is_organic==false and i!=config.iH2O and i!=config.iHp)               
              inorganion-=surrogate[i].Aaq/surrogate[i].MM/conc_org*1000.*surrogate[i].charge;                        
          
          if (inorganion<0.)
            {
              surrogate[config.iClm].Aaq+=min(-inorganion*conc_org/1000,surrogate[config.iHCl].Ag/surrogate[config.iHCl].MM)*surrogate[config.iClm].MM;
              surrogate[config.iHCl].Ag-=min(-inorganion*conc_org/1000*surrogate[config.iHCl].MM,surrogate[config.iHCl].Ag);
              
              //Recompute inorganion
              inorganion=0.0;
              for (i=0;i<n;++i)
                if (surrogate[i].is_inorganic_precursor==false and surrogate[i].is_organic==false and i!=config.iH2O and i!=config.iHp)               
                  inorganion-=surrogate[i].Aaq/surrogate[i].MM/conc_org*1000.*surrogate[i].charge;             
            }
        }

      /*
      if (inorganion<0.)
        {
          cout << "not neutral" << endl;
          exit(0);
        }*/
  
      chp=max(0.5*(inorganion+pow(pow(inorganion,2)+4*config.Ke,0.5)),1.0e-15);
      surrogate[config.iHp].Aaq=chp*conc_org/1000.;
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
	    //cout << surrogate[i].name << " " << surrogate[i].Ap << " " << surrogate[i].Aaq << endl;
          }
        else
          {	    
            if(surrogate[i].hydrophobic or all_hydrophobic)
              {        
                if (surrogate[i].kp_from_experiment)
                  surrogate[i].kpi=surrogate[i].Kp_exp_org_ssh(Temperature);
                else if (surrogate[i].kp_from_experiment==false)
                  surrogate[i].kpi=surrogate[i].Kp_eff_org_ssh(Temperature, MOW);
              } 
            if ((surrogate[i].hydrophilic)) // and surrogate[i].Aaq > 1.e-15) 
              if (config.compute_aqueous_phase_properties)
                surrogate[i].kaqi=surrogate[i].Kpart_aq_ssh(Temperature,MOW);              
              else
                {
                  if (surrogate[i].aqt==2) //diacid
                    {
                      if ((pow(gamma, 2) * chp) == 0.0)
                        {
                          cout << "gamma :" << gamma << " chp: " << chp << endl;
                          throw string("Error: division by zero for diacid.");
                        }
                      surrogate[i].kaqi=surrogate[i].Kpart_aq_ssh(Temperature,MOW)*
                        (1.0+surrogate[i].Kacidity1/(pow(gamma,2)*chp)*
                         (1.0+surrogate[i].Kacidity2/(pow(gamma,2)*chp)));
                      surrogate[i].fioni1=(surrogate[i].Kacidity1/(pow(gamma,2)*chp))/
                        (1.0+surrogate[i].Kacidity1/(pow(gamma,2)*chp)*(1.0+surrogate[i].Kacidity2/(pow(gamma,2)*chp)));
                      surrogate[i].fioni2=(surrogate[i].Kacidity1/(pow(gamma,2)*chp))*(surrogate[i].Kacidity2/(pow(gamma,2)*chp))/
                        (1.0+surrogate[i].Kacidity1/(pow(gamma,2)*chp)*(1.0+surrogate[i].Kacidity2/(pow(gamma,2)*chp)));
                    }
                  else if (surrogate[i].aqt==1) //monoacid
                    {
                      if ((pow(gamma, 2) * chp) == 0.0)
                        throw string("Error: division by zero for monoacid.");                      
                      surrogate[i].kaqi=surrogate[i].Kpart_aq_ssh(Temperature,MOW)*(1.0+surrogate[i].Kacidity1/(pow(gamma,2)*chp));
                      surrogate[i].fioni1=(surrogate[i].Kacidity1/(pow(gamma,2)*chp))/(1.0+surrogate[i].Kacidity1/(pow(gamma,2)*chp));
                    }
                  else if (surrogate[i].aqt==3) //aldehyde
                    {
                      if ((pow(10, -surrogate[i].pHref)) == 0.0)
                        throw string("Error: division by zero for aldehyde.");
                      surrogate[i].kaqi=surrogate[i].Kpart_aq_ssh(Temperature,MOW) *
                        (1.0+surrogate[i].Koligo_aq*pow(gamma*chp/pow(10,-surrogate[i].pHref),surrogate[i].beta));
                    }
                  else
                    surrogate[i].kaqi=surrogate[i].Kpart_aq_ssh(Temperature,MOW);
                }
          }
  
  if (config.compute_inorganic)
    {
      double T0=298.;
      for (i=0;i<n;i++)
        if (surrogate[i].is_inorganic_precursor)     
          if (surrogate[i].is_solid==false)    
            {
              surrogate[i].keq=surrogate[i].Kequilibrium_ssh(Temperature);
            }
          else
            surrogate[i].keq=surrogate[i].Ksol*exp(-surrogate[i].deltaH*(T0/Temperature-1.0)-surrogate[i].dCp*(1.+log(T0/Temperature)-T0/Temperature));

      Kpideal_inorganic_ssh(config, surrogate, Temperature);  
    
    }

  if (config.activity_model=="unifac")
    {
      double tval1=1.0/298.15-1.0/Temperature;
      double tval2=298.15/Temperature-1.0+log(Temperature/298.15);
      //config.Inter2_aq.resize(config.nfunc_aq,config.nfunc_aq);
      //config.Inter2_org.resize(config.nfunc_org,config.nfunc_org);
      //config.Inter2_tot.resize(config.nfunc_tot,config.nfunc_tot);
      int j,k;
      if (config.temperature_dependancy)
	for (j=0;j<config.nfunc_aq;j++)
	  for (k=0;k<config.nfunc_aq;k++)          
	    config.Inter2_aq(j,k)=exp(-config.Inter_aq(j,k)/Temperature+config.InterB_aq(j,k)*tval1+config.InterC_aq(j,k)*tval2);         
      else
	for (j=0;j<config.nfunc_aq;j++)
	  for (k=0;k<config.nfunc_aq;k++)          
	    config.Inter2_aq(j,k)=exp(-config.Inter_aq(j,k)/Temperature);

      for (i=0;i<config.nmol_aq;i++)
        for (j=0;j<config.nfunc_aq;j++)
          if (config.groups_aq(j,i)>0.0)
            {
              config.sum2mol_aq(j,i)=0.0;
              for (k=0;k<config.nfunc_aq;k++)
                if (config.groups_aq(k,i)>0.0)
                  config.sum2mol_aq(j,i)+=config.surface_fraction_molaq(k,i)*config.Inter2_aq(k,j);              
            }

      for (i=0;i<config.nmol_aq;i++)
        for (j=0;j<config.nfunc_aq;j++)
          if (config.groups_aq(j,i)>0.0)
            {
              config.group_activity_molaq(j,i)=1.0-log(config.sum2mol_aq(j,i));                     
              for (k=0;k<config.nfunc_aq;k++)      
                if (config.groups_aq(k,i)>0.0)
                  config.group_activity_molaq(j,i)-=config.surface_fraction_molaq(k,i)*config.Inter2_aq(j,k)/config.sum2mol_aq(k,i);           
            }

      if (config.temperature_dependancy)
	for (j=0;j<config.nfunc_tot;j++)
	  for (k=0;k<config.nfunc_tot;k++)          
	    config.Inter2_tot(j,k)=exp(-config.Inter_tot(j,k)/Temperature+config.InterB_tot(j,k)*tval1+config.InterC_tot(j,k)*tval2);         
      else
	for (j=0;j<config.nfunc_tot;j++)
	  for (k=0;k<config.nfunc_tot;k++)          
	    config.Inter2_tot(j,k)=exp(-config.Inter_tot(j,k)/Temperature);

      for (i=0;i<config.nmol_tot;i++)
        for (j=0;j<config.nfunc_tot;j++)
          if (config.groups_tot(j,i)>0.0)
            {
              config.sum2mol_tot(j,i)=0.0;
              for (k=0;k<config.nfunc_tot;k++)
                if (config.groups_tot(k,i)>0.0)
                  config.sum2mol_tot(j,i)+=config.surface_fraction_moltot(k,i)*config.Inter2_tot(k,j);              
            }

      for (i=0;i<config.nmol_tot;i++)
        for (j=0;j<config.nfunc_tot;j++)
          if (config.groups_tot(j,i)>0.0)
            {
              config.group_activity_moltot(j,i)=1.0-log(config.sum2mol_tot(j,i));                     
              for (k=0;k<config.nfunc_tot;k++)      
                if (config.groups_tot(k,i)>0.0)
                  config.group_activity_moltot(j,i)-=config.surface_fraction_moltot(k,i)*config.Inter2_tot(j,k)/config.sum2mol_tot(k,i);           
            }

      if (config.temperature_dependancy)
	for (j=0;j<config.nfunc_org;j++)
	  for (k=0;k<config.nfunc_org;k++)          
	    config.Inter2_org(j,k)=exp(-config.Inter_org(j,k)/Temperature+config.InterB_org(j,k)*tval1+config.InterC_org(j,k)*tval2);         
      else
	for (j=0;j<config.nfunc_aq;j++)
	  for (k=0;k<config.nfunc_aq;k++)          
	    config.Inter2_org(j,k)=exp(-config.Inter_org(j,k)/Temperature); 

      
      for (i=0;i<config.nmol_org;i++)
        for (j=0;j<config.nfunc_org;j++)
          if (config.groups_org(j,i)>0.0)
            {
              config.sum2mol_org(j,i)=0.0;
              for (k=0;k<config.nfunc_org;k++)
                if (config.groups_org(k,i)>0.0)
                  config.sum2mol_org(j,i)+=config.surface_fraction_molorg(k,i)*config.Inter2_org(k,j);              
            }

      for (i=0;i<config.nmol_org;i++)
        for (j=0;j<config.nfunc_org;j++)
          if (config.groups_org(j,i)>0.0)
            {
              config.group_activity_molorg(j,i)=1.0-log(config.sum2mol_org(j,i));                     
              for (k=0;k<config.nfunc_org;k++)      
                if (config.groups_org(k,i)>0.0)
                  config.group_activity_molorg(j,i)-=config.surface_fraction_molorg(k,i)*config.Inter2_org(j,k)/config.sum2mol_org(k,i);           
            }

    }   
}

void error_org_ssh(model_config &config, vector<species>& surrogate,double &MOinit,double &MOW,
		   double &Temperature, double &error, double &derivative, double &RH, double factor,
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
        activity_coefficients_org_ssh(config, surrogate, all_hydrophobic, Temperature, MOW);

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
                surrogate[i].Ap=(1.-factor)*surrogate[i].Ap+factor*surrogate[i].Atot*Kp*MOinit/(1+Kp*MOinit);
                MO+=surrogate[i].Ap;
                derivative+=surrogate[i].Atot*pow(Kp,2)*MOinit/(pow(1+Kp*MOinit,2))
                  -surrogate[i].Atot*Kp/(1+Kp*MOinit);
              }
            else
              {
                if (MOW == 0.0 or surrogate[i].gamma_org == 0.0)
                  throw string("Error: division by zero in error_org.");
                Kp=surrogate[i].kpi/MOW/surrogate[i].gamma_org;
                surrogate[i].Ap=(1.-factor)*surrogate[i].Ap+factor*surrogate[i].Atot*Kp*MOinit/(1+Kp*MOinit);
                MO+=surrogate[i].Ap;
                derivative+=surrogate[i].Atot*pow(Kp,2)*MOinit/(pow(1+Kp*MOinit,2))
                  -surrogate[i].Atot*Kp/(1+Kp*MOinit);
              }
  
      if (config.iH2O>=0 and config.hygroscopicity)
        if (surrogate[config.iH2O].hydrophobic) //Can H2O condense on the organic phase
          hygroscopicity_org_ssh(config, surrogate, Temperature, MOW, RH, MOinit, MO, derivative, factor);

      MO=max(MO,config.MOmin);
      
      error=MOinit-MO;
    }
  else
    {
      MO=MOinit;
      error=0.0;
    }
}

void error_ph_ssh(model_config &config, vector<species> &surrogate, double Temperature, double &chp, 
              double organion, double &error, double &derivative, double AQinit, double LWC, 
              double MMaq, double MOinit, double MOW, double conc_org)
{ 
  //This routine is used to compute the electroneutrality conditions with a 
  //method of newton raphson. The routine computes the error between two 
  //iterations and derivative of the error.
  //In this routine, organic ions are taken into account but are assumed not to
  //strongly impact the pH (the derivative of organic ions concentrations) do
  //not have to be taken into account. 
  int n=surrogate.size();
  int i;
  double cion=0.0;
  double total;

  derivative=0.0; //-1.0;

  if (conc_org == 0.0 or chp == 0.0)
    {
      cout << conc_org << " " << MOW << " " << MOinit << " " << chp << endl;
      throw string("Error: division by zero in error_ph.");
    }
  
  for (i=0;i<n;++i)
    {
      if (surrogate[i].is_organic)
        {
          if (surrogate[i].hydrophilic)
            {
              if (surrogate[i].gamma_LR == 0.0)
                throw string("Error: division by zero gamma_LR.");
              if (surrogate[i].gamma_aq == 0.0)
                throw string("Error: division by zero gamma_aq.");
              if (surrogate[i].gamma_org == 0.0)
                throw string("Error: division by zero gamma_org.");                  
               
              if (surrogate[i].aqt==1) //monoacid
                {
                  double ratio_gamma=pow(surrogate[config.iHp].gamma_LR,2.0)*surrogate[config.iHp].gamma_SRMR/surrogate[i].gamma_LR;
                  double Kaq=surrogate[i].Kpart_aq_ssh(Temperature, MMaq)*(1.0+surrogate[i].Kacidity1/(ratio_gamma*chp))/surrogate[i].gamma_aq;
                  double sumk=1.0+Kaq*AQinit;
                  if (surrogate[i].hydrophobic)                    
                    sumk+=surrogate[i].kpi/MOW/surrogate[i].gamma_org*MOinit;                                         

                  double fion1=(surrogate[i].Kacidity1/(ratio_gamma*chp))/(1.0+surrogate[i].Kacidity1/(ratio_gamma*chp));
                  double dfion1=-(surrogate[i].Kacidity1/(ratio_gamma*chp*chp))/(1.0+surrogate[i].Kacidity1/(ratio_gamma*chp))
                    +pow(surrogate[i].Kacidity1/(ratio_gamma*chp),2)/pow(1.0+surrogate[i].Kacidity1/(ratio_gamma*chp),2)/chp;
                  double dKaq=-surrogate[i].Kpart_aq_ssh(Temperature, MMaq)*(surrogate[i].Kacidity1/(ratio_gamma*chp*chp))/surrogate[i].gamma_aq;

                  double Aaq=surrogate[i].Atot*Kaq*AQinit/sumk;                            
                  double dAaq=surrogate[i].Atot*dKaq*AQinit/sumk-surrogate[i].Atot*Kaq*dKaq*pow(AQinit,2.0)/pow(sumk,2);                                                     
                  cion+=Aaq*fion1/surrogate[i].MM/conc_org*1000.0;
                  derivative+=(Aaq*dfion1+dAaq*fion1)/surrogate[i].MM/conc_org*1000.0;
                }
              else if (surrogate[i].aqt==2) //diacid
                {
                  double ratio_gamma1=pow(surrogate[config.iHp].gamma_LR,2.0)*surrogate[config.iHp].gamma_SRMR/surrogate[i].gamma_LR;
                  double ratio_gamma2=pow(surrogate[config.iHp].gamma_LR,2.0)*surrogate[config.iHp].gamma_SRMR;
                  double Kaq=surrogate[i].Kpart_aq_ssh(Temperature, MMaq)*
                    (1.0+surrogate[i].Kacidity1/(ratio_gamma1*chp)*(1.0+surrogate[i].Kacidity2/(ratio_gamma2*chp)))/surrogate[i].gamma_aq;
                  double sumk=1.0+Kaq*AQinit;
                  if (surrogate[i].hydrophobic)                    
                    sumk+=surrogate[i].kpi/MOW/surrogate[i].gamma_org*MOinit;
                  

                  double fion1=(surrogate[i].Kacidity1/(ratio_gamma1*chp))/
                    (1.0+surrogate[i].Kacidity1/(ratio_gamma1*chp)*(1.0+surrogate[i].Kacidity2/(ratio_gamma2*chp)));
                  double fion2=(surrogate[i].Kacidity1/(ratio_gamma1*chp))*(surrogate[i].Kacidity2/(ratio_gamma2*chp))/
                    (1.0+surrogate[i].Kacidity1/(ratio_gamma1*chp)*(1.0+surrogate[i].Kacidity2/(ratio_gamma2*chp)));

                  double dfion1=-(surrogate[i].Kacidity1/(ratio_gamma1*chp*chp))/
                    (1.0+surrogate[i].Kacidity1/(ratio_gamma1*chp)*(1.0+surrogate[i].Kacidity2/(ratio_gamma2*chp)))
                    +(surrogate[i].Kacidity1/(ratio_gamma1*chp))/
                    pow(1.0+surrogate[i].Kacidity1/(ratio_gamma1*chp)*(1.0+surrogate[i].Kacidity2/(ratio_gamma2*chp)),2.)
                    *(surrogate[i].Kacidity1/(ratio_gamma1*chp*chp)*(1.0+surrogate[i].Kacidity2/(ratio_gamma2*chp))+
                      surrogate[i].Kacidity1/(ratio_gamma1*chp)*surrogate[i].Kacidity2/(ratio_gamma2*chp*chp));
                  double dfion2=-2.0*(surrogate[i].Kacidity1/(ratio_gamma1*chp*chp))*(surrogate[i].Kacidity2/(ratio_gamma2*chp))/
                    (1.0+surrogate[i].Kacidity1/(ratio_gamma1*chp)*(1.0+surrogate[i].Kacidity2/(ratio_gamma2*chp)))+
                    (surrogate[i].Kacidity1/(ratio_gamma1*chp))*(surrogate[i].Kacidity2/(ratio_gamma2*chp))/
                    pow(1.0+surrogate[i].Kacidity1/(ratio_gamma1*chp)*(1.0+surrogate[i].Kacidity2/(ratio_gamma2*chp)),2.)
                    *(surrogate[i].Kacidity1/(ratio_gamma1*chp*chp)*(1.0+surrogate[i].Kacidity2/(ratio_gamma2*chp))+
                      surrogate[i].Kacidity1/(ratio_gamma1*chp)*surrogate[i].Kacidity2/(ratio_gamma2*chp*chp));

                  double Aaq=surrogate[i].Atot*Kaq*AQinit/sumk;                      
                  double dKaq=surrogate[i].Kpart_aq_ssh(Temperature, MMaq)*
                    (-surrogate[i].Kacidity1/(ratio_gamma1*chp*chp)*(1.0+surrogate[i].Kacidity2/(ratio_gamma2*chp))
                     -(1.0+surrogate[i].Kacidity1/(ratio_gamma1*chp)*surrogate[i].Kacidity2/(ratio_gamma2*chp*chp)))/surrogate[i].gamma_aq;
                  double dAaq=surrogate[i].Atot*dKaq*AQinit/sumk-surrogate[i].Atot*Kaq*dKaq*pow(AQinit,2.0)/pow(sumk,2);

                  cion+=Aaq*(fion1+2.0*fion2)/surrogate[i].MM/conc_org*1000.0;
                  derivative+=(Aaq*(dfion1+2.*dfion2)+dAaq*(fion1+2.0*fion2))/surrogate[i].MM/conc_org*1000.0;
                }
            }
        }
    }

  //adding concentrations of inorganic ions depending on pH
  //NH3:
  i=config.iNH3;
  double Kaq=surrogate[i].Kaq_inorg;
  total=surrogate[i].Ag/surrogate[i].MM+surrogate[config.iNH4p].Aaq/surrogate[config.iNH4p].MM; //gas + particle concentration
  derivative-=1000.*total*
    (Kaq/chp/(1.0+Kaq*conc_org)
     -Kaq/pow(1.0+Kaq*conc_org,2.0)*Kaq/chp*conc_org);
  cion-=1000.*total*Kaq/(1.0+Kaq*conc_org); //concentration of NH3+
 
  //HNO3:
  i=config.iHNO3;
  Kaq=surrogate[i].Kaq_inorg;
  total=surrogate[i].Ag/surrogate[i].MM+surrogate[config.iNO3m].Aaq/surrogate[config.iNO3m].MM; //gas + particle concentration
  derivative+=1000.*total*
    (-Kaq/chp/(1.0+Kaq*conc_org)
     +Kaq/pow(1.0+Kaq*conc_org,2.0)*Kaq/chp*conc_org);
  cion+=1000.*total*Kaq/(1.0+Kaq*conc_org); //concentration of NO3-

  //HCl:
  i=config.iHCl;
  Kaq=surrogate[i].Kaq_inorg;
  total=surrogate[i].Ag/surrogate[i].MM+surrogate[config.iClm].Aaq/surrogate[config.iClm].MM; //gas + particle concentration
  derivative+=1000.*total*
    (-surrogate[i].Kaq_inorg/chp/(1.0+surrogate[i].Kaq_inorg*conc_org)
     +surrogate[i].Kaq_inorg/pow(1.0+surrogate[i].Kaq_inorg*conc_org,2.0)*surrogate[i].Kaq_inorg/chp*conc_org);
  cion+=1000.*total*surrogate[i].Kaq_inorg/(1.0+surrogate[i].Kaq_inorg*conc_org); //concentration of Cl-

  //H2SO4:
  i=config.iH2SO4;
  total=surrogate[i].Ag/surrogate[i].MM+surrogate[config.iHSO4m].Aaq/surrogate[config.iHSO4m].MM
    +surrogate[config.iSO4mm].Aaq/surrogate[config.iSO4mm].MM; //gas + particle concentration
  double K=surrogate[i].keq*surrogate[config.iHSO4m].gamma_aq
    /(surrogate[config.iHp].gamma_aq*surrogate[config.iSO4mm].gamma_aq);
      
  derivative-=1000.*total/conc_org*1.0/pow(1.0+K/chp,2.0)*K/(chp*chp); //HSO4-+SO4--
  cion+=1000.*total/conc_org*(2.0-1.0/(1.0+K/chp));     //concentration of HSO4- + 2*SO4--

  //Na:
  if (config.iNa>=0) cion-=surrogate[config.iNa].Aaq/surrogate[config.iNa].MM/conc_org*1000.*surrogate[config.iNa].charge;
  if (config.iK>=0) cion-=surrogate[config.iK].Aaq/surrogate[config.iK].MM/conc_org*1000.*surrogate[config.iK].charge;
  if (config.iMg>=0) cion-=surrogate[config.iCa].Aaq/surrogate[config.iCa].MM/conc_org*1000.*surrogate[config.iCa].charge;
  if (config.iCa>=0) cion-=surrogate[config.iMg].Aaq/surrogate[config.iMg].MM/conc_org*1000.*surrogate[config.iMg].charge;
    

  //inorganion*=AQinit/conc_org;
  //cout << "inorganion: " << inorganion << endl;
  derivative=derivative*0.5*(1.0+cion/pow(pow(cion,2)+4*config.Ke,0.5))-1.0;
  
  error=0.5*(cion+pow(pow(cion,2)+4*config.Ke,0.5))-chp;
  //cout << organion+inorganion+pow(pow(organion+inorganion,2)+4*Ke,0.5) << " " << organion << " " << inorganion << " " << error << endl;
}

void error_ph_sat_ssh(model_config &config, vector<species> &surrogate, double Temperature, double &chp, 
                  double organion, double &error, double &derivative, double AQinit, double LWC, 
                  double MMaq, Array <double,1> MOinit, Array <double,1> MOW, double conc_org)
{ 
  //This routine is used to compute the electroneutrality conditions with a 
  //method of newton raphson. The routine computes the error between two 
  //iterations and derivative of the error.
  //In this routine, organic ions are taken into account but are assumed not to
  //strongly impact the pH (the derivative of organic ions concentrations) do
  //not have to be taken into account. 
  int n=surrogate.size();
  int i;
  double cion=0.0;
  double total;
  int iphase;
  int nphase=MOinit.size();

  if (conc_org == 0.0 or chp == 0.0)
    throw string("Error: division by zero in error_ph_sat.");
  
  derivative=0.0; //-1.0;
  for (i=0;i<n;++i)
    {
      if (surrogate[i].is_organic)
        {
          if (surrogate[i].hydrophilic)
            {
              if (surrogate[i].gamma_LR == 0.0)
                throw string("Error: division by zero gamma_LR in error_ph_sat.");
              if (surrogate[i].gamma_aq == 0.0)
                throw string("Error: division by zero gamma_aq in error_ph_sat.");
              if (surrogate[i].gamma_org == 0.0)
                throw string("Error: division by zero gamma_org in error_ph_sat.");                  
              
              if (surrogate[i].aqt==1) //monoacid
                {                      
                  double ratio_gamma=pow(surrogate[config.iHp].gamma_LR,2.0)*surrogate[config.iHp].gamma_SRMR/surrogate[i].gamma_LR;
                  double Kaq=surrogate[i].Kpart_aq_ssh(Temperature, MMaq)*(1.0+surrogate[i].Kacidity1/(ratio_gamma*chp))/surrogate[i].gamma_aq;
                  double sumk=1.0+Kaq*AQinit;
                  if (surrogate[i].hydrophobic)                    
                    for (iphase=0;iphase<nphase;iphase++)
                      sumk+=surrogate[i].kpi/MOW(iphase)/surrogate[i].gamma_org_sat(iphase)*MOinit(iphase);                                         

                  double fion1=(surrogate[i].Kacidity1/(ratio_gamma*chp))/(1.0+surrogate[i].Kacidity1/(ratio_gamma*chp));
                  double dfion1=-(surrogate[i].Kacidity1/(ratio_gamma*chp*chp))/(1.0+surrogate[i].Kacidity1/(ratio_gamma*chp))
                    +pow(surrogate[i].Kacidity1/(ratio_gamma*chp),2)/pow(1.0+surrogate[i].Kacidity1/(ratio_gamma*chp),2)/chp;
                  double dKaq=-surrogate[i].Kpart_aq_ssh(Temperature, MMaq)*(surrogate[i].Kacidity1/(ratio_gamma*chp*chp))/surrogate[i].gamma_aq;

                  double Aaq=surrogate[i].Atot*Kaq*AQinit/sumk;                            
                  double dAaq=surrogate[i].Atot*dKaq*AQinit/sumk-surrogate[i].Atot*Kaq*dKaq*pow(AQinit,2.0)/pow(sumk,2);                                                     
                  cion+=Aaq*fion1/surrogate[i].MM/conc_org*1000.0;
                  derivative+=(Aaq*dfion1+dAaq*fion1)/surrogate[i].MM/conc_org*1000.0;
                }
              else if (surrogate[i].aqt==2) //diacid
                {
                  double ratio_gamma1=pow(surrogate[config.iHp].gamma_LR,2.0)*surrogate[config.iHp].gamma_SRMR/surrogate[i].gamma_LR;
                  double ratio_gamma2=pow(surrogate[config.iHp].gamma_LR,2.0)*surrogate[config.iHp].gamma_SRMR;
                  double Kaq=surrogate[i].Kpart_aq_ssh(Temperature, MMaq)*
                    (1.0+surrogate[i].Kacidity1/(ratio_gamma1*chp)*(1.0+surrogate[i].Kacidity2/(ratio_gamma2*chp)))/surrogate[i].gamma_aq;
                  double sumk=1.0+Kaq*AQinit;
                  if (surrogate[i].hydrophobic)                    
                    for (iphase=0;iphase<nphase;iphase++)
                      sumk+=surrogate[i].kpi/MOW(iphase)/surrogate[i].gamma_org_sat(iphase)*MOinit(iphase);                    

                  double fion1=(surrogate[i].Kacidity1/(ratio_gamma1*chp))/
                    (1.0+surrogate[i].Kacidity1/(ratio_gamma1*chp)*(1.0+surrogate[i].Kacidity2/(ratio_gamma2*chp)));
                  double fion2=(surrogate[i].Kacidity1/(ratio_gamma1*chp))*(surrogate[i].Kacidity2/(ratio_gamma2*chp))/
                    (1.0+surrogate[i].Kacidity1/(ratio_gamma1*chp)*(1.0+surrogate[i].Kacidity2/(ratio_gamma2*chp)));

                  double dfion1=-(surrogate[i].Kacidity1/(ratio_gamma1*chp*chp))/
                    (1.0+surrogate[i].Kacidity1/(ratio_gamma1*chp)*(1.0+surrogate[i].Kacidity2/(ratio_gamma2*chp)))
                    +(surrogate[i].Kacidity1/(ratio_gamma1*chp))/
                    pow(1.0+surrogate[i].Kacidity1/(ratio_gamma1*chp)*(1.0+surrogate[i].Kacidity2/(ratio_gamma2*chp)),2.)
                    *(surrogate[i].Kacidity1/(ratio_gamma1*chp*chp)*(1.0+surrogate[i].Kacidity2/(ratio_gamma2*chp))+
                      surrogate[i].Kacidity1/(ratio_gamma1*chp)*surrogate[i].Kacidity2/(ratio_gamma2*chp*chp));
                  double dfion2=-2.0*(surrogate[i].Kacidity1/(ratio_gamma1*chp*chp))*(surrogate[i].Kacidity2/(ratio_gamma2*chp))/
                    (1.0+surrogate[i].Kacidity1/(ratio_gamma1*chp)*(1.0+surrogate[i].Kacidity2/(ratio_gamma2*chp)))+
                    (surrogate[i].Kacidity1/(ratio_gamma1*chp))*(surrogate[i].Kacidity2/(ratio_gamma2*chp))/
                    pow(1.0+surrogate[i].Kacidity1/(ratio_gamma1*chp)*(1.0+surrogate[i].Kacidity2/(ratio_gamma2*chp)),2.)
                    *(surrogate[i].Kacidity1/(ratio_gamma1*chp*chp)*(1.0+surrogate[i].Kacidity2/(ratio_gamma2*chp))+
                      surrogate[i].Kacidity1/(ratio_gamma1*chp)*surrogate[i].Kacidity2/(ratio_gamma2*chp*chp));

                  double Aaq=surrogate[i].Atot*Kaq*AQinit/sumk;                      
                  double dKaq=surrogate[i].Kpart_aq_ssh(Temperature, MMaq)*
                    (-surrogate[i].Kacidity1/(ratio_gamma1*chp*chp)*(1.0+surrogate[i].Kacidity2/(ratio_gamma2*chp))
                     -(1.0+surrogate[i].Kacidity1/(ratio_gamma1*chp)*surrogate[i].Kacidity2/(ratio_gamma2*chp*chp)))/surrogate[i].gamma_aq;
                  double dAaq=surrogate[i].Atot*dKaq*AQinit/sumk-surrogate[i].Atot*Kaq*dKaq*pow(AQinit,2.0)/pow(sumk,2);

                  cion+=Aaq*(fion1+2.0*fion2)/surrogate[i].MM/conc_org*1000.0;
                  derivative+=(Aaq*(dfion1+2.*dfion2)+dAaq*(fion1+2.0*fion2))/surrogate[i].MM/conc_org*1000.0;
                }
            }
        }
    }

  //adding concentrations of inorganic ions depending on pH
  //NH3:
  total=surrogate[config.iNH3].Ag/surrogate[config.iNH3].MM+surrogate[config.iNH4p].Aaq/surrogate[config.iNH4p].MM; //gas + particle concentration
  derivative-=1000.*total*
    (surrogate[config.iNH3].Kaq_inorg/chp/(1.0+surrogate[config.iNH3].Kaq_inorg*conc_org)
     -surrogate[config.iNH3].Kaq_inorg/pow(1.0+surrogate[config.iNH3].Kaq_inorg*conc_org,2.0)*surrogate[config.iNH3].Kaq_inorg/chp*conc_org);
  cion-=1000.*total*surrogate[config.iNH3].Kaq_inorg/(1.0+surrogate[config.iNH3].Kaq_inorg*conc_org); //concentration of NH3+
 
  //HNO3:
  total=surrogate[config.iHNO3].Ag/surrogate[config.iHNO3].MM+surrogate[config.iNO3m].Aaq/surrogate[config.iNO3m].MM; //gas + particle concentration
  derivative+=1000.*total*
    (-surrogate[config.iHNO3].Kaq_inorg/chp/(1.0+surrogate[config.iHNO3].Kaq_inorg*conc_org)
     +surrogate[config.iHNO3].Kaq_inorg/pow(1.0+surrogate[config.iHNO3].Kaq_inorg*conc_org,2.0)*surrogate[config.iHNO3].Kaq_inorg/chp*conc_org);
  cion+=1000.*total*surrogate[config.iHNO3].Kaq_inorg/(1.0+surrogate[config.iHNO3].Kaq_inorg*conc_org); //concentration of NO3-

  //HCl:
  total=surrogate[config.iHCl].Ag/surrogate[config.iHCl].MM+surrogate[config.iClm].Aaq/surrogate[config.iClm].MM; //gas + particle concentration
  derivative+=1000.*total*
    (-surrogate[config.iHCl].Kaq_inorg/chp/(1.0+surrogate[config.iHCl].Kaq_inorg*conc_org)
     +surrogate[config.iHCl].Kaq_inorg/pow(1.0+surrogate[config.iHCl].Kaq_inorg*conc_org,2.0)*surrogate[config.iHCl].Kaq_inorg/chp*conc_org);
  cion+=1000.*total*surrogate[config.iHCl].Kaq_inorg/(1.0+surrogate[config.iHCl].Kaq_inorg*conc_org); //concentration of Cl-

  //H2SO4:
  total=surrogate[config.iH2SO4].Ag/surrogate[config.iH2SO4].MM+surrogate[config.iHSO4m].Aaq/surrogate[config.iHSO4m].MM
    +surrogate[config.iSO4mm].Aaq/surrogate[config.iSO4mm].MM; //gas + particle concentration
  double K=surrogate[config.iH2SO4].keq*surrogate[config.iHSO4m].gamma_aq
    /(surrogate[config.iHp].gamma_aq*surrogate[config.iSO4mm].gamma_aq);
      
  derivative-=1000.*total/conc_org*1.0/pow(1.0+K/chp,2.0)*K/(chp*chp); //HSO4-+SO4--
  cion+=1000.*total/conc_org*(2.0-1.0/(1.0+K/chp));     //concentration of HSO4- + 2*SO4--

  //Na: 
  if (config.iNa>=0) cion-=surrogate[config.iNa].Aaq/surrogate[config.iNa].MM/conc_org*1000.*surrogate[config.iNa].charge;
  if (config.iK>=0) cion-=surrogate[config.iK].Aaq/surrogate[config.iK].MM/conc_org*1000.*surrogate[config.iK].charge;
  if (config.iMg>=0) cion-=surrogate[config.iCa].Aaq/surrogate[config.iCa].MM/conc_org*1000.*surrogate[config.iCa].charge;
  if (config.iCa>=0) cion-=surrogate[config.iMg].Aaq/surrogate[config.iMg].MM/conc_org*1000.*surrogate[config.iMg].charge;

  //inorganion*=AQinit/conc_org;
  //cout << "inorganion: " << inorganion << endl;
  derivative=derivative*0.5*(1.0+cion/pow(pow(cion,2)+4*config.Ke,0.5))-1.0;
  
  error=0.5*(cion+pow(pow(cion,2)+4*config.Ke,0.5))-chp;
  //cout << organion+inorganion+pow(pow(organion+inorganion,2)+4*Ke,0.5) << " " << organion << " " << inorganion << " " << error << endl;
}

void solidification_ssh(model_config &config, vector<species>& surrogate, double& conc_org2, double& factor)
{
  int n=surrogate.size();
  int i;
  Array <double, 1> excess,ratio;
  double Ke,prod_conc,xmol,total;     
  excess.resize(n);
  ratio.resize(n);
  excess=0.0; 
  double sum_excess=0.0;
  ratio=1.0;
  double sum_error_ap=10;
  double factor2=1.;

  double conc_org=0.;
  for (i=0;i<n;++i)
    if (surrogate[i].is_organic or i==config.iH2O)      
      conc_org+=surrogate[i].Aaq;

  conc_org=max(conc_org,1.e-5*config.MOmin); //0.00001*AQinit); //config.MOmin);
  for (i=0;i<n;++i)
    if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)      
      surrogate[i].molality=surrogate[i].Aaq/surrogate[i].MM/conc_org*1000.0;                 

  for (i=0;i<n;i++)
    {
      surrogate[i].Ag2=surrogate[i].Ag;
      surrogate[i].Ap2=surrogate[i].Ap;
      surrogate[i].Aaq2=surrogate[i].Aaq;
      surrogate[i].molality2=surrogate[i].molality;
      //if (surrogate[i].nion>2 and surrogate[i].Ap>config.MOmin)
      //  factor2=1.0; //0.5;
    }
  
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
            //surrogate[iion1].molality=surrogate[iion1].Aaq/surrogate[iion1].MM/conc_org*1000.;
            //surrogate[iion2].molality=surrogate[iion2].Aaq/surrogate[iion2].MM/conc_org*1000.;

            //cout << iion1 << " " << iion2 << " " << surrogate[iion1].molality << " " << surrogate[iion2].molality << endl;
            prod_conc=pow(surrogate[iion1].molality2,pion1)*pow(surrogate[iion2].molality2,pion2);
            xmol=-surrogate[i].Ap2/surrogate[i].MM/conc_org*1000.;            
            Ke=surrogate[i].keq/(pow(surrogate[iion1].gamma_aq,pion1)*pow(surrogate[iion2].gamma_aq,pion2));
            //cout << surrogate[i].name << " Ke: " << Ke << " " << prod_conc << " " << prod_conc2 << " " << conc_org <<  " " << factor << " " << endl;
           
            if (surrogate[i].nion>2)
              {
                //surrogate[iion3].molality=surrogate[iion3].Aaq/surrogate[iion3].MM/conc_org*1000.;
                Ke=Ke/pow(surrogate[iion3].gamma_aq,pion3);
                prod_conc*=pow(surrogate[iion3].molality2,pion3);               
              }
            else
              pion3=0.0;
            excess(i)=prod_conc-Ke;
        
            //cout << "excess: " << excess(i) << endl; 
            if ((excess(i)>0. and prod_conc>0.) or (excess(i)<0. and xmol<0.)) 
              {    

		//ratio=1.;
		if (iion1==config.iNH4p or iion2==config.iNH4p or iion3==config.iNH4p)
		  {
		    total=surrogate[config.iNH3].Ag2/surrogate[config.iNH3].MM+surrogate[config.iNH4p].Aaq2/surrogate[config.iNH4p].MM;
		    if(total>0.) ratio(config.iNH4p)=surrogate[config.iNH4p].Aaq2/surrogate[config.iNH4p].MM/total;
		  }

                if  (iion1==config.iNO3m or iion2==config.iNO3m or iion3==config.iNO3m)
		  {
		    total=surrogate[config.iHNO3].Ag2/surrogate[config.iHNO3].MM+surrogate[config.iNO3m].Aaq2/surrogate[config.iNO3m].MM;
		    if(total>0.) ratio(config.iNO3m)=surrogate[config.iNO3m].Aaq2/surrogate[config.iNO3m].MM/total;
		  }


                if (iion1==config.iClm or iion2==config.iClm or iion3==config.iClm)
		  {
		    total=surrogate[config.iHCl].Ag2/surrogate[config.iHCl].MM+surrogate[config.iClm].Aaq2/surrogate[config.iClm].MM;
		    if(total>0.) ratio(config.iClm)=surrogate[config.iClm].Aaq2/surrogate[config.iClm].MM/total;
		  }

		if (iion1==config.iSO4mm or iion1==config.iHSO4m or iion2==config.iSO4mm or iion2==config.iHSO4m or iion3==config.iSO4mm or iion3==config.iHSO4m)
		  {
		    total=surrogate[config.iSO4mm].Aaq2/surrogate[config.iSO4mm].MM+surrogate[config.iHSO4m].Aaq2/surrogate[config.iHSO4m].MM;
		    if(total>0.) ratio(config.iSO4mm)=surrogate[config.iSO4mm].Aaq2/surrogate[config.iSO4mm].MM/total;
		    if(total>0.) ratio(config.iHSO4m)=1.0-ratio(config.iSO4mm);
		  }

		//ratio=1.;
		/*
                if (iion1==config.iNH4p or iion2==config.iNH4p or iion3==config.iNH4p)
		  {
		    total=surrogate[config.iNH3].Ag2/surrogate[config.iNH3].MM+surrogate[config.iNH4p].Aaq2/surrogate[config.iNH4p].MM;
		    if(total>0.) ratio(config.iNH4p)=surrogate[config.iNH4p].Aaq2/surrogate[config.iNH4p].MM/total;
		  }

                if  (iion1==config.iNO3m or iion2==config.iNO3m or iion3==config.iNO3m)
		  {
		    total=surrogate[config.iHNO3].Ag2/surrogate[config.iHNO3].MM+surrogate[config.iNO3m].Aaq2/surrogate[config.iNO3m].MM;
		    if(total>0.) ratio(config.iNO3m)=surrogate[config.iNO3m].Aaq2/surrogate[config.iNO3m].MM/total;
		  }

                if (iion1==config.iClm or iion2==config.iClm or iion3==config.iClm)
		  {
		    total=surrogate[config.iHCl].Ag2/surrogate[config.iHCl].MM+surrogate[config.iClm].Aaq2/surrogate[config.iClm].MM;
		    if(total>0.) ratio(config.iClm)=surrogate[config.iClm].Aaq2/surrogate[config.iClm].MM/total;
		  }

                if (iion1==config.iSO4mm or iion1==config.iHSO4m or iion2==config.iSO4mm or iion2==config.iHSO4m or iion3==config.iSO4mm or iion3==config.iHSO4m)
		  {
		    total=surrogate[config.iSO4mm].Aaq2/surrogate[config.iSO4mm].MM+surrogate[config.iHSO4m].Aaq2/surrogate[config.iHSO4m].MM;
		    if(total>0.) ratio(config.iSO4mm)=surrogate[config.iSO4mm].Aaq2/surrogate[config.iSO4mm].MM/total;
		    if(total>0.) ratio(config.iHSO4m)=1.0-ratio(config.iSO4mm);
		  }*/

                //double xmol_save=xmol;        
		/*
		  if (excess(i)>0. or excess(i)<0.) //prod_conc2>Ke)
                  {*/

		if (excess(i)<0.) // and iion3<1)
		  xmol=-surrogate[i].Ap2/surrogate[i].MM/conc_org*1000.; //0.; //-pow(Ke,1./3);
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
		    //cout << iter << endl;
                           
		    error=pow(m1,pion1)*pow(m2,pion2);
		    derror=-ratio(iion1)*pion1*pow(m1,pion1-1)*pion1*pow(m2,pion2)*pow(m3,pion3)-ratio(iion2)*pion2*pion2*pow(m2,pion2-1)*pow(m1,pion1)*pow(m3,pion3);
		    if (iion3>0) 
		      {
			derror-=ratio(iion3)*pion3*pion3*pow(m3,pion3-1)*pow(m2,pion2)*pow(m1,pion1);                   
			error*=pow(m3,pion3);
		      }
		    error=error-Ke;                  
		    if (abs(derror)>0.)                      
		      xmol=xmol-error/derror;
		    /*
		      if (surrogate[iion1].molality-ratio(iion1)*xmol*pion1<=0.0)                                             
		      xmol=0.99*min(surrogate[iion1].molality/ratio(iion1)/pion1,xmol_save);
		      if (surrogate[iion2].molality-ratio(iion2)*xmol*pion2<=0.0)                                             
		      xmol=0.99*min(surrogate[iion2].molality/ratio(iion2)/pion2,xmol_save);*/
     
		    iter++;
		    //cout << "error: " << error << " " << error/Ke << endl;
		    //xmol_save=xmol;
		    //cout << "X: " << iter << " " << xmol << " " << error << " " << derror << " " << Ke << " " << pow(m1,pion1)*pow(m2,pion2) << endl;                                                                       
                             
		  }
                /*
                if (i==46 and excess(i)<0.)
                  {
                    cout << xmol << " " << error << " " << Ke << endl;
                    //exit(0);
                  }
                
                if (iter==2000)
                  {
                    cout << "problem " << endl;
                    exit(0);
                    }*/
		//if (xmol<0.) cout << surrogate[i].name << " " << iter << " " << xmol << " " << -surrogate[i].Ap/conc_org*1000./surrogate[i].MM << endl;
		//xmol=max(xmol,xmol_save);
		//xmol=min(xmol,surrogate[iion1].molality/ratio(iion1)/pion1);
		//xmol=min(xmol,surrogate[iion2].molality/ratio(iion2)/pion2);
		//if (iion3>0) xmol=min(xmol,surrogate[iion3].molality/ratio(iion3)/pion3);
		excess(iion1)-=xmol*ratio(iion1)*pion1;
		excess(iion2)-=xmol*ratio(iion2)*pion2;
		//  }
                excess(i)=xmol;
                //cout << xmol << endl;
              }             

            if (prod_conc==0.0) excess(i)=0.0;             
            if (xmol>0.) sum_excess+=xmol;

            //if (xmol<0.) excess(i)=xmol_save;
            //cout << surrogate[i].name << " " << excess(i) << " " << xmol << " " << surrogate[i].Ap << " " << surrogate[i].keq << endl;

	    if (iion1==config.iSO4mm or iion1==config.iHSO4m or iion2==config.iSO4mm or iion2==config.iHSO4m or iion3==config.iSO4mm or iion3==config.iHSO4m)
	      {
		total=surrogate[config.iSO4mm].Aaq2/surrogate[config.iSO4mm].MM+surrogate[config.iHSO4m].Aaq2/surrogate[config.iHSO4m].MM;
		if(total>0.) ratio(config.iSO4mm)=surrogate[config.iSO4mm].Aaq2/surrogate[config.iSO4mm].MM/total;
		if(total>0.) ratio(config.iHSO4m)=1.0-ratio(config.iSO4mm);       
	      }

            //excess*=factor2;
            if (excess(i)!=excess(i)) exit(0);
        
            if (excess(i)<0.) //dissolution
              {
                excess(i)=factor2*max(-surrogate[i].Ap2/surrogate[i].MM/conc_org*1000,excess(i));
                //cout << "dissolution " << surrogate[i].Ap << " " << factor << " " << surrogate[config.iH2O].Aaq << " " << surrogate[i].gamma_aq << " "<< excess(i) << endl;
                //surrogate[i].Ap=factor*excess(i)*surrogate[i].MM*conc_org/1000.+surrogate[i].Ap;
                if (surrogate[i].Ap2>0.) sum_error_ap+=abs(excess(i)*surrogate[i].MM*conc_org/1000./surrogate[i].Ap2);
                surrogate[i].Ap2=max(excess(i)*surrogate[i].MM*conc_org/1000.+surrogate[i].Ap2,0.);
                //cout << surrogate[i].Ap << endl;            
                //int iion1=surrogate[i].iion1;
                //int pion1=surrogate[i].pion1;
                surrogate[iion1].Aaq2=-excess(i)*surrogate[iion1].MM*conc_org/1000.*pion1+surrogate[iion1].Aaq2;
                //int iion2=surrogate[i].iion2;
                //int pion2=surrogate[i].pion2;
                surrogate[iion2].Aaq2=-excess(i)*surrogate[iion2].MM*conc_org/1000.*pion2+surrogate[iion2].Aaq2;
                if (surrogate[i].nion>2)
                  {
                    //int iion3=surrogate[i].iion3;
                    //int pion3=surrogate[i].pion3;
                    surrogate[iion3].Aaq2=-excess(i)*surrogate[iion3].MM*conc_org/1000.*pion3+surrogate[iion3].Aaq2;
                  }

		if (iion1==config.iNH4p or iion2==config.iNH4p or iion3==config.iNH4p)
		  {
		    total=surrogate[config.iNH3].Ag2/surrogate[config.iNH3].MM+surrogate[config.iNH4p].Aaq2/surrogate[config.iNH4p].MM;
		    surrogate[config.iNH3].Ag2=total*surrogate[config.iNH3].MM*(1.0-ratio(config.iNH4p));
		    surrogate[config.iNH4p].Aaq2=total*surrogate[config.iNH4p].MM*ratio(config.iNH4p);
		  }

		if (iion1==config.iNO3m or iion2==config.iNO3m or iion3==config.iNO3m)
		  {
		    total=surrogate[config.iHNO3].Ag2/surrogate[config.iHNO3].MM+surrogate[config.iNO3m].Aaq2/surrogate[config.iNO3m].MM;
		    surrogate[config.iHNO3].Ag2=total*surrogate[config.iHNO3].MM*(1.0-ratio(config.iNO3m));
		    surrogate[config.iNO3m].Aaq2=total*surrogate[config.iNO3m].MM*ratio(config.iNO3m);
		  }


                if (iion1==config.iClm or iion2==config.iClm or iion3==config.iClm)
		  {
		    total=surrogate[config.iHCl].Ag2/surrogate[config.iHCl].MM+surrogate[config.iClm].Aaq2/surrogate[config.iClm].MM;
		    surrogate[config.iHCl].Ag2=total*surrogate[config.iHCl].MM*(1.0-ratio(config.iClm));
		    surrogate[config.iClm].Aaq2=total*surrogate[config.iClm].MM*ratio(config.iClm);
		  }

		if (iion1==config.iSO4mm or iion1==config.iHSO4m or iion2==config.iSO4mm or iion2==config.iHSO4m or iion3==config.iSO4mm or iion3==config.iHSO4m)
		  {
		    total=surrogate[config.iSO4mm].Aaq2/surrogate[config.iSO4mm].MM+surrogate[config.iHSO4m].Aaq2/surrogate[config.iHSO4m].MM;
		    surrogate[config.iSO4mm].Aaq2=total*surrogate[config.iSO4mm].MM*ratio(config.iSO4mm);
		    surrogate[config.iHSO4m].Aaq2=total*surrogate[config.iHSO4m].MM*(1.-ratio(config.iSO4mm));
		  }
              }

  
            //cout << "avant " << min_excess << " " << max_excess << " " << 0.225*surrogate[42].Ap << " " << 0.225*surrogate[42].Ap+0.27272727*surrogate[43].Ap+surrogate[config.iNH3].Ag*18./17.+surrogate[config.iNH4p].Aaq << endl;

            //cout << "prendant " << surrogate[config.iNa].Aaq << " " << surrogate[44].Ap*2*23./surrogate[44].MM+surrogate[config.iNa].Aaq << endl;
        
            if (excess(i)>0.) //solidification
              {
                //cout << "solid " << surrogate[i].Ap << " " << factor << " " << surrogate[config.iH2O].Aaq << " " << surrogate[config.iH2O].gamma_aq << " " << config.iH2O << endl;
                //surrogate[i].Ap=max(factor*excess(i)*surrogate[i].MM*conc_org/1000.+surrogate[i].Ap,0.)
                if (surrogate[i].Ap2>0.) sum_error_ap+=abs(excess(i)*surrogate[i].MM*conc_org/1000./surrogate[i].Ap2);
                excess(i)*=factor2;
                surrogate[i].Ap2=excess(i)*surrogate[i].MM*conc_org/1000.+surrogate[i].Ap2;
                //int iion1=surrogate[i].iion1;
                //int pion1=surrogate[i].pion1;
                surrogate[iion1].Aaq2=-excess(i)*surrogate[iion1].MM*conc_org/1000.*pion1+surrogate[iion1].Aaq2;
                //int iion2=surrogate[i].iion2;
                //int pion2=surrogate[i].pion2;
                surrogate[iion2].Aaq2=-excess(i)*surrogate[iion2].MM*conc_org/1000.*pion2+surrogate[iion2].Aaq2;
                if (surrogate[i].nion>2)
                  {
                    //int iion3=surrogate[i].iion3;
                    //int pion3=surrogate[i].pion3;
                    surrogate[iion3].Aaq2=-excess(i)*surrogate[iion3].MM*conc_org/1000.*pion3+surrogate[iion3].Aaq2;
                  }

                //cout << "pendant " << min_excess << " " << max_excess << " " << 0.27272727*surrogate[43].Ap << " " << 0.27272727*surrogate[43].Ap+surrogate[config.iNH3].Ag*18./17.+surrogate[config.iNH4p].Aaq << endl;
                double total=0.0;
                //cout << surrogate[config.iNH3].Ag/surrogate[config.iNH3].MM+surrogate[config.iNH4p].Aaq/surrogate[config.iNH4p].MM << endl;
		if (iion1==config.iNH4p or iion2==config.iNH4p or iion3==config.iNH4p)
		  {
		    total=surrogate[config.iNH3].Ag2/surrogate[config.iNH3].MM+surrogate[config.iNH4p].Aaq2/surrogate[config.iNH4p].MM;
		    surrogate[config.iNH3].Ag2=max(total*surrogate[config.iNH3].MM*(1.0-ratio(config.iNH4p)),0.);
		    surrogate[config.iNH4p].Aaq2=total*surrogate[config.iNH4p].MM*ratio(config.iNH4p);
		  }
                //cout << surrogate[config.iNH3].Ag/surrogate[config.iNH3].MM << " " << surrogate[config.iNH4p].Aaq/surrogate[config.iNH4p].MM << endl;

                //cout << surrogate[config.iNH3].Ag/surrogate[config.iNH3].MM+surrogate[config.iNH4p].Aaq/surrogate[config.iNH4p].MM << endl;
		if (iion1==config.iNO3m or iion2==config.iNO3m or iion3==config.iNO3m)
		  {
		    total=surrogate[config.iHNO3].Ag2/surrogate[config.iHNO3].MM+surrogate[config.iNO3m].Aaq2/surrogate[config.iNO3m].MM;
		    surrogate[config.iHNO3].Ag2=max(total*surrogate[config.iHNO3].MM*(1.0-ratio(config.iNO3m)),0.);
		    surrogate[config.iNO3m].Aaq2=total*surrogate[config.iNO3m].MM*ratio(config.iNO3m);
		  }


                if (iion1==config.iClm or iion2==config.iClm or iion3==config.iClm)
		  {
		    total=surrogate[config.iHCl].Ag2/surrogate[config.iHCl].MM+surrogate[config.iClm].Aaq2/surrogate[config.iClm].MM;
		    surrogate[config.iHCl].Ag2=max(total*surrogate[config.iHCl].MM*(1.0-ratio(config.iClm)),0.);
		    surrogate[config.iClm].Aaq2=total*surrogate[config.iClm].MM*ratio(config.iClm);
		  }

		if (iion1==config.iSO4mm or iion1==config.iHSO4m or iion2==config.iSO4mm or iion2==config.iHSO4m or iion3==config.iSO4mm or iion3==config.iHSO4m)
		  {
		    total=surrogate[config.iSO4mm].Aaq2/surrogate[config.iSO4mm].MM+surrogate[config.iHSO4m].Aaq2/surrogate[config.iHSO4m].MM;
		    surrogate[config.iSO4mm].Aaq2=total*surrogate[config.iSO4mm].MM*ratio(config.iSO4mm);
		    surrogate[config.iHSO4m].Aaq2=total*surrogate[config.iHSO4m].MM*ratio(config.iHSO4m);
		  }
              }
   
            //cout << "fin " << min_excess << " " << max_excess << " " << 0.27272727*surrogate[43].Ap << " " << 0.27272727*surrogate[43].Ap+surrogate[config.iNH3].Ag*18./17.+surrogate[config.iNH4p].Aaq << endl;
            //cout << surrogate[i].Ap << endl;
            //cout << "la " << surrogate[i].Ap << " " << factor << " " << surrogate[config.iH2O].Aaq << " " << surrogate[config.iH2O].gamma_aq << " " << config.iH2O << endl;
        
            if (surrogate[i].Ap2<0. and surrogate[config.iH2O].Aaq2>1.0e4) exit(0);
            int j;
            for (j=0;j<n;++j)
              if (surrogate[j].is_organic==false and j!=config.iH2O and surrogate[j].is_inorganic_precursor==false)
                if (j!=config.iHp)
                  {

                    //if (surrogate[j].Aaq>0.) cout << surrogate[j].name << " " << surrogate[j].Aaq << endl;
                    surrogate[j].molality2=surrogate[j].Aaq2/surrogate[j].MM/conc_org*1000.;
                    //if (surrogate[j].molality<0.) cout << surrogate[j].name << " " << surrogate[j].molality << " " << excess(i) << endl;
                    //if (surrogate[j].molality<0.) exit(0);
                  }

            //cout << "apres " << surrogate[config.iNa].Aaq << " " << surrogate[44].Ap*2*23./surrogate[44].MM+surrogate[config.iNa].Aaq << endl;

          }
      //cout << sum_error_ap/factor2 << endl;
      //if ((sum_error_ap_save-sum_error_ap)<1.e-4*sum_error_ap) factor2=max(factor2/2,0.5);
      k++;
    }
  
  for (i=0;i<n;i++)
    {
      surrogate[i].Ag=factor*surrogate[i].Ag2+(1.-factor)*surrogate[i].Ag;
      surrogate[i].Ap=factor*surrogate[i].Ap2+(1.-factor)*surrogate[i].Ap;
      surrogate[i].Aaq=factor*surrogate[i].Aaq2+(1.-factor)*surrogate[i].Aaq;
      surrogate[i].molality=factor*surrogate[i].molality2+(1.-factor)*surrogate[i].molality;
    }
}

void error_inorg_aq_ssh(model_config &config, vector<species>& surrogate,
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
  double chp2=chp;
  double MOinit=0.0;
  double MOW=0.0;
  //MMaq: Mean molar mass of the aqueous phase
  
  AQinit=max(AQinit,config.MOmin);

  //if (config.compute_rho_aqueous)
  //compute_density_aqueous_phase_ssh(config, surrogate, LWC, Temperature);  
  double conc_org=LWC;
  for (i=0;i<n;++i)
    if (surrogate[i].is_organic or i==config.iH2O)      
      conc_org+=surrogate[i].Aaq;

  for (i=0;i<n;++i)
    if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)          
      conc_org=max(conc_org,surrogate[i].Aaq/surrogate[i].MM/config.molalmax*1000.0);
  
  conc_org=max(conc_org,1.e-5*config.MOmin); //0.00001*AQinit); //config.MOmin);
  
  compute_ionic_strenght2_ssh(config, surrogate, Temperature, AQinit, conc_inorganic, ionic, chp, organion, ionic_organic, conc_org,factor);

  //initialize AQ
  AQ=0.0;
  
  if (config.solids)
    {
      solidification_ssh(config,surrogate,conc_org,factor);
      conc_org=LWC;
      for (i=0;i<n;++i)
        if (surrogate[i].is_organic or i==config.iH2O)      
          conc_org+=surrogate[i].Aaq;
      
      for (i=0;i<n;++i)
        if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)          
	  conc_org=max(conc_org,surrogate[i].Aaq/surrogate[i].MM/config.molalmax*1000.0);
      conc_org=max(conc_org,1.e-5*config.MOmin); //0.00001*AQinit); //config.MOmin);
    }
      
  //compute acitivity coefficients and MMaq  
  if (compute_activity_coefficients)
    {            
      activity_coefficients_aq_ssh(config, surrogate, Temperature, LWC, MMaq, XH2O, conc_org);     
      if (config.compute_long_and_medium_range_interactions)
        activity_coefficients_LR_MR_ssh(config, surrogate, Temperature, LWC, ionic);     
    }
  
  for (i=0;i<n;++i)
    if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)
      {
        surrogate[i].gamma_aq=max(surrogate[i].gamma_aq,0.1*surrogate[i].gamma_aq_old);
        surrogate[i].gamma_aq=min(surrogate[i].gamma_aq,10.*surrogate[i].gamma_aq_old);
      }
  
  //cout << surrogate[config.iSO4mm].gamma_aq << " " << surrogate[config.iSO4mm].gamma_aq_old << endl;

  /*
    cout << surrogate[config.iH2O].Aaq << " " << surrogate[config.iH2O].gamma_aq << endl;
    for (i=0;i<n;i++)
    if (surrogate[i].Aaq>0.)
    cout << surrogate[i].name << " " << surrogate[i].Aaq << " " << surrogate[i].gamma_aq << endl;*/
  
  //pH computation 
  
  if (config.compute_aqueous_phase_properties)
    {
      //If inorganic ion concentrations are computed by SOAP, used a method of
      //newton raphson

      if (config.compute_inorganic)
        {
          double error_h=1000.0;
          double derivative_h;
          chp2=max(chp,1.0e-15);
          int index=0;          
          while(abs(error_h/chp2)>1.0e-3 and index<20)
            {
              Kpreal_inorganic_ssh(config, surrogate, chp2);
              error_ph_ssh(config, surrogate, Temperature, chp2, organion, error_h, derivative_h,AQinit,LWC,MMaq,MOinit,MOW,conc_org);
              if (chp2-error_h/derivative_h>0.0 and derivative_h!=0.0)
                chp2=chp2-error_h/derivative_h;
              else
                chp2=chp2+error_h;

              //Too high ph may cause instability
              chp2=max(chp2,1.0e-8);	                    
              index++;
            }              
          chp2=min(10*chp,max(0.1*chp,chp2));
          //Too low ph may cause instability
          chp2=min(chp2,30./surrogate[config.iHp].gamma_aq);
          chp=max(factor*chp2+(1.0-factor)*chp,1.e-15);  
          Kpreal_inorganic_ssh(config, surrogate, chp);
          /*
          if (config.solids)
          solidification_ssh(config,surrogate,conc_org,factor);*/

        }
      else
        {
          double inorganion=0.0;
          for (i=0;i<n;++i)
            if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)
              if (i!=config.iHp)
                inorganion-=surrogate[i].molality*surrogate[i].charge;
      
          chp=max(factor*0.5*(organion+inorganion+pow(pow(organion+inorganion,2)+4*config.Ke,0.5))+(1.0-factor)*chp,1.e-15);
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
              Kp=surrogate[i].Kp_eff_aqreal_ssh(config, Temperature, ionic, chp,surrogate[config.iHp].gamma_LR,
                                            surrogate[config.iHp].gamma_SRMR,MMaq,fion1,fion2)
                /surrogate[i].gamma_aq;
              surrogate[i].Aaq=factor*surrogate[i].Atot*Kp*AQinit/(1+Kp*AQinit)+
                (1.0-factor)*surrogate[i].Aaq;
              AQ+=surrogate[i].Aaq;
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

  

  //compute the concentrations of inorganics
  double total;
  //Kpreal_inorganic_ssh(config, surrogate, chp);

  //sulfate	    
  i=config.iH2SO4;
  total=surrogate[i].Ag+surrogate[config.iHSO4m].Aaq/surrogate[config.iHSO4m].MM*surrogate[i].MM
    +surrogate[config.iSO4mm].Aaq/surrogate[config.iSO4mm].MM*surrogate[i].MM; //gas+particle concentrations
  double Keq=surrogate[i].keq/chp*surrogate[config.iHSO4m].gamma_aq/surrogate[config.iHp].gamma_aq/surrogate[config.iSO4mm].gamma_aq;
  surrogate[i].Ag=(1.0-factor)*surrogate[i].Ag;
  surrogate[config.iHSO4m].Aaq=factor*total*surrogate[config.iHSO4m].MM/surrogate[i].MM/(1.0+Keq)
    +(1.0-factor)*surrogate[config.iHSO4m].Aaq; //HSO4-
  surrogate[config.iSO4mm].Aaq=factor*total*surrogate[config.iSO4mm].MM/surrogate[i].MM*Keq/
    (1.0+Keq)+(1.0-factor)*surrogate[config.iSO4mm].Aaq; //SO4--               
            
  //NH3
  i=config.iNH3;
  total=surrogate[i].Ag+surrogate[config.iNH4p].Aaq/surrogate[config.iNH4p].MM*surrogate[i].MM; //gas+particle concentrations
  surrogate[i].Ag=factor*total/(1.0+surrogate[i].Kaq_inorg*conc_org)+(1.0-factor)*surrogate[i].Ag; 
  surrogate[config.iNH4p].Aaq=factor*total*surrogate[i].Kaq_inorg*conc_org/(1.0+surrogate[i].Kaq_inorg*conc_org)*surrogate[config.iNH4p].MM/surrogate[i].MM
    +(1.0-factor)*surrogate[config.iNH4p].Aaq; //NH4+                           

  //HNO3
  i=config.iHNO3;
  total=surrogate[i].Ag+surrogate[config.iNO3m].Aaq/surrogate[config.iNO3m].MM*surrogate[i].MM; //gas+particle concentrations
  surrogate[i].Ag=factor*total/(1.0+surrogate[i].Kaq_inorg*conc_org)+(1.0-factor)*surrogate[i].Ag;
  //cout << "NO3 " << surrogate[i].Kaq_inorg << endl;
  surrogate[config.iNO3m].Aaq=factor*total*surrogate[i].Kaq_inorg*conc_org/(1.0+surrogate[i].Kaq_inorg*conc_org)*surrogate[config.iNO3m].MM/surrogate[i].MM
    +(1.0-factor)*surrogate[config.iNO3m].Aaq; //NO3-
      
  //HCl
  i=config.iHCl;
  total=surrogate[i].Ag+surrogate[config.iClm].Aaq/surrogate[config.iClm].MM*surrogate[i].MM; //gas+particle concentrations
  surrogate[i].Ag=factor*total/(1.0+surrogate[i].Kaq_inorg*conc_org)+(1.0-factor)*surrogate[i].Ag;
  surrogate[config.iClm].Aaq=factor*total*surrogate[i].Kaq_inorg*conc_org/(1.0+surrogate[i].Kaq_inorg*conc_org)*surrogate[config.iClm].MM/surrogate[i].MM
    +(1.0-factor)*surrogate[config.iClm].Aaq; //Cl-

  conc_inorganic=0.0;
  for (i=0;i<n;i++)
    if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)
      conc_inorganic+=surrogate[i].Aaq; //compute total inorganic concentration

  AQ+=conc_inorganic;  
  if (config.iH2O>=0)
    //compute hygroscopicity
    {      
      LWC=0.0; //In the case where SOAP compute inorganics, LWC is only used for the initialization
      hygroscopicity_tot_ssh(config, surrogate, Temperature, RH, AQinit, conc_inorganic,MMaq,AQ,derivative,factor);
    }



  AQ=max(AQ,config.MOmin);
  error=AQinit-AQ;
  ionic_organic=ionic_organic_tmp;
  organion=organion_tmp;  
}


void error_aq_ssh(model_config &config, vector<species>& surrogate,
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
  double chp2=chp;
  double MOinit=0.0;
  double MOW=0.0;
  //MMaq: Mean molar mass of the aqueous phase
  
  AQinit=max(AQinit,config.MOmin);

  //if (config.compute_rho_aqueous)
  //compute_density_aqueous_phase_ssh(config, surrogate, LWC, Temperature);  
  double conc_org=LWC;
  for (i=0;i<n;++i)
    if (surrogate[i].is_organic or i==config.iH2O)      
      conc_org+=surrogate[i].Aaq;
  
  for (i=0;i<n;++i)
    if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)          
      conc_org=max(conc_org,surrogate[i].Aaq/surrogate[i].MM/config.molalmax*1000.0);    
  conc_org=max(conc_org,1.e-5*config.MOmin); //0.00001*AQinit); //config.MOmin);
  
  compute_ionic_strenght2_ssh(config, surrogate, Temperature, AQinit, conc_inorganic, ionic, chp, organion, ionic_organic, conc_org,factor);

  //initialize AQ
  AQ=0.0;
  
  if (config.solids)
    {
      solidification_ssh(config,surrogate,conc_org,factor);
      conc_org=LWC;
      for (i=0;i<n;++i)
        if (surrogate[i].is_organic or i==config.iH2O)      
          conc_org+=surrogate[i].Aaq;

      for (i=0;i<n;++i)
        if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)          
	  conc_org=max(conc_org,surrogate[i].Aaq/surrogate[i].MM/config.molalmax*1000.0);
      conc_org=max(conc_org,1.e-5*config.MOmin); //0.00001*AQinit); //config.MOmin);
    }
      
  //compute acitivity coefficients and MMaq  
  if (compute_activity_coefficients)
    {            
      activity_coefficients_aq_ssh(config, surrogate, Temperature, LWC, MMaq, XH2O, conc_org);     
      if (config.compute_long_and_medium_range_interactions)
        activity_coefficients_LR_MR_ssh(config, surrogate, Temperature, LWC, ionic);     
    }

  if (config.compute_inorganic)
    for (i=0;i<n;++i)
      if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)
        {
          surrogate[i].gamma_aq=max(surrogate[i].gamma_aq,0.1*surrogate[i].gamma_aq_old);
          surrogate[i].gamma_aq=min(surrogate[i].gamma_aq,10.*surrogate[i].gamma_aq_old);
        }
  
  //cout << surrogate[config.iSO4mm].gamma_aq << " " << surrogate[config.iSO4mm].gamma_aq_old << endl;

  /*
    cout << surrogate[config.iH2O].Aaq << " " << surrogate[config.iH2O].gamma_aq << endl;
    for (i=0;i<n;i++)
    if (surrogate[i].Aaq>0.)
    cout << surrogate[i].name << " " << surrogate[i].Aaq << " " << surrogate[i].gamma_aq << endl;*/
  
  //pH computation 
  
  if (config.compute_aqueous_phase_properties)
    {
      //If inorganic ion concentrations are computed by SOAP, used a method of
      //newton raphson

      if (config.compute_inorganic)
        {
          double error_h=1000.0;
          double derivative_h;
          chp2=max(chp,1.0e-15);
          int index=0;          
          while(abs(error_h/chp2)>1.0e-3 and index<20)
            {
              Kpreal_inorganic_ssh(config, surrogate, chp2);
              error_ph_ssh(config, surrogate, Temperature, chp2, organion, error_h, derivative_h,AQinit,LWC,MMaq,MOinit,MOW,conc_org);
              if (chp2-error_h/derivative_h>0.0 and derivative_h!=0.0)
                chp2=chp2-error_h/derivative_h;
              else
                chp2=chp2+error_h;

              //Too high ph may cause instability
              chp2=max(chp2,1.0e-8);	                    
              index++;
            }              
          chp2=min(10*chp,max(0.1*chp,chp2));
          //Too low ph may cause instability
          if (surrogate[config.iHp].gamma_aq == 0.0)
            throw string("Error: zero division in error_aq.");
          chp2=min(chp2,30./surrogate[config.iHp].gamma_aq);
          chp=max(factor*chp2+(1.0-factor)*chp,1.e-15);  
          Kpreal_inorganic_ssh(config, surrogate, chp);
          /*
          if (config.solids)
          solidification_ssh(config,surrogate,conc_org,factor);*/

        }
      else
        {
          double inorganion=0.0;
          for (i=0;i<n;++i)
            if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)
              if (i!=config.iHp)
                inorganion-=surrogate[i].molality*surrogate[i].charge;
      
          chp=max(factor*0.5*(organion+inorganion+pow(pow(organion+inorganion,2)+4*config.Ke,0.5))+(1.0-factor)*chp,1.0e-15);
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

              if (surrogate[i].gamma_aq == 0.0)
                throw string("Error: zero division gamma_aq in error_aq.");
              Kp=surrogate[i].Kp_eff_aqreal_ssh(config, Temperature, ionic, chp,surrogate[config.iHp].gamma_LR,
                                            surrogate[config.iHp].gamma_SRMR,MMaq,fion1,fion2)
                /surrogate[i].gamma_aq;

              surrogate[i].Aaq=factor*surrogate[i].Atot*Kp*AQinit/(1+Kp*AQinit)+
                (1.0-factor)*surrogate[i].Aaq;
              AQ+=surrogate[i].Aaq;
              derivative+=surrogate[i].Atot*pow(Kp,2)*AQinit/(pow(1+Kp*AQinit,2))
                -surrogate[i].Atot*Kp/(1+Kp*AQinit);
              //molality1: molality of ions HA- or A-
              if (conc_org == 0.0)
                throw string("Error: division by zero conc_org.");
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
      //Kpreal_inorganic_ssh(config, surrogate, chp);

      //sulfate	    
      i=config.iH2SO4;
      total=surrogate[i].Ag+surrogate[config.iHSO4m].Aaq/surrogate[config.iHSO4m].MM*surrogate[i].MM
        +surrogate[config.iSO4mm].Aaq/surrogate[config.iSO4mm].MM*surrogate[i].MM; //gas+particle concentrations
      double Keq=surrogate[i].keq/chp*surrogate[config.iHSO4m].gamma_aq/surrogate[config.iHp].gamma_aq/surrogate[config.iSO4mm].gamma_aq;
      surrogate[i].Ag=(1.0-factor)*surrogate[i].Ag;
      surrogate[config.iHSO4m].Aaq=factor*total*surrogate[config.iHSO4m].MM/surrogate[i].MM/(1.0+Keq)
        +(1.0-factor)*surrogate[config.iHSO4m].Aaq; //HSO4-
      surrogate[config.iSO4mm].Aaq=factor*total*surrogate[config.iSO4mm].MM/surrogate[i].MM*Keq/
        (1.0+Keq)+(1.0-factor)*surrogate[config.iSO4mm].Aaq; //SO4--               
            
      //NH3
      i=config.iNH3;
      total=surrogate[i].Ag+surrogate[config.iNH4p].Aaq/surrogate[config.iNH4p].MM*surrogate[i].MM; //gas+particle concentrations
      surrogate[i].Ag=factor*total/(1.0+surrogate[i].Kaq_inorg*conc_org)+(1.0-factor)*surrogate[i].Ag; 
      surrogate[config.iNH4p].Aaq=factor*total*surrogate[i].Kaq_inorg*conc_org/(1.0+surrogate[i].Kaq_inorg*conc_org)*surrogate[config.iNH4p].MM/surrogate[i].MM
        +(1.0-factor)*surrogate[config.iNH4p].Aaq; //NH4+                
      derivative+=total*surrogate[config.iNH4p].MM/surrogate[i].MM*pow(surrogate[i].Kaq_inorg,2)*conc_org/(pow(1.+surrogate[i].Kaq_inorg*conc_org,2))
        -total*surrogate[config.iNH4p].MM/surrogate[i].MM*surrogate[i].Kaq_inorg/(1.+surrogate[i].Kaq_inorg*conc_org);                

      //HNO3
      i=config.iHNO3;
      total=surrogate[i].Ag+surrogate[config.iNO3m].Aaq/surrogate[config.iNO3m].MM*surrogate[i].MM; //gas+particle concentrations
      surrogate[i].Ag=factor*total/(1.0+surrogate[i].Kaq_inorg*conc_org)+(1.0-factor)*surrogate[i].Ag;
      //cout << "NO3 " << surrogate[i].Kaq_inorg << endl;
      surrogate[config.iNO3m].Aaq=factor*total*surrogate[i].Kaq_inorg*conc_org/(1.0+surrogate[i].Kaq_inorg*conc_org)*surrogate[config.iNO3m].MM/surrogate[i].MM
        +(1.0-factor)*surrogate[config.iNO3m].Aaq; //NO3-
      derivative+=total*surrogate[config.iNO3m].MM/surrogate[i].MM*pow(surrogate[i].Kaq_inorg,2)*conc_org/(pow(1.+surrogate[i].Kaq_inorg*conc_org,2))
        -total*surrogate[config.iNO3m].MM/surrogate[i].MM*surrogate[i].Kaq_inorg/(1.+surrogate[i].Kaq_inorg*conc_org);

      //HCl
      i=config.iHCl;
      total=surrogate[i].Ag+surrogate[config.iClm].Aaq/surrogate[config.iClm].MM*surrogate[i].MM; //gas+particle concentrations
      surrogate[i].Ag=factor*total/(1.0+surrogate[i].Kaq_inorg*conc_org)+(1.0-factor)*surrogate[i].Ag;
      surrogate[config.iClm].Aaq=factor*total*surrogate[i].Kaq_inorg*conc_org/(1.0+surrogate[i].Kaq_inorg*conc_org)*surrogate[config.iClm].MM/surrogate[i].MM
        +(1.0-factor)*surrogate[config.iClm].Aaq; //Cl-
      derivative+=total*surrogate[config.iClm].MM/surrogate[i].MM*pow(surrogate[i].Kaq_inorg,2)*conc_org/(pow(1.+surrogate[i].Kaq_inorg*conc_org,2))
        -total*surrogate[config.iClm].MM/surrogate[i].MM*surrogate[i].Kaq_inorg/(1.+surrogate[i].Kaq_inorg*conc_org);

      conc_inorganic=0.0;
      for (i=0;i<n;i++)
        if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)
          conc_inorganic+=surrogate[i].Aaq; //compute total inorganic concentration

      AQ+=conc_inorganic;
      if (config.iH2O>=0)
        //compute hygroscopicity
        {
          LWC=0.0; //In the case where SOAP compute inorganics, LWC is only used for the initialization
          hygroscopicity_tot_ssh(config, surrogate, Temperature, RH, AQinit, conc_inorganic,MMaq,AQ,derivative,factor);
        }

    }
  else
    {
      AQ+=LWC+conc_inorganic;  
      if (config.iH2O>=0) //compute hygroscopicity
        if (surrogate[config.iH2O].hydrophilic and config.hygroscopicity) //Can H2O condense on the organic phase
          hygroscopicity_aq_ssh(config, surrogate, Temperature, RH, AQinit, LWC,
                            conc_inorganic,MMaq,AQ,derivative,factor);
    }
  
  AQ=max(AQ,config.MOmin);
  error=AQinit-AQ;
  ionic_organic=ionic_organic_tmp;
  organion=organion_tmp;  
}

void error_coupled_ssh(model_config &config, vector<species>& surrogate,
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

  //if (config.compute_rho_aqueous)
  //compute_density_aqueous_phase_ssh(config, surrogate, LWC, Temperature);

  double conc_org=LWC;
  for (i=0;i<n;++i)
    if (surrogate[i].is_organic or i==config.iH2O)      
      conc_org+=surrogate[i].Aaq;

  for (i=0;i<n;++i)
    if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)          
      conc_org=max(conc_org,surrogate[i].Aaq/surrogate[i].MM/config.molalmax*1000.0);
  conc_org=max(conc_org,1.e-5*config.MOmin);
  //conc_org=max(conc_org,config.MOmin);
  
  compute_ionic_strenght2_ssh(config, surrogate, Temperature, AQinit, conc_inorganic, ionic, chp2,
                          organion, ionic_organic, conc_org, factor);

  if (config.solids)
    {
      solidification_ssh(config,surrogate,conc_org,factor);
      conc_org=LWC;
      for (i=0;i<n;++i)
        if (surrogate[i].is_organic or i==config.iH2O)      
          conc_org+=surrogate[i].Aaq;

      for (i=0;i<n;++i)
        if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)          
	  conc_org=max(conc_org,surrogate[i].Aaq/surrogate[i].MM/config.molalmax*1000.0);      
      conc_org=max(conc_org,1.e-5*config.MOmin); //0.00001*AQinit); //config.MOmin);
    }
  
  //initialize AQ and MO
  AQ=0.0;
  MO=0.0;
  deriv_error1_MO=1.0;
  deriv_error1_AQ=0.0;
  deriv_error2_MO=0.0;
  deriv_error2_AQ=1.0;

  AQinit=max(AQinit,config.MOmin);
  MOinit=max(MOinit,config.MOmin);

  //cout << " ici conc_org " << conc_org << endl;

  //compute acitivity coefficients, MOW and MMaq
  if (compute_activity_coefficients)
    {
      if (config.compute_organic)
        activity_coefficients_org_ssh(config, surrogate, false, Temperature, MOW);
      activity_coefficients_aq_ssh(config, surrogate, Temperature, LWC, MMaq, XH2O,conc_org);
      if (config.compute_long_and_medium_range_interactions)
        activity_coefficients_LR_MR_ssh(config, surrogate, Temperature, LWC, ionic);
    }

  //Prevent strong variations of activity coefficients for inorganic ions
  /*
  for (i=0;i<n;++i)
    if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)
      {
        surrogate[i].gamma_aq=max(surrogate[i].gamma_aq,0.1*surrogate[i].gamma_aq_old);
        surrogate[i].gamma_aq=min(surrogate[i].gamma_aq,10.*surrogate[i].gamma_aq_old);
      }*/
  
  /* 
     if (MMaq==0.0)
     cout << "MMaq is zero" << endl;*/
  
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
          while(abs(error_h/chp2)>1.0e-3 and index<20)
            {
              Kpreal_inorganic_ssh(config, surrogate, chp2);
              error_ph_ssh(config, surrogate, Temperature, chp2, organion, error_h, derivative_h,AQinit,LWC,MMaq,MOinit,MOW,conc_org);
              if (chp2-error_h/derivative_h>0.0 and derivative_h!=0.0)
                chp2=chp2-error_h/derivative_h;
              else
                chp2=chp2+error_h;
              
              chp2=max(chp2,1.0e-8);	                    
              index++;
            }          
          chp2=min(10*chp,max(0.1*chp,chp2));
          //Too low ph may cause instability
          chp2=min(chp2,30./surrogate[config.iHp].gamma_aq);
          chp=max(factor*chp2+(1.0-factor)*chp,1.0e-15);          
          Kpreal_inorganic_ssh(config, surrogate, chp);
        }
      else
        {
          double inorganion=0.0;
          for (i=0;i<n;++i)
            if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)
              if (i!=config.iHp)
                inorganion-=surrogate[i].molality*surrogate[i].charge;
      
          chp=max(factor*0.5*(organion+inorganion+pow(pow(organion+inorganion,2)+4*config.Ke,0.5))+(1.0-factor)*chp,1.0e-15);
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
              Kp=surrogate[i].Kp_eff_aqreal_ssh(config, Temperature, ionic, chp,surrogate[config.iHp].gamma_LR,
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
                  (1.0+MOinit/AQinit*MMaq/MOW*surrogate[i].gamma_aq/surrogate[i].gamma_org)+(1.0-factor)*surrogate[i].Aaq;
                surrogate[i].Ap=surrogate[i].Atot-surrogate[i].Aaq;
                MO+=surrogate[i].Ap;
                AQ+=surrogate[i].Aaq;
                deriv_error1_MO+=-surrogate[i].Atot/
                  (pow(1+MOinit/AQinit*MMaq/MOW*surrogate[i].gamma_aq/surrogate[i].gamma_org,2))
                  *(MMaq/MOW*surrogate[i].gamma_aq/surrogate[i].gamma_org)/AQinit;
                deriv_error1_AQ+=surrogate[i].Atot/
                  (pow(1+MOinit/AQinit*MMaq/MOW*surrogate[i].gamma_aq/surrogate[i].gamma_org,2))
                  *(MMaq/MOW*surrogate[i].gamma_aq/surrogate[i].gamma_org)*MOinit/(pow(AQinit,2));
                deriv_error2_AQ+=-surrogate[i].Atot/
                  (pow(1+MOinit/AQinit*MMaq/MOW*surrogate[i].gamma_aq/surrogate[i].gamma_org,2))
                  *(MMaq/MOW*surrogate[i].gamma_aq/surrogate[i].gamma_org)*MOinit/(pow(AQinit,2));
                deriv_error2_MO+=surrogate[i].Atot/
                  (pow(1+MOinit/AQinit*MMaq/MOW*surrogate[i].gamma_aq/surrogate[i].gamma_org,2))
                  *(MMaq/MOW*surrogate[i].gamma_aq/surrogate[i].gamma_org)/AQinit;
              }
          else
            {
              fion1=0.0;
              fion2=0.0;
              Kp_aq=surrogate[i].Kp_eff_aqreal_ssh(config, Temperature, ionic, chp,
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

      i=config.iH2SO4; //sulfate
      total=surrogate[i].Ag+surrogate[config.iHSO4m].Aaq/surrogate[config.iHSO4m].MM*surrogate[i].MM
        +surrogate[config.iSO4mm].Aaq/surrogate[config.iSO4mm].MM*surrogate[i].MM; //gas+particle concentrations
      double Keq=surrogate[i].keq/chp*surrogate[config.iHSO4m].gamma_aq/surrogate[config.iHp].gamma_aq/surrogate[config.iSO4mm].gamma_aq;
      surrogate[i].Ag=(1.0-factor)*surrogate[i].Ag;
      surrogate[config.iHSO4m].Aaq=factor*total*surrogate[config.iHSO4m].MM/surrogate[i].MM/(1.0+Keq)
        +(1.0-factor)*surrogate[config.iHSO4m].Aaq; //HSO4-
      surrogate[config.iSO4mm].Aaq=factor*total*surrogate[config.iSO4mm].MM/surrogate[i].MM*Keq/
        (1.0+Keq)+(1.0-factor)*surrogate[config.iSO4mm].Aaq; //SO4--
 
      i=config.iNH3;//ammoniac	      
      total=surrogate[i].Ag+surrogate[config.iNH4p].Aaq/surrogate[config.iNH4p].MM*surrogate[i].MM; //gas+particle concentrations
      surrogate[i].Ag=factor*total/(1.0+surrogate[i].Kaq_inorg*conc_org)+(1.0-factor)*surrogate[i].Ag; 
      surrogate[config.iNH4p].Aaq=factor*total*surrogate[i].Kaq_inorg*conc_org/(1.0+surrogate[i].Kaq_inorg*conc_org)*surrogate[config.iNH4p].MM/surrogate[i].MM
        +(1.0-factor)*surrogate[config.iNH4p].Aaq; //NH4+                
      deriv_error2_AQ+=total*surrogate[config.iNH4p].MM/surrogate[i].MM*pow(surrogate[i].Kaq_inorg,2)*conc_org/(pow(1.+surrogate[i].Kaq_inorg*conc_org,2))
        -total*surrogate[config.iNH4p].MM/surrogate[i].MM*surrogate[i].Kaq_inorg/(1.+surrogate[i].Kaq_inorg*conc_org);

      i=config.iHNO3; //nitrate	      
      total=surrogate[i].Ag+surrogate[config.iNO3m].Aaq/surrogate[config.iNO3m].MM*surrogate[i].MM; //gas+particle concentrations
      surrogate[i].Ag=factor*total/(1.0+surrogate[i].Kaq_inorg*conc_org)+(1.0-factor)*surrogate[i].Ag;
      //cout << "NO3 " << surrogate[i].Kaq_inorg << endl;
      surrogate[config.iNO3m].Aaq=factor*total*surrogate[i].Kaq_inorg*conc_org/(1.0+surrogate[i].Kaq_inorg*conc_org)*surrogate[config.iNO3m].MM/surrogate[i].MM
        +(1.0-factor)*surrogate[config.iNO3m].Aaq; //NO3-
      deriv_error2_AQ+=total*surrogate[config.iNO3m].MM/surrogate[i].MM*pow(surrogate[i].Kaq_inorg,2)*conc_org/(pow(1.+surrogate[i].Kaq_inorg*conc_org,2))
        -total*surrogate[config.iNO3m].MM/surrogate[i].MM*surrogate[i].Kaq_inorg/(1.+surrogate[i].Kaq_inorg*conc_org);

      i=config.iHCl; //chloride	      
      total=surrogate[i].Ag+surrogate[config.iClm].Aaq/surrogate[config.iClm].MM*surrogate[i].MM; //gas+particle concentrations
      surrogate[i].Ag=factor*total/(1.0+surrogate[i].Kaq_inorg*conc_org)+(1.0-factor)*surrogate[i].Ag;
      surrogate[config.iClm].Aaq=factor*total*surrogate[i].Kaq_inorg*conc_org/(1.0+surrogate[i].Kaq_inorg*conc_org)*surrogate[config.iClm].MM/surrogate[i].MM
        +(1.0-factor)*surrogate[config.iClm].Aaq; //Cl-
      deriv_error2_AQ+=total*surrogate[config.iClm].MM/surrogate[i].MM*pow(surrogate[i].Kaq_inorg,2)*conc_org/(pow(1.+surrogate[i].Kaq_inorg*conc_org,2))
        -total*surrogate[config.iClm].MM/surrogate[i].MM*surrogate[i].Kaq_inorg/(1.+surrogate[i].Kaq_inorg*conc_org);

      for (i=0;i<n;i++)
        if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)
          conc_inorganic+=surrogate[i].Aaq; //compute total inorganic concentrations

      AQ+=conc_inorganic;
      if (config.iH2O>=0)
        {
          //compute hygroscopicity
          LWC=0.0; //In the case where SOAP compute inorganics, LWC is only used for the initialization
          if (config.hygroscopicity)
	    hygroscopicity_coupled_ssh(config, surrogate, Temperature, MOW, MMaq, RH, LWC, conc_inorganic,
				       MOinit, MO, AQinit, AQ, deriv_error1_MO, 
				       deriv_error1_AQ, deriv_error2_MO, deriv_error2_AQ, factor);
	  
          //hygroscopicity_tot_ssh(config, surrogate, Temperature, RH, AQinit, conc_inorganic,MMaq,AQ,deriv_error2_AQ,factor);
          //if (config.hygroscopicity and config.compute_organic)
          //  if (surrogate[config.iH2O].hydrophobic) //Can H2O condense on the organic phase
          //  hygroscopicity_org_ssh(config, surrogate, MOW, RH, MOinit, MO, deriv_error1_MO,factor);
        }
    }
  else
    {
      AQ+=LWC+conc_inorganic;  
      if (config.iH2O>=0 and config.hygroscopicity)        
        if (config.hygroscopicity)
	  {
	    if (surrogate[config.iH2O].hydrophobic) //Can H2O condense on the organic phase
	      hygroscopicity_org_ssh(config, surrogate, Temperature, MOW, RH, MOinit, MO, deriv_error1_MO, factor);
	    if (surrogate[config.iH2O].hydrophilic) //Can H2O condense on the organic phase
	      hygroscopicity_aq_ssh(config, surrogate, Temperature, RH, AQinit, LWC,
				    conc_inorganic,MMaq,AQ,deriv_error2_AQ,factor);
	    /*
	    hygroscopicity_coupled_tot_ssh(config, surrogate, Temperature, MOW, MMaq, RH, LWC, conc_inorganic,
					   MOinit, MO, AQinit, AQ, deriv_error1_MO, 
					   deriv_error1_AQ, deriv_error2_MO, deriv_error2_AQ, factor);       */
	  }
    }  
  /*
    for (i=0;i<n;i++)
    cout << surrogate[i].name << " " << surrogate[i].Aaq << endl;*/
  
  AQ=max(AQ,config.MOmin);
  MO=max(MO,config.MOmin);
  
  
  error1=MOinit-MO;
  error2=AQinit-AQ;
  ionic_organic=ionic_organic_tmp;
  organion=organion_tmp;

}

void error_coupled_inorg_ssh(model_config &config, vector<species>& surrogate,
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

  //if (config.compute_rho_aqueous)
  //compute_density_aqueous_phase(config, surrogate, LWC, Temperature);
  double conc_org=LWC;
  for (i=0;i<n;++i)
    if (surrogate[i].is_organic or i==config.iH2O)      
      conc_org+=surrogate[i].Aaq;

  for (i=0;i<n;++i)
    if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)          
      conc_org=max(conc_org,surrogate[i].Aaq/surrogate[i].MM/config.molalmax*1000.0);

  conc_org=max(conc_org,1.e-5*config.MOmin);
  //conc_org=max(conc_org,config.MOmin);
  
  compute_ionic_strenght2_ssh(config, surrogate, Temperature, AQinit, conc_inorganic, ionic, chp2,
                          organion, ionic_organic, conc_org, factor);
  if (config.solids)
    {
      solidification_ssh(config,surrogate,conc_org,factor);
      conc_org=LWC;
      for (i=0;i<n;++i)
        if (surrogate[i].is_organic or i==config.iH2O)      
          conc_org+=surrogate[i].Aaq;

      for (i=0;i<n;++i)
        if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)          
	  conc_org=max(conc_org,surrogate[i].Aaq/surrogate[i].MM/config.molalmax*1000.0);
      conc_org=max(conc_org,1.e-5*config.MOmin); //0.00001*AQinit); //config.MOmin);
    }
      

  //initialize AQ and MO
  AQ=0.0;
  MO=0.0;
  deriv_error1_MO=1.0;
  deriv_error1_AQ=0.0;
  deriv_error2_MO=0.0;
  deriv_error2_AQ=1.0;

  AQinit=max(AQinit,config.MOmin);
  MOinit=max(MOinit,config.MOmin);
 
  //cout << " ici conc_org " << conc_org << endl;

  //compute acitivity coefficients, MOW and MMaq
  if (compute_activity_coefficients)
    {
      if (config.compute_organic)
        activity_coefficients_org_ssh(config, surrogate, false, Temperature, MOW);
      activity_coefficients_aq_ssh(config, surrogate, Temperature, LWC, MMaq, XH2O, conc_org);
      if (config.compute_long_and_medium_range_interactions)
        activity_coefficients_LR_MR_ssh(config, surrogate, Temperature, LWC, ionic);
    }
  
  for (i=0;i<n;++i)
    if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)
      {
        surrogate[i].gamma_aq=max(surrogate[i].gamma_aq,0.1*surrogate[i].gamma_aq_old);
        surrogate[i].gamma_aq=min(surrogate[i].gamma_aq,10.*surrogate[i].gamma_aq_old);
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
          while(abs(error_h/chp2)>1.0e-3 and index<20)
            {
              Kpreal_inorganic_ssh(config, surrogate, chp2);
              error_ph_ssh(config, surrogate, Temperature, chp2, organion, error_h, derivative_h,AQinit,LWC,MMaq,MOinit,MOW,conc_org);
              if (chp2-error_h/derivative_h>0.0 and derivative_h!=0.0)
                chp2=chp2-error_h/derivative_h;
              else
                chp2=chp2+error_h;

              //Too high ph may cause instability
              chp2=max(chp2,1.0e-8);	 
              //chp2=max(chp2,1.0e-20);	                    
              index++;
            }          
          chp2=min(10*chp,max(0.1*chp,chp2));
          //Too low ph may cause instability
          chp2=min(chp2,30./surrogate[config.iHp].gamma_aq);
          chp=max(factor*chp2+(1.0-factor)*chp,1.0e-15);          
          Kpreal_inorganic_ssh(config, surrogate, chp);
        }
      else
        {
          double inorganion=0.0;
          for (i=0;i<n;++i)
            if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)
              if (i!=config.iHp)
                inorganion-=surrogate[i].molality*surrogate[i].charge;
      
          chp=max(factor*0.5*(organion+inorganion+pow(pow(organion+inorganion,2)+4*config.Ke,0.5))+(1.0-factor)*chp,1.0e-15);
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
              }
            else
              {
                Kp=surrogate[i].kpi/MOW/surrogate[i].gamma_org;
                surrogate[i].Aaq=0.0;
                surrogate[i].Ap=factor*surrogate[i].Atot*Kp*MOinit/(1+Kp*MOinit)+(1.0-factor)*surrogate[i].Ap;
                MO+=surrogate[i].Ap;
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
              Kp=surrogate[i].Kp_eff_aqreal_ssh(config, Temperature, ionic, chp,surrogate[config.iHp].gamma_LR,
                                            surrogate[config.iHp].gamma_SRMR,MMaq, fion1, fion2)
                /surrogate[i].gamma_aq;
              surrogate[i].Aaq=factor*surrogate[i].Atot*Kp*AQinit/(1+Kp*AQinit)+(1.0-factor)*surrogate[i].Aaq;
              surrogate[i].Ap=0.0;
              AQ+=surrogate[i].Aaq;
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
                  (1.0+MOinit/AQinit*MMaq/MOW*surrogate[i].gamma_aq/surrogate[i].gamma_org)+(1.0-factor)*surrogate[i].Aaq;
                surrogate[i].Ap=surrogate[i].Atot-surrogate[i].Aaq;
                MO+=surrogate[i].Ap;
                AQ+=surrogate[i].Aaq;
              }
          else
            {
              fion1=0.0;
              fion2=0.0;
              Kp_aq=surrogate[i].Kp_eff_aqreal_ssh(config, Temperature, ionic, chp,
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

      i=config.iH2SO4; //sulfate
      total=surrogate[i].Ag+surrogate[config.iHSO4m].Aaq/surrogate[config.iHSO4m].MM*surrogate[i].MM
        +surrogate[config.iSO4mm].Aaq/surrogate[config.iSO4mm].MM*surrogate[i].MM; //gas+particle concentrations
      double Keq=surrogate[i].keq/chp*surrogate[config.iHSO4m].gamma_aq/surrogate[config.iHp].gamma_aq/surrogate[config.iSO4mm].gamma_aq;
      surrogate[i].Ag=(1.0-factor)*surrogate[i].Ag;
      surrogate[config.iHSO4m].Aaq=factor*total*surrogate[config.iHSO4m].MM/surrogate[i].MM/(1.0+Keq)
        +(1.0-factor)*surrogate[config.iHSO4m].Aaq; //HSO4-
      surrogate[config.iSO4mm].Aaq=factor*total*surrogate[config.iSO4mm].MM/surrogate[i].MM*Keq/
        (1.0+Keq)+(1.0-factor)*surrogate[config.iSO4mm].Aaq; //SO4--
	      
      i=config.iNH3; //ammoniac	      
      total=surrogate[i].Ag+surrogate[config.iNH4p].Aaq/surrogate[config.iNH4p].MM*surrogate[i].MM; //gas+particle concentrations
      surrogate[i].Ag=factor*total/(1.0+surrogate[i].Kaq_inorg*conc_org)+(1.0-factor)*surrogate[i].Ag; 
      surrogate[config.iNH4p].Aaq=factor*total*surrogate[i].Kaq_inorg*conc_org/(1.0+surrogate[i].Kaq_inorg*conc_org)*surrogate[config.iNH4p].MM/surrogate[i].MM
        +(1.0-factor)*surrogate[config.iNH4p].Aaq; //NH4+
	      
      i=config.iHNO3; //nitrate	      
      total=surrogate[i].Ag+surrogate[config.iNO3m].Aaq/surrogate[config.iNO3m].MM*surrogate[i].MM; //gas+particle concentrations
      surrogate[i].Ag=factor*total/(1.0+surrogate[i].Kaq_inorg*conc_org)+(1.0-factor)*surrogate[i].Ag;
      //cout << "NO3 " << surrogate[i].Kaq_inorg << endl;
      surrogate[config.iNO3m].Aaq=factor*total*surrogate[i].Kaq_inorg*conc_org/(1.0+surrogate[i].Kaq_inorg*conc_org)*surrogate[config.iNO3m].MM/surrogate[i].MM
        +(1.0-factor)*surrogate[config.iNO3m].Aaq; //NO3-
	      
      i=config.iHCl; //chloride	      
      total=surrogate[i].Ag+surrogate[config.iClm].Aaq/surrogate[config.iClm].MM*surrogate[i].MM; //gas+particle concentrations
      surrogate[i].Ag=factor*total/(1.0+surrogate[i].Kaq_inorg*conc_org)+(1.0-factor)*surrogate[i].Ag;
      surrogate[config.iClm].Aaq=factor*total*surrogate[i].Kaq_inorg*conc_org/(1.0+surrogate[i].Kaq_inorg*conc_org)*surrogate[config.iClm].MM/surrogate[i].MM
        +(1.0-factor)*surrogate[config.iClm].Aaq; //Cl-	      	  

      for (i=0;i<n;i++)
        if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)
          conc_inorganic+=surrogate[i].Aaq; //compute total inorganic concentrations

      AQ+=conc_inorganic;
      if (config.iH2O>=0)
        {
          //compute hygroscopicity
          LWC=0.0; //In the case where SOAP compute inorganics, LWC is only used for the initialization
          if (config.hygroscopicity)
	    {
	      /*
	      hygroscopicity_tot_ssh(config, surrogate, Temperature, RH, AQinit, conc_inorganic,MMaq,AQ,deriv_error2_AQ,factor);
	      hygroscopicity_org_ssh(config, surrogate, Temperature, MOW, RH, MOinit, MO, deriv_error1_MO, factor);
	      */
	      hygroscopicity_coupled_tot_ssh(config, surrogate, Temperature, MOW, MMaq, RH, LWC, conc_inorganic,
					     MOinit, MO, AQinit, AQ, deriv_error1_MO, 
					     deriv_error1_AQ, deriv_error2_MO, deriv_error2_AQ, factor);
	    }
	  
          //hygroscopicity_tot_ssh(config, surrogate, Temperature, RH, AQinit, conc_inorganic,MMaq,AQ,deriv_error2_AQ,factor);
          //if (config.hygroscopicity and config.compute_organic)
          //  if (surrogate[config.iH2O].hydrophobic) //Can H2O condense on the organic phase
          //  hygroscopicity_org_ssh(config, surrogate, MOW, RH, MOinit, MO, deriv_error1_MO,factor);
        }
    }
  else
    {
      AQ+=LWC+conc_inorganic;  
      if (config.iH2O>=0 and config.hygroscopicity)        
        if (config.hygroscopicity)
          hygroscopicity_coupled_ssh(config, surrogate, Temperature, MOW, MMaq, RH, LWC, conc_inorganic,
                                 MOinit, MO, AQinit, AQ, deriv_error1_MO, 
                                 deriv_error1_AQ, deriv_error2_MO, deriv_error2_AQ, factor);        
    }  
  /*
    for (i=0;i<n;i++)
    cout << surrogate[i].name << " " << surrogate[i].Aaq << endl;*/
  
  AQ=max(AQ,config.MOmin);
  MO=max(MO,config.MOmin);
  
  
  error1=MOinit-MO;
  error2=AQinit-AQ;
  ionic_organic=ionic_organic_tmp;
  organion=organion_tmp;

}

void newton_raphson_coupled_ssh(double &MO, double &AQ,
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

void init_saturation2_ssh(model_config &config, vector<species>& surrogate,
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


void init_saturation_ssh(model_config &config, vector<species>& surrogate,
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
	  
      if (sumX == 0.0)
        throw string("Error: division by zero sumX in equilibrium.cxx.");
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

      if (sumX == 0.0)
        throw string("Error: division by zero sumX #2 in equilibrium.cxx.");      
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
            surrogate[i].Ap_sat(j)=0.0; //surrogate[i].Ap_sat_old(j);
	  
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
  /*    else
        AQinit+=surrogate[i].Aaq;*/
  
}

void error_saturation_ssh(model_config &config, vector<species>& surrogate,
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
  
  for (i=0;i<n;++i)
    if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)          
      conc_org=max(conc_org,surrogate[i].Aaq/surrogate[i].MM/config.molalmax*1000.0);
  
  conc_org=max(conc_org,1.e-5*config.MOmin);
  //conc_org=max(conc_org,config.MOmin);

  //if (config.compute_rho_aqueous)
  //compute_density_aqueous_phase(config, surrogate, LWC, Temperature);
  compute_ionic_strenght2_ssh(config, surrogate, Temperature, AQinit, conc_inorganic, ionic, chp, organion, ionic_organic, conc_org, factor);  
  
  for (i=0;i<nphase+1;++i)
    {
      for (j=0;j<nphase+1;++j)
        Jacobian(i,j)=0.0;
      Jacobian(i,i)=1.0;
    }

  //compute acitivity coefficients, MOW and MMaq
  if (compute_activity_coefficients)
    {
      activity_coefficients_saturation_ssh(config, surrogate, false, Temperature, MOW);
      activity_coefficients_aq_ssh(config, surrogate, Temperature, LWC, MMaq, XH2O, conc_org);
      if (config.compute_long_and_medium_range_interactions)
        activity_coefficients_LR_MR_ssh(config, surrogate, Temperature, LWC, ionic);
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
              Kpreal_inorganic_ssh(config, surrogate, chp2);
              error_ph_sat_ssh(config, surrogate, Temperature, chp2, organion, error_h, derivative_h,AQinit,LWC,MMaq,MOinit,MOW,conc_org);
              if (chp2-error_h/derivative_h>0.0 and derivative_h!=0.0)
                chp2=chp2-error_h/derivative_h;
              else
                chp2=chp2+error_h;    
            }
          chp=max(factor*chp2+(1.0-factor)*chp,1.0e-15);
        }
      else
        {
          double inorganion=0.0;
          for (i=0;i<n;++i)
            if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)
              if (i!=config.iHp)
                inorganion-=surrogate[i].molality*surrogate[i].charge;
      
          chp=max(factor*0.5*(organion+inorganion+pow(pow(organion+inorganion,2)+4*config.Ke,0.5))+(1.0-factor)*chp,1.0e-15);
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
            sum+=MOinit(j)/MOW(j)/surrogate[i].gamma_org_sat(j);
		  
          surrogate[i].Aaq=0.0;
          for (j=0;j<nphase;++j)
            {
              surrogate[i].Ap_sat(j)=factor*surrogate[i].Atot*
                (MOinit(j)/MOW(j)/surrogate[i].gamma_org_sat(j))/sum+(1.0-factor)*surrogate[i].Ap_sat(j);
              MO(j)+=surrogate[i].Ap_sat(j);
              for (k=0;k<nphase;++k)
                if (j==k)
                  Jacobian(j,k)-=surrogate[i].Atot/(sum*surrogate[i].gamma_org_sat(j)*MOW(j))
                    -surrogate[i].Atot*MOinit(j)/(pow(sum*surrogate[i].gamma_org_sat(j)*MOW(j),2));
                else
                  Jacobian(j,k)+=surrogate[i].Atot*MOinit(j)/
                    (pow(sum,2)*surrogate[i].gamma_org_sat(j)*MOW(j)*surrogate[i].gamma_org_sat(k)*MOW(k));	  
            }
        }
      else
        if (surrogate[i].kp_from_experiment)
          {
            Kp=surrogate[i].kpi;
            for (j=0;j<nphase;++j)
              if (j!=surrogate[i].jmain_phase)
                surrogate[i].Ap_sat(j)=0.0; //(1.0-factor)*surrogate[i].Ap_sat(j);
			
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
          Kp=surrogate[i].Kp_eff_aqreal_ssh(config, Temperature, ionic, chp,surrogate[config.iHp].gamma_LR,
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
          sum=AQinit/surrogate[i].gamma_aq/MMaq;
          for (j=0;j<nphase;++j)
            sum+=MOinit(j)/surrogate[i].gamma_org_sat(j)/MOW(j);
		  
          surrogate[i].Aaq=factor*surrogate[i].Atot*AQinit/(surrogate[i].gamma_aq*MMaq*sum)
            +(1.0-factor)*surrogate[i].Aaq;
          AQ+=surrogate[i].Aaq;
		  
          for (j=0;j<nphase;++j)
            {
              surrogate[i].Ap_sat(j)=(1.0-factor)*surrogate[i].Ap_sat(j)+factor*surrogate[i].Atot*
                (MOinit(j)/surrogate[i].gamma_org_sat(j)/MOW(j))/sum;
              MO(j)+=surrogate[i].Ap_sat(j);
              for (k=0;k<nphase;++k)
                if (j==k)
                  Jacobian(j,k)-=surrogate[i].Atot/(sum*surrogate[i].gamma_org_sat(j)*MOW(j))
                    -surrogate[i].Atot*MOinit(j)/(pow(sum*surrogate[i].gamma_org_sat(j)*MOW(j),2));
                else
                  Jacobian(j,k)+=surrogate[i].Atot*MOinit(j)/
                    (pow(sum,2)*surrogate[i].gamma_org_sat(j)*MOW(j)*surrogate[i].gamma_org_sat(k)*MOW(k));
              Jacobian(j,nphase)+=surrogate[i].Atot*MOinit(j)/
                (pow(sum,2)*surrogate[i].gamma_org_sat(j)*surrogate[i].gamma_aq*MMaq);
            }
		  
          for (k=0;k<nphase;++k)
            Jacobian(nphase,k)+=surrogate[i].Atot*AQinit/
              (surrogate[i].gamma_aq*MMaq*pow(sum,2)*surrogate[i].gamma_org_sat(k)*MOW(k));

          Jacobian(nphase,nphase)-=surrogate[i].Atot/(sum*surrogate[i].gamma_aq*MMaq)
            -surrogate[i].Atot*AQinit/(pow(sum*MMaq*surrogate[i].gamma_aq,2));
        }
      else
        {
          fion1=0.0;
          fion2=0.0;
          Kp_aq=surrogate[i].Kp_eff_aqreal_ssh(config, Temperature, ionic, chp,
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
      Kpreal_inorganic_ssh(config, surrogate, chp);

      i=config.iH2SO4; //sulfate
      total=surrogate[i].Ag+surrogate[config.iHSO4m].Aaq/surrogate[config.iHSO4m].MM*surrogate[i].MM
        +surrogate[config.iSO4mm].Aaq/surrogate[config.iSO4mm].MM*surrogate[i].MM;
      double Keq=surrogate[i].keq/chp*surrogate[config.iHSO4m].gamma_aq/surrogate[config.iHp].gamma_aq/surrogate[config.iSO4mm].gamma_aq;
      surrogate[i].Ag=(1.0-factor)*surrogate[i].Ag;
      surrogate[config.iHSO4m].Aaq=factor*total*surrogate[config.iHSO4m].MM/surrogate[i].MM/(1.0+Keq)
        +(1.0-factor)*surrogate[config.iHSO4m].Aaq;
      surrogate[config.iSO4mm].Aaq=factor*total*surrogate[config.iSO4mm].MM/surrogate[i].MM*Keq/
        (1.0+Keq)+(1.0-factor)*surrogate[config.iSO4mm].Aaq;

      i=config.iNH3; //ammonium      
      total=surrogate[i].Ag+surrogate[config.iNH4p].Aaq/surrogate[config.iNH4p].MM*surrogate[i].MM;
      surrogate[i].Ag=factor*total/(1.0+surrogate[i].Kaq_inorg*conc_org)+(1.0-factor)*surrogate[i].Ag;
      surrogate[config.iNH4p].Aaq=factor*total*surrogate[i].Kaq_inorg*conc_org/(1.0+surrogate[i].Kaq_inorg*conc_org)*surrogate[config.iNH4p].MM/surrogate[i].MM
        +(1.0-factor)*surrogate[config.iNH4p].Aaq;
      Jacobian(nphase,nphase)+=total*surrogate[config.iNH4p].MM/surrogate[i].MM*pow(surrogate[i].Kaq_inorg,2)*conc_org/(pow(1.+surrogate[i].Kaq_inorg*conc_org,2))
        -total*surrogate[config.iNH4p].MM/surrogate[i].MM*surrogate[i].Kaq_inorg/(1.+surrogate[i].Kaq_inorg*conc_org);
            
      i=config.iHNO3;  //nitrate
      total=surrogate[i].Ag+surrogate[config.iNO3m].Aaq/surrogate[config.iNO3m].MM*surrogate[i].MM;
      surrogate[i].Ag=factor*total/(1.0+surrogate[i].Kaq_inorg*conc_org)+(1.0-factor)*surrogate[i].Ag;
      surrogate[config.iNO3m].Aaq=factor*total*surrogate[i].Kaq_inorg*conc_org/(1.0+surrogate[i].Kaq_inorg*conc_org)*surrogate[config.iNO3m].MM/surrogate[i].MM
        +(1.0-factor)*surrogate[config.iNO3m].Aaq;
      Jacobian(nphase,nphase)+=total*surrogate[config.iNO3m].MM/surrogate[i].MM*pow(surrogate[i].Kaq_inorg,2)*conc_org/(pow(1.+surrogate[i].Kaq_inorg*conc_org,2))
        -total*surrogate[config.iNO3m].MM/surrogate[i].MM*surrogate[i].Kaq_inorg/(1.+surrogate[i].Kaq_inorg*conc_org);
            
      i=config.iHCl;  //chloride            
      total=surrogate[i].Ag+surrogate[config.iClm].Aaq/surrogate[config.iClm].MM*surrogate[i].MM;
      surrogate[i].Ag=factor*total/(1.0+surrogate[i].Kaq_inorg*conc_org)+(1.0-factor)*surrogate[i].Ag;
      surrogate[config.iClm].Aaq=factor*total*surrogate[i].Kaq_inorg*conc_org/(1.0+surrogate[i].Kaq_inorg*conc_org)*surrogate[config.iClm].MM/surrogate[i].MM
        +(1.0-factor)*surrogate[config.iClm].Aaq;
      Jacobian(nphase,nphase)+=total*surrogate[config.iClm].MM/surrogate[i].MM*pow(surrogate[i].Kaq_inorg,2)*conc_org/(pow(1.+surrogate[i].Kaq_inorg*conc_org,2))
        -total*surrogate[config.iClm].MM/surrogate[i].MM*surrogate[i].Kaq_inorg/(1.+surrogate[i].Kaq_inorg*conc_org);           

      for (i=0;i<n;i++)
        if (surrogate[i].is_organic==false and i!=config.iH2O and surrogate[i].is_inorganic_precursor==false)
          conc_inorganic+=surrogate[i].Aaq; //compute total inorganic concentrations
  
      AQ+=conc_inorganic;
      if (config.iH2O>=0)
        {
          LWC=0.0;
          /*
            hygroscopicity_tot_ssh(config, surrogate, Temperature, RH, AQinit, conc_inorganic,MMaq,AQ,Jacobian(nphase,nphase),factor);
            if (config.hygroscopicity)
            if (surrogate[config.iH2O].hydrophobic) //Can H2O condense on the organic phase
            for (j=0;j<nphase;++j)
            hygroscopicity_org_sat_ssh(config, surrogate,
            j, MOW(j), RH, MOinit(j), MO(j),Jacobian(j,j),factor);  */
	  
          if (surrogate[config.iH2O].hydrophilic)
            hygroscopicity_org_sat_ssh(config, surrogate, Temperature, MOW, MMaq, RH, LWC, conc_inorganic, MOinit, MO, AQinit, AQ, Jacobian,
                                   true, factor);
          else
            hygroscopicity_org_sat_ssh(config, surrogate, Temperature, MOW, MMaq, RH, LWC, conc_inorganic, MOinit, MO, AQinit, AQ, Jacobian,
                                   false, factor);	      
        }
      
    }
  else
    {	
      AQ+=LWC+conc_inorganic;  
      if (config.iH2O>=0 and config.hygroscopicity)
        if (surrogate[config.iH2O].hydrophilic)
          hygroscopicity_org_sat_ssh(config, surrogate, Temperature, MOW, MMaq, RH, LWC, conc_inorganic, MOinit, MO, AQinit, AQ, Jacobian,
                                 true, factor);
        else
          hygroscopicity_org_sat_ssh(config, surrogate, Temperature, MOW, MMaq, RH, LWC, conc_inorganic, MOinit, MO, AQinit, AQ, Jacobian,
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

void error_saturation_hydrophobic_ssh(model_config &config, vector<species>& surrogate,
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
      MOinit(i)=MOinit(i); //max(MOinit(i),config.MOmin);
    }
  
  for (i=0;i<nphase+1;++i)
    {
      for (j=0;j<nphase+1;++j)
        Jacobian(i,j)=0.0;
      Jacobian(i,i)=1.0;
    }
  
  //compute acitivity coefficients and MOW
  if (compute_activity_coefficients)
    activity_coefficients_saturation_ssh(config, surrogate, true, Temperature, MOW);
  
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
            sum+=MOinit(j)/surrogate[i].gamma_org_sat(j)/MOW(j);
		  
          surrogate[i].Aaq=0.0;
          for (j=0;j<nphase;++j)
            {
              surrogate[i].Ap_sat(j)=factor*surrogate[i].Atot*
                (MOinit(j)/surrogate[i].gamma_org_sat(j)/MOW(j))/sum+(1.0-factor)*surrogate[i].Ap_sat(j);
              MO(j)+=surrogate[i].Ap_sat(j);
              for (k=0;k<nphase;++k)
                if (j==k)
                  Jacobian(j,k)-=surrogate[i].Atot/(sum*surrogate[i].gamma_org_sat(j)*MOW(j))
                    -surrogate[i].Atot*MOinit(j)/(pow(sum*surrogate[i].gamma_org_sat(j)*MOW(j),2));
                else
                  Jacobian(j,k)+=surrogate[i].Atot*MOinit(j)/
                    (pow(sum,2)*surrogate[i].gamma_org_sat(j)*MOW(j)*surrogate[i].gamma_org_sat(k)*MOW(k));	  
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
              Kp_org(j)=surrogate[i].kpi/MOW(j)/surrogate[i].gamma_org_sat(j);

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
          hygroscopicity_org_sat_ssh(config, surrogate, Temperature, MOW, MMaq, RH, 0.0, 0.0, MOinit, MO, AQinit, AQ, Jacobian,
                                 false, factor);	  
        }

  for (j=0;j<nphase;++j)
    {
      MO(j)=max(MO(j),config.MOmin);
      error(j)=MOinit(j)-MO(j);      
    }
  error(nphase)=0.0;
  
}

void stability_ssh(model_config &config, vector<species>& surrogate, double &Temperature,
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
    else if (surrogate[i].is_inorganic_precursor and surrogate[i].is_solid==false)
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
              log(surrogate[i].partial_pressure/surrogate[i].Psat_ssh(Temperature));

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
              log(surrogate[i].partial_pressure_old/surrogate[i].Psat_ssh(Temperature));

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


void mineur_ssh(int &n, int &i, int &j, Array<double, 2> &Jacobian, Array<double, 2> &com_Jacobian)
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


double determinant_ssh(int n, Array<double, 2> &Jacobian)
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
                  mineur_ssh(m, i,j, com_Jacobian_old, com_Jacobian);
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

void compute_comatrice_ssh(int &n, Array<double, 2> &Jacobian, Array<double, 2> &comatrice)
{
  int i,j;
  Array<double, 2> com_Jacobian;
  comatrice.resize(n,n);
  for (i=0;i<n;++i)
    for (j=0 ; j<n; ++j)
      {
        mineur_ssh(n,i,j,Jacobian,com_Jacobian);
        comatrice(i,j)=pow(-1.0,i+j)*determinant_ssh(n-1,com_Jacobian);
      }
}

void inversion_matrice_ssh(int &n, Array<double, 2> &Jacobian, Array<double, 2> &inverse)
{
  int i,j;
  Array<double, 2> comatrice;
  compute_comatrice_ssh(n,Jacobian,comatrice);
  double det=determinant_ssh(n,Jacobian);
  inverse.resize(n,n);
  
  for (i=0;i<n;++i)
    for (j=0;j<n;++j)
      inverse(i,j)=comatrice(j,i)/det;
}

void newton_raphson_sat_ssh(model_config &config, Array<double, 1> &MO_sat, double &AQ_sat,
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
  
  if (determinant_ssh(n,Jacobian)!=0.0 and point_fixe==false)
    {
      inversion_matrice_ssh(n, Jacobian, inverse);
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


void saturation_ssh(model_config &config, vector<species>& surrogate,
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
          for (j=0;j<int(surrogate[0].Ap_sat.size());++j)
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
      //double AQ_sat_old=0.0;
      //double AQ_sat_old2=0.0;
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
      //bool point_fixe=false;
      if (config.initialized_saturation and surrogate[0].Ap_sat_save.size()>1)
        {	  
          for (i=0;i<n;++i)
            {
              for (j=0;j<nphase;++j)
                surrogate[i].Ap_sat(j)=surrogate[i].Ap_sat_save(j);
              if (nphase<int(surrogate[0].Ap_sat_save.size()))
                for (j=nphase;j<int(surrogate[0].Ap_sat_save.size());++j)
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
        init_saturation_ssh(config,surrogate,MO_sat,AQ_sat,LWC,conc_inorganic,all_hydrophobic);      
      
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
            error_saturation_hydrophobic_ssh(config,surrogate,MO_sat,MOW,
                                         Temperature,RH,error_sat,Jacobian,1.0/nh, compute_activity_coefficients);
          else
            error_saturation_ssh(config,surrogate,MO_sat,MOW,MMaq,AQ_sat,LWC,conc_inorganic,
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
      stability_ssh(config,surrogate,Temperature,stable_system,all_hydrophobic,MO_sat);
      /*if (stable_system)
        cout << index_iter << " stable " << nphase << endl;
        else
        cout << index_iter << " instable " << nphase << endl;*/
      
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
