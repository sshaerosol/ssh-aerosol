//!!-----------------------------------------------------------------------
//!!     Copyright (C) 2019 CEREA (ENPC) - INERIS
//!!     SSH-aerosol is distributed under the GNU General Public License v3
//!!-----------------------------------------------------------------------

void kinetic_ssh(model_config &config, vector<species>& surrogate,
                double &MOinit,double &MOW,double &MMaq,
                double &LWC, double &AQinit, double &ionic, double &chp,
                double &Temperature, double &RH, double &deltat, int index)
{
  int n=surrogate.size();
  int i,j,jion,jmol;
  double gamma=1.+0.5*pow(2,0.5);  
  double tiny=0.; 
  for (i=0;i<n;++i)        
    surrogate[i].flux_chem_tot(index)=0.0;          
      
  double Xmono=0.0;          
  double Xmonoaq=0.0;  
  for (i=0;i<n;i++)
    if (surrogate[i].is_monomer and surrogate[i].is_organic)
      {
        if (surrogate[i].hydrophobic or LWC<=config.LWClimit)
          Xmono+=surrogate[i].gamma_org*surrogate[i].Ap*MOW/surrogate[i].MM/max(MOinit,1.0e-10);              
        if (surrogate[i].hydrophilic)
          Xmonoaq+=surrogate[i].gamma_aq*surrogate[i].GAMMAinf*surrogate[i].Aaq*MMaq/surrogate[i].MM/max(AQinit,1.0e-10);                      
      }  

  double XH2O=RH;
  double xmin=1.0e-10;
  double fion1,fion2; 
  double conc_org=LWC;
  for (i=0;i<n;++i)
    if (surrogate[i].is_organic or i==config.iH2O)      
      conc_org+=surrogate[i].Aaq;


  for (i=0;i<n;++i)
    if (surrogate[i].is_organic)
      {
        if (surrogate[i].is_monomer and surrogate[i].Ap+surrogate[i].Aaq>0.0)
          {
            j=surrogate[i].ioligo;  
            double Keq2=surrogate[j].MM/surrogate[i].MM*pow(surrogate[i].Keq_oligo,surrogate[i].moligo-1)*pow(max(Xmono,xmin),surrogate[i].moligo-2)/pow(max(XH2O,xmin),surrogate[i].moligo-1);                      
            double Kaq2=surrogate[j].MM/surrogate[i].MM*pow(surrogate[i].Keq_oligo,surrogate[i].moligo-1)*pow(max(Xmonoaq,xmin),surrogate[i].moligo-2)/pow(max(XH2O,xmin),surrogate[i].moligo-1);  
            //double K2=(surrogate[i].gamma_org*surrogate[i].Ap*Keq2+surrogate[i].GAMMAinf*surrogate[i].gamma_aq*surrogate[i].Aaq*Kaq2)/
            //  (surrogate[i].gamma_org*surrogate[i].Ap+surrogate[i].GAMMAinf*surrogate[i].gamma_aq*surrogate[i].Aaq);
            double fraci=0.0;
            double fracj=0.0;
            double fraciaq=0.0;
            double fracjaq=0.0;
            if (surrogate[i].nonvolatile)
              {
                if (surrogate[i].hydrophilic and surrogate[i].hydrophobic and LWC>config.LWClimit)
                  {                
                    double ratio=AQinit/max(MOinit,1.0e-10)*surrogate[i].gamma_org/(surrogate[i].gamma_aq*surrogate[i].GAMMAinf)*MOW/MMaq;
                    fraci=surrogate[i].gamma_org*1.0/(1.0+ratio)*Xmono; //+surrogate[i].GAMMAinf*surrogate[i].gamma_aq*ratio/(1.0+ratio)*Xmonoaq;
                    fraciaq=surrogate[i].GAMMAinf*surrogate[i].gamma_aq*ratio/(1.0+ratio)*Xmonoaq;
                    //cout << ratio << " " << surrogate[i].gamma_org << " " << surrogate[i].gamma_aq*surrogate[i].GAMMAinf << " " << surrogate[i].Ap << " " << surrogate[i].Aaq << endl;
                  }
                else if (surrogate[i].hydrophobic or LWC<=config.LWClimit)                
                  fraci=surrogate[i].gamma_org*Xmono; 
                else if (surrogate[i].hydrophilic)                
                  fraciaq=surrogate[i].gamma_aq*surrogate[i].GAMMAinf*Xmonoaq;                  
              }
            else
              {
                double Kaq=0.0;
                double Kp=0.0;
                double sum=1.0;
                if (surrogate[i].hydrophobic or LWC<=config.LWClimit)
                  {
		    if (surrogate[i].kp_from_experiment)
		      Kp=surrogate[i].kpi;
		    else
		      Kp=surrogate[i].kpi/MOW/surrogate[i].gamma_org;
                    sum+=Kp*MOinit;
                  }
 
                if (surrogate[i].hydrophilic and LWC>config.LWClimit) 
                  {
                    Kaq=surrogate[i].Kp_eff_aqreal_ssh(config, Temperature, ionic, chp,surrogate[config.iHp].gamma_LR,
                                                   surrogate[config.iHp].gamma_SRMR,MMaq,fion1,fion2)
                      /surrogate[i].gamma_aq;
                    sum+=Kaq*AQinit;
                  }             

                if (surrogate[i].hydrophilic and surrogate[i].hydrophobic and LWC>config.LWClimit) 
                  {               
                    fraciaq=surrogate[i].GAMMAinf*surrogate[i].gamma_aq*Kaq*AQinit/sum*Xmonoaq; //+surrogate[i].gamma_org*Kp*MOinit/sum*Xmono;                
                    fraci=surrogate[i].gamma_org*Kp*MOinit/sum*Xmono;
                  }
                else if (surrogate[i].hydrophobic or LWC<=config.LWClimit)                
                  fraci=surrogate[i].gamma_org*Kp*MOinit/sum*Xmono; 
                else if (surrogate[i].hydrophilic)                
                  fraciaq=surrogate[i].GAMMAinf*surrogate[i].gamma_aq*Kaq*AQinit/sum*Xmonoaq;                
                             
              }

            if (surrogate[j].nonvolatile)
              {
                if (surrogate[j].hydrophilic and surrogate[j].hydrophobic and LWC>config.LWClimit)
                  {
                    double ratio=AQinit/max(MOinit,1.0e-10)*surrogate[j].gamma_org/(surrogate[j].GAMMAinf*surrogate[j].gamma_aq);
                    fracj=surrogate[j].gamma_org*1.0/(1.0+ratio)*Xmono; //+surrogate[j].GAMMAinf*surrogate[j].gamma_aq*ratio/(1.0+ratio)*Xmonoaq;
                    fracjaq=surrogate[j].GAMMAinf*surrogate[j].gamma_aq*ratio/(1.0+ratio)*Xmonoaq;                   
                  }
                else if (surrogate[j].hydrophobic or LWC<=config.LWClimit)                
                  fracj=surrogate[j].gamma_org*Xmono;   
                else if (surrogate[j].hydrophilic)                
                  fracjaq=surrogate[j].gamma_aq*surrogate[j].GAMMAinf*Xmonoaq;                

              }
            else
              {
                double Kaq=0.0;
                double Kp=0.0;
                double sum=1.0;
                if (surrogate[j].hydrophobic or LWC<=config.LWClimit)
                  {
		    if (surrogate[j].kp_from_experiment)
		      Kp=surrogate[j].kpi;
		    else
		      Kp=surrogate[j].kpi/MOW/surrogate[j].gamma_org;
                    sum+=Kp*MOinit;
                  }

                if (surrogate[j].hydrophilic and LWC>config.LWClimit) 
                  {
                    Kaq=surrogate[j].Kp_eff_aqreal_ssh(config, Temperature, ionic, chp,surrogate[config.iHp].gamma_LR,
                                                   surrogate[config.iHp].gamma_SRMR,MMaq,fion1,fion2)
                      /surrogate[j].gamma_aq;
                    sum+=Kaq*AQinit;
                  }             

                if (surrogate[j].hydrophilic and surrogate[j].hydrophobic and LWC>config.LWClimit)                
                  {
                    fracj=surrogate[j].gamma_org*Kp*MOinit/sum*Xmono;                
                    fracjaq=surrogate[j].GAMMAinf*surrogate[j].gamma_aq*Kaq*AQinit/sum*Xmonoaq; 
                  }
                else if (surrogate[j].hydrophobic or LWC<=config.LWClimit)  
                  fracj=surrogate[j].gamma_org*Kp*MOinit/sum*Xmono;   
                else if (surrogate[j].hydrophilic)                
                  fracjaq=surrogate[j].GAMMAinf*surrogate[j].gamma_aq*Kaq*AQinit/sum*Xmonoaq;                                      
              }

	    double flux;
	    if (surrogate[i].catalyzed_ph and config.isorropia_ph==false)
	      flux=surrogate[i].koligo*(surrogate[i].Atot*fraci-surrogate[j].Atot*fracj/Keq2)*deltat+
              surrogate[i].koligo*(surrogate[i].Atot*fraciaq-surrogate[j].Atot*fracjaq/Kaq2)*deltat*chp*surrogate[config.iHp].gamma_aq/config.chp_org_ref;
	    else if (surrogate[i].catalyzed_ph and config.isorropia_ph==true)
	      flux=surrogate[i].koligo*(surrogate[i].Atot*fraci-surrogate[j].Atot*fracj/Keq2)*deltat+
              surrogate[i].koligo*(surrogate[i].Atot*fraciaq-surrogate[j].Atot*fracjaq/Kaq2)*deltat*chp/config.chp_org_ref;
	    else
	      flux=surrogate[i].koligo*(surrogate[i].Atot*fraci-surrogate[j].Atot*fracj/Keq2)*deltat+
		surrogate[i].koligo*(surrogate[i].Atot*fraciaq-surrogate[j].Atot*fracjaq/Kaq2)*deltat;                                                                 
            surrogate[i].flux_chem_tot(index)+=-flux;                        
            surrogate[j].flux_chem_tot(index)+=flux;			                                          
          }

        if (surrogate[i].rion and surrogate[i].Aaq>0.0) // and config.compute_inorganic)
          for (jion=0;jion<surrogate[i].nion_chem;jion++)
            {            
              double molality=surrogate[surrogate[i].iion(jion)].gamma_aq*surrogate[surrogate[i].iion(jion)].Aaq/surrogate[i].MM/conc_org*1000.0;
              double flux=surrogate[i].kion[jion]*molality*deltat*surrogate[i].GAMMAinf*surrogate[i].gamma_aq*surrogate[i].Aaq;

	      if (surrogate[i].rion_ph_catalyzed[jion] and config.isorropia_ph==false)
		flux=flux*chp*surrogate[config.iHp].gamma_aq;
	      else if (surrogate[i].rion_ph_catalyzed[jion] and config.isorropia_ph==true)
		flux=flux*chp;
	      if (surrogate[i].rion_water_catalyzed[jion])
		flux=flux*RH;

              surrogate[i].flux_chem_tot(index)+=-flux;                    
              if (surrogate[i].rion_catalyzed[jion]==false and molality>0.)
                {
                  if (surrogate[i].iion(jion)==config.iHSO4m or surrogate[i].iion(jion)==config.iSO4mm)
                    {
                      double totsulf=surrogate[config.iHSO4m].Aaq/surrogate[config.iHSO4m].MM+surrogate[config.iSO4mm].Aaq/surrogate[config.iSO4mm].MM;
                      surrogate[config.iHSO4m].flux_chem_tot(index)-=surrogate[config.iHSO4m].Aaq/surrogate[config.iHSO4m].MM/totsulf*flux/surrogate[i].MM*surrogate[config.iHSO4m].MM;
                      surrogate[config.iSO4mm].flux_chem_tot(index)-=surrogate[config.iSO4mm].Aaq/surrogate[config.iSO4mm].MM/totsulf*flux/surrogate[i].MM*surrogate[config.iSO4mm].MM;
                    }
                  if (surrogate[i].iion(jion)==config.iNO3m) // and config.compute_inorganic)
                    {
                      double totno3=surrogate[config.iNO3m].Aaq+surrogate[config.iHNO3].Ag/surrogate[config.iHNO3].MM*surrogate[config.iNO3m].MM;
                      surrogate[config.iNO3m].flux_chem_tot(index)-=surrogate[config.iNO3m].Aaq/totno3*flux/surrogate[i].MM*surrogate[config.iNO3m].MM;
                      surrogate[config.iHNO3].flux_chem_tot(index)-=surrogate[config.iHNO3].Ag/surrogate[config.iHNO3].MM*surrogate[config.iNO3m].MM/totno3*flux/surrogate[i].MM*surrogate[config.iHNO3].MM;
                    }
		  /*
		  else if (surrogate[i].iion(jion)==config.iNO3m and config.inorganic_chemistry)
		    {
                      surrogate[config.iNO3m].flux_chem_tot(index)-=flux/surrogate[i].MM*surrogate[config.iNO3m].MM;
		      }*/
                  if (surrogate[i].iion(jion)==config.iClm) // and config.compute_organic)
                    {
                      double totcl=surrogate[config.iClm].Aaq+surrogate[config.iHCl].Ag/surrogate[config.iHCl].MM*surrogate[config.iClm].MM;
                      surrogate[config.iClm].flux_chem_tot(index)-=surrogate[config.iClm].Aaq/totcl*flux/surrogate[i].MM*surrogate[config.iClm].MM;
                      surrogate[config.iHCl].flux_chem_tot(index)-=surrogate[config.iHCl].Ag/surrogate[config.iHCl].MM/totcl*flux/surrogate[i].MM*surrogate[config.iHCl].MM;
                    }
		  /*
		  else if (surrogate[i].iion(jion)==config.iClm and config.inorganic_chemistry)
		    {
                      surrogate[config.iClm].flux_chem_tot(index)-=flux/surrogate[i].MM*surrogate[config.iClm].MM;
                      //flux2+=flux/surrogate[i].MM*surrogate[config.iClm].MM;
		      }*/
                  if (surrogate[i].iion(jion)==config.iNH4p) // and config.compute_organic)
                    {
                      double totnh4=surrogate[config.iNH4p].Aaq+surrogate[config.iNH3].Ag/surrogate[config.iNH3].MM*surrogate[config.iNH4p].MM;
                      surrogate[config.iNH4p].flux_chem_tot(index)-=surrogate[config.iNH4p].Aaq/totnh4*flux/surrogate[i].MM*surrogate[config.iNH4p].MM;
                      surrogate[config.iNH3].flux_chem_tot(index)-=surrogate[config.iNH3].Ag/totnh4*flux/surrogate[config.iNH3].MM/surrogate[i].MM*surrogate[config.iNH3].MM;
                      //flux2+=flux/surrogate[i].MM*surrogate[config.iNH4p].MM;
                    }
		  /*
		  else if (surrogate[i].iion(jion)==config.iNH4p and config.inorganic_chemistry)
		    {
                      surrogate[config.iNH4p].flux_chem_tot(index)-=flux/surrogate[i].MM*surrogate[config.iNH4p].MM;
                      //flux2+=flux/surrogate[i].MM*surrogate[config.iNH4p].MM;
		      }*/

                }
              surrogate[surrogate[i].iproduct(jion)].flux_chem_tot(index)+=flux*surrogate[surrogate[i].iproduct(jion)].MM/surrogate[i].MM;
            }

	for (jmol=0;jmol<surrogate[i].n_irreversible;jmol++)
	  if (surrogate[i].i_irreversible[jmol]>=0)
	    {
	      double flux=0.;
	      double fluxloc=0.;
	      if (surrogate[i].hydrophilic)
		{		
		  fluxloc=surrogate[i].k_irreversible[jmol]*deltat*surrogate[i].GAMMAinf*surrogate[i].gamma_aq*surrogate[i].Aaq;
		  if (surrogate[i].irr_catalyzed_water[jmol])
		    fluxloc=fluxloc*RH;
		  if (surrogate[i].irr_catalyzed_pH[jmol] and config.isorropia_ph==false)
		    fluxloc=fluxloc*chp*surrogate[config.iHp].gamma_aq;
		  else if (surrogate[i].irr_catalyzed_pH[jmol] and config.isorropia_ph==true)
		    fluxloc=fluxloc*chp;
		  flux+=fluxloc;
		}
	      if (surrogate[i].hydrophobic)
		{
		  fluxloc=surrogate[i].k_irreversible[jmol]*deltat*surrogate[i].gamma_org*surrogate[i].Ap;
		  if (surrogate[i].irr_catalyzed_water[jmol])
		    fluxloc=fluxloc*RH;
		  if (surrogate[i].irr_catalyzed_pH[jmol])
		    fluxloc=fluxloc*config.chp_org_ref;
		  flux+=fluxloc;
		}
	     
	      surrogate[i].flux_chem_tot(index)+=-flux;
	      if (surrogate[i].irr_mass_conserving[jmol])
		surrogate[surrogate[i].i_irreversible[jmol]].flux_chem_tot(index)+=flux;
	      else
		surrogate[surrogate[i].i_irreversible[jmol]].flux_chem_tot(index)+=flux*surrogate[surrogate[i].i_irreversible[jmol]].MM/surrogate[i].MM;

	      if (surrogate[i].iother_irreversible[jmol]>=0)
		if (surrogate[i].iother_irreversible[jmol]==config.iNO3m)
		  {
		    double totno3=surrogate[config.iNO3m].Aaq+surrogate[config.iHNO3].Ag/surrogate[config.iHNO3].MM*surrogate[config.iNO3m].MM;
		    if (totno3>0.)
		      {
			surrogate[surrogate[i].iother_irreversible[jmol]].flux_chem_tot(index)+=surrogate[config.iNO3m].Aaq/totno3*flux*surrogate[surrogate[i].iother_irreversible[jmol]].MM/surrogate[i].MM;
			surrogate[config.iHNO3].flux_chem_tot(index)+=surrogate[config.iHNO3].Ag/surrogate[config.iHNO3].MM*surrogate[config.iNO3m].MM/totno3*flux*surrogate[config.iHNO3].MM/surrogate[i].MM;
		      }
		  }
		else if (surrogate[i].iother_irreversible[jmol]==config.iClm)
		  {
		    double totcl=surrogate[config.iClm].Aaq+surrogate[config.iHCl].Ag/surrogate[config.iHCl].MM*surrogate[config.iClm].MM;
		    if (totcl>0.)
		      {
			surrogate[surrogate[i].iother_irreversible[jmol]].flux_chem_tot(index)+=surrogate[config.iClm].Aaq/totcl*flux*surrogate[surrogate[i].iother_irreversible[jmol]].MM/surrogate[i].MM;
			surrogate[config.iHCl].flux_chem_tot(index)+=surrogate[config.iHCl].Ag/surrogate[config.iHCl].MM*surrogate[config.iClm].MM/totcl*flux*surrogate[config.iHCl].MM/surrogate[i].MM;
		      }
		  }
		else if (surrogate[i].iother_irreversible[jmol]==config.iNH4p)
		  {
		    double totnh4=surrogate[config.iNH4p].Aaq+surrogate[config.iNH3].Ag/surrogate[config.iNH3].MM*surrogate[config.iNH4p].MM;
		    if (totnh4>0.)
		      {
			surrogate[surrogate[i].iother_irreversible[jmol]].flux_chem_tot(index)+=surrogate[config.iNH4p].Aaq/totnh4*flux*surrogate[surrogate[i].iother_irreversible[jmol]].MM/surrogate[i].MM;
			surrogate[config.iNH3].flux_chem_tot(index)+=surrogate[config.iNH3].Ag/surrogate[config.iNH3].MM*surrogate[config.iNH4p].MM/totnh4*flux*surrogate[config.iNH3].MM/surrogate[i].MM;
		      }
		  }
		else if (surrogate[i].iother_irreversible[jmol]==config.iSO4mm or surrogate[i].iother_irreversible[jmol]==config.iHSO4m)
		  {
		    double totsulf=surrogate[config.iHSO4m].Aaq/surrogate[config.iHSO4m].MM+surrogate[config.iSO4mm].Aaq/surrogate[config.iSO4mm].MM;
		    if (totsulf>0.)
		      {
			surrogate[config.iSO4mm].flux_chem_tot(index)+=surrogate[config.iSO4mm].Aaq/surrogate[config.iSO4mm].MM/totsulf*flux*surrogate[config.iSO4mm].MM/surrogate[i].MM;
			surrogate[config.iHSO4m].flux_chem_tot(index)+=surrogate[config.iHSO4m].Aaq/surrogate[config.iHSO4m].MM/totsulf*flux*surrogate[config.iHSO4m].MM/surrogate[i].MM;
		      }
		  }
		else
		  surrogate[surrogate[i].iother_irreversible[jmol]].flux_chem_tot(index)+=flux*surrogate[surrogate[i].iother_irreversible[jmol]].MM/surrogate[i].MM;
	    }
      }
   
  for (i=0;i<n;++i)
    {           
      if (surrogate[i].is_organic)
        if (surrogate[i].Atot>tiny)
          surrogate[i].Jdn_gas(index)=min(surrogate[i].flux_chem_tot(index)/surrogate[i].Atot/deltat,0.0);      
        else if (i==config.iHSO4m or i==config.iSO4mm)
          {
            double totsulf=surrogate[config.iHSO4m].Aaq/surrogate[config.iHSO4m].MM+surrogate[config.iSO4mm].Aaq/surrogate[config.iSO4mm].MM;
            if (totsulf>0.0)
              surrogate[i].Jdn_gas(index)=min((surrogate[config.iHSO4m].flux_chem_tot(index)/surrogate[config.iHSO4m].MM
                                               +surrogate[config.iSO4mm].flux_chem_tot(index)/surrogate[config.iSO4mm].MM)/totsulf/deltat,0.0);
            
          }
        else if (i==config.iNO3m)
          {
	    if (config.compute_inorganic or config.inorganic_chemistry)
	      {
		double totno3=surrogate[config.iNO3m].Aaq+surrogate[config.iHNO3].Ag/surrogate[config.iHNO3].MM*surrogate[config.iNO3m].MM;
		if (totno3>0.0)
		  surrogate[i].Jdn_gas(index)=min((surrogate[config.iNO3m].flux_chem_tot(index)+surrogate[config.iHNO3].flux_chem_tot(index)/surrogate[config.iHNO3].MM*surrogate[config.iNO3m].MM)/
						  totno3/deltat,0.0);
	      }
          }
        else if (i==config.iClm)
          {
	    if (config.compute_inorganic or config.inorganic_chemistry)
	      {
		double totcl=surrogate[config.iClm].Aaq+surrogate[config.iHCl].Ag/surrogate[config.iHCl].MM*surrogate[config.iClm].MM;
		if (totcl>0.0)
		  surrogate[i].Jdn_gas(index)=min((surrogate[config.iClm].flux_chem_tot(index)+surrogate[config.iHCl].flux_chem_tot(index)/surrogate[config.iHCl].MM*surrogate[config.iClm].MM)/
						  totcl/deltat,0.0);
	      }		
          }
        else if (i==config.iNH4p)
          {
	    if (config.compute_inorganic or config.inorganic_chemistry)
	      { 
		double totnh4=surrogate[config.iNH4p].Aaq+surrogate[config.iNH3].Ag/surrogate[config.iNH3].MM*surrogate[config.iNH4p].MM;
		if (totnh4>0.0)
		  surrogate[i].Jdn_gas(index)=min((surrogate[config.iNH4p].flux_chem_tot(index)+surrogate[config.iNH3].flux_chem_tot(index)/surrogate[config.iNH3].MM*surrogate[config.iNH4p].MM)/
                                              totnh4/deltat,0.0);
	      }
		
          }
    }

  for (i=0;i<n;++i)    
    surrogate[i].flux_chem_tot(index)=0.0;   
  
  for (i=0;i<n;++i)
    if (surrogate[i].is_organic)
      {
        if (surrogate[i].is_monomer and surrogate[i].Ap+surrogate[i].Aaq>0.0)
          {
            j=surrogate[i].ioligo;  
            double Keq2=surrogate[j].MM/surrogate[i].MM*pow(surrogate[i].Keq_oligo,surrogate[i].moligo-1)*pow(max(Xmono,xmin),surrogate[i].moligo-2)/pow(max(XH2O,xmin),surrogate[i].moligo-1);                      
            double Kaq2=surrogate[j].MM/surrogate[i].MM*pow(surrogate[i].Keq_oligo,surrogate[i].moligo-1)*pow(max(Xmonoaq,xmin),surrogate[i].moligo-2)/pow(max(XH2O,xmin),surrogate[i].moligo-1);  
            //double K2=(surrogate[i].gamma_org*surrogate[i].Ap*Keq2+surrogate[i].GAMMAinf*surrogate[i].gamma_aq*surrogate[i].Aaq*Kaq2)/
            //  (surrogate[i].gamma_org*surrogate[i].Ap+surrogate[i].GAMMAinf*surrogate[i].gamma_aq*surrogate[i].Aaq);
            double fraci=0.0;
            double fracj=0.0;
            double fraciaq=0.0;
            double fracjaq=0.0;
            if (surrogate[i].nonvolatile)
              {
                if (surrogate[i].hydrophilic and surrogate[i].hydrophobic and LWC>config.LWClimit)
                  {                
                    double ratio=AQinit/max(MOinit,1.0e-10)*surrogate[i].gamma_org/(surrogate[i].gamma_aq*surrogate[i].GAMMAinf)*MOW/MMaq;
                    fraci=surrogate[i].gamma_org*1.0/(1.0+ratio)*Xmono; //+surrogate[i].GAMMAinf*surrogate[i].gamma_aq*ratio/(1.0+ratio)*Xmonoaq;
                    fraciaq=surrogate[i].GAMMAinf*surrogate[i].gamma_aq*ratio/(1.0+ratio)*Xmonoaq;
                  }
                else if (surrogate[i].hydrophobic or LWC<=config.LWClimit)                
                  fraci=surrogate[i].gamma_org*Xmono; 
                else if (surrogate[i].hydrophilic)                
                  fraciaq=surrogate[i].gamma_aq*surrogate[i].GAMMAinf*Xmonoaq;                  
              }
            else
              {
                double Kaq=0.0;
                double Kp=0.0;
                double sum=1.0;
                if (surrogate[i].hydrophobic or LWC<=config.LWClimit)
                  {
		    if (surrogate[i].kp_from_experiment)
		      Kp=surrogate[i].kpi;
		    else
		      Kp=surrogate[i].kpi/MOW/surrogate[i].gamma_org;
                    sum+=Kp*MOinit;
                  }

                if (surrogate[i].hydrophilic and LWC>config.LWClimit) 
                  {
                    Kaq=surrogate[i].Kp_eff_aqreal_ssh(config, Temperature, ionic, chp,surrogate[config.iHp].gamma_LR,
                                                   surrogate[config.iHp].gamma_SRMR,MMaq,fion1,fion2)
                      /surrogate[i].gamma_aq;
                    sum+=Kaq*AQinit;
                  }             

                if (surrogate[i].hydrophilic and surrogate[i].hydrophobic and LWC>config.LWClimit) 
                  {               
                    fraciaq=surrogate[i].GAMMAinf*surrogate[i].gamma_aq*Kaq*AQinit/sum*Xmonoaq; //+surrogate[i].gamma_org*Kp*MOinit/sum*Xmono;                
                    fraci=surrogate[i].gamma_org*Kp*MOinit/sum*Xmono;
                  }
                else if (surrogate[i].hydrophobic or LWC<=config.LWClimit)                
                  fraci=surrogate[i].gamma_org*Kp*MOinit/sum*Xmono; 
                else if (surrogate[i].hydrophilic)                
                  fraciaq=surrogate[i].GAMMAinf*surrogate[i].gamma_aq*Kaq*AQinit/sum*Xmonoaq;                
                             
              }

            if (surrogate[j].nonvolatile)
              {
                if (surrogate[j].hydrophilic and surrogate[j].hydrophobic and LWC>config.LWClimit)
                  {
                    double ratio=AQinit/max(MOinit,1.0e-10)*surrogate[j].gamma_org/(surrogate[j].GAMMAinf*surrogate[j].gamma_aq);
                    fracj=surrogate[j].gamma_org*1.0/(1.0+ratio)*Xmono; //+surrogate[j].GAMMAinf*surrogate[j].gamma_aq*ratio/(1.0+ratio)*Xmonoaq;
                    fracjaq=surrogate[j].GAMMAinf*surrogate[j].gamma_aq*ratio/(1.0+ratio)*Xmonoaq;
                  }
                else if (surrogate[j].hydrophobic or LWC<=config.LWClimit)                
                  fracj=surrogate[j].gamma_org*Xmono;   
                else if (surrogate[j].hydrophilic)                
                  fracjaq=surrogate[j].gamma_aq*surrogate[j].GAMMAinf*Xmonoaq;                

              }
            else
              {
                double Kaq=0.0;
                double Kp=0.0;
                double sum=1.0;
                if (surrogate[j].hydrophobic or LWC<=config.LWClimit)
                  {
		    if (surrogate[j].kp_from_experiment)
		      Kp=surrogate[j].kpi;
		    else
		      Kp=surrogate[j].kpi/MOW/surrogate[j].gamma_org;
                    sum+=Kp*MOinit;
                  }

                if (surrogate[j].hydrophilic and LWC>config.LWClimit) 
                  { 
                    Kaq=surrogate[j].Kp_eff_aqreal_ssh(config, Temperature, ionic, chp,surrogate[config.iHp].gamma_LR,
                                                   surrogate[config.iHp].gamma_SRMR,MMaq,fion1,fion2)
                      /surrogate[j].gamma_aq;
                    sum+=Kaq*AQinit;
                  }             

                if (surrogate[j].hydrophilic and surrogate[j].hydrophobic and LWC>config.LWClimit)                
                  {
                    fracj=surrogate[j].gamma_org*Kp*MOinit/sum*Xmono;                
                    fracjaq=surrogate[j].GAMMAinf*surrogate[j].gamma_aq*Kaq*AQinit/sum*Xmonoaq; 
                  }
                else if (surrogate[j].hydrophobic or LWC<=config.LWClimit)  
                  fracj=surrogate[j].gamma_org*Kp*MOinit/sum*Xmono;   
                else if (surrogate[j].hydrophilic)                
                  fracjaq=surrogate[j].GAMMAinf*surrogate[j].gamma_aq*Kaq*AQinit/sum*Xmonoaq;                                      
              }

	    double flux;
	    if (surrogate[i].catalyzed_ph and config.isorropia_ph==false)
	      flux=surrogate[i].koligo*(surrogate[i].Atot*fraci-surrogate[j].Atot*fracj/Keq2)*deltat+
              surrogate[i].koligo*(surrogate[i].Atot*fraciaq-surrogate[j].Atot*fracjaq/Kaq2)*deltat*chp*surrogate[config.iHp].gamma_aq/config.chp_org_ref;
	    else if (surrogate[i].catalyzed_ph and config.isorropia_ph==true)
	      flux=surrogate[i].koligo*(surrogate[i].Atot*fraci-surrogate[j].Atot*fracj/Keq2)*deltat+
              surrogate[i].koligo*(surrogate[i].Atot*fraciaq-surrogate[j].Atot*fracjaq/Kaq2)*deltat*chp/config.chp_org_ref;	    
            else
	      flux=surrogate[i].koligo*(surrogate[i].Atot*fraci-surrogate[j].Atot*fracj/Keq2)*deltat+
              surrogate[i].koligo*(surrogate[i].Atot*fraciaq-surrogate[j].Atot*fracjaq/Kaq2)*deltat;            

            if (flux>0.0)                        
              flux=flux/(1.0-gamma*surrogate[i].Jdn_gas(index)*deltat);                                                  
            else
              flux=flux/(1.0-gamma*surrogate[j].Jdn_gas(index)*deltat);

            surrogate[i].flux_chem_tot(index)+=-flux;                        
            surrogate[j].flux_chem_tot(index)+=flux;			                                          
          }    
        if (surrogate[i].rion and surrogate[i].Aaq>0.0) // and config.compute_inorganic)
          for (jion=0;jion<surrogate[i].nion_chem;jion++)
            {              
              double molality=surrogate[surrogate[i].iion(jion)].gamma_aq*surrogate[surrogate[i].iion(jion)].Aaq/surrogate[i].MM/conc_org*1000.0;
              double flux=surrogate[i].kion[jion]*molality*deltat*surrogate[i].GAMMAinf*surrogate[i].gamma_aq*surrogate[i].Aaq;
	      if (surrogate[i].rion_ph_catalyzed[jion] and config.isorropia_ph==false)
		flux=flux*chp*surrogate[config.iHp].gamma_aq;
	      else if (surrogate[i].rion_ph_catalyzed[jion] and config.isorropia_ph==true)
		flux=flux*chp;
	      if (surrogate[i].rion_water_catalyzed[jion])
		flux=flux*RH;
             
              if (surrogate[i].rion_catalyzed[jion]==false and molality>0.)
                {
		  
                  if (surrogate[i].iion(jion)==config.iHSO4m or surrogate[i].iion(jion)==config.iSO4mm)
                    {
                      flux=flux/(1.0-gamma*surrogate[i].Jdn_gas(index)*deltat
                                 -gamma*surrogate[config.iHSO4m].Jdn_gas(index)/surrogate[config.iHSO4m].MM*surrogate[i].MM*deltat); 
                      double totsulf=surrogate[config.iHSO4m].Aaq/surrogate[config.iHSO4m].MM+surrogate[config.iSO4mm].Aaq/surrogate[config.iSO4mm].MM;
                      surrogate[config.iHSO4m].flux_chem_tot(index)-=surrogate[config.iHSO4m].Aaq/surrogate[config.iHSO4m].MM/totsulf*flux/surrogate[i].MM*surrogate[config.iHSO4m].MM;
                      surrogate[config.iSO4mm].flux_chem_tot(index)-=surrogate[config.iSO4mm].Aaq/surrogate[config.iSO4mm].MM/totsulf*flux/surrogate[i].MM*surrogate[config.iSO4mm].MM;
                    }
                  else if (surrogate[i].iion(jion)==config.iNO3m) // and config.compute_inorganic)
                    {
                      flux=flux/(1.0-gamma*surrogate[i].Jdn_gas(index)*deltat
                                 -gamma*surrogate[config.iNO3m].Jdn_gas(index)/surrogate[config.iNO3m].MM*surrogate[i].MM*deltat); 
                      double totno3=surrogate[config.iNO3m].Aaq+surrogate[config.iHNO3].Ag/surrogate[config.iHNO3].MM*surrogate[config.iNO3m].MM;
                      surrogate[config.iNO3m].flux_chem_tot(index)-=surrogate[config.iNO3m].Aaq/totno3*flux/surrogate[i].MM*surrogate[config.iNO3m].MM;
                      surrogate[config.iHNO3].flux_chem_tot(index)-=surrogate[config.iHNO3].Ag/surrogate[config.iHNO3].MM*surrogate[config.iNO3m].MM/totno3*flux/surrogate[i].MM*surrogate[config.iHNO3].MM;
                    }
		  /*
		  else if (surrogate[i].iion(jion)==config.iNO3m and config.inorganic_chemistry)
                    {
                      flux=flux/(1.0-gamma*surrogate[i].Jdn_gas(index)*deltat
                                 -gamma*surrogate[config.iNO3m].Jdn_gas(index)/surrogate[config.iNO3m].MM*surrogate[i].MM*deltat); 
                      surrogate[config.iNO3m].flux_chem_tot(index)-=flux/surrogate[i].MM*surrogate[config.iNO3m].MM;
		      }*/
                  else if (surrogate[i].iion(jion)==config.iClm) // and config.compute_inorganic)
                    {
                      flux=flux/(1.0-gamma*surrogate[i].Jdn_gas(index)*deltat
                                 -gamma*surrogate[config.iClm].Jdn_gas(index)/surrogate[config.iClm].MM*surrogate[i].MM*deltat); 
                      double totcl=surrogate[config.iClm].Aaq+surrogate[config.iHCl].Ag/surrogate[config.iHCl].MM*surrogate[config.iClm].MM;
                      surrogate[config.iClm].flux_chem_tot(index)-=surrogate[config.iClm].Aaq/totcl*flux/surrogate[i].MM*surrogate[config.iClm].MM;
                      surrogate[config.iHCl].flux_chem_tot(index)-=surrogate[config.iHCl].Ag/surrogate[config.iHCl].MM*surrogate[config.iClm].MM/totcl*flux/surrogate[i].MM*surrogate[config.iHCl].MM;
                    }
		  /*
		  else if (surrogate[i].iion(jion)==config.iClm and config.inorganic_chemistry)
                    {
                      flux=flux/(1.0-gamma*surrogate[i].Jdn_gas(index)*deltat
                                 -gamma*surrogate[config.iClm].Jdn_gas(index)/surrogate[config.iClm].MM*surrogate[i].MM*deltat); 
                      surrogate[config.iClm].flux_chem_tot(index)-=flux/surrogate[i].MM*surrogate[config.iClm].MM;
		      }*/
                  else if (surrogate[i].iion(jion)==config.iNH4p) // and config.compute_inorganic)
                    {
                      flux=flux/(1.0-gamma*surrogate[i].Jdn_gas(index)*deltat
                                 -gamma*surrogate[config.iNH4p].Jdn_gas(index)/surrogate[config.iNH4p].MM*surrogate[i].MM*deltat); 
                      double totnh4=surrogate[config.iNH4p].Aaq+surrogate[config.iNH3].Ag/surrogate[config.iNH3].MM*surrogate[config.iNH4p].MM;
                      surrogate[config.iNH4p].flux_chem_tot(index)-=surrogate[config.iNH4p].Aaq/totnh4*flux/surrogate[i].MM*surrogate[config.iNH4p].MM;
                      surrogate[config.iNH3].flux_chem_tot(index)-=surrogate[config.iNH3].Ag/surrogate[config.iNH3].MM*surrogate[config.iNH4p].MM/totnh4*flux/surrogate[i].MM*surrogate[config.iNH3].MM;
                    }
		  /*
		  else if (surrogate[i].iion(jion)==config.iNH4p) and config.inorganic_chemistry)
                    {
                      flux=flux/(1.0-gamma*surrogate[i].Jdn_gas(index)*deltat
                                 -gamma*surrogate[config.iNH4p].Jdn_gas(index)/surrogate[config.iNH4p].MM*surrogate[i].MM*deltat); 
                      surrogate[config.iNH4p].flux_chem_tot(index)-=flux/surrogate[i].MM*surrogate[config.iNH4p].MM;
		      }*/

                }
              else
                flux=flux/(1.0-gamma*surrogate[i].Jdn_gas(index)*deltat);

              surrogate[i].flux_chem_tot(index)+=-flux; 
              surrogate[surrogate[i].iproduct(jion)].flux_chem_tot(index)+=flux*surrogate[surrogate[i].iproduct(jion)].MM/surrogate[i].MM;
            }

	for (jmol=0;jmol<surrogate[i].n_irreversible;jmol++)
	  if (surrogate[i].i_irreversible[jmol]>=0)
	    {
	      double flux=0.;
	      double fluxloc=0.;
	      if (surrogate[i].hydrophilic)
		{		
		  fluxloc=surrogate[i].k_irreversible[jmol]*deltat*surrogate[i].GAMMAinf*surrogate[i].gamma_aq*surrogate[i].Aaq;
		  if (surrogate[i].irr_catalyzed_water[jmol])
		    fluxloc=fluxloc*RH;
		  if (surrogate[i].irr_catalyzed_pH[jmol] and config.isorropia_ph==false)
		    fluxloc=fluxloc*chp*surrogate[config.iHp].gamma_aq;
		  else if (surrogate[i].irr_catalyzed_pH[jmol] and config.isorropia_ph==true)
		    fluxloc=fluxloc*chp;
		  flux+=fluxloc;
		}
	      
	      if (surrogate[i].hydrophobic)
		{
		  fluxloc=surrogate[i].k_irreversible[jmol]*deltat*surrogate[i].gamma_org*surrogate[i].Ap;
		  if (surrogate[i].irr_catalyzed_water[jmol])
		    fluxloc=fluxloc*RH;
		  if (surrogate[i].irr_catalyzed_pH[jmol])
		    fluxloc=fluxloc*config.chp_org_ref;
		  flux+=fluxloc;
		}

	      flux=flux/(1.0-gamma*surrogate[i].Jdn_gas(index)*deltat);
	    
	      surrogate[i].flux_chem_tot(index)+=-flux;
	      if (surrogate[i].irr_mass_conserving[jmol])
		surrogate[surrogate[i].i_irreversible[jmol]].flux_chem_tot(index)+=flux;
	      else
		surrogate[surrogate[i].i_irreversible[jmol]].flux_chem_tot(index)+=flux*surrogate[surrogate[i].i_irreversible[jmol]].MM/surrogate[i].MM;

	      if (surrogate[i].iother_irreversible[jmol]>=0)	      
		{
		  if (surrogate[i].iother_irreversible[jmol]==config.iNO3m)
		    {
		      double totno3=surrogate[config.iNO3m].Aaq+surrogate[config.iHNO3].Ag/surrogate[config.iHNO3].MM*surrogate[config.iNO3m].MM;
		      if (totno3>0.)
			{
			  surrogate[surrogate[i].iother_irreversible[jmol]].flux_chem_tot(index)+=surrogate[config.iNO3m].Aaq/totno3*flux*surrogate[surrogate[i].iother_irreversible[jmol]].MM/surrogate[i].MM;
			  surrogate[config.iHNO3].flux_chem_tot(index)+=surrogate[config.iHNO3].Ag/surrogate[config.iHNO3].MM*surrogate[config.iNO3m].MM/totno3*flux*surrogate[config.iHNO3].MM/surrogate[i].MM;
			}
		    }
		  else if (surrogate[i].iother_irreversible[jmol]==config.iClm)
		    {
		      double totcl=surrogate[config.iClm].Aaq+surrogate[config.iHCl].Ag/surrogate[config.iHCl].MM*surrogate[config.iClm].MM;
		      if (totcl>0.)
			{
			  surrogate[surrogate[i].iother_irreversible[jmol]].flux_chem_tot(index)+=surrogate[config.iClm].Aaq/totcl*flux*surrogate[surrogate[i].iother_irreversible[jmol]].MM/surrogate[i].MM;
			  surrogate[config.iHCl].flux_chem_tot(index)+=surrogate[config.iHCl].Ag/surrogate[config.iHCl].MM*surrogate[config.iClm].MM/totcl*flux*surrogate[config.iHCl].MM/surrogate[i].MM;
			}
		    }
		  else if (surrogate[i].iother_irreversible[jmol]==config.iNH4p)
		    {
		      double totnh4=surrogate[config.iNH4p].Aaq+surrogate[config.iNH3].Ag/surrogate[config.iNH3].MM*surrogate[config.iNH4p].MM;
		      if (totnh4>0.)
			{
			  surrogate[surrogate[i].iother_irreversible[jmol]].flux_chem_tot(index)+=surrogate[config.iNH4p].Aaq/totnh4*flux*surrogate[surrogate[i].iother_irreversible[jmol]].MM/surrogate[i].MM;
			  surrogate[config.iNH3].flux_chem_tot(index)+=surrogate[config.iNH3].Ag/surrogate[config.iNH3].MM*surrogate[config.iNH4p].MM/totnh4*flux*surrogate[config.iNH3].MM/surrogate[i].MM;
			}
		    }
		  else if (surrogate[i].iother_irreversible[jmol]==config.iSO4mm or surrogate[i].iother_irreversible[jmol]==config.iHSO4m)
		    {
		      double totsulf=surrogate[config.iHSO4m].Aaq/surrogate[config.iHSO4m].MM+surrogate[config.iSO4mm].Aaq/surrogate[config.iSO4mm].MM;
		      if (totsulf>0.)
			{
			  surrogate[config.iSO4mm].flux_chem_tot(index)+=surrogate[config.iSO4mm].Aaq/surrogate[config.iSO4mm].MM/totsulf*flux*surrogate[config.iSO4mm].MM/surrogate[i].MM;
			  surrogate[config.iHSO4m].flux_chem_tot(index)+=surrogate[config.iHSO4m].Aaq/surrogate[config.iHSO4m].MM/totsulf*flux*surrogate[config.iHSO4m].MM/surrogate[i].MM;
			}
		    }		
		  else
		    surrogate[surrogate[i].iother_irreversible[jmol]].flux_chem_tot(index)+=flux*surrogate[surrogate[i].iother_irreversible[jmol]].MM/surrogate[i].MM;
		}
	
	    }
	
        
      }
} 


void kinetic_sat_ssh(model_config &config, vector<species>& surrogate,
                 Array <double, 1> &MOinit, Array <double, 1> &MOW,double &MMaq,
                 double &LWC, double &AQinit, double &ionic, double &chp,
                 double &Temperature, double &RH, double &deltat, int index)
{
  int n=surrogate.size();
  int i,iphase,j,jion,jmol;
  double gamma=1.+0.5*pow(2,0.5);  
  double tiny=0.;  
  for (i=0;i<n;++i)        
    surrogate[i].flux_chem_tot(index)=0.0;          
      
  int nphase=surrogate[0].Ap_sat.size();
  Array <double, 1> Xmono;          
  Xmono.resize(nphase);
  double Xmonoaq=0.0;  
  for (i=0;i<n;i++)
    if (surrogate[i].is_monomer and surrogate[i].is_organic)
      {
        if (surrogate[i].hydrophobic or LWC<=config.LWClimit)
          for (iphase=0;iphase<nphase;iphase++)
            Xmono(iphase)+=surrogate[i].gamma_org_sat(iphase)*surrogate[i].Ap_sat(iphase)*MOW(iphase)/surrogate[i].MM/max(MOinit(iphase),1.0e-10);              

        if (surrogate[i].hydrophilic)
          Xmonoaq+=surrogate[i].gamma_aq*surrogate[i].GAMMAinf*surrogate[i].Aaq*MMaq/surrogate[i].MM/max(AQinit,1.0e-10);                      
      }  

  double XH2O=RH;
  double xmin=1.0e-10;
  double fion1,fion2; 
  double conc_org=LWC;
  for (i=0;i<n;++i)
    if (surrogate[i].is_organic or i==config.iH2O)      
      conc_org+=surrogate[i].Aaq;      

  for (i=0;i<n;++i)
    if (surrogate[i].is_organic)
      {
        if (surrogate[i].is_monomer and surrogate[i].Atot>0.0)
          {
            j=surrogate[i].ioligo; 
            Array <double, 1> Keq2;
            Keq2.resize(nphase);
            for (iphase=0;iphase<nphase;iphase++)
              Keq2(iphase)=surrogate[j].MM/surrogate[i].MM*pow(surrogate[i].Keq_oligo,surrogate[i].moligo-1)*pow(max(Xmono(iphase),xmin),surrogate[i].moligo-2)/pow(max(XH2O,xmin),surrogate[i].moligo-1);                      
            double Kaq2=surrogate[j].MM/surrogate[i].MM*pow(surrogate[i].Keq_oligo,surrogate[i].moligo-1)*pow(max(Xmonoaq,xmin),surrogate[i].moligo-2)/pow(max(XH2O,xmin),surrogate[i].moligo-1);  
            Array <double, 1> fraci,fracj;
            fraci.resize(nphase);
            fracj.resize(nphase);
            fraci=0.0;
            fracj=0.0;            
            double fraciaq=0.0;
            double fracjaq=0.0;
            if (surrogate[i].nonvolatile)
              {
                double sum=0.0;
                if (surrogate[i].hydrophobic or LWC<=config.LWClimit)
                  for (iphase=0;iphase<nphase;iphase++)                                   
                    sum+=MOinit(iphase)/MOW(iphase)/surrogate[i].gamma_org_sat(iphase);                    
 
                if (surrogate[i].hydrophilic and LWC>config.LWClimit) 
                  sum+=AQinit/(MMaq*surrogate[i].gamma_aq*surrogate[i].GAMMAinf);

                if (sum>0.0)
                  {
                    if (surrogate[i].hydrophobic or LWC<=config.LWClimit)
                      for (iphase=0;iphase<nphase;iphase++)                                   
                        fraci(iphase)=MOinit(iphase)/MOW(iphase)/sum*Xmono(iphase);                    
                    
                    if (surrogate[i].hydrophilic and LWC>config.LWClimit) 
                      fraciaq=AQinit/MMaq/sum*Xmonoaq;                 
                  }
              }
            else
              {
                double Kaq=0.0;
                Array<double, 1> Kp;
                Kp.resize(nphase);                
                double sum=1.0;
                
                if (surrogate[i].hydrophobic or LWC<=config.LWClimit)
                  for (iphase=0;iphase<nphase;iphase++) 
                    {
		      if (surrogate[i].kp_from_experiment)
			{
			  if (iphase==surrogate[i].jmain_phase)
			    Kp(iphase)=surrogate[i].kpi;
			  else
			    Kp(iphase)=0.;
			}
		      else
			Kp(iphase)=surrogate[i].kpi/MOW(iphase)/surrogate[i].gamma_org_sat(iphase);
                      sum+=Kp(iphase)*MOinit(iphase);
                    }
 
                if (surrogate[i].hydrophilic and LWC>config.LWClimit) 
                  {
                    Kaq=surrogate[i].Kp_eff_aqreal_ssh(config, Temperature, ionic, chp,surrogate[config.iHp].gamma_LR,
                                                   surrogate[config.iHp].gamma_SRMR,MMaq,fion1,fion2)
                      /surrogate[i].gamma_aq;
                    sum+=Kaq*AQinit;
                  }             

                if (surrogate[i].hydrophilic and surrogate[i].hydrophobic and LWC>config.LWClimit) 
                  {               
                    fraciaq=surrogate[i].GAMMAinf*surrogate[i].gamma_aq*Kaq*AQinit/sum*Xmonoaq; //+surrogate[i].gamma_org*Kp*MOinit/sum*Xmono;                
                    for (iphase=0;iphase<nphase;iphase++) 
                      fraci(iphase)=surrogate[i].gamma_org_sat(iphase)*Kp(iphase)*MOinit(iphase)/sum*Xmono(iphase);
                  }
                else if (surrogate[i].hydrophobic or LWC<=config.LWClimit)                
                  for (iphase=0;iphase<nphase;iphase++) 
                    fraci(iphase)=surrogate[i].gamma_org_sat(iphase)*Kp(iphase)*MOinit(iphase)/sum*Xmono(iphase);                  
                else if (surrogate[i].hydrophilic)  
                  fraciaq=surrogate[i].GAMMAinf*surrogate[i].gamma_aq*Kaq*AQinit/sum*Xmonoaq;                                             
              }

            if (surrogate[j].nonvolatile)
              {
                double sum=0.0;
                if (surrogate[j].hydrophobic or LWC<=config.LWClimit)
                  for (iphase=0;iphase<nphase;iphase++)                                   
                    sum+=MOinit(iphase)/MOW(iphase)/surrogate[j].gamma_org_sat(iphase);                    
 
                if (surrogate[j].hydrophilic and LWC>config.LWClimit) 
                  sum+=AQinit/(MMaq*surrogate[j].gamma_aq*surrogate[j].GAMMAinf);

                if (sum>0.0)
                  {
                    if (surrogate[j].hydrophobic or LWC<=config.LWClimit)
                      for (iphase=0;iphase<nphase;iphase++)                                   
                        fraci(iphase)=MOinit(iphase)/MOW(iphase)/sum*Xmono(iphase);                    
                    
                    if (surrogate[j].hydrophilic and LWC>config.LWClimit) 
                      fraciaq=AQinit/MMaq/sum*Xmonoaq;                 
                  }
              }
            else
              {
                double Kaq=0.0;
                Array<double, 1> Kp;
                Kp.resize(nphase);                
                double sum=1.0;
                
                if (surrogate[j].hydrophobic or LWC<=config.LWClimit)
                  for (iphase=0;iphase<nphase;iphase++) 
                    {
		      if (surrogate[i].kp_from_experiment)
			{
			  if (iphase==surrogate[i].jmain_phase)
			    Kp(iphase)=surrogate[j].kpi;
			  else
			    Kp(iphase)=0.;
			}
		      else
			Kp(iphase)=surrogate[j].kpi/MOW(iphase)/surrogate[j].gamma_org_sat(iphase);
                      sum+=Kp(iphase)*MOinit(iphase);
                    }
 
                if (surrogate[j].hydrophilic and LWC>config.LWClimit) 
                  {
                    Kaq=surrogate[j].Kp_eff_aqreal_ssh(config, Temperature, ionic, chp,surrogate[config.iHp].gamma_LR,
                                                   surrogate[config.iHp].gamma_SRMR,MMaq,fion1,fion2)
                      /surrogate[j].gamma_aq;
                    sum+=Kaq*AQinit;
                  }             

                if (surrogate[j].hydrophilic and surrogate[j].hydrophobic and LWC>config.LWClimit) 
                  {               
                    fraciaq=surrogate[j].GAMMAinf*surrogate[j].gamma_aq*Kaq*AQinit/sum*Xmonoaq; //+surrogate[j].gamma_org*Kp*MOinit/sum*Xmono;                
                    for (iphase=0;iphase<nphase;iphase++) 
                      fraci(iphase)=surrogate[j].gamma_org_sat(iphase)*Kp(iphase)*MOinit(iphase)/sum*Xmono(iphase);
                  }
                else if (surrogate[j].hydrophobic or LWC<=config.LWClimit)                
                  for (iphase=0;iphase<nphase;iphase++) 
                    fraci(iphase)=surrogate[j].gamma_org_sat(iphase)*Kp(iphase)*MOinit(iphase)/sum*Xmono(iphase);                  
                else if (surrogate[j].hydrophilic)  
                  fraciaq=surrogate[j].GAMMAinf*surrogate[j].gamma_aq*Kaq*AQinit/sum*Xmonoaq;                                             
              }
            
            double flux=surrogate[i].koligo*(surrogate[i].Atot*fraciaq-surrogate[j].Atot*fracjaq/Kaq2)*deltat;
	    if (surrogate[i].catalyzed_ph and config.isorropia_ph==false)
	      flux=flux*chp*surrogate[config.iHp].gamma_aq/config.chp_org_ref;
	    else if (surrogate[i].catalyzed_ph and config.isorropia_ph==true)
	      flux=flux*chp/config.chp_org_ref;
	    
            for (iphase=0;iphase<nphase;iphase++) 
              flux+=surrogate[i].koligo*(surrogate[i].Atot*fraci(iphase)-surrogate[j].Atot*fracj(iphase)/Keq2(iphase))*deltat;

            surrogate[i].flux_chem_tot(index)+=-flux;                        
            surrogate[j].flux_chem_tot(index)+=flux;			                                          
          }

        if (surrogate[i].rion and surrogate[i].Aaq>0.0) // and config.compute_inorganic)
          for (jion=0;jion<surrogate[i].nion_chem;jion++)
            {            
              double molality=surrogate[surrogate[i].iion(jion)].gamma_aq*surrogate[surrogate[i].iion(jion)].Aaq/surrogate[i].MM/conc_org*1000.0;
              double flux=surrogate[i].kion[jion]*molality*deltat*surrogate[i].GAMMAinf*surrogate[i].gamma_aq*surrogate[i].Aaq;	      
	      if (surrogate[i].rion_ph_catalyzed[jion] and config.isorropia_ph==false)
		flux=flux*chp*surrogate[config.iHp].gamma_aq;
	      else if (surrogate[i].rion_ph_catalyzed[jion] and config.isorropia_ph==true)
		flux=flux*chp;
	      if (surrogate[i].rion_water_catalyzed[jion])
		flux=flux*RH;
            
              surrogate[i].flux_chem_tot(index)+=-flux;
	      if (surrogate[i].rion_catalyzed[jion]==false and molality>0.)
                {
		  if (surrogate[i].iion(jion)==config.iHSO4m or surrogate[i].iion(jion)==config.iSO4mm)
                    {
                      double totsulf=surrogate[config.iHSO4m].Aaq/surrogate[config.iHSO4m].MM+surrogate[config.iSO4mm].Aaq/surrogate[config.iSO4mm].MM;
                      surrogate[config.iHSO4m].flux_chem_tot(index)-=surrogate[config.iHSO4m].Aaq/surrogate[config.iHSO4m].MM/totsulf*flux/surrogate[i].MM;
                      surrogate[config.iSO4mm].flux_chem_tot(index)-=surrogate[config.iSO4mm].Aaq/surrogate[config.iSO4mm].MM/totsulf*flux/surrogate[i].MM;
                    }
                  if (surrogate[i].iion(jion)==config.iNO3m) // and config.compute_inorganic)
                    {
                      double totno3=surrogate[config.iNO3m].Aaq+surrogate[config.iHNO3].Ag/surrogate[config.iHNO3].MM*surrogate[config.iNO3m].MM;
                      surrogate[config.iNO3m].flux_chem_tot(index)-=surrogate[config.iNO3m].Aaq/totno3*flux/surrogate[i].MM*surrogate[config.iNO3m].MM;
                      surrogate[config.iHNO3].flux_chem_tot(index)-=surrogate[config.iHNO3].Ag/surrogate[config.iHNO3].MM*surrogate[config.iNO3m].MM/totno3*flux/surrogate[i].MM*surrogate[config.iHNO3].MM;
                    }
		  /*
		  else if (surrogate[i].iion(jion)==config.iNO3m and config.inorganic_chemistry)
		    {
                      surrogate[config.iNO3m].flux_chem_tot(index)-=flux/surrogate[i].MM*surrogate[config.iNO3m].MM;
		      }*/
                  if (surrogate[i].iion(jion)==config.iClm) // and config.compute_organic)
                    {
                      double totcl=surrogate[config.iClm].Aaq+surrogate[config.iHCl].Ag/surrogate[config.iHCl].MM*surrogate[config.iClm].MM;
                      surrogate[config.iClm].flux_chem_tot(index)-=surrogate[config.iClm].Aaq/totcl*flux/surrogate[i].MM*surrogate[config.iClm].MM;
                      surrogate[config.iHCl].flux_chem_tot(index)-=surrogate[config.iHCl].Ag/surrogate[config.iHCl].MM/totcl*flux/surrogate[i].MM*surrogate[config.iHCl].MM;
                    }
		  /*
		  else if (surrogate[i].iion(jion)==config.iClm and config.inorganic_chemistry)
		    {
                      surrogate[config.iClm].flux_chem_tot(index)-=flux/surrogate[i].MM*surrogate[config.iClm].MM;
                      //flux2+=flux/surrogate[i].MM*surrogate[config.iClm].MM;
		      }*/
                  if (surrogate[i].iion(jion)==config.iNH4p) // and config.compute_organic)
                    {
                      double totnh4=surrogate[config.iNH4p].Aaq+surrogate[config.iNH3].Ag/surrogate[config.iNH3].MM*surrogate[config.iNH4p].MM;
                      surrogate[config.iNH4p].flux_chem_tot(index)-=surrogate[config.iNH4p].Aaq/totnh4*flux/surrogate[i].MM*surrogate[config.iNH4p].MM;
                      surrogate[config.iNH3].flux_chem_tot(index)-=surrogate[config.iNH3].Ag/totnh4*flux/surrogate[config.iNH3].MM/surrogate[i].MM*surrogate[config.iNH3].MM;
                      //flux2+=flux/surrogate[i].MM*surrogate[config.iNH4p].MM;
                    }
		}

              surrogate[surrogate[i].iproduct(jion)].flux_chem_tot(index)+=flux*surrogate[surrogate[i].iproduct(jion)].MM/surrogate[i].MM;
            }

	for (jmol=0;jmol<surrogate[i].n_irreversible;jmol++)
	  if (surrogate[i].i_irreversible[jmol]>=0)
	    {
	      double flux=0.;
	      double fluxloc=0.;
	      if (surrogate[i].hydrophilic)
		{		
		  fluxloc=surrogate[i].k_irreversible[jmol]*deltat*surrogate[i].GAMMAinf*surrogate[i].gamma_aq*surrogate[i].Aaq;
		  if (surrogate[i].irr_catalyzed_water[jmol])
		    fluxloc=fluxloc*RH;
		  if (surrogate[i].irr_catalyzed_pH[jmol] and config.isorropia_ph==false)
		    fluxloc=fluxloc*chp*surrogate[config.iHp].gamma_aq;
		  else if (surrogate[i].irr_catalyzed_pH[jmol] and config.isorropia_ph==true)
		    fluxloc=fluxloc*chp;
		  flux+=fluxloc;
		}
	      if (surrogate[i].hydrophobic)
		{
		  fluxloc=0.;
		  for (iphase=0;iphase<nphase;iphase++)
		    fluxloc+=surrogate[i].k_irreversible[jmol]*deltat*surrogate[i].gamma_org*surrogate[i].Ap_sat(iphase);

		  if (surrogate[i].irr_catalyzed_water[jmol])
		    fluxloc=fluxloc*RH;
		  if (surrogate[i].irr_catalyzed_pH[jmol])
		    fluxloc=fluxloc*config.chp_org_ref;
		  flux+=fluxloc;
		}

	      surrogate[i].flux_chem_tot(index)+=-flux;
	      if (surrogate[i].irr_mass_conserving[jmol])
		surrogate[surrogate[i].i_irreversible[jmol]].flux_chem_tot(index)+=flux;
	      else
		surrogate[surrogate[i].i_irreversible[jmol]].flux_chem_tot(index)+=flux*surrogate[surrogate[i].i_irreversible[jmol]].MM/surrogate[i].MM;

	      if (surrogate[i].iother_irreversible[jmol]>=0)
		if (surrogate[i].iother_irreversible[jmol]==config.iNO3m)
		  {
		    double totno3=surrogate[config.iNO3m].Aaq+surrogate[config.iHNO3].Ag/surrogate[config.iHNO3].MM*surrogate[config.iNO3m].MM;
		    if (totno3>0.)
		      {
			surrogate[surrogate[i].iother_irreversible[jmol]].flux_chem_tot(index)+=surrogate[config.iNO3m].Aaq/totno3*flux*surrogate[surrogate[i].iother_irreversible[jmol]].MM/surrogate[i].MM;
			surrogate[config.iHNO3].flux_chem_tot(index)+=surrogate[config.iHNO3].Ag/surrogate[config.iHNO3].MM*surrogate[config.iNO3m].MM/totno3*flux*surrogate[config.iHNO3].MM/surrogate[i].MM;
		      }
		  }
		else if (surrogate[i].iother_irreversible[jmol]==config.iClm)
		  {
		    double totcl=surrogate[config.iClm].Aaq+surrogate[config.iHCl].Ag/surrogate[config.iHCl].MM*surrogate[config.iClm].MM;
		    if (totcl>0.)
		      {
			surrogate[surrogate[i].iother_irreversible[jmol]].flux_chem_tot(index)+=surrogate[config.iClm].Aaq/totcl*flux*surrogate[surrogate[i].iother_irreversible[jmol]].MM/surrogate[i].MM;
			surrogate[config.iHCl].flux_chem_tot(index)+=surrogate[config.iHCl].Ag/surrogate[config.iHCl].MM*surrogate[config.iClm].MM/totcl*flux*surrogate[config.iHCl].MM/surrogate[i].MM;
		      }
		  }
		else if (surrogate[i].iother_irreversible[jmol]==config.iNH4p)
		  {
		    double totnh4=surrogate[config.iNH4p].Aaq+surrogate[config.iNH3].Ag/surrogate[config.iNH3].MM*surrogate[config.iNH4p].MM;
		    if (totnh4>0.)
		      {
			surrogate[surrogate[i].iother_irreversible[jmol]].flux_chem_tot(index)+=surrogate[config.iNH4p].Aaq/totnh4*flux*surrogate[surrogate[i].iother_irreversible[jmol]].MM/surrogate[i].MM;
			surrogate[config.iNH3].flux_chem_tot(index)+=surrogate[config.iNH3].Ag/surrogate[config.iNH3].MM*surrogate[config.iNH4p].MM/totnh4*flux*surrogate[config.iNH3].MM/surrogate[i].MM;
		      }
		  }
		else if (surrogate[i].iother_irreversible[jmol]==config.iSO4mm or surrogate[i].iother_irreversible[jmol]==config.iHSO4m)
		  {
		    double totsulf=surrogate[config.iHSO4m].Aaq/surrogate[config.iHSO4m].MM+surrogate[config.iSO4mm].Aaq/surrogate[config.iSO4mm].MM;
		    if (totsulf>0.)
		      {
			surrogate[config.iSO4mm].flux_chem_tot(index)+=surrogate[config.iSO4mm].Aaq/surrogate[config.iSO4mm].MM/totsulf*flux*surrogate[config.iSO4mm].MM/surrogate[i].MM;
			surrogate[config.iHSO4m].flux_chem_tot(index)+=surrogate[config.iHSO4m].Aaq/surrogate[config.iHSO4m].MM/totsulf*flux*surrogate[config.iHSO4m].MM/surrogate[i].MM;
		      }
		  }
		else
		  surrogate[surrogate[i].iother_irreversible[jmol]].flux_chem_tot(index)+=flux*surrogate[surrogate[i].iother_irreversible[jmol]].MM/surrogate[i].MM;
	    }
	
      }
   
  for (i=0;i<n;++i)
    {           
      if (surrogate[i].is_organic)
        if (surrogate[i].Atot>tiny)
          surrogate[i].Jdn_gas(index)=min(surrogate[i].flux_chem_tot(index)/surrogate[i].Atot/deltat,0.0);     
        else if (i==config.iHSO4m)
          {
            double totsulf=surrogate[config.iHSO4m].Aaq/surrogate[config.iHSO4m].MM+surrogate[config.iSO4mm].Aaq/surrogate[config.iSO4mm].MM;
            if (totsulf>0.0)
              surrogate[i].Jdn_gas(index)=min((surrogate[config.iHSO4m].flux_chem_tot(index)/surrogate[config.iHSO4m].MM
                                               +surrogate[config.iSO4mm].flux_chem_tot(index)/surrogate[config.iSO4mm].MM)/totsulf/deltat,0.0);
            
          }
        else if (i==config.iNO3m)
          {
            double totno3=surrogate[config.iNO3m].Aaq+surrogate[config.iHNO3].Ag/surrogate[config.iHNO3].MM*surrogate[config.iNO3m].MM;
            if (totno3>0.0)
              surrogate[i].Jdn_gas(index)=min((surrogate[config.iNO3m].flux_chem_tot(index)+surrogate[config.iHNO3].flux_chem_tot(index)/surrogate[config.iHNO3].MM*surrogate[config.iNO3m].MM)/
                                              totno3/deltat,0.0);            
          }
        else if (i==config.iClm)
          {
            double totcl=surrogate[config.iClm].Aaq+surrogate[config.iHCl].Ag/surrogate[config.iHCl].MM*surrogate[config.iClm].MM;
            if (totcl>0.0)
              surrogate[i].Jdn_gas(index)=min((surrogate[config.iClm].flux_chem_tot(index)+surrogate[config.iHCl].flux_chem_tot(index)/surrogate[config.iHCl].MM*surrogate[config.iClm].MM)/
                                              totcl/deltat,0.0);            
          }
        else if (i==config.iNH4p)
          {
            double totnh4=surrogate[config.iNH4p].Aaq+surrogate[config.iNH3].Ag/surrogate[config.iNH3].MM*surrogate[config.iNH4p].MM;
            if (totnh4>0.0)
              surrogate[i].Jdn_gas(index)=min((surrogate[config.iNH4p].flux_chem_tot(index)+surrogate[config.iNH3].flux_chem_tot(index)/surrogate[config.iNH3].MM*surrogate[config.iNH4p].MM)/
                                              totnh4/deltat,0.0);            
          }
    }

  for (i=0;i<n;++i)    
    surrogate[i].flux_chem_tot(index)=0.0;   

  for (i=0;i<n;++i)
    if (surrogate[i].is_organic)
      {
        if (surrogate[i].is_monomer and surrogate[i].Atot>0.0)
          {
            j=surrogate[i].ioligo; 
            Array <double, 1> Keq2;
            Keq2.resize(nphase);
            for (iphase=0;iphase<nphase;iphase++)
              Keq2(iphase)=surrogate[j].MM/surrogate[i].MM*pow(surrogate[i].Keq_oligo,surrogate[i].moligo-1)*pow(max(Xmono(iphase),xmin),surrogate[i].moligo-2)/pow(max(XH2O,xmin),surrogate[i].moligo-1);                      
            double Kaq2=surrogate[j].MM/surrogate[i].MM*pow(surrogate[i].Keq_oligo,surrogate[i].moligo-1)*pow(max(Xmonoaq,xmin),surrogate[i].moligo-2)/pow(max(XH2O,xmin),surrogate[i].moligo-1);  
            Array <double, 1> fraci,fracj;
            fraci.resize(nphase);
            fracj.resize(nphase);
            fraci=0.0;
            fracj=0.0;            
            double fraciaq=0.0;
            double fracjaq=0.0;
            if (surrogate[i].nonvolatile)
              {
                double sum=0.0;
                if (surrogate[i].hydrophobic or LWC<=config.LWClimit)
                  for (iphase=0;iphase<nphase;iphase++)                                   
                    sum+=MOinit(iphase)/MOW(iphase)/surrogate[i].gamma_org_sat(iphase);                    
 
                if (surrogate[i].hydrophilic and LWC>config.LWClimit) 
                  sum+=AQinit/(MMaq*surrogate[i].gamma_aq*surrogate[i].GAMMAinf);

                if (sum>0.0)
                  {
                    if (surrogate[i].hydrophobic or LWC<=config.LWClimit)
                      for (iphase=0;iphase<nphase;iphase++)                                   
                        fraci(iphase)=MOinit(iphase)/MOW(iphase)/sum*Xmono(iphase);                    
                    
                    if (surrogate[i].hydrophilic and LWC>config.LWClimit) 
                      fraciaq=AQinit/MMaq/sum*Xmonoaq;                 
                  }
              }
            else
              {
                double Kaq=0.0;
                Array<double, 1> Kp;
                Kp.resize(nphase);                
                double sum=1.0;
                
                if (surrogate[i].hydrophobic or LWC<=config.LWClimit)
                  for (iphase=0;iphase<nphase;iphase++) 
                    {
		      if (surrogate[i].kp_from_experiment)
			{
			  if (iphase==surrogate[i].jmain_phase)
			    Kp(iphase)=surrogate[i].kpi;
			  else
			    Kp(iphase)=0.;
			}
		      else
			Kp(iphase)=surrogate[i].kpi/MOW(iphase)/surrogate[i].gamma_org_sat(iphase);
                      sum+=Kp(iphase)*MOinit(iphase);
                    }
 
                if (surrogate[i].hydrophilic and LWC>config.LWClimit) 
                  {
                    Kaq=surrogate[i].Kp_eff_aqreal_ssh(config, Temperature, ionic, chp,surrogate[config.iHp].gamma_LR,
                                                   surrogate[config.iHp].gamma_SRMR,MMaq,fion1,fion2)
                      /surrogate[i].gamma_aq;
                    sum+=Kaq*AQinit;
                  }             

                if (surrogate[i].hydrophilic and surrogate[i].hydrophobic and LWC>config.LWClimit) 
                  {               
                    fraciaq=surrogate[i].GAMMAinf*surrogate[i].gamma_aq*Kaq*AQinit/sum*Xmonoaq; //+surrogate[i].gamma_org*Kp*MOinit/sum*Xmono;                
                    for (iphase=0;iphase<nphase;iphase++) 
                      fraci(iphase)=surrogate[i].gamma_org_sat(iphase)*Kp(iphase)*MOinit(iphase)/sum*Xmono(iphase);
                  }
                else if (surrogate[i].hydrophobic or LWC<=config.LWClimit)                
                  for (iphase=0;iphase<nphase;iphase++) 
                    fraci(iphase)=surrogate[i].gamma_org_sat(iphase)*Kp(iphase)*MOinit(iphase)/sum*Xmono(iphase);                  
                else if (surrogate[i].hydrophilic)  
                  fraciaq=surrogate[i].GAMMAinf*surrogate[i].gamma_aq*Kaq*AQinit/sum*Xmonoaq;                                             
              }

            if (surrogate[j].nonvolatile)
              {
                double sum=0.0;
                if (surrogate[j].hydrophobic or LWC<=config.LWClimit)
                  for (iphase=0;iphase<nphase;iphase++)                                   
                    sum+=MOinit(iphase)/MOW(iphase)/surrogate[j].gamma_org_sat(iphase);                    
 
                if (surrogate[j].hydrophilic and LWC>config.LWClimit) 
                  sum+=AQinit/(MMaq*surrogate[j].gamma_aq*surrogate[j].GAMMAinf);

                if (sum>0.0)
                  {
                    if (surrogate[j].hydrophobic or LWC<=config.LWClimit)
                      for (iphase=0;iphase<nphase;iphase++)                                   
                        fraci(iphase)=MOinit(iphase)/MOW(iphase)/sum*Xmono(iphase);                    
                    
                    if (surrogate[j].hydrophilic and LWC>config.LWClimit) 
                      fraciaq=AQinit/MMaq/sum*Xmonoaq;                 
                  }
              }
            else
              {
                double Kaq=0.0;
                Array<double, 1> Kp;
                Kp.resize(nphase);                
                double sum=1.0;
                
                if (surrogate[j].hydrophobic or LWC<=config.LWClimit)
                  for (iphase=0;iphase<nphase;iphase++) 
                    {
		      if (surrogate[i].kp_from_experiment)
			{
			  if (iphase==surrogate[i].jmain_phase)
			    Kp(iphase)=surrogate[j].kpi;
			  else
			    Kp(iphase)=0.;
			}
		      else
			Kp(iphase)=surrogate[j].kpi/MOW(iphase)/surrogate[j].gamma_org_sat(iphase);
                      sum+=Kp(iphase)*MOinit(iphase);
                    }
 
                if (surrogate[j].hydrophilic and LWC>config.LWClimit) 
                  {
                    Kaq=surrogate[j].Kp_eff_aqreal_ssh(config, Temperature, ionic, chp,surrogate[config.iHp].gamma_LR,
                                                   surrogate[config.iHp].gamma_SRMR,MMaq,fion1,fion2)
                      /surrogate[j].gamma_aq;
                    sum+=Kaq*AQinit;
                  }             

                if (surrogate[j].hydrophilic and surrogate[j].hydrophobic and LWC>config.LWClimit) 
                  {               
                    fraciaq=surrogate[j].GAMMAinf*surrogate[j].gamma_aq*Kaq*AQinit/sum*Xmonoaq; //+surrogate[j].gamma_org*Kp*MOinit/sum*Xmono;                
                    for (iphase=0;iphase<nphase;iphase++) 
                      fraci(iphase)=surrogate[j].gamma_org_sat(iphase)*Kp(iphase)*MOinit(iphase)/sum*Xmono(iphase);
                  }
                else if (surrogate[j].hydrophobic or LWC<=config.LWClimit)                
                  for (iphase=0;iphase<nphase;iphase++) 
                    fraci(iphase)=surrogate[j].gamma_org_sat(iphase)*Kp(iphase)*MOinit(iphase)/sum*Xmono(iphase);                  
                else if (surrogate[j].hydrophilic)  
                  fraciaq=surrogate[j].GAMMAinf*surrogate[j].gamma_aq*Kaq*AQinit/sum*Xmonoaq;                                             
              }
            
            double flux=surrogate[i].koligo*(surrogate[i].Atot*fraciaq-surrogate[j].Atot*fracjaq/Kaq2)*deltat;
	    if (surrogate[i].catalyzed_ph and config.isorropia_ph==false)
	      flux=flux*chp*surrogate[config.iHp].gamma_aq/config.chp_org_ref;
	    else if (surrogate[i].catalyzed_ph and config.isorropia_ph==true)
	      flux=flux*chp/config.chp_org_ref;
	    
            for (iphase=0;iphase<nphase;iphase++) 
              flux+=surrogate[i].koligo*(surrogate[i].Atot*fraci(iphase)-surrogate[j].Atot*fracj(iphase)/Keq2(iphase))*deltat;

            if (flux>0.0)                        
              flux=flux/(1.0-gamma*surrogate[i].Jdn_gas(index)*deltat);                                                  
            else
              flux=flux/(1.0-gamma*surrogate[j].Jdn_gas(index)*deltat);

            surrogate[i].flux_chem_tot(index)+=-flux;                        
            surrogate[j].flux_chem_tot(index)+=flux;			                                          
          }

        if (surrogate[i].rion and surrogate[i].Aaq>0.0) // and config.compute_inorganic)
          for (jion=0;jion<surrogate[i].nion_chem;jion++)
            {              
              double molality=surrogate[surrogate[i].iion(jion)].gamma_aq*surrogate[surrogate[i].iion(jion)].Aaq/surrogate[i].MM/conc_org*1000.0;
              double flux=surrogate[i].kion[jion]*molality*deltat*surrogate[i].GAMMAinf*surrogate[i].gamma_aq*surrogate[i].Aaq;	      
	      if (surrogate[i].rion_ph_catalyzed[jion] and config.isorropia_ph==false)
		flux=flux*chp*surrogate[config.iHp].gamma_aq;
	      else if (surrogate[i].rion_ph_catalyzed[jion] and config.isorropia_ph==true)
		flux=flux*chp;
	      
	      if (surrogate[i].rion_water_catalyzed[jion])
		flux=flux*RH;
if (surrogate[i].rion_catalyzed[jion]==false and molality>0.)
                {
		  
                  if (surrogate[i].iion(jion)==config.iHSO4m or surrogate[i].iion(jion)==config.iSO4mm)
                    {
                      flux=flux/(1.0-gamma*surrogate[i].Jdn_gas(index)*deltat
                                 -gamma*surrogate[config.iHSO4m].Jdn_gas(index)/surrogate[config.iHSO4m].MM*surrogate[i].MM*deltat); 
                      double totsulf=surrogate[config.iHSO4m].Aaq/surrogate[config.iHSO4m].MM+surrogate[config.iSO4mm].Aaq/surrogate[config.iSO4mm].MM;
                      surrogate[config.iHSO4m].flux_chem_tot(index)-=surrogate[config.iHSO4m].Aaq/surrogate[config.iHSO4m].MM/totsulf*flux/surrogate[i].MM;
                      surrogate[config.iSO4mm].flux_chem_tot(index)-=surrogate[config.iSO4mm].Aaq/surrogate[config.iSO4mm].MM/totsulf*flux/surrogate[i].MM;
                    }
                  else if (surrogate[i].iion(jion)==config.iNO3m) // and config.compute_inorganic)
                    {
                      flux=flux/(1.0-gamma*surrogate[i].Jdn_gas(index)*deltat
                                 -gamma*surrogate[config.iNO3m].Jdn_gas(index)/surrogate[config.iNO3m].MM*surrogate[i].MM*deltat); 
                      double totno3=surrogate[config.iNO3m].Aaq+surrogate[config.iHNO3].Ag/surrogate[config.iHNO3].MM*surrogate[config.iNO3m].MM;
                      surrogate[config.iNO3m].flux_chem_tot(index)-=surrogate[config.iNO3m].Aaq/totno3*flux/surrogate[i].MM*surrogate[config.iNO3m].MM;
                      surrogate[config.iHNO3].flux_chem_tot(index)-=surrogate[config.iHNO3].Ag/surrogate[config.iHNO3].MM*surrogate[config.iNO3m].MM/totno3*flux/surrogate[i].MM*surrogate[config.iHNO3].MM;
                    }
		  /*
		  else if (surrogate[i].iion(jion)==config.iNO3m and config.inorganic_chemistry)
                    {
                      flux=flux/(1.0-gamma*surrogate[i].Jdn_gas(index)*deltat
                                 -gamma*surrogate[config.iNO3m].Jdn_gas(index)/surrogate[config.iNO3m].MM*surrogate[i].MM*deltat); 
                      surrogate[config.iNO3m].flux_chem_tot(index)-=flux/surrogate[i].MM*surrogate[config.iNO3m].MM;
		      }*/
                  else if (surrogate[i].iion(jion)==config.iClm) // and config.compute_inorganic)
                    {
                      flux=flux/(1.0-gamma*surrogate[i].Jdn_gas(index)*deltat
                                 -gamma*surrogate[config.iClm].Jdn_gas(index)/surrogate[config.iClm].MM*surrogate[i].MM*deltat); 
                      double totcl=surrogate[config.iClm].Aaq+surrogate[config.iHCl].Ag/surrogate[config.iHCl].MM*surrogate[config.iClm].MM;
                      surrogate[config.iClm].flux_chem_tot(index)-=surrogate[config.iClm].Aaq/totcl*flux/surrogate[i].MM*surrogate[config.iClm].MM;
                      surrogate[config.iHCl].flux_chem_tot(index)-=surrogate[config.iHCl].Ag/surrogate[config.iHCl].MM*surrogate[config.iClm].MM/totcl*flux/surrogate[i].MM*surrogate[config.iHCl].MM;
                    }
		  /*
		  else if (surrogate[i].iion(jion)==config.iClm and config.inorganic_chemistry)
                    {
                      flux=flux/(1.0-gamma*surrogate[i].Jdn_gas(index)*deltat
                                 -gamma*surrogate[config.iClm].Jdn_gas(index)/surrogate[config.iClm].MM*surrogate[i].MM*deltat); 
                      surrogate[config.iClm].flux_chem_tot(index)-=flux/surrogate[i].MM*surrogate[config.iClm].MM;
		      }*/
                  else if (surrogate[i].iion(jion)==config.iNH4p) // and config.compute_inorganic)
                    {
                      flux=flux/(1.0-gamma*surrogate[i].Jdn_gas(index)*deltat
                                 -gamma*surrogate[config.iNH4p].Jdn_gas(index)/surrogate[config.iNH4p].MM*surrogate[i].MM*deltat); 
                      double totnh4=surrogate[config.iNH4p].Aaq+surrogate[config.iNH3].Ag/surrogate[config.iNH3].MM*surrogate[config.iNH4p].MM;
                      surrogate[config.iNH4p].flux_chem_tot(index)-=surrogate[config.iNH4p].Aaq/totnh4*flux/surrogate[i].MM*surrogate[config.iNH4p].MM;
                      surrogate[config.iNH3].flux_chem_tot(index)-=surrogate[config.iNH3].Ag/surrogate[config.iNH3].MM*surrogate[config.iNH4p].MM/totnh4*flux/surrogate[i].MM*surrogate[config.iNH3].MM;
                    }
		  /*
		  else if (surrogate[i].iion(jion)==config.iNH4p) and config.inorganic_chemistry)
                    {
                      flux=flux/(1.0-gamma*surrogate[i].Jdn_gas(index)*deltat
                                 -gamma*surrogate[config.iNH4p].Jdn_gas(index)/surrogate[config.iNH4p].MM*surrogate[i].MM*deltat); 
                      surrogate[config.iNH4p].flux_chem_tot(index)-=flux/surrogate[i].MM*surrogate[config.iNH4p].MM;
		      }*/

                }
              else
                flux=flux/(1.0-gamma*surrogate[i].Jdn_gas(index)*deltat);

              surrogate[i].flux_chem_tot(index)+=-flux; 
              surrogate[surrogate[i].iproduct(jion)].flux_chem_tot(index)+=flux*surrogate[surrogate[i].iproduct(jion)].MM/surrogate[i].MM;
            }

	for (jmol=0;jmol<surrogate[i].n_irreversible;jmol++)
	  if (surrogate[i].i_irreversible[jmol]>=0)
	    {
	      double flux=0.;
	      double fluxloc=0.;
	      if (surrogate[i].hydrophilic)
		{		
		  fluxloc=surrogate[i].k_irreversible[jmol]*deltat*surrogate[i].GAMMAinf*surrogate[i].gamma_aq*surrogate[i].Aaq;
		  if (surrogate[i].irr_catalyzed_water[jmol])
		    fluxloc=fluxloc*RH;
		  if (surrogate[i].irr_catalyzed_pH[jmol] and config.isorropia_ph==false)
		    fluxloc=fluxloc*chp*surrogate[config.iHp].gamma_aq;
		  else if (surrogate[i].irr_catalyzed_pH[jmol] and config.isorropia_ph==true)
		    fluxloc=fluxloc*chp;
		  flux+=fluxloc;
		}
	      if (surrogate[i].hydrophobic)
		{
		  fluxloc=0.;
		  for (iphase=0;iphase<nphase;iphase++)
		    fluxloc+=surrogate[i].k_irreversible[jmol]*deltat*surrogate[i].gamma_org*surrogate[i].Ap_sat(iphase);

		  if (surrogate[i].irr_catalyzed_water[jmol])
		    fluxloc=fluxloc*RH;
		  if (surrogate[i].irr_catalyzed_pH[jmol])
		    fluxloc=fluxloc*config.chp_org_ref;
		  flux+=fluxloc;
		}
	     
	      flux=flux/(1.0-gamma*surrogate[i].Jdn_gas(index)*deltat);         
	      surrogate[i].flux_chem_tot(index)+=-flux;
	      if (surrogate[i].irr_mass_conserving[jmol])
		surrogate[surrogate[i].i_irreversible[jmol]].flux_chem_tot(index)+=flux;
	      else
		surrogate[surrogate[i].i_irreversible[jmol]].flux_chem_tot(index)+=flux*surrogate[surrogate[i].i_irreversible[jmol]].MM/surrogate[i].MM;

	      if (surrogate[i].iother_irreversible[jmol]>=0)
		{
		  if (surrogate[i].iother_irreversible[jmol]==config.iNO3m)
		    {
		      double totno3=surrogate[config.iNO3m].Aaq+surrogate[config.iHNO3].Ag/surrogate[config.iHNO3].MM*surrogate[config.iNO3m].MM;
		      if (totno3>0.)
			{
			  surrogate[surrogate[i].iother_irreversible[jmol]].flux_chem_tot(index)+=surrogate[config.iNO3m].Aaq/totno3*flux*surrogate[surrogate[i].iother_irreversible[jmol]].MM/surrogate[i].MM;
			  surrogate[config.iHNO3].flux_chem_tot(index)+=surrogate[config.iHNO3].Ag/surrogate[config.iHNO3].MM*surrogate[config.iNO3m].MM/totno3*flux*surrogate[config.iHNO3].MM/surrogate[i].MM;
			}
		    }
		  else if (surrogate[i].iother_irreversible[jmol]==config.iClm)
		    {
		      double totcl=surrogate[config.iClm].Aaq+surrogate[config.iHCl].Ag/surrogate[config.iHCl].MM*surrogate[config.iClm].MM;
		      if (totcl>0.)
			{
			  surrogate[surrogate[i].iother_irreversible[jmol]].flux_chem_tot(index)+=surrogate[config.iClm].Aaq/totcl*flux*surrogate[surrogate[i].iother_irreversible[jmol]].MM/surrogate[i].MM;
			  surrogate[config.iHCl].flux_chem_tot(index)+=surrogate[config.iHCl].Ag/surrogate[config.iHCl].MM*surrogate[config.iClm].MM/totcl*flux*surrogate[config.iHCl].MM/surrogate[i].MM;
			}
		    }
		  else if (surrogate[i].iother_irreversible[jmol]==config.iNH4p)
		    {
		      double totnh4=surrogate[config.iNH4p].Aaq+surrogate[config.iNH3].Ag/surrogate[config.iNH3].MM*surrogate[config.iNH4p].MM;
		      if (totnh4>0.)
			{
			  surrogate[surrogate[i].iother_irreversible[jmol]].flux_chem_tot(index)+=surrogate[config.iNH4p].Aaq/totnh4*flux*surrogate[surrogate[i].iother_irreversible[jmol]].MM/surrogate[i].MM;
			  surrogate[config.iNH3].flux_chem_tot(index)+=surrogate[config.iNH3].Ag/surrogate[config.iNH3].MM*surrogate[config.iNH4p].MM/totnh4*flux*surrogate[config.iNH3].MM/surrogate[i].MM;
			}
		    }
		  else if (surrogate[i].iother_irreversible[jmol]==config.iSO4mm or surrogate[i].iother_irreversible[jmol]==config.iHSO4m)
		    {
		      double totsulf=surrogate[config.iHSO4m].Aaq/surrogate[config.iHSO4m].MM+surrogate[config.iSO4mm].Aaq/surrogate[config.iSO4mm].MM;
		      if (totsulf>0.)
			{
			  surrogate[config.iSO4mm].flux_chem_tot(index)+=surrogate[config.iSO4mm].Aaq/surrogate[config.iSO4mm].MM/totsulf*flux*surrogate[config.iSO4mm].MM/surrogate[i].MM;
			  surrogate[config.iHSO4m].flux_chem_tot(index)+=surrogate[config.iHSO4m].Aaq/surrogate[config.iHSO4m].MM/totsulf*flux*surrogate[config.iHSO4m].MM/surrogate[i].MM;
			}
		    }	       
		  else
		    surrogate[surrogate[i].iother_irreversible[jmol]].flux_chem_tot(index)+=flux*surrogate[surrogate[i].iother_irreversible[jmol]].MM/surrogate[i].MM;
		}

	    }
        
      }
} 

void integer_chem_ssh(model_config &config, vector<species>& surrogate,
                  double &MOinit,double &MOW, double &MMaq,
                  double &LWC, double &AQinit, double &ionic, double &chp,
                  double &Temperature, double &RH, double &deltat, bool &compute_activity_coefficients)
{
  //compute the dynamic evolution of the organic-phase and the aqueous-phase concentrations with
  //         the roschem algorythm
  //FOR THE COUPLED SYSTEM
  //compute_activity_coefficients: should the activity coefficients be recomputed?
  //tequilibrium: criterium to distinguish cases under equilibrium (characteristic time lower
  //                  than tequilibrium)
  int n=surrogate.size();
  int i; 
  
  //double tiny=0.; //1.0e-26;
  double conc_available=0.0;
  double sum_rates;
  //double gamma=1.0+pow(2.,0.5)/2.; //1.7071;
  //double Kaq, Kp;      
  double redmax=0.01;  
  bool all_hydrophobic=false;
  if (LWC<=config.LWClimit)
    all_hydrophobic=true;
  
  //     gamma > .25     ->  A-stable  ( x(n+1)-x(n))*f(x(n)) > 0
  //     gamma < .25     ->  monotonic behaviour (x(n+1)-x*)(x(n)-x*) > 0
  //     gamma =  1+-sqrt(1/2) ->  L-stability
          

  MOinit=0.0;
  AQinit=LWC;
  double fion1,fion2;
  double XH2O,viscosity;

  for (i=0;i<n;++i)
    {
      if (surrogate[i].hydrophobic or (all_hydrophobic and (surrogate[i].is_organic or i==config.iH2O)))        
        MOinit+=surrogate[i].Ap;        
      if (surrogate[i].hydrophilic)        
        AQinit+=surrogate[i].Aaq;          
    }

  double conc_org=LWC;
  for (i=0;i<n;++i)
    if (surrogate[i].is_organic or i==config.iH2O)     
      conc_org+=surrogate[i].Aaq;
      
  conc_org=max(conc_org,1.0e-5*config.MOmin);

  if (compute_activity_coefficients)
    {
      if (config.compute_organic)
        if (LWC>config.LWClimit)
          activity_coefficients_org_ssh(config, surrogate, false, Temperature, MOW, viscosity);
        else
          activity_coefficients_org_ssh(config, surrogate, true, Temperature, MOW, viscosity);          
      activity_coefficients_aq_ssh(config, surrogate, Temperature, LWC, MMaq, XH2O, conc_org);
      if (config.compute_long_and_medium_range_interactions)
        activity_coefficients_LR_MR_ssh(config, surrogate, Temperature, LWC, ionic);
    }
  
  /*
    double sumX=0.0;
    double sumXaq=0.0;
    for (i=0;i<n;++i)
    {
    surrogate[i].X=0.0;
    surrogate[i].Xaq=0.0;
    if (surrogate[i].hydrophobic or all_hydrophobic)
    {
    MOinit+=surrogate[i].Ap;
    surrogate[i].X=surrogate[i].Ap/surrogate[i].MM;
    sumX+=surrogate[i].X;
    }
    if (surrogate[i].hydrophilic)
    {
    AQinit+=surrogate[i].Aaq;
    surrogate[i].Xaq=surrogate[i].Aaq/surrogate[i].MM;
    sumXaq+=surrogate[i].Xaq;
    }
    }

    MOW=0.0;
    if (sumX>0.0)    
    {
    for (i=0;i<n;++i)    
    //if (surrogate[i].hydrophobic or all_hydrophobic)                       
    {
    surrogate[i].X/=sumX;        
    MOW+=surrogate[i].X*surrogate[i].MM;
    }
    }
    else
    MOW=280;
    
    MMaq=0.0;
    if (sumXaq>0.0)
    {
    for (i=0;i<n;++i)                   
    //if (surrogate[i].hydrophilic)
    {          
    surrogate[i].Xaq/=sumXaq;
    MMaq+=surrogate[i].Xaq*surrogate[i].MM;
    }        
    }*/

  //save initial concentrations 
  for (i=0;i<n;++i)    
    {
      surrogate[i].Atot0=surrogate[i].Atot;    
      surrogate[i].Aaq0=surrogate[i].Aaq;
      surrogate[i].Ag0=surrogate[i].Ag;
    }

  kinetic_ssh(config, surrogate, MOinit, MOW, MMaq, LWC, AQinit, ionic, chp, Temperature, RH, deltat, 0);  
  //compute the first evaluation of concentrations
  for (i=0;i<n;++i)
    if(surrogate[i].is_organic) 
      {
	conc_available=0.0;
	sum_rates=0.0;
	
        surrogate[i].Atot=max(redmax*surrogate[i].Atot,surrogate[i].Atot+surrogate[i].flux_chem_tot(0));	     
        surrogate[i].Atot1=surrogate[i].Atot;
        if (surrogate[i].nonvolatile)
          {
            if (surrogate[i].hydrophilic and surrogate[i].hydrophobic and LWC>config.LWClimit)
              {                                  
                double ratio=AQinit/max(MOinit,1.0e-10)*surrogate[i].gamma_org/(surrogate[i].GAMMAinf*surrogate[i].gamma_aq)*MOW/MMaq;
                surrogate[i].Ap=1.0/(1.0+ratio)*surrogate[i].Atot;
                surrogate[i].Aaq=ratio/(1.0+ratio)*surrogate[i].Atot;
              }
            else if (surrogate[i].hydrophobic or LWC<=config.LWClimit)                
              surrogate[i].Ap=surrogate[i].Atot;
            else if (surrogate[i].hydrophilic)      
              surrogate[i].Aaq=surrogate[i].Atot;            
          }
        else
          {
            double Kaq=0.0;
            double Kp=0.0;
            double sum=1.0;
            if (surrogate[i].hydrophobic or LWC<=config.LWClimit)
              {
		if (surrogate[i].kp_from_experiment)
		  Kp=surrogate[i].kpi;
		else
		  Kp=surrogate[i].kpi/MOW/surrogate[i].gamma_org;
                sum+=Kp*MOinit;
              }

            if (surrogate[i].hydrophilic and LWC>config.LWClimit) 
              {
                Kaq=surrogate[i].Kp_eff_aqreal_ssh(config, Temperature, ionic, chp,surrogate[config.iHp].gamma_LR,
                                               surrogate[config.iHp].gamma_SRMR,MMaq,fion1,fion2)
                  /surrogate[i].gamma_aq;
                sum+=Kaq*AQinit;
              }             
           
            if (surrogate[i].hydrophilic and surrogate[i].hydrophobic and LWC>config.LWClimit)                
              {
                surrogate[i].Ap=surrogate[i].Atot*Kp*MOinit/sum; 
                surrogate[i].Aaq=surrogate[i].Atot*Kaq*AQinit/sum;                 
              }                
            else if (surrogate[i].hydrophobic or LWC<=config.LWClimit)                
              surrogate[i].Ap=surrogate[i].Atot*Kp*MOinit/sum; 
            else if (surrogate[i].hydrophilic) 
              surrogate[i].Aaq=surrogate[i].Atot*Kaq*AQinit/sum;                                             
          }                                 
      }
    else if (i==config.iHSO4m or i==config.iSO4mm or i==config.iNO3m or i==config.iClm or i==config.iNH4p)
      {
	surrogate[i].Aaq=max(redmax*surrogate[i].Aaq,surrogate[i].Aaq+surrogate[i].flux_chem_tot(0));
	surrogate[i].Aaq1=surrogate[i].Aaq; 
	//cout << "computed " << surrogate[i].name << " " << surrogate[i].Aaq << " " << surrogate[i].Aaq0 << endl;
      }
    else if (i==config.iHNO3 or i==config.iHCl or i==config.iNH3)
      {
	surrogate[i].Ag=max(redmax*surrogate[i].Ag,surrogate[i].Ag+surrogate[i].flux_chem_tot(0));      
	surrogate[i].Ag1=surrogate[i].Ag;
      }

  MOinit=0.0;
  AQinit=LWC;

  for (i=0;i<n;++i)
    {
      if (surrogate[i].hydrophobic or (all_hydrophobic and (surrogate[i].is_organic or i==config.iH2O)))        
        MOinit+=surrogate[i].Ap;        
      if (surrogate[i].hydrophilic)        
        AQinit+=surrogate[i].Aaq;          
    }

  conc_org=LWC;
  for (i=0;i<n;++i)
    if (surrogate[i].is_organic or i==config.iH2O)     
      conc_org+=surrogate[i].Aaq;
      
  conc_org=max(conc_org,1.0e-5*config.MOmin);

  if (compute_activity_coefficients)
    {
      if (config.compute_organic)
        activity_coefficients_org_ssh(config, surrogate, false, Temperature, MOW, viscosity);
      activity_coefficients_aq_ssh(config, surrogate, Temperature, LWC, MMaq, XH2O, conc_org);
      if (config.compute_long_and_medium_range_interactions)
        activity_coefficients_LR_MR_ssh(config, surrogate, Temperature, LWC, ionic);
    }
  
  
  kinetic_ssh(config, surrogate, MOinit, MOW, MMaq, LWC, AQinit, ionic, chp, Temperature, RH, deltat, 1);

  for (i=0;i<n;++i)
    if(surrogate[i].is_organic) 
      {
	conc_available=0.0;
	sum_rates=0.0;	      	     
        surrogate[i].Atot=max(redmax*surrogate[i].Atot0,surrogate[i].Atot0+0.5*(surrogate[i].flux_chem_tot(0)+surrogate[i].flux_chem_tot(1))); 
        if (surrogate[i].nonvolatile)
          {
            if (surrogate[i].hydrophilic and surrogate[i].hydrophobic and LWC>config.LWClimit)
              {                                  
                double ratio=AQinit/max(MOinit,1.0e-10)*surrogate[i].gamma_org/(surrogate[i].GAMMAinf*surrogate[i].gamma_aq)*MOW/MMaq;
                surrogate[i].Ap=1.0/(1.0+ratio)*surrogate[i].Atot;
                surrogate[i].Aaq=ratio/(1.0+ratio)*surrogate[i].Atot;
              }
            else if (surrogate[i].hydrophobic or LWC<=config.LWClimit)                
              surrogate[i].Ap=surrogate[i].Atot;
            else if (surrogate[i].hydrophilic)      
              surrogate[i].Aaq=surrogate[i].Atot;            
          }
        else
          {
            double Kaq=0.0;
            double Kp=0.0;
            double sum=1.0;
            if (surrogate[i].hydrophobic or LWC<=config.LWClimit)
              {
		if (surrogate[i].kp_from_experiment)
		  Kp=surrogate[i].kpi;
		else
		  Kp=surrogate[i].kpi/MOW/surrogate[i].gamma_org;
                sum+=Kp*MOinit;
              }

            if (surrogate[i].hydrophilic and LWC>config.LWClimit) 
              {
                Kaq=surrogate[i].Kp_eff_aqreal_ssh(config, Temperature, ionic, chp,surrogate[config.iHp].gamma_LR,
                                               surrogate[config.iHp].gamma_SRMR,MMaq,fion1,fion2)
                  /surrogate[i].gamma_aq;
                sum+=Kaq*AQinit;
              }             
           
            if (surrogate[i].hydrophilic and surrogate[i].hydrophobic and LWC>config.LWClimit)                
              {
                surrogate[i].Ap=surrogate[i].Atot*Kp*MOinit/sum; 
                surrogate[i].Aaq=surrogate[i].Atot*Kaq*AQinit/sum;                 
              }                
            else if (surrogate[i].hydrophobic or LWC<=config.LWClimit)                
              surrogate[i].Ap=surrogate[i].Atot*Kp*MOinit/sum; 
            else if (surrogate[i].hydrophilic) 
              surrogate[i].Aaq=surrogate[i].Atot*Kaq*AQinit/sum;                
                             
          }        
      }
    else if (i==config.iHSO4m or i==config.iSO4mm or i==config.iSO4mm or i==config.iNO3m or i==config.iClm or i==config.iNH4p)
      {
	surrogate[i].Aaq=max(redmax*surrogate[i].Aaq0,surrogate[i].Aaq0+0.5*(surrogate[i].flux_chem_tot(0)+surrogate[i].flux_chem_tot(1)));
      }
    else if (i==config.iHNO3 or i==config.iHCl or i==config.iNH3)      
      surrogate[i].Ag=max(redmax*surrogate[i].Ag0,surrogate[i].Ag0+0.5*(surrogate[i].flux_chem_tot(0)+surrogate[i].flux_chem_tot(1)));      
}


void integer_chem_sat_ssh(model_config &config, vector<species>& surrogate,
                      Array <double, 1> &MOinit, Array <double, 1> &MOW, double &MMaq,
                      double &LWC, double &AQinit, double &ionic, double &chp,
                      double &Temperature, double &RH, double &deltat, bool &compute_activity_coefficients)
{
  //compute the dynamic evolution of the organic-phase and the aqueous-phase concentrations with
  //         the roschem algorythm
  //FOR THE COUPLED SYSTEM
  //compute_activity_coefficients: should the activity coefficients be recomputed?
  //tequilibrium: criterium to distinguish cases under equilibrium (characteristic time lower
  //                  than tequilibrium)
  int n=surrogate.size();
  int i,iphase; 
  int nphase=surrogate[0].Ap_sat.size();
  //double tiny=0.; //1.0e-26;
  double conc_available=0.0;
  double sum_rates;
  //double gamma=1.0+pow(2.,0.5)/2.; //1.7071;
  //double Kaq, Kp;      
  double redmax=0.01;  
  bool all_hydrophobic=false;
  if (LWC<=config.LWClimit)
    all_hydrophobic=true;
  
  //     gamma > .25     ->  A-stable  ( x(n+1)-x(n))*f(x(n)) > 0
  //     gamma < .25     ->  monotonic behaviour (x(n+1)-x*)(x(n)-x*) > 0
  //     gamma =  1+-sqrt(1/2) ->  L-stability
          

  MOinit=0.0;
  AQinit=LWC;
  double fion1,fion2;
  double XH2O;
  double viscosity;

  for (i=0;i<n;++i)
    {
      if (surrogate[i].hydrophobic or (all_hydrophobic and (surrogate[i].is_organic or i==config.iH2O)))        
        for (iphase=0;iphase<nphase;iphase++)
          MOinit(iphase)+=surrogate[i].Ap_sat(iphase);        
      if (surrogate[i].hydrophilic)        
        AQinit+=surrogate[i].Aaq;          
    }

  double conc_org=LWC;
  for (i=0;i<n;++i)
    if (surrogate[i].is_organic or i==config.iH2O)     
      conc_org+=surrogate[i].Aaq;
      
  conc_org=max(conc_org,1.0e-5*config.MOmin);

  if (compute_activity_coefficients)
    {
      if (config.compute_organic)
        if (LWC>config.LWClimit)
          activity_coefficients_saturation_ssh(config, surrogate, false, Temperature, MOW, viscosity);       
        else
          activity_coefficients_saturation_ssh(config, surrogate, true, Temperature, MOW, viscosity); 
      activity_coefficients_aq_ssh(config, surrogate, Temperature, LWC, MMaq, XH2O, conc_org);
      if (config.compute_long_and_medium_range_interactions)
        activity_coefficients_LR_MR_ssh(config, surrogate, Temperature, LWC, ionic);
    }  

  //save initial concentrations 
  for (i=0;i<n;++i)    
    {
      surrogate[i].Atot0=surrogate[i].Atot;    
      surrogate[i].Aaq0=surrogate[i].Aaq;
      surrogate[i].Ag0=surrogate[i].Ag;
    }

  kinetic_sat_ssh(config, surrogate, MOinit, MOW, MMaq, LWC, AQinit, ionic, chp, Temperature, RH, deltat, 0);  
  //compute the first evaluation of concentrations
  for (i=0;i<n;++i)
    if(surrogate[i].is_organic) 
      {
	conc_available=0.0;
	sum_rates=0.0;	      
        surrogate[i].Atot=max(redmax*surrogate[i].Atot,surrogate[i].Atot+surrogate[i].flux_chem_tot(0));	     
        surrogate[i].Atot1=surrogate[i].Atot;

        if (surrogate[i].nonvolatile)
          {
            double sum=0.0;
            if (surrogate[i].hydrophobic or LWC<=config.LWClimit)
              for (iphase=0;iphase<nphase;iphase++)                                   
                sum+=MOinit(iphase)/MOW(iphase)/surrogate[i].gamma_org_sat(iphase);                    
              
            if (surrogate[i].hydrophilic and LWC>config.LWClimit) 
              sum+=AQinit/(MMaq*surrogate[i].gamma_aq*surrogate[i].GAMMAinf);
              
            if (sum>0.0)
              {
                if (surrogate[i].hydrophobic or LWC<=config.LWClimit)
                  for (iphase=0;iphase<nphase;iphase++)                                   
                    surrogate[i].Ap_sat(iphase)=surrogate[i].Atot*MOinit(iphase)/MOW(iphase)/surrogate[i].gamma_org_sat(iphase)/sum;                    
                  
                if (surrogate[i].hydrophilic and LWC>config.LWClimit) 
                  surrogate[i].Aaq=surrogate[i].Atot*AQinit/(MMaq*surrogate[i].gamma_aq*surrogate[i].GAMMAinf)/sum;
              }        
          }
        else
          {
            double Kaq=0.0;
            Array <double, 1> Kp;
            Kp.resize(nphase);
            double sum=1.0;
            if (surrogate[i].hydrophobic or LWC<=config.LWClimit)
              {
		if (surrogate[i].kp_from_experiment)
		  {
		    if (iphase==surrogate[i].jmain_phase)
		      Kp(iphase)=surrogate[i].kpi;
		    else
		      Kp(iphase)=0;
		  }
		else
		  Kp(iphase)=surrogate[i].kpi/MOW(iphase)/surrogate[i].gamma_org_sat(iphase);
                sum+=Kp(iphase)*MOinit(iphase);
              }

            if (surrogate[i].hydrophilic and LWC>config.LWClimit) 
              {
                Kaq=surrogate[i].Kp_eff_aqreal_ssh(config, Temperature, ionic, chp,surrogate[config.iHp].gamma_LR,
                                               surrogate[config.iHp].gamma_SRMR,MMaq,fion1,fion2)
                  /surrogate[i].gamma_aq;
                sum+=Kaq*AQinit;
              }             

            if (surrogate[i].hydrophilic and surrogate[i].hydrophobic and LWC>config.LWClimit)                
              {
                for (iphase=0;iphase<nphase;iphase++)
                  surrogate[i].Ap_sat(iphase)=surrogate[i].Atot*Kp(iphase)*MOinit(iphase)/sum; 
                surrogate[i].Aaq=surrogate[i].Atot*Kaq*AQinit/sum; 
              }                
            else if (surrogate[i].hydrophobic or LWC<=config.LWClimit)                
              for (iphase=0;iphase<nphase;iphase++)
                surrogate[i].Ap_sat(iphase)=surrogate[i].Atot*Kp(iphase)*MOinit(iphase)/sum; 
            else if (surrogate[i].hydrophilic) 
              surrogate[i].Aaq=surrogate[i].Atot*Kaq*AQinit/sum;                                             
          }        
      }
    else if (i==config.iHSO4m or i==config.iSO4mm or i==config.iNO3m or i==config.iClm or i==config.iNH4p)
      {
	surrogate[i].Aaq=max(redmax*surrogate[i].Aaq,surrogate[i].Aaq+surrogate[i].flux_chem_tot(0));
	surrogate[i].Aaq1=surrogate[i].Aaq; 
	//cout << "computed " << surrogate[i].name << " " << surrogate[i].Aaq << " " << surrogate[i].Aaq0 << endl;
      }
    else if (i==config.iHNO3 or i==config.iHCl or i==config.iNH3)
      {
	surrogate[i].Ag=max(redmax*surrogate[i].Ag,surrogate[i].Ag+surrogate[i].flux_chem_tot(0));      
	surrogate[i].Ag1=surrogate[i].Ag;
      }
  /*
    else if (i==config.iHSO4m or i==config.iSO4mm or i==config.iNO3m or i==config.iClm or i==config.iNH4p)      
      surrogate[i].Aaq=max(redmax*surrogate[i].Aaq,surrogate[i].Aaq+surrogate[i].flux_chem_tot(0)*surrogate[i].aqratio);
    else if (i==config.iHNO3 or i==config.iHCl or i==config.iNH3)      
    surrogate[i].Ag=max(redmax*surrogate[i].Ag,surrogate[i].Ag+surrogate[i].flux_chem_tot(0)*surrogate[i].aqratio);      */
  
  
  MOinit=0.0;
  AQinit=LWC;
  
  for (i=0;i<n;++i)
    {
      if (surrogate[i].hydrophobic or (all_hydrophobic and (surrogate[i].is_organic or i==config.iH2O)))        
        for (iphase=0;iphase<nphase;iphase++)
          MOinit(iphase)+=surrogate[i].Ap_sat(iphase);        
      if (surrogate[i].hydrophilic)        
        AQinit+=surrogate[i].Aaq;          
    }

  conc_org=LWC;
  for (i=0;i<n;++i)
    if (surrogate[i].is_organic or i==config.iH2O)     
      conc_org+=surrogate[i].Aaq;
      
  conc_org=max(conc_org,1.0e-5*config.MOmin);
  
  if (compute_activity_coefficients)
    {
      if (config.compute_organic)
        activity_coefficients_saturation_ssh(config, surrogate, false, Temperature, MOW, viscosity);
      activity_coefficients_aq_ssh(config, surrogate, Temperature, LWC, MMaq, XH2O, conc_org);
      if (config.compute_long_and_medium_range_interactions)
        activity_coefficients_LR_MR_ssh(config, surrogate, Temperature, LWC, ionic);
    }
  
  
  kinetic_sat_ssh(config, surrogate, MOinit, MOW, MMaq, LWC, AQinit, ionic, chp, Temperature, RH, deltat, 1);

  for (i=0;i<n;++i)
    if(surrogate[i].is_organic) 
      {
        conc_available=0.0;
        sum_rates=0.0;	      	     
        surrogate[i].Atot=max(redmax*surrogate[i].Atot0,surrogate[i].Atot0+0.5*(surrogate[i].flux_chem_tot(0)+surrogate[i].flux_chem_tot(1))); 
        if (surrogate[i].nonvolatile)
          {
            double sum=0.0;
            if (surrogate[i].hydrophobic or LWC<=config.LWClimit)
              for (iphase=0;iphase<nphase;iphase++)                                   
                sum+=MOinit(iphase)/MOW(iphase)/surrogate[i].gamma_org_sat(iphase);                    
              
            if (surrogate[i].hydrophilic and LWC>config.LWClimit) 
              sum+=AQinit/(MMaq*surrogate[i].gamma_aq*surrogate[i].GAMMAinf);
              
            if (sum>0.0)
              {
                if (surrogate[i].hydrophobic or LWC<=config.LWClimit)
                  for (iphase=0;iphase<nphase;iphase++)                                   
                    surrogate[i].Ap_sat(iphase)=surrogate[i].Atot*MOinit(iphase)/MOW(iphase)/surrogate[i].gamma_org_sat(iphase)/sum;                    
                  
                if (surrogate[i].hydrophilic and LWC>config.LWClimit) 
                  surrogate[i].Aaq=surrogate[i].Atot*AQinit/(MMaq*surrogate[i].gamma_aq*surrogate[i].GAMMAinf)/sum;
              }        
          }
        else
          {
            double Kaq=0.0;
            Array <double, 1> Kp;
            Kp.resize(nphase);
            double sum=1.0;
            if (surrogate[i].hydrophobic or LWC<=config.LWClimit)
              {
		if (surrogate[i].kp_from_experiment)
		  {
		    if (iphase==surrogate[i].jmain_phase)
		      Kp(iphase)=surrogate[i].kpi;
		    else
		      Kp(iphase)=0;
		  }
		else
		  Kp(iphase)=surrogate[i].kpi/MOW(iphase)/surrogate[i].gamma_org_sat(iphase);
                sum+=Kp(iphase)*MOinit(iphase);
              }

            if (surrogate[i].hydrophilic and LWC>config.LWClimit) 
              {
                Kaq=surrogate[i].Kp_eff_aqreal_ssh(config, Temperature, ionic, chp,surrogate[config.iHp].gamma_LR,
                                               surrogate[config.iHp].gamma_SRMR,MMaq,fion1,fion2)
                  /surrogate[i].gamma_aq;
                sum+=Kaq*AQinit;
              }             

            if (surrogate[i].hydrophilic and surrogate[i].hydrophobic and LWC>config.LWClimit)                
              {
                for (iphase=0;iphase<nphase;iphase++)
                  surrogate[i].Ap_sat(iphase)=surrogate[i].Atot*Kp(iphase)*MOinit(iphase)/sum; 
                surrogate[i].Aaq=surrogate[i].Atot*Kaq*AQinit/sum; 
              }                
            else if (surrogate[i].hydrophobic or LWC<=config.LWClimit)                
              for (iphase=0;iphase<nphase;iphase++)
                surrogate[i].Ap_sat(iphase)=surrogate[i].Atot*Kp(iphase)*MOinit(iphase)/sum; 
            else if (surrogate[i].hydrophilic) 
              surrogate[i].Aaq=surrogate[i].Atot*Kaq*AQinit/sum;                                             
          }        
      }        
    else if (i==config.iHSO4m or i==config.iSO4mm or i==config.iSO4mm or i==config.iNO3m or i==config.iClm or i==config.iNH4p)         
      surrogate[i].Aaq=max(redmax*surrogate[i].Aaq0,surrogate[i].Aaq0+0.5*(surrogate[i].flux_chem_tot(0)+surrogate[i].flux_chem_tot(1)));               
    else if (i==config.iHNO3 or i==config.iHCl or i==config.iNH3)      
      surrogate[i].Ag=max(redmax*surrogate[i].Ag0,surrogate[i].Ag0+0.5*(surrogate[i].flux_chem_tot(0)+surrogate[i].flux_chem_tot(1)));      
}


void adapstep_chem_ssh(model_config &config, vector<species>& surrogate, double &deltat1, double &t, double &tend, double &deltatmin)
{
  int n=surrogate.size();
  int i;
  double n2err=0.0;
  int m=0;
  double EPSER=config.EPSER;
  double R=10.;
  double tinym=1.0e-5;
  for (i=0;i<n;i++)
    if (surrogate[i].is_organic)
      {        
        if (surrogate[i].Atot1 > tinym or surrogate[i].Atot > tinym)
          {
            ++m;
            n2err+=pow((surrogate[i].Atot-surrogate[i].Atot1)/(surrogate[i].Atot1+tinym),2);    
          }
      }
    else if (i==config.iHSO4m or i==config.iSO4mm or i==config.iNO3m or i==config.iClm or i==config.iNH4p)
      {
	if (surrogate[i].Atot1 > tinym or surrogate[i].Atot > tinym)
          {
            ++m;
            n2err+=pow((surrogate[i].Aaq-surrogate[i].Aaq1)/(surrogate[i].Aaq1+tinym),2);    
          }
      }
    else if (i==config.iHNO3 or i==config.iHCl or i==config.iNH3)
      {
	if (surrogate[i].Ag1 > tinym or surrogate[i].Ag > tinym)
          {
            ++m;
            n2err+=pow((surrogate[i].Ag-surrogate[i].Ag1)/(surrogate[i].Ag1+tinym),2);    
          }
      }


  n2err=min(n2err,EPSER*R*R);
  n2err=max(n2err,EPSER/(R*R));
  deltat1=deltat1*pow(EPSER/n2err,0.5);
  deltat1=max(deltatmin,deltat1);
  deltat1=min(tend-t,deltat1);
  //deltat1=min(deltat1,0.1);
}

/*
void solve_chemistry_ssh(model_config &config, vector<species>& surrogate,
                     double &MOinit,double &MOW,
                     double &LWC, double &AQinit, double &ionic, double &chp,
                     double &Temperature, double &RH, double deltat, bool &compute_activity_coefficients)
{
  double t=0;
  double dt2=config.dtchem_min;
  double dt1=dt2;
  double MMaq;
  int nphase=1;
  if (config.compute_saturation)        
    nphase=surrogate[0].Ap_sat.size();                  

  if (config.compute_saturation and nphase>1)
    {
      Array <double, 1> MOinit_sat,MOW_sat;
      MOinit_sat.resize(nphase);
      MOW_sat.resize(nphase);      
      while (t<deltat)
        {      
          dt1=min(deltat-t,dt1);	 
          dt2=dt1;            
          integer_chem_sat_ssh(config, surrogate, MOinit_sat, MOW_sat, MMaq, LWC, AQinit, ionic, chp, Temperature, RH, dt1, compute_activity_coefficients);
	  
          //compute the new time step so that changes are small
          adapstep_chem_ssh(config,surrogate,dt1,t,deltat,config.dtchem_min);      
          t+=dt2;     
        }
    }
  else
    {
      while (t<deltat)
        {      
          dt1=min(deltat-t,dt1);	 
          dt2=dt1;            
          integer_chem_ssh(config, surrogate, MOinit, MOW, MMaq, LWC, AQinit, ionic, chp, Temperature, RH, dt1, compute_activity_coefficients);

          //compute the new time step so that changes are small
          adapstep_chem_ssh(config,surrogate,dt1,t,deltat,config.dtchem_min);
	  
          t+=dt2;     
	}
    }
    }*/
