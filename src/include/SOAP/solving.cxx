//!!-----------------------------------------------------------------------
//!!     Copyright (C) 2019-2024 CEREA (ENPC) - INERIS
//!!     SSH-aerosol is distributed under the GNU General Public License v3
//!!-----------------------------------------------------------------------

#include "properties.cxx"
#include "equilibrium.cxx"
#include "dynamic.cxx"
#include "chemistry.cxx"
using namespace ssh_soap;

void solve_equilibrium_ssh(model_config &config, vector<species>& surrogate,
                           double &MOinit,double &MOW,
                           double &LWC, double &AQinit, double &ionic, double &chp,
                           double &Temperature, double &RH)
{  
  double organion=0.0;
  double conc_inorganic=0.0;
  double ionic_organic=0.0;
  double MMaq=18.0;
  int i; //,it;
  int n=surrogate.size();
  for (i=0;i<n;++i)
    if (surrogate[i].is_organic==false and i!=config.iH2O)
      conc_inorganic+=surrogate[i].Aaq;
  bool compute_activity_coefficients=true;
    
  //if (config.compute_inorganic)
  //  config.coupled_phases=true;
  
  if (config.coupled_phases or 
      (RH>=config.RHcoupling and surrogate[config.iH2O].hydrophilic and surrogate[config.iH2O].hydrophobic))
    {      
      // If the syst is coupled particulate hydrophilic compounds concentrations
      //must be computed simultaneously with the particulate hydrophobic compounds concentrations
      //cout << "The system is coupled." << endl;
      if (LWC>config.LWClimit)
        {      
          bool all_hydrophobic=false;
          double error1=1000.0;
          double error2=1000.0;
          double deriv_error1_MO,deriv_error1_AQ,deriv_error2_MO,deriv_error2_AQ;
          //double derivative=0.0;
          int index_iter=0;
          double MO=MOinit;
          double AQ=AQinit;
          // minimize the functions errors with a method of newton raphson
          // error1 = MOinit - sum of concentrations of organic compounds in the organic phase
          //           - hygroscopicity of the organic phase  
          // error2 = AQinit - sum of concentrations of organic compounds in the aqueous phase
          //           - hygroscopicity of organic compounds
          //           - LWC - sum of inorganic ions

          if (config.compute_inorganic==true)
            {
              int nh=config.nh_inorg_init;
              //Array<double, 1> vec_error_aq,vec_error_chp,vec_error_compaq;
	      Array<double, 1> error_spec,error_spec_old;
              //vec_error_aq.resize(config.max_iter);
	      //vec_error_chp.resize(config.max_iter);
	      //vec_error_compaq.resize(config.max_iter);
              error_spec.resize(n);
              error_spec_old.resize(n);
              //bool non_convergence;
              //int iiter=0;
	      double chp_old=chp;
              double error3=1000.0;
	      double error4=1000.0;
              double var_rej=0.2;
              double var_max=0.12;
              double var_min=0.08;
              double factor_min=0.01;
              double factor_max=1; //1.0/pow(2.,nh);
              if (config.solids)
                {
                  if (RH<0.8)
                    factor_max=0.5;
                }
              else
                if (RH<0.2)
                  factor_max=0.5;

              double factor=1.0/pow(2.0,nh);
              //double factor=1.0/nh;
              double relprec=1.0e-3;
              //int maxiiter=30;
              //cout << nh << endl;
              double error1_old=0.0;
              double error2_old=0.0;
              double error3_old=0.0;
              double error4_old=0.0;
              double factor_old=factor;
              double factor_old2=factor;
              double chp_old2=chp_old;
              double ionic_old;
              for (i=0;i<n;i++)
                if (surrogate[i].hydrophilic)
                  {                              
                    surrogate[i].gamma_aq=1.0;
                  }
              error_spec=1000.;

	      //double error4_min=1000;
	      //int nerr4=0;
              int ntoo_big=0;             
	      //int ntoo_low=0;
              //cout << "ici" << endl;
              //if (RH<0.4 or RH>0.9)
              //config.first_evaluation_activity_coefficients=true;                            
              //config.compute_organic=false; 

              double RHsave=RH;
              //RH=max(RH,0.2);
              double error_hyd=0.0;
	      
              while ((index_iter < config.max_iter) and (abs(error1)/factor_old > config.precision or abs(error2)/factor_old > config.precision or
                                                         ((abs(error3)/factor_old>relprec  or abs(error4)/factor_old>relprec) and AQ>100*config.MOmin) or
                                                         config.first_evaluation_activity_coefficients==true or RH>RHsave))
                {
		  /*
                    if (config.first_evaluation_activity_coefficients==true)
                    {
                    if (abs(error2)/factor_old <= config.precision and abs(error3)/factor_old<=relprec and abs(error4)/factor_old<=relprec and index_iter>0)
                    {
                    config.first_evaluation_activity_coefficients=false;                          
                    RH=RHsave;
                    }
                    }  */          
                    
                  if (config.first_evaluation_activity_coefficients==false or index_iter==0)
                    compute_activity_coefficients=true;
                  else
                    compute_activity_coefficients=false;

                  double Asol=0.;
		  for (i=0;i<n;i++)
                    {
		      surrogate[i].Ag_old=surrogate[i].Ag;
                      if (surrogate[i].hydrophilic and surrogate[i].is_inorganic_precursor==false)
                        {                          
                          surrogate[i].Aaq_old=surrogate[i].Aaq;
                          surrogate[i].gamma_aq_old=surrogate[i].gamma_aq;
                        }
                      if (surrogate[i].hydrophobic)
                        {
                          surrogate[i].Ap_old=surrogate[i].Ap;
                          surrogate[i].gamma_org_old=surrogate[i].gamma_org;
                        }
                      else if (surrogate[i].is_solid)
                        {
                          Asol+=surrogate[i].Ap;
                          surrogate[i].Ap_old=surrogate[i].Ap;
                        }
                    }

                  if (AQ<=config.MOmin)
                    error3_old=0.;
                  else
                    error3_old=(chp-chp_old)/max(chp_old,1.0e-20)/factor_old;                                
                    
                  error4_old=error4/factor_old;
                  chp_old2=chp_old;
                  chp_old=chp;		  
                  double AQsave=max(AQ,config.MOmin);
                  double MOsave=max(MO,config.MOmin);
                  double negligeable=max(config.MOmin,0.01*(AQsave+Asol));
		  double negligeable_org=max(config.MOmin,0.01*MO);
                  
		  error1_old=error1/factor_old;
                  error2_old=error2/factor_old;                  
                  factor_old2=factor_old;
                  factor_old=factor;
                  error_spec_old=error_spec;
                  ionic_old=ionic;
                  
                  error_coupled_inorg_ssh(config,surrogate,MO,MOW,MMaq,AQ,LWC,conc_inorganic,ionic,ionic_organic,
                                          chp, organion,
                                          Temperature,RH,
                                          error1,deriv_error1_MO,deriv_error1_AQ,
                                          error2,deriv_error2_MO,deriv_error2_AQ,factor, compute_activity_coefficients);
		  hydratation(config,surrogate,RH,error_hyd,false,false);
		  
                                            
                  double var=0.;                    
		  error4=error_hyd;
                  error_spec=0.0;
		  for (i=0;i<n;i++)
		    {
		      if (surrogate[i].hydrophilic  and surrogate[i].is_inorganic_precursor==false)
			{
			  if (surrogate[i].Aaq_old>negligeable and surrogate[i].Aaq>negligeable)                                              
			    {
			      error_spec(i)=(surrogate[i].Aaq-surrogate[i].Aaq_old)/surrogate[i].Aaq_old;                          
			      error4=max(error4,abs(error_spec(i)));
			      //cout << surrogate[i].name << endl;
			    }
			}
		      else if (surrogate[i].is_solid and (surrogate[i].Ap>negligeable and surrogate[i].Ap_old>negligeable))
			{
			  error_spec(i)=(surrogate[i].Ap-surrogate[i].Ap_old)/surrogate[i].Ap_old;                          
			  error4=max(error4,abs(error_spec(i)));
			  //cout << surrogate[i].name << " " << surrogate[i].Ap << " " << surrogate[i].Ap_old << endl;
			}

		      if (surrogate[i].hydrophobic)
			{
			  if (surrogate[i].Ap_old>negligeable_org and surrogate[i].Ap>negligeable_org)                                              
			    {
			      error_spec(i)=max(error_spec(i),(surrogate[i].Ap-surrogate[i].Ap_old)/surrogate[i].Ap_old);                          
			      error4=max(error4,abs(error_spec(i)));
			      //cout << surrogate[i].name << endl;
			    }
			}
		    }

                  int m=0;
                  for (i=0;i<n;i++)
		    if (surrogate[i].hydrophilic and surrogate[i].is_inorganic_precursor==false)
                      //if (error_spec(i)/factor_old>relprec)
                      {
                        if (0.5*surrogate[i].gamma_aq+0.5*surrogate[i].gamma_aq_old>1.0e-6 and (surrogate[i].Aaq>negligeable or surrogate[i].Aaq_old>negligeable))
                          {
                            m++;
                            var+=abs((surrogate[i].gamma_aq_old-surrogate[i].gamma_aq)/(0.5*surrogate[i].gamma_aq+0.5*surrogate[i].gamma_aq_old));
                          }
                      }
                  
                  for (i=0;i<n;i++)
		    if (surrogate[i].hydrophobic)
                      //if (error_spec(i)/factor_old>relprec)
                      {
                        if (0.5*surrogate[i].gamma_org+0.5*surrogate[i].gamma_org_old>1.0e-6 and (surrogate[i].Ap>negligeable_org or surrogate[i].Ap_old>negligeable_org))
                          {
                            m++;
                            var+=abs((surrogate[i].gamma_org_old-surrogate[i].gamma_org)/(0.5*surrogate[i].gamma_org+0.5*surrogate[i].gamma_org_old));
                          }
                      }

                  if (m>0)
                    var=var/m; 

                  if (AQ<=config.MOmin)
                    error3=0.;
                  else
                    error3=(chp-chp_old)/max(chp_old,1.0e-15);
                  
                  
                  if (var<=var_rej or factor==factor_min)
                    {
                      AQ=AQ-error2;
                      MO=MO-error1;
                      AQ=max(config.MOmin,max(0.1*AQsave,min(AQ,10.0*AQsave)));
                      MO=max(config.MOmin,max(0.1*MOsave,min(MO,10.0*MOsave)));
                    }                                     

                  if (var>var_max and index_iter>0 and factor>factor_min)
                    {                      
                      factor=max(factor_old*max(var_max/var,0.1),factor_min);            
                      if (var>var_rej)
                        for (i=0;i<n;i++)
                          {
			    surrogate[i].Ag=surrogate[i].Ag_old;
			    if (surrogate[i].hydrophilic and surrogate[i].is_inorganic_precursor==false)
                              surrogate[i].Aaq=surrogate[i].Aaq_old;                        
                            if (surrogate[i].hydrophobic or surrogate[i].is_solid)
			      surrogate[i].Ap=surrogate[i].Ap_old; 
                          }
                    }
                  else if (var<var_min and index_iter>0 and factor<1.)
                    {
                      if (var>0.)
                        factor=min(factor_old*min(var_max/var,10.),factor_max);
                      else
                        factor=min(factor_old*10.,1.0);
                    }                    

                  if (index_iter>1000)
                    factor_max=min(factor_max,0.5);

                  if (index_iter>0)
                    {             
                      if (error2*error2_old<0.0 and abs(error2)/factor_old>config.precision and abs(error1)/factor_old<config.precision)                        
                        ntoo_big=ntoo_big+10;
		      else if (error1*error1_old<0.0 and abs(error1)/factor_old>config.precision and abs(error2)/factor_old<config.precision)                        
                        ntoo_big=ntoo_big+10;
		      else if (error2*error2_old<0.0 and abs(error2)/factor_old>config.precision and error1*error1_old<0.0 and abs(error1)/factor_old>config.precision)
			ntoo_big=ntoo_big+10;
		      
                      int a=0;
                      //cout << error3/factor_old << " " << error3_old << endl;
                      if (error3*error3_old<0.0 and abs(error3)/factor_old>relprec)
                        {
                          //if (max(abs(error3),abs(error3_old))/(abs(error3)+abs(error3_old))<0.6)
                            
                          factor=min(factor,factor_old*(max(abs(error3)/factor_old,abs(error3_old))/(abs(error3)/factor_old+abs(error3_old)))); 
                          ntoo_big=ntoo_big+1;
                        }                     
                      
                      if (error4/factor_old>relprec)
                        for (i=0;i<n;i++)
                          if (error_spec(i)*error_spec_old(i)<0.0 and abs(error_spec(i))/factor_old>relprec)     
                            {                       
                              //factor=min(factor,0.5*factor_old); //
                              //factor=min(factor,factor_old*(max(abs(error_spec(i)),abs(error_spec_old(i)))/(abs(error_spec(i))+abs(error_spec_old(i)))));                              
                              //ntoo_big++;
                              a=a+1;
                            }

                      //factpr=min(factor,factor_old*(ionic

                      if (error4>var_rej)
                        {
                          //factor=min(factor,factor_old*var_rej/error4);
                          //a=a+1;
                          ntoo_big++; //(error4/factor_old/var_rej);
                        }               
                      if (a>0)
                        ntoo_big++;
                      factor=max(factor,factor_min);
                    }                  
                                
                  if (ntoo_big>200 and factor_max>factor_min)
                    {
                      factor_max=max(factor_max/2,factor_min);
                      
                      ntoo_big=0;                  
                    }          
                  
                  factor=min(factor,factor_max);                   
                  ++index_iter;
                  
                }	
	   
              if (index_iter>=config.max_iter)
		{

		  config.iiter=0;
		  cout << "inorg tot " << index_iter << " " << error1/factor_old << " " << error2/factor_old << " " << error3/factor_old  << " " << error4/factor_old << " " << factor_old << " " << RH << " " << Temperature << " " << MO << " " << AQ << " " << surrogate[config.iH2O].Ap << " " << surrogate[config.iH2O].Aaq << " " << factor_max << " " << factor << endl;
		  //throw string("Stop.");
		}
              else
                config.iiter=1;
            }
          else
	    {
	      int max_iter2=min(config.max_iter,50); //config.max_iter);
	      //index_iter=0;
	      double factor=1;
	      double error1_save=1000;
	      double error1_save2=1000;
	      double error2_save=1000;
	      double error2_save2=1000;
	      double error_hyd=0.0;
	      while ((index_iter < max_iter2) and ((abs(error1)/factor > config.precision)
						   or (abs(error2)/factor > config.precision) or error_hyd>1.0e-3))
		{
		  if (config.first_evaluation_activity_coefficients==false or index_iter==0)
		    compute_activity_coefficients=true;
		  else
		    compute_activity_coefficients=false;
		  
		  error_coupled_ssh(config,surrogate,MO,MOW,MMaq,AQ,LWC,conc_inorganic,ionic,ionic_organic,
				    chp, organion,
				    Temperature,RH,
				    error1,deriv_error1_MO,deriv_error1_AQ,
				    error2,deriv_error2_MO,deriv_error2_AQ,factor, compute_activity_coefficients);
		  //solve the system with a method of newton raphson
		  if ((abs(error1_save2-error1)<0.01*abs(error1) and error1_save*error1<0. and error2/factor<config.precision) or
		      (abs(error2_save2-error2)<0.01*abs(error2) and error2_save*error2<0. and error1/factor<config.precision) or
		      (abs(error1_save2-error1)<0.01*abs(error1) and error1_save*error1<0. and abs(error2_save2-error2)<0.01*abs(error2) and error2_save*error2<0))
		    {
		      error1_save=1000;
		      error1_save2=1000;
		      error2_save=1000;
		      error2_save2=1000;
		      factor=max(factor/2,0.01);
		      //cout << "non convergence" << endl;
		    }
		  
		  newton_raphson_coupled_ssh(MO,AQ,error1,deriv_error1_MO,deriv_error1_AQ,
					     error2,deriv_error2_MO,deriv_error2_AQ);

		  hydratation(config,surrogate,RH,error_hyd,false,false);
		  
		  error1_save2=error1_save;
		  error1_save=error1;
		  error2_save2=error2_save;
		  error2_save=error2;
		  ++index_iter;
		  //cout << index_iter << " " << MO << " " << AQ << " " << error1 << " " << error2 << " " << surrogate[config.iH2O].Ap << " " << surrogate[config.iH2O].Aaq << endl;
		}

	      // If newton raphson did not solve the system, changed of resolution
	      if ((index_iter < config.max_iter) and ((abs(error1) > config.precision)
						      or (abs(error2) > config.precision)))
		{
		  int nh=config.nh_inorg_init;              
		  //Array<double, 1> vec_error_aq,vec_error_compaq;
		  Array<double, 1> error_spec,error_spec_old;
		  //vec_error_aq.resize(config.max_iter);
		  //vec_error_compaq.resize(config.max_iter);
		  error_spec.resize(n);
		  error_spec_old.resize(n);
		  //bool non_convergence;
		  //int iiter=0;
		  double error4=1000.0;
		  double var_rej=0.2;
		  double var_max=0.12;
		  double var_min=0.08;
		  double factor_min=0.01;
		  double factor_max=1; //1.0/pow(2.,nh);
		  if (config.solids)
		    {
		      if (RH<0.8)
			factor_max=0.5;
		    }
		  else
		    if (RH<0.2)
		      factor_max=0.5;

		  if (pow(2.0, nh) == 0.0)
		    throw string("Error: division by zero in solving_equilibrium.");
		  double factor=1.0/pow(2.0,nh);
		  //double factor=1.0/nh;
		  double relprec=1.0e-3;
		  //int maxiiter=30;
		  //cout << nh << endl;
		  double error1_old=0.0;
		  double error2_old=0.0;
		  double error4_old=0.0;
		  double factor_old=factor;
		  double factor_old2=factor;
		  double ionic_old;
		  for (i=0;i<n;i++)
		    if (surrogate[i].hydrophilic)
		      {                              
			surrogate[i].gamma_aq=1.0;
		      }
		  error_spec=1000.;

		  int ntoo_big=0;             

		  double RHsave=RH;
		  //RH=max(RH,0.2);              

		  while ((index_iter < config.max_iter) and (abs(error1)/factor_old>config.precision or abs(error2)/factor_old > config.precision or abs(error4)/factor_old>relprec or config.first_evaluation_activity_coefficients==true or RH>RHsave))
		    {                                            
		      if (config.first_evaluation_activity_coefficients==false)
			compute_activity_coefficients=true;
		      else
			compute_activity_coefficients=false;

		      double Asol=0.;
		      for (i=0;i<n;i++)
			{
			  if (surrogate[i].hydrophilic)
			    {
			      surrogate[i].Ag_old=surrogate[i].Ag;
			      surrogate[i].Aaq_old=surrogate[i].Aaq;
			      surrogate[i].gamma_aq_old=surrogate[i].gamma_aq;
			    }
			  if (surrogate[i].hydrophobic)
			    {
			      surrogate[i].Ag_old=surrogate[i].Ag;
			      surrogate[i].Ap_old=surrogate[i].Ap;
			      surrogate[i].gamma_org_old=surrogate[i].gamma_org;
			    }
			}
                    
		      error4_old=error4/factor_old;		  
		      double AQsave=max(AQ,config.MOmin);
		      double negligeable=max(config.MOmin,0.01*(AQsave+Asol));
		      double negligeable_org=max(config.MOmin,0.01*MO);
                  
		      error1_old=error1/factor_old;
		      error2_old=error2/factor_old;                  
		      factor_old2=factor_old;
		      factor_old=factor;
		      error_spec_old=error_spec;
		      ionic_old=ionic;

		      error_coupled_ssh(config,surrogate,MO,MOW,MMaq,AQ,LWC,conc_inorganic,ionic,ionic_organic,
					chp, organion,
					Temperature,RH,
					error1,deriv_error1_MO,deriv_error1_AQ,
					error2,deriv_error2_MO,deriv_error2_AQ,factor, compute_activity_coefficients);
		      hydratation(config,surrogate,RH,error_hyd,false,false);
                                                
		      double var=0.;                    
		      int m=0;
		      for (i=0;i<n;i++)
			if (surrogate[i].hydrophilic)
			  {
			    if (0.5*surrogate[i].gamma_aq+0.5*surrogate[i].gamma_aq_old>1.0e-6 and (surrogate[i].Aaq>negligeable or surrogate[i].Aaq_old>negligeable))
			      {
				m++;
				var+=abs((surrogate[i].gamma_aq_old-surrogate[i].gamma_aq)/(0.5*surrogate[i].gamma_aq+0.5*surrogate[i].gamma_aq_old));
			      }
			  }
		      for (i=0;i<n;i++)
			if (surrogate[i].hydrophobic)
			  {
			    if (0.5*surrogate[i].gamma_org+0.5*surrogate[i].gamma_org_old>1.0e-6 and (surrogate[i].Ap>negligeable_org or surrogate[i].Ap_old>negligeable_org))
			      {
				m++;
				var+=abs((surrogate[i].gamma_org_old-surrogate[i].gamma_org)/(0.5*surrogate[i].gamma_org+0.5*surrogate[i].gamma_org_old));
			      }
			  }

		      error4=error_hyd;
		      error_spec=0.0;
		      for (i=0;i<n;i++)
			{
			  if (surrogate[i].hydrophilic)
			    {
			      if (surrogate[i].Aaq_old>negligeable and surrogate[i].Aaq>negligeable)                                              
				{
				  error_spec(i)=(surrogate[i].Aaq-surrogate[i].Aaq_old)/surrogate[i].Aaq_old;                          
				  error4=max(error4,abs(error_spec(i)));
				  //cout << surrogate[i].name << endl;
				}
			    }
			  else if (surrogate[i].is_solid and (surrogate[i].Ap>negligeable and surrogate[i].Ap_old>negligeable))
			    {
			      error_spec(i)=(surrogate[i].Ap-surrogate[i].Ap_old)/surrogate[i].Ap_old;                          
			      error4=max(error4,abs(error_spec(i)));
			      //cout << surrogate[i].name << " " << surrogate[i].Ap << " " << surrogate[i].Ap_old << endl;
			    }
			  
			  if (surrogate[i].hydrophobic)
			    {
			      if (surrogate[i].Ap_old>negligeable_org and surrogate[i].Ap>negligeable_org)                                              
				{
				  error_spec(i)=max(error_spec(i),(surrogate[i].Ap-surrogate[i].Ap_old)/surrogate[i].Ap_old);                          
				  error4=max(error4,abs(error_spec(i)));
				  //cout << surrogate[i].name << endl;
				}
			    }
			}

		      //If takes too much iterations make sure that factor_max is not above 0.5
		      if (index_iter>1000)
			factor_max=min(factor_max,0.5);

		      if (m>0)
			var=var/m; 
                  
		      if (var<=var_rej or factor==factor_min)
			{
			  AQ=AQ-error2;
			  MO=MO-error1;
			  AQ=max(config.MOmin,max(0.1*AQsave,min(AQ,10.0*AQsave)));
			}                                     

		      if (var>var_max and index_iter>0 and factor>factor_min)
			{                      
			  factor=max(factor_old*max(var_max/var,0.1),factor_min);            
			  if (var>var_rej)
			    for (i=0;i<n;i++)
			      {
				surrogate[i].Ag=surrogate[i].Ag_old;
				if (surrogate[i].hydrophobic)
				  {                              
				    surrogate[i].Ap=surrogate[i].Ap_old;                                                        
				  }
                            
				if (surrogate[i].hydrophilic)
				  {                              
				    surrogate[i].Aaq=surrogate[i].Aaq_old;                                                        
				  }
				else if (surrogate[i].is_solid)
				  surrogate[i].Ap=surrogate[i].Ap_old;
			      }
			}
		      else if (var<var_min and index_iter>0 and factor<1.)
			{
			  if (var>0.)
			    factor=min(factor_old*min(var_max/var,10.),factor_max);
			  else
			    factor=min(factor_old*10.,1.0);
			}                    
                    

		      if (index_iter>0)
			{                                         
			  if (error2*error2_old<0.0 and abs(error2)/factor_old>config.precision)                        
			    ntoo_big=ntoo_big+10;

			  if (error1*error1_old<0.0 and abs(error1)/factor_old>config.precision)                        
			    ntoo_big=ntoo_big+10;

			  int a=0;               
			  if (error4/factor_old>relprec)
			    for (i=0;i<n;i++)
			      if (error_spec(i)*error_spec_old(i)<0.0 and abs(error_spec(i))/factor_old>relprec)     
				{                       
				  a=a+1;
				}

			  if (error4>var_rej)
			    {
			      ntoo_big++; //(error4/factor_old/var_rej);
			    }               
			  if (a>0)
			    ntoo_big++;
			  factor=max(factor,factor_min);
			}                  
                                
		      if (ntoo_big>200 and factor_max>factor_min)
			{
			  factor_max=max(factor_max/2,factor_min);
			  ntoo_big=0;                  
			}                            
                  
		      factor=min(factor,factor_max);
                  
		      //vec_error_aq(index_iter)=abs(error2);                      
		      //vec_error_compaq(index_iter)=abs(error4);
		      ++index_iter;
		    }
		

		  if (index_iter==config.max_iter)
		    {
		      cout << "coupled tot2 " << config.max_iter << " " << error1/factor_old << " " << error2/factor_old << " " << abs(error4)/factor_old << " " << factor_old << " " << RH << " " << Temperature << " " << MO << " " << AQ << " " << surrogate[config.iH2O].Ap << " " << surrogate[config.iH2O].Aaq << endl;
		      //throw string("Stop.");
		    }
		}
	    }                         

          if (config.compute_saturation and MO > 0.0 and config.compute_organic)
            saturation_ssh(config,surrogate,all_hydrophobic,LWC,ionic,conc_inorganic,ionic_organic,
                           chp,organion,Temperature,RH);
	}
      else
        {
          //if no water assumes that all compounds are hydrophobic
          bool all_hydrophobic=true;
          double error3=1000.0;
          double derivative=0.0;
          int index_iter=0;
          double MO=MOinit;
          // minimize the function error with a method of newton raphson
          // error3 = MOinit - sum of concentrations of organic compounds in the organic phase
          //           - hygroscopicity of the organic phase
	  double error_hyd=0.;
          while ((index_iter < config.max_iter) and (abs(error3) > config.precision or error_hyd>1.0e-3))
            {
              if (config.first_evaluation_activity_coefficients==false or index_iter==0)
                compute_activity_coefficients=true;
              else
                compute_activity_coefficients=false;
              
              error_org_ssh(config,surrogate,MO,MOW,Temperature,error3,derivative,RH,1.0,
                            all_hydrophobic, compute_activity_coefficients);
	      hydratation(config,surrogate,RH,error_hyd,true,false);
	      
              //solve the system with a method of newton raphson
              if (derivative != 0.0 and MO-error3/derivative >= 0.0)
                MO=MO-error3/derivative;
              else
                MO=MO-error3;

              ++index_iter;
            }         

	  if (index_iter==config.max_iter)
	    cout << "all_hygroscopic " << config.max_iter << endl;
       		  
          for (i=0;i<n;++i)
	    if (surrogate[i].is_organic)
	      surrogate[i].Aaq=0.0;

          if (config.compute_saturation and MO > 0.0 and config.compute_organic)
            saturation_ssh(config,surrogate,all_hydrophobic,LWC,ionic,conc_inorganic,ionic_organic,
                           chp,organion,Temperature,RH);
        }
    }
  else
    {
      // If the syst is not coupled particulate hydrophilic compounds concentrations
      //can be computed independently of the particulate hydrophobic compounds concentrations
      //cout << "The system is not coupled." << endl;
      bool all_hydrophobic;
      double MO=MOinit;
      if (LWC>config.LWClimit)
        {
          all_hydrophobic=false;
          double error=1000.0;
          double derivative=0.0;
          int index_iter=0;
	  //int max_iter2=50;
		  
          // minimize the function error with a method of newton raphson
          // error = MOinit - sum of concentrations of organic compounds in the organic phase
          //           - hygroscopicity of the organic phase
	  int max_iter2=min(config.max_iter,50);
	  double error_save2=1000.;
	  double error_save=1000.;
	  error=1000;
	  double factor=1.;
	  double min_error=error_save;
	  int count_iter=0;
	  //int index_iter=0;
	  double error_hyd=0;
	  while ((index_iter < max_iter2) and (abs(error)/factor > config.precision or error_hyd>1.0e-3))
	    {
	      if (config.first_evaluation_activity_coefficients==false or index_iter==0)
		compute_activity_coefficients=true;
	      else
		compute_activity_coefficients=false;
	      
	      error_org_ssh(config,surrogate,MO,MOW,Temperature,error,derivative,RH,factor,
			    all_hydrophobic, compute_activity_coefficients);
	      hydratation(config,surrogate,RH,error_hyd,false,false);
	      //solve the system with a method of newton raphson
	      
	      if (abs(error_save2-error)<0.01*abs(error) and error_save*error<0.)
		{
		  error_save=1000;
		  error_save2=1000;
		  factor=max(factor/2,0.01);
		  count_iter=-1;
		  //cout << "non convergence" << endl;
		}
	      min_error=min(abs(error),min_error);
	      if (abs(error)>abs(error_save) and abs(error)>10*min_error and count_iter>=100)
		{
		  surrogate[config.iH2O].Ap=0.;
		  factor=max(factor/2,0.01);
		  error_save=1000;
		  error_save2=1000;
		  min_error=1000;
		  count_iter=-1;
		}
	      
	      if (derivative != 0.0 and MO-error/derivative >= 0.0)
		MO=MO-error/derivative;
	      else
		MO=MO-error;
	      
	      error_save2=error_save;
	      error_save=error;
	      
	      //cout << MO << " " << error << " " << surrogate[config.iH2O].Ap << " " << factor << endl;
	      
	      ++index_iter;
	      ++count_iter;
	    }

	  if (index_iter < config.max_iter and abs(error) > config.precision)
	    {
	      int nh=config.nh_inorg_init;              
	      //Array<double, 1> vec_error_aq,vec_error_compaq;
	      Array<double, 1> error_spec,error_spec_old;
	      //vec_error_aq.resize(config.max_iter);
	      //vec_error_compaq.resize(config.max_iter);
	      error_spec.resize(n);
	      error_spec_old.resize(n);
	      //bool non_convergence;
	      //int iiter=0;
	      double error4=1000.0;
	      double var_rej=0.2;
	      double var_max=0.12;
	      double var_min=0.08;
	      double factor_min=0.01;
	      double factor_max=1; //1.0/pow(2.,nh);
	      double factor=1.0/pow(2.0,nh);
	      //double factor=1.0/nh;
	      double relprec=1.0e-3;
	      //int maxiiter=30;
	      //cout << nh << endl;
	      double error1_old=0.0;
	      double error4_old=0.0;
	      double factor_old=factor;
	      double factor_old2=factor;
	      error_spec=1000.;

	      int ntoo_big=0;             

	      double RHsave=RH;
	      //RH=max(RH,0.2);              

	      while ((index_iter < config.max_iter) and (abs(error)/factor_old>config.precision or abs(error4)/factor_old>relprec))
		{                                            
		  if (config.first_evaluation_activity_coefficients==false)
		    compute_activity_coefficients=true;
		  else
		    compute_activity_coefficients=false;

		  double Asol=0.;
		  for (i=0;i<n;i++)
		    {
		      if (surrogate[i].hydrophobic)
			{
			  surrogate[i].Ag_old=surrogate[i].Ag;
			  surrogate[i].Ap_old=surrogate[i].Ap;
			  surrogate[i].gamma_org_old=surrogate[i].gamma_org;
			}
		    }
                    
		  error4_old=error4/factor_old;		  
		  double negligeable_org=max(config.MOmin,0.01*MO);
                  
		  error1_old=error/factor_old;                
		  factor_old2=factor_old;
		  factor_old=factor;
		  error_spec_old=error_spec;

		  error_org_ssh(config,surrogate,MO,MOW,Temperature,error,derivative,RH,factor,
				all_hydrophobic, compute_activity_coefficients);
		  hydratation(config,surrogate,RH,error_hyd,false,false);
                                                
		  double var=0.;                    
		  int m=0;
                  for (i=0;i<n;i++)
		    if (surrogate[i].hydrophobic)
                      {
                        if (0.5*surrogate[i].gamma_org+0.5*surrogate[i].gamma_org_old>1.0e-6 and (surrogate[i].Ap>negligeable_org or surrogate[i].Ap_old>negligeable_org))
                          {
                            m++;
                            var+=abs((surrogate[i].gamma_org_old-surrogate[i].gamma_org)/(0.5*surrogate[i].gamma_org+0.5*surrogate[i].gamma_org_old));
                          }
                      }

		  error4=error_hyd;
		  error_spec=0.0;
		  for (i=0;i<n;i++)
		    if (surrogate[i].hydrophobic)
		      {
			if (surrogate[i].Ap_old>negligeable_org and surrogate[i].Ap>negligeable_org)                                              
			  {
			    error_spec(i)=max(error_spec(i),(surrogate[i].Ap-surrogate[i].Ap_old)/surrogate[i].Ap_old);                          
			    error4=max(error4,abs(error_spec(i)));
			    //cout << surrogate[i].name << endl;
			  }
		      }
		  
		  if (m>0)
		    var=var/m; 
                  
		  if (var<=var_rej or factor==factor_min)
		    MO=MO-error;                                   

		  if (var>var_max and index_iter>0 and factor>factor_min)
		    {                      
		      factor=max(factor_old*max(var_max/var,0.1),factor_min);            
		      if (var>var_rej)
			for (i=0;i<n;i++)
			  {
			    surrogate[i].Ag=surrogate[i].Ag_old;
			    if (surrogate[i].hydrophobic)                           
			      surrogate[i].Ap=surrogate[i].Ap_old;    
			  }                                                    
			      
		    }
		  else if (var<var_min and index_iter>0 and factor<1.)
		    {
		      if (var>0.)
			factor=min(factor_old*min(var_max/var,10.),factor_max);
		      else
			factor=min(factor_old*10.,1.0);
		    }                    
                    

		  if (index_iter>0)
		    {                                     
		      if (error*error1_old<0.0 and abs(error)/factor_old>config.precision)                        
			ntoo_big=ntoo_big+10;

		      int a=0;               
		      if (error4/factor_old>relprec)
			for (i=0;i<n;i++)
			  if (error_spec(i)*error_spec_old(i)<0.0 and abs(error_spec(i))/factor_old>relprec)                            
			    a=a+1;
			    
		      if (error4>var_rej)
			{
			  ntoo_big++; //(error4/factor_old/var_rej);
			}               
		      if (a>0)
			ntoo_big++;
		      factor=max(factor,factor_min);
		    }                  
                                
		  if (ntoo_big>200 and factor_max>factor_min)
		    {
		      factor_max=max(factor_max/2,factor_min);
		      ntoo_big=0;                  
		    }                            
                  
		  factor=min(factor,factor_max);
                  
		  //vec_error_aq(index_iter)=abs(error2);                      
		  //vec_error_compaq(index_iter)=abs(error4);
		  ++index_iter;
		}
		

	      if (index_iter==config.max_iter)
		{
		  cout << "org " << config.max_iter << " " << error/factor_old << " " << abs(error4)/factor_old << " " << factor_old << " " << RH << " " << Temperature << " " << MO << " " << surrogate[config.iH2O].Ap << endl;
		  //throw string("Stop.");
		}
	    }
		  
	  double error2=1000.0;
	  double AQ=AQinit;
	  index_iter=0;
	  // error2 = AQinit - sum of concentrations of organic compounds in the aqueous phase
	  //           - hygroscopicity of organic compounds
	  //           - LWC - sum of inorganic ions
	  if (config.compute_inorganic)
	    {
	      int nh=config.nh_inorg_init;
	      //Array<double, 1> vec_error_aq,vec_error_chp,vec_error_compaq;
	      Array<double, 1> error_spec,error_spec_old;
	      //vec_error_aq.resize(config.max_iter);
	      //vec_error_chp.resize(config.max_iter);
	      //vec_error_compaq.resize(config.max_iter);
	      error_spec.resize(n);
	      error_spec_old.resize(n);
	      //bool non_convergence;
	      //int iiter=0;
	      double chp_old=chp;
	      double error3=1000.0;
	      double error4=1000.0;
	      double var_rej=0.2;
	      double var_max=0.12;
	      double var_min=0.08;
	      double factor_min=0.01;
	      double factor_max=1; //1.0/pow(2.,nh);
	      if (config.solids)
		{
		  if (RH<0.8)
		    factor_max=0.5;
		}
	      else
		if (RH<0.2)
		  factor_max=0.5;

	      double factor=1.0/pow(2.0,nh);             
	      double relprec=1.0e-3;              
	      double error2_old=0.0;
	      double error3_old=0.0;
	      double error4_old=0.0;
	      double factor_old=factor;
	      double factor_old2=factor;
	      double chp_old2=chp_old;
	      double ionic_old;
	      for (i=0;i<n;i++)
		if (surrogate[i].hydrophilic)
		  {                              
		    surrogate[i].gamma_aq=1.0;
		  }
	      error_spec=1000.;

	      int ntoo_big=0;             
         
	      //if (RH<0.4 or RH>0.9)
	      config.first_evaluation_activity_coefficients=true;                            
	      //config.compute_organic=false; 

	      double RHsave=RH;
	      double error_hyd=0.;
	      //RH=max(RH,0.2);
              
	      while ((index_iter < config.max_iter) and (abs(error2)/factor_old > config.precision or
							 ((abs(error3)/factor_old>relprec  or abs(error4)/factor_old>relprec) and AQ>100*config.MOmin)
							 or config.first_evaluation_activity_coefficients==true or RH>RHsave))
		{                                    
		  if (config.first_evaluation_activity_coefficients==true)
		    {
		      if (abs(error2)/factor_old <= config.precision and abs(error3)/factor_old<=relprec and abs(error4)/factor_old<=relprec and index_iter>0)
			{
			  config.first_evaluation_activity_coefficients=false;                          
			  RH=RHsave;
			}
		    }            
                    
		  if (config.first_evaluation_activity_coefficients==false or index_iter==0)
		    compute_activity_coefficients=true;
		  else
		    compute_activity_coefficients=false;

		  double Asol=0.;
		  for (i=0;i<n;i++)
		    if (surrogate[i].hydrophilic)
		      {
			surrogate[i].Ag_old=surrogate[i].Ag;
			surrogate[i].Aaq_old=surrogate[i].Aaq;
			surrogate[i].gamma_aq_old=surrogate[i].gamma_aq;
		      }
		    else if (surrogate[i].is_solid)
		      {
			Asol+=surrogate[i].Ap;
			surrogate[i].Ap_old=surrogate[i].Ap;
		      }

		  if (AQ<=config.MOmin)
		    error3_old=0.;
		  else
		    error3_old=(chp-chp_old)/max(chp_old,1.0e-20)/factor_old;                                
                    
		  error4_old=error4/factor_old;
		  chp_old2=chp_old;
		  chp_old=chp;		  
		  double AQsave=max(AQ,config.MOmin);
		  double negligeable=max(config.MOmin,0.01*(AQsave+Asol));
                  
		  error2_old=error2/factor_old;                  
		  factor_old2=factor_old;
		  factor_old=factor;
		  error_spec_old=error_spec;
		  ionic_old=ionic;

		  error_inorg_aq_ssh(config,surrogate,AQ,LWC,conc_inorganic,ionic,chp,MMaq,
				     Temperature,error2,derivative,RH,
				     organion,ionic_organic,factor, compute_activity_coefficients);
		  hydratation(config,surrogate,RH,error_hyd,false,false);
                                            
		  double var=0.;                    
		  int m=0;
		  for (i=0;i<n;i++)
		    if (surrogate[i].hydrophilic)
		      {
			if (0.5*surrogate[i].gamma_aq+0.5*surrogate[i].gamma_aq_old>1.0e-6 and (surrogate[i].Aaq>negligeable or surrogate[i].Aaq_old>negligeable))
			  {
			    m++;
			    var+=abs((surrogate[i].gamma_aq_old-surrogate[i].gamma_aq)/(0.5*surrogate[i].gamma_aq+0.5*surrogate[i].gamma_aq_old));
			  }
		      }

		  error4=error_hyd;
		  error_spec=0.0;
		  for (i=0;i<n;i++)
		    if (surrogate[i].hydrophilic)
		      {
			if (surrogate[i].Aaq_old>negligeable and surrogate[i].Aaq>negligeable)                                              
			  {
			    error_spec(i)=(surrogate[i].Aaq-surrogate[i].Aaq_old)/surrogate[i].Aaq_old;                          
			    error4=max(error4,abs(error_spec(i)));                              
			  }
		      }
		    else if (surrogate[i].is_solid and (surrogate[i].Ap>negligeable and surrogate[i].Ap_old>negligeable))
		      {
			error_spec(i)=(surrogate[i].Ap-surrogate[i].Ap_old)/surrogate[i].Ap_old;                          
			error4=max(error4,abs(error_spec(i)));                          
		      }

		  if (m>0)
		    var=var/m; 

		  if (AQ<=config.MOmin)
		    error3=0.;
		  else
		    error3=(chp-chp_old)/max(chp_old,1.0e-15);
                  
                  
		  if (var<=var_rej or factor==factor_min)
		    {
		      AQ=AQ-error2;                      
		      AQ=max(config.MOmin,max(0.1*AQsave,min(AQ,10.0*AQsave)));
		    }                                     

		  if (var>var_max and index_iter>0 and factor>factor_min)
		    {                      
		      factor=max(factor_old*max(var_max/var,0.1),factor_min);            
		      if (var>var_rej)
			for (i=0;i<n;i++)
			  if (surrogate[i].hydrophilic)
			    {                              
			      surrogate[i].Aaq=surrogate[i].Aaq_old;
			      surrogate[i].Ag=surrogate[i].Ag_old;                           
			    }
			  else if (surrogate[i].is_solid)
			    surrogate[i].Ap=surrogate[i].Ap_old;
		    }
		  else if (var<var_min and index_iter>0 and factor<1.)
		    {
		      if (var>0.)
			factor=min(factor_old*min(var_max/var,10.),factor_max);
		      else
			factor=min(factor_old*10.,1.0);
		    }                    
                    

		  if (index_iter>0)
		    {                                         
		      if (error2*error2_old<0.0 and abs(error2)/factor_old>config.precision)                        
			ntoo_big=ntoo_big+10;                        

		      int a=0;    
		      if (error3*error3_old<0.0 and abs(error3)/factor_old>relprec)
			{
			  factor=min(factor,factor_old*(max(abs(error3)/factor_old,abs(error3_old))/(abs(error3)/factor_old+abs(error3_old)))); 
			  ntoo_big=ntoo_big+1;
			}                     
                      
		      if (error4/factor_old>relprec)
			for (i=0;i<n;i++)
			  if (error_spec(i)*error_spec_old(i)<0.0 and abs(error_spec(i))/factor_old>relprec)     
			    a=a+1;                              

		      if (error4>var_rej)
			ntoo_big++; 

		      if (a>0)
			ntoo_big++;
		      factor=max(factor,factor_min);
		    }                  
                                
		  if (ntoo_big>200 and factor_max>factor_min)
		    {
		      factor_max=max(factor_max/2,factor_min);                      
		      ntoo_big=0;                  
		    }                            
                  
		  factor=min(factor,factor_max);
       
		  //vec_error_aq(index_iter)=abs(error2);
		  //vec_error_chp(index_iter)=abs(error3);                      
		  //vec_error_compaq(index_iter)=abs(error4);

		  ++index_iter;
                  
		}	

	      if (index_iter>=config.max_iter)
		config.iiter=0;
	      else
		config.iiter=1;

             
	      if (index_iter>=config.max_iter)
		cout << " non convergence " << abs(error2)/factor_old << " " << abs(error3)/factor_old << " " << abs(error4)/factor_old << " "<< factor_old << " " << ionic << " " << AQ << " " << surrogate[config.iH2O].Aaq << " " << factor << endl;
	    }
	  else
	    {
	      int max_iter2=min(config.max_iter,50);
	      double error2_save=1000;
	      double error2_save2=1000;
	      double factor=1.0;
	      double error_min=1000;
	      double error_hyd=0.;
	      //int nmin=0;
	      //double AQmin=AQ;
	      while ((index_iter < max_iter2) and (abs(error2)/factor > config.precision or error_hyd>1.0e-3))
		{
		  if (config.first_evaluation_activity_coefficients==false or index_iter==0)
		    compute_activity_coefficients=true;
		  else
		    compute_activity_coefficients=false;

		  error_aq_ssh(config,surrogate,AQ,LWC,conc_inorganic,ionic,chp,MMaq,
			       Temperature,error2,derivative,RH,
			       organion,ionic_organic,factor,compute_activity_coefficients);
		  hydratation(config,surrogate,RH,error_hyd,false,false);

		  if (abs(error2_save2-error2)<0.01*abs(error2) and error2_save*error2<0.)
		    {
		      error2_save=1000;
		      error2_save2=1000;
		      factor=max(factor/2,0.01);
		      //cout << "non convergence" << endl;
		    }
	     
		  if (derivative != 0.0 and AQ-error2/derivative >= 0.0)
		    AQ=AQ-error2/derivative;
		  else
		    AQ=AQ-factor*error2;
		  
		  error2_save2=error2_save;
		  error2_save=error2;

		  ++index_iter;
		}


	      
	      // If newton raphson did not solve the system, changed of resolution
	      if ((index_iter < config.max_iter) and (abs(error2) > config.precision))
		{      
		  int nh=config.nh_inorg_init;              
		  //Array<double, 1> vec_error_aq,vec_error_compaq;
		  Array<double, 1> error_spec,error_spec_old;
		  //vec_error_aq.resize(config.max_iter);
		  //vec_error_compaq.resize(config.max_iter);
		  error_spec.resize(n);
		  error_spec_old.resize(n);
		  //bool non_convergence;
		  //int iiter=0;
		  double error4=1000.0;
		  double var_rej=0.2;
		  double var_max=0.12;
		  double var_min=0.08;
		  double factor_min=0.01;
		  double factor_max=1; //1.0/pow(2.,nh);
		  if (config.solids)
		    {
		      if (RH<0.8)
			factor_max=0.5;
		    }
		  else
		    if (RH<0.2)
		      factor_max=0.5;

		  double factor=1.0/pow(2.0,nh);
		  //double factor=1.0/nh;
		  double relprec=1.0e-3;
		  //int maxiiter=30;
		  //cout << nh << endl;
		  double error1=0.0;
		  double error2_old=0.0;
		  double error4_old=0.0;
		  double factor_old=factor;
		  double factor_old2=factor;
		  double deriv_error1_MO,deriv_error1_AQ,deriv_error2_MO,deriv_error2_AQ;
		  double ionic_old;
		  for (i=0;i<n;i++)
		    if (surrogate[i].hydrophilic)		                                   
		      surrogate[i].gamma_aq=1.0;
		      
		  error_spec=1000.;

		  int ntoo_big=0;             
		 
		  //if (RH<0.4 or RH>0.9)
		  //config.first_evaluation_activity_coefficients=true;                            
		  //config.compute_organic=false; 

		  double RHsave=RH;
		  while ((index_iter < config.max_iter) and (abs(error2)/factor_old > config.precision or abs(error4)/factor_old>relprec or
							     config.first_evaluation_activity_coefficients==true or RH>RHsave))
		    {                                    
                    
		      if (config.first_evaluation_activity_coefficients==false or index_iter==0)
			compute_activity_coefficients=true;
		      else
			compute_activity_coefficients=false;

		      double Asol=0.;
		      for (i=0;i<n;i++)
			{
			  if (surrogate[i].hydrophilic)
			    {
			      surrogate[i].Ag_old=surrogate[i].Ag;
			      surrogate[i].Aaq_old=surrogate[i].Aaq;
			      surrogate[i].gamma_aq_old=surrogate[i].gamma_aq;
			    }
			}
                    
		      error4_old=error4/factor_old;		  
		      double AQsave=max(AQ,config.MOmin);
		      double negligeable=max(config.MOmin,0.01*(AQsave+Asol));
                  
		      error2_old=error2/factor_old;                  
		      factor_old2=factor_old;
		      factor_old=factor;
		      error_spec_old=error_spec;
		      ionic_old=ionic;

		      error_aq_ssh(config,surrogate,AQ,LWC,conc_inorganic,ionic,chp,MMaq,
				   Temperature,error2,derivative,RH,
				   organion,ionic_organic,factor,compute_activity_coefficients);
		      hydratation(config,surrogate,RH,error_hyd,false,false);
                                                
		      double var=0.;                    
		      int m=0;
		      for (i=0;i<n;i++)
			if (surrogate[i].hydrophilic)
			  {
			    if (0.5*surrogate[i].gamma_aq+0.5*surrogate[i].gamma_aq_old>1.0e-6 and (surrogate[i].Aaq>negligeable or surrogate[i].Aaq_old>negligeable))
			      {
				m++;
				var+=abs((surrogate[i].gamma_aq_old-surrogate[i].gamma_aq)/(0.5*surrogate[i].gamma_aq+0.5*surrogate[i].gamma_aq_old));
			      }
			  }

		      error4=error_hyd;
		      error_spec=0.0;
		      for (i=0;i<n;i++)
			if (surrogate[i].hydrophilic)
			  {
			    if (surrogate[i].Aaq_old>negligeable and surrogate[i].Aaq>negligeable)                                              
			      {
				error_spec(i)=(surrogate[i].Aaq-surrogate[i].Aaq_old)/surrogate[i].Aaq_old;                          
				error4=max(error4,abs(error_spec(i)));
				//cout << surrogate[i].name << endl;
			      }
			  }
			else if (surrogate[i].is_solid and (surrogate[i].Ap>negligeable and surrogate[i].Ap_old>negligeable))
			  {
			    error_spec(i)=(surrogate[i].Ap-surrogate[i].Ap_old)/surrogate[i].Ap_old;                          
			    error4=max(error4,abs(error_spec(i)));
			    //cout << surrogate[i].name << " " << surrogate[i].Ap << " " << surrogate[i].Ap_old << endl;
			  }

		      if (m>0)
			var=var/m; 
                  
		      if (var<=var_rej or factor==factor_min)
			{
			  AQ=AQ-error2;		       
			  AQ=max(config.MOmin,max(0.1*AQsave,min(AQ,10.0*AQsave)));
			}                                     

		      if (var>var_max and index_iter>0 and factor>factor_min)
			{                      
			  factor=max(factor_old*max(var_max/var,0.1),factor_min);            
			  if (var>var_rej)
			    for (i=0;i<n;i++)
			      {
				surrogate[i].Ag=surrogate[i].Ag_old;				                            
				if (surrogate[i].hydrophilic)				                                
				  surrogate[i].Aaq=surrogate[i].Aaq_old;                                                        
				  
				else if (surrogate[i].is_solid)
				  surrogate[i].Ap=surrogate[i].Ap_old;
			      }
			}
		      else if (var<var_min and index_iter>0 and factor<1.)
			{
			  if (var>0.)
			    factor=min(factor_old*min(var_max/var,10.),factor_max);
			  else
			    factor=min(factor_old*10.,1.0);
			}                    
                    

		      if (index_iter>0)
			{                                         
			  if (error2*error2_old<0.0 and abs(error2)/factor_old>config.precision)                        
			    ntoo_big=ntoo_big+10;
			  int a=0;               
			  if (error4/factor_old>relprec)
			    for (i=0;i<n;i++)
			      if (error_spec(i)*error_spec_old(i)<0.0 and abs(error_spec(i))/factor_old>relprec)     
				{                       
				  a=a+1;
				}

			  if (error4>var_rej)
			    {
			      ntoo_big++; //(error4/factor_old/var_rej);
			    }               
			  if (a>0)
			    ntoo_big++;
			  factor=max(factor,factor_min);
			}                  
                                
		      if (ntoo_big>200 and factor_max>factor_min)
			{
			  factor_max=max(factor_max/2,factor_min);
			  ntoo_big=0;                  
			}                            
                  
		      factor=min(factor,factor_max);
                  
		      //vec_error_aq(index_iter)=abs(error2);                      
		      //vec_error_compaq(index_iter)=abs(error4);
		      ++index_iter;		    
		    }	      
		}   

	      if (index_iter==config.max_iter)
		{
		  cout << "aq " << index_iter << " " << error2 << " " << derivative << " " << RH << " " << Temperature << " " << AQ << " " ;
		  for (i=0;i<n;i++)
		    if (surrogate[i].hydrophilic)
		      cout << surrogate[i].name << " " << surrogate[i].Aaq << " " << surrogate[i].Atot << " ";
		  cout << config.precision << " " << endl; 
		}
	      //cout << index_iter << endl;
	    }
	}
      else
	{
	  //if low concentrations of water assumes that all compounds are hydrophobic
	  all_hydrophobic=true;
	  double error3=1000.0;
	  double derivative=0.0;
	  int index_iter=0;
	  double error_hyd=0.;
	  // minimize the function error with a method of newton raphson
	  // error3 = MOinit - sum of concentrations of organic compounds in the organic phase
	  //           - hygroscopicity of the organic phase  
	  while ((index_iter < config.max_iter) and (abs(error3) > config.precision or error_hyd>1.0e-3))
	    {
	      if (config.first_evaluation_activity_coefficients==false or index_iter==0)
		compute_activity_coefficients=true;
	      else
		compute_activity_coefficients=false;
              
	      error_org_ssh(config,surrogate,MO,MOW,Temperature,error3,derivative,RH,1.0,
			    all_hydrophobic, compute_activity_coefficients);
	      hydratation(config,surrogate,RH,error_hyd,true,false);
	      if (derivative != 0.0 and MO-error3/derivative >= 0.0)
		MO=MO-error3/derivative;
	      else
		MO=MO-error3;

	      ++index_iter;
	    }
	  if (index_iter==config.max_iter)
	    {
	      cout << "tot " << index_iter << " " << RH << " " << Temperature << " " << MO << " " ;
	      for (i=0;i<n;i++)
		if (surrogate[i].is_organic)
		  cout << surrogate[i].name << " " << surrogate[i].Atot << " " ;
	      cout << endl;
	    }
          
	  for (i=0;i<n;++i)
	    if (surrogate[i].is_organic)
	      surrogate[i].Aaq=0.0;
	}
	  
      if (config.compute_saturation and MO > 0.0 and config.compute_organic)
	saturation_ssh(config,surrogate,all_hydrophobic,LWC,ionic,conc_inorganic,ionic_organic,
		       chp,organion,Temperature,RH);
    }

  for(i=0;i<n;i++)
    if (surrogate[i].is_organic)
      surrogate[i].Ag=max(surrogate[i].Atot-surrogate[i].Aaq
			  -surrogate[i].Ap, 0.0);
}

void solve_system_ssh(model_config &config, vector<species>& surrogate,
                      double &MOinit,double &MOW,
                      double &LWC, double &AQinit, double &ionic, double &chp,
                      double &Temperature, double &RH, double &deltat)
{
  //double organion=0.0;
  double conc_inorganic=0.0;
  //double ionic_organic=0.0;
  //double MMaq;  
  int i,it;
  
  
  int n=surrogate.size();
  for (i=0;i<n;++i)
    if (surrogate[i].is_organic==false and i!=config.iH2O)
      conc_inorganic+=surrogate[i].Aaq;
  bool compute_activity_coefficients=true;
   
  if (LWC>config.LWClimit)
    initialisation_eq_ssh(config,surrogate,Temperature,RH,ionic, chp, AQinit,LWC,false);
  else
    initialisation_eq_ssh(config,surrogate,Temperature,RH,ionic,chp,AQinit,LWC,true);   
  
  if (config.chemistry)
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
		solve_equilibrium_ssh(config, surrogate, MOinit,MOW, LWC, AQinit, ionic, chp,
                              Temperature, RH);
		
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
		solve_equilibrium_ssh(config, surrogate, MOinit,MOW, LWC, AQinit, ionic, chp,
                              Temperature, RH);
		
		dt1=min(deltat-t,dt1);	 
		dt2=dt1;
		//cout << dt1 << endl;
		integer_chem_ssh(config, surrogate, MOinit, MOW, MMaq, LWC, AQinit, ionic, chp, Temperature, RH, dt1, compute_activity_coefficients);

		//compute the new time step so that changes are small
		adapstep_chem_ssh(config,surrogate,dt1,t,deltat,config.dtchem_min);
	  
		t+=dt2;     
	      }
	  }

	/*
    for (it=0;it<config.nt;it++)
      { 
        solve_equilibrium_ssh(config, surrogate, MOinit,MOW, LWC, AQinit, ionic, chp,
                              Temperature, RH);
        solve_chemistry_ssh(config, surrogate, MOinit,MOW, LWC, AQinit, ionic, chp,
                            Temperature, RH, deltat/config.nt, compute_activity_coefficients);
      }*/
    }
 
  solve_equilibrium_ssh(config, surrogate, MOinit,MOW, LWC, AQinit, ionic, chp,
                        Temperature, RH);
}

void global_equilibrium_ssh(model_config &config, vector<species>& surrogate,
                            double &MOinit,double &MOW,
                            double &LWC, double &AQinit, double &ionic, double &chp,
                            double &Temperature, double &RH, double &deltat)
{
  int icycle;
  if (config.compute_inorganic)
    {
      LWC=1.0;
      config.LWClimit=-1.0;
      if (config.compute_long_and_medium_range_interactions)
	config.compute_aqueous_phase_properties=true;
      else
	config.compute_aqueous_phase_properties=false;   
    }
  else
    config.compute_aqueous_phase_properties=false;

  if (config.coupling_organic_inorganic or config.compute_organic==false 
      or config.compute_inorganic==false)
    {

      //config.solids=false;
      /*solve_system_ssh(config, surrogate, MOinit, MOW, LWC, AQinit, ionic, 
        chp, Temperature, RH);*/
      //config.solids=true;
      //config.compute_organic=false;
      solve_system_ssh(config, surrogate, MOinit, MOW, LWC, AQinit, ionic, 
                       chp, Temperature, RH, deltat);   
    }
  else
    for (icycle=0;icycle<config.number_of_org_inorg_cycles;icycle++)
      {
        config.compute_organic=false;
        solve_system_ssh(config, surrogate, MOinit, MOW, LWC, AQinit, ionic, 
                         chp, Temperature, RH, deltat);
        config.compute_inorganic=false;
	

        config.compute_organic=true;
        solve_system_ssh(config, surrogate, MOinit, MOW, LWC, AQinit, ionic, 
                         chp, Temperature, RH, deltat);
        config.compute_inorganic=true;
      }
}

void solve_local_equilibriums_uncoupled_ssh(model_config config, vector<species> &surrogate,
                                            Array<double, 3> &MOinit, Array<double, 3> &MOW, Array<double, 1> &number,
                                            Array<double, 1> &Vsol,
                                            Array<double, 1> &LWC, Array<double, 1> &AQinit, Array<double, 1> &ionic,
                                            Array<double, 1> &chp,
                                            double &Temperature, double &RH,
                                            Array<double, 1> &AQ, Array<double, 3> &MO,
                                            Array<double, 1> &conc_inorganic, Array<double, 1> &ionic_organic, Array<double, 1> &organion, Array<double, 1> &MMaq)
{
  //The method of newton raphson cannot be used due to the high number of variables
  //(nbins*(1+nlayer*nphase)
  //
  //Method to reach equilibrium
  // - initialisation of MOinit and AQinit
  // - initialisation of MO and AQ
  // - compute error_tot = maximum of relative errors
  //   if error_tot < 0.1% the system has converged
  // - if error_tot = error_old2 (error from iteration - 2) the system cannot converge
  //   due to strong variation of the composition.
  //   in that case a method is used so that variations are smaller:
  //     Ap(iter+1)=Ap(iter)*(1-1.0/nh)+1.0/nh*Ap(computed)
  //     1.0/nh weight factor
  //     nh is increased to ensure convergence
  //     error_tot must be inferior to factor * 0.1% to reach equilibrium
  
  int b,ilayer,iphase,i;
  int n=surrogate.size();
  double error_org=10.0;
  double error_aq=10.0;
  int index=0;
  Array<double, 1> vec_error_org,vec_error_aq;
  vec_error_org.resize(config.max_iter);
  vec_error_aq.resize(config.max_iter);
  bool non_convergence;
  double LWCtot=0.0;
  for (b=0;b<config.nbins;++b)
    LWCtot+=LWC(b);
  Array<double, 1> AQsave;
  Array<double, 3> MOsave;
  AQsave.resize(config.nbins);
  MOsave.resize(config.nbins,config.nlayer,config.max_number_of_phases);

  int nh_aq;
  if (config.compute_inorganic and config.compute_organic==false)
    nh_aq=config.nh_inorg_init;
  else if (config.compute_inorganic and config.compute_organic)
    nh_aq=max(config.nh_aq_init,config.nh_inorg_init);
  else
    nh_aq=config.nh_aq_init;
  int nh_org=config.nh_org_init;

  for (i=0;i<n;i++)
    {
      surrogate[i].Ag0=surrogate[i].Ag;
      for(b=0;b<config.nbins;b++)
        {
          if (surrogate[i].hydrophilic)
            surrogate[i].Aaq_bins_init0(b)=surrogate[i].Aaq_bins_init(b);
          if (surrogate[i].hydrophobic)
            for (ilayer=0;ilayer<config.nlayer;++ilayer)
              for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                surrogate[i].Ap_layer_init0(b,ilayer,iphase)=surrogate[i].Ap_layer_init(b,ilayer,iphase);          
        }
    }

  for(b=0;b<config.nbins;b++)
    {
      AQsave(b)=AQinit(b);
      for (ilayer=0;ilayer<config.nlayer;++ilayer)
        for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
          MOsave(b,ilayer,iphase)=MOinit(b,ilayer,iphase); 
    }

  int iorg=0;
  int iaq=0;

  while ((error_org>1.0/nh_org*config.relative_precision
          or error_aq>1.0/nh_aq*config.relative_precision)
         and index < config.max_iter)
    {
      if (index>2)
        {
          non_convergence=false;
          if (iorg>20)
            for (i=max(index-20,2);i<index-1;i++)
              if (((abs(vec_error_org(i)-vec_error_org(index-1))/vec_error_org(index-1)*nh_org<1.0e-4 and vec_error_org(index-1)*nh_org>config.precision) and
                   (abs(vec_error_org(i-1)-vec_error_org(index-2))/vec_error_org(index-2)*nh_org<1.0e-4) and (abs(vec_error_org(i-2)-vec_error_org(index-3))/vec_error_org(index-3)*nh_org<1.0e-4))
                  or vec_error_org(index-1)*nh_org>10.0)
                non_convergence=true;
              
          if (non_convergence and nh_org<config.nh_max)
            { 
              if (vec_error_org(index-1)*nh_org>100.0)
                {
                  for (i=0;i<n;i++)
                    if (surrogate[i].hydrophobic)
                      {
                        surrogate[i].Ag=surrogate[i].Ag0;
                        for(b=0;b<config.nbins;b++)
                          for (ilayer=0;ilayer<config.nlayer;++ilayer)
                            for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                              surrogate[i].Ap_layer_init(b,ilayer,iphase)=surrogate[i].Ap_layer_init0(b,ilayer,iphase);                       
                      }
                }
              ++nh_org;
              for (i=max(index-50,0);i<index;i++)
                vec_error_org(i)=-1.0;
            }
          
          non_convergence=false;
          if (iaq>20)
            for (i=max(index-20,2);i<index-1;i++)
              if (((abs(vec_error_aq(i)-vec_error_aq(index-1))/vec_error_aq(index-1)*nh_aq<1.0e-4 and vec_error_aq(index-1)*nh_aq>config.precision) and
                   (abs(vec_error_aq(i-1)-vec_error_aq(index-2))/vec_error_aq(index-2)*nh_aq<1.0e-4) and (abs(vec_error_aq(i-2)-vec_error_aq(index-3))/vec_error_aq(index-3)*nh_aq<1.0e-4))
                  or vec_error_aq(index-1)*nh_aq>10.0)
                non_convergence=true;
              
          if (non_convergence and nh_aq<config.nh_max)
            { 
              if (vec_error_aq(index-1)*nh_aq>100.0)
                {
                  for (i=0;i<n;i++)
                    if (surrogate[i].hydrophilic)
                      {
                        surrogate[i].Ag=surrogate[i].Ag0;
                        for(b=0;b<config.nbins;b++)
                          surrogate[i].Aaq_bins_init(b)=surrogate[i].Aaq_bins_init0(b);                       
                      }
                }
              ++nh_aq;
              for (i=max(index-50,0);i<index;i++)
                vec_error_aq(i)=-1.0;
            }
        }
      iaq++;
      iorg++;

      for (b=0;b<config.nbins;++b)
        for (ilayer=0;ilayer<config.nlayer;++ilayer)
          for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
            for (i=0;i<n;++i)
              surrogate[i].Ap_layer(b,ilayer,iphase)=
                surrogate[i].Ap_layer_init(b,ilayer,iphase);

      for (b=0;b<config.nbins;++b)
        for (i=0;i<n;i++)
          surrogate[i].Aaq_bins(b)=surrogate[i].Aaq_bins_init(b);
		  
      //if the system has not converged for organic concentrations
      if (error_org>1.0/nh_org*config.relative_precision)
        {
          //compute concentrations at equilibrium for the organic phase
          //if (config.first_evaluation_activity_coefficients==false)
          if (config.first_evaluation_activity_coefficients==false)
            {
              equilibrium_org_ssh(config, surrogate, config.tequilibrium, MOinit,MO,
                                  Temperature, MOW, true,1.0/nh_org);
              if (config.compute_saturation and config.compute_organic)
                phase_repartition_ssh(config,surrogate,Temperature,MOinit,MO,MOW,1.0/nh_org);
            }
          else
            {
              equilibrium_org_ssh(config, surrogate, config.tequilibrium, MOinit,MO,
                                  Temperature, MOW, false,1.0/nh_org);
              if (config.compute_saturation and config.compute_organic)
                phase_repartition_ssh(config,surrogate,Temperature,MOinit,MO,MOW,1.0/nh_org);
            }

          //redistribute concentrations to ensure that the volume of layers are constant
          redistribution_ssh(config, surrogate,MOinit,MO);
        }

      if (error_aq>1.0/nh_aq*config.relative_precision and LWCtot>config.LWClimit)
        density_aqueous_phase_ssh(config, surrogate, LWC, Temperature);
		  
      //if the system has not converged for aqueous concentrations

      if (error_aq>1.0/nh_aq*config.relative_precision and LWCtot>config.LWClimit)
        if (config.first_evaluation_activity_coefficients==false)
          equilibrium_aq_ssh(config, surrogate, config.tequilibrium, AQinit, AQ,
                             MOinit,conc_inorganic, ionic, ionic_organic, organion,
                             chp, LWC, Temperature, RH, MMaq, true,1.0/nh_aq);
        else
          equilibrium_aq_ssh(config, surrogate, config.tequilibrium, AQinit, AQ,
                             MOinit,conc_inorganic, ionic, ionic_organic, organion,
                             chp, LWC, Temperature, RH, MMaq, false,1.0/nh_aq);

      water_concentration_ssh(config, surrogate, Temperature, RH);

      //compute the new diameters of particles
      compute_diameters_ssh(config, surrogate, Vsol, number, LWC, LWCtot);		  		  
      if (config.explicit_representation)
	compute_morphology_ssh(config, Vsol, number);
      //Computation of error_aq and error_tot
      error_aq=0.0;
      error_org=0.0;          
          
      for (b=0;b<config.nbins;++b)
        {
          for (ilayer=0;ilayer<config.nlayer;++ilayer)
            for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
              {
                if (MO(b,ilayer,iphase)>1.0e-5)
                  error_org=max(error_org,
                                abs(MO(b,ilayer,iphase)-MOinit(b,ilayer,iphase))
                                /MO(b,ilayer,iphase));
                MOinit(b,ilayer,iphase)=MO(b,ilayer,iphase);

                for (i=0;i<n;i++)
                  if (surrogate[i].hydrophobic and
                      surrogate[i].Ap_layer_init(b,ilayer,iphase)>1.0e-5 and surrogate[i].Ap_layer(b,ilayer,iphase)>1.0e-5)
                    error_org=max(error_org,
                                  abs(surrogate[i].Ap_layer_init(b,ilayer,iphase)
                                      -surrogate[i].Ap_layer(b,ilayer,iphase))
                                  /surrogate[i].Ap_layer(b,ilayer,iphase));   
              }

          if (LWCtot>config.LWClimit)
            {
              if (AQ(b)>1.0e-5)
                error_aq=max(error_aq,abs(AQ(b)-AQinit(b))/AQ(b));
              AQinit(b)=AQ(b);

              for (i=0;i<n;i++)
                if (surrogate[i].hydrophilic and
                    surrogate[i].Aaq_bins_init(b)>1.0e-5 and surrogate[i].Aaq_bins(b)>1.0e-5 and
                    i!=config.iHp and i!=config.iOHm)
                  error_aq=max(error_aq,abs(surrogate[i].Aaq_bins_init(b)
                                            -surrogate[i].Aaq_bins(b))
                               /surrogate[i].Aaq_bins(b));

            }
        }

      vec_error_org(index)=error_org;
      vec_error_aq(index)=error_aq;
		  
      //Computation of characteristic times to reach equilibrium
      tau_dif_ssh(config, surrogate, number, Vsol, Temperature);
      tau_kmt_ssh(config, surrogate, Temperature, number);
      /*      if (config.explicit_representation)
              compute_morphology_ssh(config, Vsol);*/

      characteristic_time_ssh(config, surrogate, MOinit, AQinit, LWCtot); 
      if (LWCtot>config.LWClimit)
        characteristic_time_aq_ssh(config, surrogate, Temperature, chp, LWC, AQinit, MOinit, MMaq, ionic);  

      ++index;
    }


  if (config.compute_saturation and config.first_evaluation_of_saturation==false and config.compute_organic)
    {
      number_org_phases_ssh(config,surrogate,Temperature,MOinit,MOW);
      compute_kp_org_ssh(config, surrogate, MOinit, Temperature, MOW);  
      tau_dif_ssh(config, surrogate, number, Vsol, Temperature);
      tau_kmt_ssh(config, surrogate, Temperature, number);
      characteristic_time_ssh(config, surrogate, MOinit, AQinit, LWCtot);
      if (LWCtot>config.LWClimit)
        characteristic_time_aq_ssh(config, surrogate, Temperature, chp, LWC, AQinit, MOinit, MMaq, ionic);
    }

  if (error_org>config.relative_precision/nh_org or error_aq>config.relative_precision/nh_aq)
    {
      cout << "The model did not converged..." << endl;
      cout << "uncoupled" << endl;
      cout << config.nphase << endl;
      cout << "error org" << error_org*nh_org << " " << nh_aq << endl;
      cout << "error aq" << error_aq*nh_aq << " " << nh_org << endl;
      cout << "avant" << MOinit << endl;
      cout << "avant" << AQinit << endl;
      cout << "apres" << MO << endl;
      cout << "apres" << AQ << endl;
    }
}

void solve_local_equilibriums_coupled_ssh(model_config config, vector<species> &surrogate,
                                          Array<double, 3> &MOinit, Array<double, 3> &MOW, Array<double, 1> &number,
                                          Array<double, 1> &Vsol,
                                          Array<double, 1> &LWC, Array<double, 1> &AQinit, Array<double, 1> &ionic,
                                          Array<double, 1> &chp,
                                          double &Temperature, double &RH,
                                          Array<double, 1> &AQ, Array<double, 3> &MO,
                                          Array<double, 1> &conc_inorganic, Array<double, 1> &ionic_organic, Array<double, 1> &organion, Array<double, 1> &MMaq)
{
  int b,ilayer,iphase,i;
  int n=surrogate.size();
  double error_tot=10.0;
  int index=0;
  Array<double, 1> vec_error_org,vec_error_aq;
  vec_error_org.resize(config.max_iter);
  vec_error_aq.resize(config.max_iter);
  bool non_convergence;
  double LWCtot=0.0;
  for (b=0;b<config.nbins;++b)
    LWCtot+=LWC(b);
  Array<double, 1> AQsave;
  Array<double, 1> chp_save;
  Array<double, 3> MOsave;
  AQsave.resize(config.nbins);
  chp_save.resize(config.nbins);
  MOsave.resize(config.nbins,config.nlayer,config.max_number_of_phases);

  //The method of newton raphson cannot be used due to the high number of variables
  //(nbins*(1+nlayer*nphase)
  //
  //Method to reach equilibrium
  // - initialisation of MOinit and AQinit
  // - initialisation of MO and AQ
  // - compute error_tot = maximum of relative errors
  //   if error_tot < 0.1% the system has converged
  // - if error_tot = error_old2 (error from iteration - 2) the system cannot converge
  //   due to strong variation of the composition.
  //   in that case a method is used so that variations are smaller:
  //     Ap(iter+1)=Ap(iter)*(1-1.0/nh)+1.0/nh*Ap(computed)
  //     1.0/nh weight factor
  //     nh is increased to ensure convergence
  //     error_tot must be inferior to factor * 0.1% to reach equilibrium
  int iiter=0;
  int nh;

  if (config.compute_inorganic and config.compute_organic==false)
    nh=config.nh_inorg_init;
  else if (config.compute_inorganic and config.compute_organic)
    nh=max(config.nh_inorg_init,max(config.nh_aq_init,config.nh_org_init));
  else
    nh=max(config.nh_aq_init,config.nh_org_init);

  for (i=0;i<n;i++)
    {
      surrogate[i].Ag0=surrogate[i].Ag;
      for(b=0;b<config.nbins;b++)
        {
          if (surrogate[i].hydrophilic)
            surrogate[i].Aaq_bins_init0(b)=surrogate[i].Aaq_bins_init(b);
          if (surrogate[i].hydrophobic)
            for (ilayer=0;ilayer<config.nlayer;++ilayer)
              for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                surrogate[i].Ap_layer_init0(b,ilayer,iphase)=surrogate[i].Ap_layer_init(b,ilayer,iphase); 
        }
    }

  for(b=0;b<config.nbins;b++)
    {
      AQsave(b)=AQinit(b);
      for (ilayer=0;ilayer<config.nlayer;++ilayer)
        for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
          MOsave(b,ilayer,iphase)=MOinit(b,ilayer,iphase); 
    }

  while (error_tot>config.relative_precision/nh and index < config.max_iter) 
    {       
      if (index>2)
        {
          //ensure that the system can converge
          non_convergence=false;   
          if (iiter>20)
            for (i=max(index-20,2);i<index-1;i++)
              if (((abs(vec_error_org(i)-vec_error_org(index-1))/vec_error_org(index-1)*nh<1.0e-4 and vec_error_org(index-1)*nh>config.precision) and
                   (abs(vec_error_org(i-1)-vec_error_org(index-2))/vec_error_org(index-2)*nh<1.0e-4) and (abs(vec_error_org(i-2)-vec_error_org(index-3))/vec_error_org(index-3)*nh<1.0e-4))
                  or (abs(vec_error_aq(i)-vec_error_aq(index-1))/vec_error_aq(index-1)*nh<1.0e-4 and vec_error_aq(index-1)*nh>config.precision and 
                      (abs(vec_error_aq(i-1)-vec_error_aq(index-2))/vec_error_aq(index-2)*nh<1.0e-4) and (abs(vec_error_aq(i-2)-vec_error_aq(index-3))/vec_error_aq(index-3)*nh<1.0e-4))                  
                  or vec_error_org(index-1)*nh>10.0 or vec_error_aq(index-1)*nh>10.0)
                non_convergence=true;

          if (iiter>100)
            if (vec_error_aq(index-1)*nh>config.precision)
              {
                double a=1;
                double b=1;
                for (i=index-100;i<index-90;i++)
                  a*=vec_error_aq(i);
              
                
                for (i=index-11;i<index-1;i++)
                  b*=vec_error_aq(i);
                  
                if (b>a*0.999)
                  non_convergence=true;
              }

          if (iiter>100)
            if (vec_error_org(index-1)*nh>config.precision)
              {
                double a=1;
                double b=1;
                for (i=index-100;i<index-90;i++)
                  a*=vec_error_org(i);
                
                for (i=index-11;i<index-1;i++)
                  b*=vec_error_org(i);
                  
                if (b>a*0.999)
                  non_convergence=true;
              }
          
              
          if (non_convergence and nh<config.nh_max)
            { 
              if (vec_error_org(index-1)*nh>100.0 or vec_error_aq(index-1)*nh>100.0)
                {
                  for (i=0;i<n;i++)
                    {
                      surrogate[i].Ag=surrogate[i].Ag0;
                      if (surrogate[i].hydrophobic)
                        for(b=0;b<config.nbins;b++)
                          for (ilayer=0;ilayer<config.nlayer;++ilayer)
                            for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                              surrogate[i].Ap_layer_init(b,ilayer,iphase)=surrogate[i].Ap_layer_init0(b,ilayer,iphase);          
                      if (surrogate[i].hydrophilic)
                        for(b=0;b<config.nbins;b++)
                          surrogate[i].Aaq_bins_init(b)=surrogate[i].Aaq_bins_init0(b);             
                    }
                }
              ++nh;
              for (i=max(index-50,0);i<index;i++)
                {
                  vec_error_org(i)=-1.0;
                  vec_error_aq(i)=-1.0;
                }
              iiter=0;

              for(b=0;b<config.nbins;b++)
                {
                  AQinit(b)=AQsave(b);
                  for (ilayer=0;ilayer<config.nlayer;++ilayer)
                    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                      MOinit(b,ilayer,iphase)=MOsave(b,ilayer,iphase); 
                }
            }
        }

      for (b=0;b<config.nbins;++b)
        for (ilayer=0;ilayer<config.nlayer;++ilayer)
          for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
            for (i=0;i<n;++i)
              surrogate[i].Ap_layer(b,ilayer,iphase)=
                surrogate[i].Ap_layer_init(b,ilayer,iphase);

      for (b=0;b<config.nbins;++b)
        for (i=0;i<n;i++)
          surrogate[i].Aaq_bins(b)=surrogate[i].Aaq_bins_init(b);

      iiter++;   
      if (LWCtot>config.LWClimit)
        density_aqueous_phase_ssh(config, surrogate, LWC, Temperature);
          
      for (b=0;b<config.nbins;b++)
	chp_save(b)=chp(b);
      
      if (config.first_evaluation_activity_coefficients==false)
        {	  
          equilibrium_tot_ssh(config, surrogate, config.tequilibrium, AQinit, AQ, conc_inorganic,
                              ionic, ionic_organic, organion, chp, LWC, MOinit, MO,
                              MOW, Temperature, RH, MMaq, true, 1.0/nh);       
          if (config.compute_saturation and config.compute_organic)
            phase_repartition_ssh(config,surrogate,Temperature,MOinit,MO,MOW,1.0/nh);	  
        }
      else
        {
          equilibrium_tot_ssh(config, surrogate, config.tequilibrium, AQinit, AQ, conc_inorganic,
                              ionic, ionic_organic, organion, chp, LWC, MOinit, MO,
                              MOW, Temperature, RH, MMaq, false, 1.0/nh);
          if (config.compute_saturation and config.compute_organic)
            phase_repartition_ssh(config,surrogate,Temperature,MOinit,MO,MOW,1.0/nh);
	}      

      //redistribute concentrations to ensure that the volume of layers are constant      
      redistribution_ssh(config, surrogate,MOinit,MO);

      water_concentration_ssh(config, surrogate, Temperature, RH);

      //compute the new diameters of particle due to the growth of particles by condensation
      compute_diameters_ssh(config, surrogate, Vsol, number, LWC, LWCtot);
      if (config.explicit_representation)
	compute_morphology_ssh(config, Vsol, number);
		  
      //Computation of error_tot
      vec_error_org(index)=0.0;
      vec_error_aq(index)=0.0;
      for (b=0;b<config.nbins;++b)
        {
          AQ(b)=max(AQ(b),config.MOmin);
          for (ilayer=0;ilayer<config.nlayer;++ilayer)
            for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
              {
                if (MO(b,ilayer,iphase)>1.0e-5)
                  vec_error_org(index)=max(vec_error_org(index),
                                           abs(MO(b,ilayer,iphase)-MOinit(b,ilayer,iphase))
                                           /MO(b,ilayer,iphase));
                
                for (i=0;i<n;i++)
                  if (surrogate[i].hydrophobic and
                      surrogate[i].Ap_layer(b,ilayer,iphase)>1.0e-5)
                    vec_error_org(index)=max(vec_error_org(index),
                                             abs(surrogate[i].Ap_layer_init(b,ilayer,iphase)
                                                 -surrogate[i].Ap_layer(b,ilayer,iphase))
                                             /surrogate[i].Ap_layer(b,ilayer,iphase));   

                MOinit(b,ilayer,iphase)=max(MO(b,ilayer,iphase),config.MOmin*config.Vlayer(ilayer));
              }

          if (AQ(b)>1.0e-5)
            vec_error_aq(index)=max(vec_error_aq(index),abs(AQ(b)-AQinit(b))/AQ(b));

	  if (chp(b)>0.)
	    vec_error_aq(index)=max(vec_error_aq(index),abs(chp(b)-chp_save(b))/chp(b));
       
          for (i=0;i<n;i++)
            if (surrogate[i].hydrophilic and i!=config.iHp and i!=config.iOHm)
	      if (surrogate[i].Aaq_bins(b)>1.0e-5)
		vec_error_aq(index)=max(vec_error_aq(index),abs(surrogate[i].Aaq_bins_init(b)
								-surrogate[i].Aaq_bins(b))
                                        /surrogate[i].Aaq_bins(b));
			  
          AQinit(b)=AQ(b); //max(AQ(b),config.MOmin);

	  
        }

      vec_error_org(index)=max(vec_error_org(index),0.001*config.precision/nh);
      vec_error_aq(index)=max(vec_error_aq(index),0.001*config.precision/nh);

      error_tot=max(vec_error_org(index),vec_error_aq(index));
      //cout << error_tot << " " << surrogate[config.iH2O].gamma_aq_bins << " " << chp << endl;
      //cout << config.gamma_MR_ions << endl;
      //Computation of characteristic times to reach equilibrium
      ; 
      tau_dif_ssh(config, surrogate, number, Vsol, Temperature);
      tau_kmt_ssh(config, surrogate, Temperature, number);
      characteristic_time_ssh(config, surrogate, MOinit, AQinit, LWCtot); 
      if (LWCtot>config.LWClimit)
        characteristic_time_aq_ssh(config, surrogate, Temperature, chp, LWC, AQinit, MOinit, MMaq, ionic);      

      ++index;  
    }

  //cout << "niter: " << index << " " << chp(0) << " " << chp_save(0) << endl; 
  
  if (config.compute_saturation and config.first_evaluation_of_saturation==false and config.compute_organic)
    {
      number_org_phases_ssh(config,surrogate,Temperature,MOinit,MOW);
      compute_kp_org_ssh(config, surrogate, MOinit, Temperature, MOW);
      tau_dif_ssh(config, surrogate, number, Vsol, Temperature);
      tau_kmt_ssh(config, surrogate, Temperature, number);
      characteristic_time_ssh(config, surrogate, MOinit, AQinit, LWCtot);
      if (LWCtot>config.LWClimit)
        characteristic_time_aq_ssh(config, surrogate, Temperature, chp, LWC, AQinit, MOinit, MMaq, ionic);
    }

  if (error_tot>config.relative_precision/nh)
    {
      cout << "The model did not converged..." << endl;
      cout << "coupled" << endl;
      for (b=0;b<config.nbins;b++)
	cout << b << " " << abs(chp(b)-chp_save(b))/chp(b)*nh << endl;
      cout << config.nphase << " " << nh << endl;
      cout << surrogate[config.iH2O].time_aq << endl;
      cout << "error" << error_tot << endl;
      cout << "avant" << MOinit << endl;
      cout << "avant" << AQinit << endl;
      cout << "apres" << MO << endl;
      cout << "apres" << AQ << endl;
      for (i=0;i<n;i++)
	cout << surrogate[i].name << " " << sum(surrogate[i].Aaq_bins_init) << " " << surrogate[i].gamma_aq_bins << endl;
      cout << chp << endl;
    }
  
}


void solve_implicit_water_coupled_ssh(model_config config, vector<species> &surrogate,
                                      Array<double, 3> &MOinit, Array<double, 3> &MOW, Array<double, 1> &number,
                                      Array<double, 1> &Vsol,
                                      Array<double, 1> &LWC, Array<double, 1> &AQinit, Array<double, 1> &ionic,
                                      Array<double, 1> &chp,
                                      double &Temperature, double &RH,
                                      Array<double, 1> &AQ, Array<double, 3> &MO,
                                      Array<double, 1> &conc_inorganic, Array<double, 1> &ionic_organic, Array<double, 1> &organion, Array<double, 1> &MMaq,
                                      double t, double deltat)
{

  int b,ilayer,iphase,i;
  int n=surrogate.size();
  double error_tot=10.0;
  int index=0;
  Array<double, 1> vec_error_org,vec_error_aq,vec_error_gas,vec_error_chp; 
  vec_error_org.resize(config.max_iter);
  vec_error_aq.resize(config.max_iter);
  vec_error_gas.resize(config.max_iter);
  vec_error_chp.resize(config.max_iter);
  vec_error_org=0.;
  vec_error_aq=0.;
  vec_error_gas=0.;
  vec_error_chp=0.;
  bool non_convergence;
  double LWCtot=0.0;
  for (b=0;b<config.nbins;++b)
    LWCtot+=LWC(b);
  Array<double, 1> AQsave;
  Array<double, 1> chp_save;
  Array<double, 3> MOsave;
  AQsave.resize(config.nbins);
  chp_save.resize(config.nbins);
  MOsave.resize(config.nbins,config.nlayer,config.max_number_of_phases);

  for (i=0;i<n;i++)
    //if (i!=config.iH2O)
    {
      surrogate[i].Atot0=surrogate[i].Atot;
      surrogate[i].Ag0=surrogate[i].Ag;
      for(b=0;b<config.nbins;b++)
        {
          if (surrogate[i].hydrophilic)
            surrogate[i].Aaq_bins_init0(b)=surrogate[i].Aaq_bins_init(b);
          if (surrogate[i].hydrophobic)
            for (ilayer=0;ilayer<config.nlayer;++ilayer)
              for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                surrogate[i].Ap_layer_init0(b,ilayer,iphase)=surrogate[i].Ap_layer_init(b,ilayer,iphase); 
        }
      
      if (surrogate[i].is_solid)
	for(b=0;b<config.nbins;b++)
	  surrogate[i].Asol_bins_init0(b)=surrogate[i].Asol_bins_init(b);
    }

  /*
    if (sum(surrogate[config.iH2O].Ap_layer_init)+sum(surrogate[config.iH2O].Aaq_bins_init)>100.)
    {      
    cout << "intermediaire " << sum(surrogate[config.iH2O].Ap_layer_init) << " " << sum(surrogate[config.iH2O].Aaq_bins_init) << " " << chp << " " << RH << " " << sum(surrogate[config.iH2O].Ap_layer_init0) << " " << sum(surrogate[config.iH2O].Aaq_bins_init0) <<endl;
    }*/
  
  index=0;
  int iiter=0;
  error_tot=1000;
  bool compute_inorganic_save=config.compute_inorganic;
  bool compute_organic_save=config.compute_organic;
  config.compute_inorganic=false;
  config.compute_organic=false;
  int nh;
  if (config.compute_inorganic and config.compute_organic==false)
    nh=config.nh_inorg_init;
  else if (config.compute_inorganic and config.compute_organic)
    nh=max(config.nh_inorg_init,max(config.nh_aq_init,config.nh_org_init));
  else
    nh=max(config.nh_aq_init,config.nh_org_init);
  config.nh_max=5;
  
  double factor=pow(0.5,nh-1);
  double factor_min=pow(0.5,config.nh_max-1);
  vec_error_org=-1;
  vec_error_aq=-1;
  vec_error_gas=-1;
  vec_error_chp=-1;
  double maxaq=0.;

  for(b=0;b<config.nbins;b++)
    {
      AQsave(b)=max(AQinit(b)-surrogate[config.iH2O].Aaq_bins_init(b)+surrogate[config.iH2O].Aaq_bins_init0(b),config.MOmin);
      for (ilayer=0;ilayer<config.nlayer;++ilayer)
        for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
          MOsave(b,ilayer,iphase)=max(MOinit(b,ilayer,iphase)-surrogate[config.iH2O].Ap_layer_init(b)+surrogate[config.iH2O].Ap_layer_init0(b),config.MOmin); 
    }

  while ((error_tot>config.relative_precision*factor or maxaq>0.01) and index < config.max_iter) 
    {     
      water_concentration_ssh(config, surrogate, Temperature, RH);
      if (index>2)
	{
	  //ensure that the system can converge
	  non_convergence=false;   
	  if (iiter>20)
	    for (i=max(index-20,2);i<index-1;i++)
	      if (((abs((vec_error_org(i)-vec_error_org(index-1))/vec_error_org(index-1)/factor)<1.0e-4 and abs(vec_error_org(index-1))/factor>config.relative_precision) and
		   (abs((vec_error_org(i-1)-vec_error_org(index-2))/vec_error_org(index-2))/factor<1.0e-4) and (abs((vec_error_org(i-2)-vec_error_org(index-3))/vec_error_org(index-3))/factor<1.0e-4))
		  or (abs((vec_error_aq(i)-vec_error_aq(index-1))/vec_error_aq(index-1))/factor<1.0e-4 and abs(vec_error_aq(index-1))/factor>config.relative_precision and 
		      (abs((vec_error_aq(i-1)-vec_error_aq(index-2))/vec_error_aq(index-2))/factor<1.0e-4) and (abs((vec_error_aq(i-2)-vec_error_aq(index-3))/vec_error_aq(index-3))/factor<1.0e-4))
		  or (abs((vec_error_gas(i)-vec_error_gas(index-1))/vec_error_gas(index-1))/factor<1.0e-4 and abs(vec_error_gas(index-1))/factor>config.relative_precision and 
		      (abs((vec_error_gas(i-1)-vec_error_gas(index-2))/vec_error_gas(index-2))/factor<1.0e-4) and (abs((vec_error_gas(i-2)-vec_error_gas(index-3))/vec_error_gas(index-3))/factor<1.0e-4)) 
		  or abs(vec_error_org(index-1))/factor>10.0 or abs(vec_error_aq(index-1))/factor>10.0 or abs(vec_error_gas(index-1))/factor>10.0)
		non_convergence=true;

          if (iiter>20)
            if (((abs((vec_error_org(index-3)-vec_error_org(index-1))/vec_error_org(index-1))/factor<1.0e-2 and abs(vec_error_org(index-1))/factor>config.relative_precision) and
                 (abs((vec_error_org(index-4)-vec_error_org(index-2))/vec_error_org(index-2))/factor<1.0e-2) and vec_error_org(index-1)*vec_error_org(index-2)<0.) or
                ((abs((vec_error_aq(index-3)-vec_error_aq(index-1))/vec_error_aq(index-1))/factor<1.0e-2 and abs(vec_error_aq(index-1))/factor>config.relative_precision) and
                 (abs((vec_error_aq(index-4)-vec_error_aq(index-2))/vec_error_aq(index-2))/factor<1.0e-2) and vec_error_aq(index-1)*vec_error_aq(index-2)<0.) or
                ((abs((vec_error_gas(index-3)-vec_error_gas(index-1))/vec_error_gas(index-1))/factor<1.0e-2 and abs(vec_error_gas(index-1))/factor>config.relative_precision) and
                 (abs((vec_error_gas(index-4)-vec_error_gas(index-2))/vec_error_gas(index-2))/factor<1.0e-2) and vec_error_gas(index-1)*vec_error_gas(index-2)<0.))
              non_convergence=true;

	  if (iiter>100)
	    if (abs(vec_error_aq(index-1))/factor>config.relative_precision)
	      {
		double a=1;
		double b=1;
		for (i=index-100;i<index-90;i++)
		  a*=abs(vec_error_aq(i));
              
                
		for (i=index-11;i<index-1;i++)
		  b*=abs(vec_error_aq(i));
                  
		if (b>a*0.999)
		  non_convergence=true;
	      }

	  if (iiter>100)
	    if (abs(vec_error_org(index-1))/factor>config.relative_precision)
	      {
		double a=1;
		double b=1;
		for (i=index-100;i<index-90;i++)
		  a*=abs(vec_error_org(i));
                
		for (i=index-11;i<index-1;i++)
		  b*=abs(vec_error_org(i));
                  
		if (b>a*0.999)
		  non_convergence=true;
	      }
	  
	  if (iiter>100)
	    if (vec_error_gas(index-1)/factor>config.relative_precision)
	      {
		double a=1;
		double b=1;
		for (i=index-100;i<index-90;i++)
		  a*=abs(vec_error_gas(i));
                
		for (i=index-11;i<index-1;i++)
		  b*=abs(vec_error_gas(i));
                  
		if (b>a*0.999)
		  non_convergence=true;
	      }

	  if (non_convergence and nh<config.nh_max)
	    { 
	      if (abs(vec_error_org(index-1))/factor>100.0 or abs(vec_error_aq(index-1))/factor>100.0 or abs(vec_error_gas(index-1))/factor>100.0)
		{
		  i=config.iH2O;
		  surrogate[i].Atot=surrogate[i].Atot1;
		  surrogate[i].Ag=surrogate[i].Ag1;  
		  if (surrogate[i].hydrophobic)
		    surrogate[i].Ap_layer_init=surrogate[i].Ap_layer;        
			 
		  if (surrogate[i].hydrophilic)
		    for(b=0;b<config.nbins;b++)
		      surrogate[i].Aaq_bins_init(b)=surrogate[i].Aaq_bins(b);

		  if (config.solids and surrogate[i].is_solid)
		    for(b=0;b<config.nbins;b++)
		      surrogate[i].Asol_bins_init(b)=surrogate[i].Asol_bins(b);
		}
	      ++nh;
	      factor=pow(0.5,nh-1);
	      for (i=max(index-50,0);i<index;i++)
		{
		  vec_error_org(i)=-1.0;
		  vec_error_aq(i)=-1.0;
		  vec_error_gas(i)=-1.0;
		  vec_error_chp(i)=-1.0;
		}
	      iiter=0;

	      for(b=0;b<config.nbins;b++)
		{
		  AQinit(b)=AQsave(b);
		  for (ilayer=0;ilayer<config.nlayer;++ilayer)
		    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		      MOinit(b,ilayer,iphase)=MOsave(b,ilayer,iphase); 
		}
	    }
	}

      //i=config.iH2O;
      for (i=0;i<n;i++)
        {
          for (b=0;b<config.nbins;++b)
            for (ilayer=0;ilayer<config.nlayer;++ilayer)
              for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                surrogate[i].Ap_layer(b,ilayer,iphase)=
                  surrogate[i].Ap_layer_init(b,ilayer,iphase);
          
          for (b=0;b<config.nbins;++b)
            surrogate[i].Aaq_bins(b)=surrogate[i].Aaq_bins_init(b);

	  for (b=0;b<config.nbins;++b)
            surrogate[i].Asol_bins(b)=surrogate[i].Asol_bins_init(b);
          
          surrogate[i].Ag1=surrogate[i].Ag;
          surrogate[i].Atot1=surrogate[i].Atot;
        }

      iiter++;   
      //if (LWCtot>config.LWClimit)
      //  density_aqueous_phase_ssh(config, surrogate, LWC, Temperature);
          
      for (b=0;b<config.nbins;b++)
	chp_save(b)=chp(b);
      
      if (config.first_evaluation_activity_coefficients==false)
	{	  
	  twostep_tot_ssh(config, surrogate, config.tequilibrium, AQinit, AQ, conc_inorganic,
			  ionic, ionic_organic, organion, chp, LWC, MOinit, MO,
			  MOW, Temperature, RH, MMaq, true, factor, t, deltat, index);       
	  if (config.compute_saturation and config.compute_organic)
	    phase_repartition_ssh(config,surrogate,Temperature,MOinit,MO,MOW,factor);	  
	}
      else
	{
	  twostep_tot_ssh(config, surrogate, config.tequilibrium, AQinit, AQ, conc_inorganic,
			  ionic, ionic_organic, organion, chp, LWC, MOinit, MO,
			  MOW, Temperature, RH, MMaq, false, factor, t, deltat, index);
	  if (config.compute_saturation and config.compute_organic)
	    phase_repartition_ssh(config,surrogate,Temperature,MOinit,MO,MOW,factor);
	}

      //redistribute concentrations to ensure that the volume of layers are constant      
      redistribution_ssh(config, surrogate,MOinit,MO);

      water_concentration_ssh(config, surrogate, Temperature, RH);
      
      tau_dif_ssh(config, surrogate, number, Vsol, Temperature);
      tau_kmt_ssh(config, surrogate, Temperature, number);
      

      //compute the new diameters of particle due to the growth of particles by condensation
      compute_diameters_ssh(config, surrogate, Vsol, number, LWC, LWCtot);
      if (config.explicit_representation)
	compute_morphology_ssh(config, Vsol, number);
		  
      //Computation of error_tot
      vec_error_org(index)=0.0;
      vec_error_aq(index)=0.0;
      vec_error_gas(index)=0.0;
      vec_error_chp(index)=0.0;
      i=config.iH2O;
      if (surrogate[i].is_organic or (surrogate[i].is_inorganic_precursor and surrogate[i].is_solid==false))
	if (surrogate[i].Atot1>1.0e-10)
          {
            double errloc=(surrogate[i].Atot-surrogate[i].Atot1)/surrogate[i].Atot1;
            if (abs(errloc)>abs(vec_error_gas(index)))
              vec_error_gas(index)=errloc;	
          }


      for (i=0;i<n;i++)
	if (surrogate[i].is_ion)
	  {
	    double btot=sum(surrogate[i].Aaq_bins_init);
	    double btot1=sum(surrogate[i].Aaq_bins);
	    if (i==config.iNO3m)
	      {
		btot+=surrogate[config.iHNO3].Ag/surrogate[config.iHNO3].MM*surrogate[config.iNO3m].MM;
		btot1+=surrogate[config.iHNO3].Ag1/surrogate[config.iHNO3].MM*surrogate[config.iNO3m].MM;
	      }
	    else if (i==config.iClm)
	      {
		btot+=surrogate[config.iHCl].Ag/surrogate[config.iHCl].MM*surrogate[i].MM;
		btot1+=surrogate[config.iHCl].Ag1/surrogate[config.iHCl].MM*surrogate[i].MM;
	      }
	    else if (i==config.iNH4p)
	      {
		btot+=surrogate[config.iNH3].Ag/surrogate[config.iNH3].MM*surrogate[i].MM;
		btot1+=surrogate[config.iNH3].Ag1/surrogate[config.iNH3].MM*surrogate[i].MM;
	      }
	    else if (i==config.iSO4mm)
	      {
		btot+=sum(surrogate[config.iHSO4m].Aaq_bins_init)/surrogate[config.iHSO4m].MM*surrogate[i].MM;
		btot1+=sum(surrogate[config.iHSO4m].Aaq_bins)/surrogate[config.iHSO4m].MM*surrogate[i].MM;		
	      }
	    else
	      btot1=0.;

	    
	    if (btot1>1.0e-10)
              {
		//Add strict convergence criteria on the total mass of inorganics when chemistry is accounted for to avoid numerical creation of loss of matter
		// max(60./deltat,1.0) create a stricter parameter for lower time step
                double errloc=(btot-btot1)/(1.0e-5)*max(60./deltat,1.0);
                if (abs(errloc)>abs(vec_error_gas(index)))
                  vec_error_gas(index)=errloc;
              }
	  }
	

      
      maxaq=0.;
      for (b=0;b<config.nbins;++b)
	{
	  AQ(b)=max(AQ(b),config.MOmin);
	  for (ilayer=0;ilayer<config.nlayer;++ilayer)
	    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
	      {
                MOinit(b,ilayer,iphase)=max(MO(b,ilayer,iphase),config.MOmin*config.Vlayer(ilayer));
		if (MO(b,ilayer,iphase)>config.MOmin)
                  {
                    double errloc=(MO(b,ilayer,iphase)-MOinit(b,ilayer,iphase))/MO(b,ilayer,iphase);
                    if (abs(errloc)>abs(vec_error_org(index)))
                      vec_error_org(index)=errloc;
                  }
		if (surrogate[i].hydrophobic and
		    surrogate[i].Ap_layer(b,ilayer,iphase)>0 and (surrogate[i].Ap_layer(b,ilayer,iphase)>config.MOmin or surrogate[i].Ap_layer_init(b,ilayer,iphase)>config.MOmin))
		  {
                    double errloc=(surrogate[i].Ap_layer_init(b,ilayer,iphase)-surrogate[i].Ap_layer(b,ilayer,iphase))
                      /surrogate[i].Ap_layer(b,ilayer,iphase);
                    if (abs(errloc)>abs(vec_error_org(index)))
                      vec_error_org(index)=errloc;
		  }
	      }

	  maxaq=max(maxaq,abs(AQinit(b)-AQ(b)));
          if (abs(AQinit(b)-AQ(b))>config.precision)
            {
              if (AQ(b)>1.0e-5)
                {
                  double errloc=(AQ(b)-AQinit(b))/AQ(b);
                  if (abs(errloc)>abs(vec_error_aq(index)))
                    {
                      vec_error_aq(index)=errloc;
                      vec_error_chp(index)=surrogate[i].Aaq_bins_init(b);
                    }
                  /*
                    if (chp(b)>0.)
                    {
                    vec_error_chp(index)=max(vec_error_chp(index),abs(chp(b)-chp_save(b))/chp(b));
                    }*/
                }
       
              if (surrogate[i].hydrophilic and i!=config.iHp and i!=config.iOHm)
                if (surrogate[i].Aaq_bins(b)>0 and (surrogate[i].Aaq_bins(b)>config.MOmin or surrogate[i].Aaq_bins_init(b)>config.MOmin))
                  {
                    double errloc=(surrogate[i].Aaq_bins_init(b)-surrogate[i].Aaq_bins(b))
                      /surrogate[i].Aaq_bins(b);
                    if (abs(errloc)>abs(vec_error_aq(index)))
                      {
                        vec_error_aq(index)=errloc;
                        vec_error_chp(index)=surrogate[i].Aaq_bins_init(b);
                      }
                  }
            }
			  
	  AQinit(b)=max(AQ(b),config.MOmin);
	  
	}

      if (iiter>200 and factor>factor_min)
        {
          factor=max(factor/2,factor_min);
          iiter=0;
        }
      
      error_tot=0.;
      for (b=0;b<config.nbins;b++)
        if (surrogate[config.iH2O].Aaq_bins_init(b)>0.)
          error_tot=max(error_tot,abs(surrogate[config.iH2O].Aaq_bins_init(b)-surrogate[config.iH2O].Aaq_bins(b))/surrogate[config.iH2O].Aaq_bins_init(b));

      //cout << surrogate[config.iH2O].Aaq_bins_init(0)/AQinit(0) << " " << surrogate[config.iH2O].Atot << " " << surrogate[config.iH2O].Ag << endl;
      //else
      //error_tot=0;
      for (b=0;b<config.nbins;b++)
        if (surrogate[config.iH2O].Ap_layer_init(b)>0.)
          error_tot=max(abs(surrogate[config.iH2O].Ap_layer_init(b)-surrogate[config.iH2O].Ap_layer(b))/surrogate[config.iH2O].Ap_layer_init(b),error_tot);

      if (vec_error_org(index)>0.)
        vec_error_org(index)=max(vec_error_org(index),0.001*config.relative_precision*factor);
      else
        vec_error_org(index)=min(vec_error_org(index),-0.001*config.relative_precision*factor);

      if (vec_error_aq(index)>0.)
        vec_error_aq(index)=max(vec_error_aq(index),0.001*config.relative_precision*factor);
      else
        vec_error_aq(index)=min(vec_error_aq(index),-0.001*config.relative_precision*factor);

      if (vec_error_gas(index)>0.)
        vec_error_gas(index)=max(vec_error_gas(index),0.001*config.relative_precision*factor);
      else
        vec_error_gas(index)=min(vec_error_gas(index),-0.001*config.relative_precision*factor);

      //cout << surrogate[config.iH2O].Ap_layer_init << " " << surrogate[config.iH2O].Aaq_bins_init << endl;
      
      //error_tot=vec_error_gas(index);
      //error_tot=max(error_tot,vec_error_gas(index));
      error_tot=max(max(abs(vec_error_org(index)),abs(vec_error_aq(index))),abs(vec_error_gas(index)));
      /*cout << index << " " << error_tot/factor << endl;
        if ((error_tot>config.relative_precision*factor or maxaq>0.01) and index < config.max_iter)
        cout << "true1" << endl;
        if (error_tot>config.relative_precision*factor or maxaq>0.01)
        cout << "true2" << endl;
        if (error_tot>config.relative_precision*factor)
        cout << "true3" << endl;
        cout << AQinit << endl;
        cout << AQ << endl;
        cout << maxaq << endl;
        cout << surrogate[config.iH2O].gamma_aq_bins << endl;*/
      ++index;  
    }

  config.compute_inorganic=compute_inorganic_save;
  config.compute_organic=compute_organic_save;

  if (index==config.max_iter)
    {
      cout << "water " << index << " " << factor << " " << error_tot/factor << " " << RH << " "<< endl;
      cout << vec_error_org(index-1)/factor << " " << vec_error_aq(index-1)/factor << " " << vec_error_gas(index-1)/factor << endl;
      cout << surrogate[config.iH2O].Ap_layer_init << endl;
      cout << surrogate[config.iH2O].Ap_layer << endl;
      
      //for(index=0;index<config.max_iter;index++)
      //  cout << index << " " << vec_error_aq(index) << " " << vec_error_chp(index) << endl;

      for (i=0;i<n;i++)
        for (b=0;b<config.nbins;b++)
          cout << b << " " << surrogate[i].name << " " << surrogate[i].Aaq_bins_init(b) << " " << surrogate[i].Aaq_bins(b) << endl;
      
      error_tot=0.;
      for (b=0;b<config.nbins;b++)
        if (sum(surrogate[config.iH2O].Aaq_bins_init)>0.)
          error_tot=abs(surrogate[config.iH2O].Aaq_bins_init(b)-surrogate[config.iH2O].Aaq_bins(b))/surrogate[config.iH2O].Aaq_bins_init(b);
      cout << error_tot << endl;
      cout << "AQ: " << AQinit << endl;
      cout << "MO: " << MOinit << endl;
      cout << "chp: " << chp << endl;
      cout << "chp_save: " << chp_save << endl;
      throw string("stop");
    }

}



void solve_implicit_coupled_ssh(model_config config, vector<species> &surrogate,
                                Array<double, 3> &MOinit, Array<double, 3> &MOW, Array<double, 1> &number,
                                Array<double, 1> &Vsol,
                                Array<double, 1> &LWC, Array<double, 1> &AQinit, Array<double, 1> &ionic,
                                Array<double, 1> &chp,
                                double &Temperature, double &RH,
                                Array<double, 1> &AQ, Array<double, 3> &MO,
                                Array<double, 1> &conc_inorganic, Array<double, 1> &ionic_organic, Array<double, 1> &organion, Array<double, 1> &MMaq,
                                double t, double deltat, int &index, bool &reject_step, Array<double, 1> &chp20)
{
  int b,ilayer,iphase,i;
  int n=surrogate.size();
  double error_tot=10.0;
  int index_max=index;
  index=0;
  Array<double, 1> vec_error_org,vec_error_aq,vec_error_gas,vec_error_chp,vec_error_compo; 
  vec_error_org.resize(config.max_iter);
  vec_error_aq.resize(config.max_iter);
  vec_error_gas.resize(config.max_iter);
  vec_error_chp.resize(config.max_iter);
  vec_error_compo.resize(config.max_iter);
  vec_error_org=0.;
  vec_error_aq=0.;
  vec_error_gas=0.;
  vec_error_chp=0.;
  vec_error_compo=0.;
  bool non_convergence;
  double LWCtot=0.0;
  for (b=0;b<config.nbins;++b)
    LWCtot+=LWC(b);
  Array<double, 1> AQsave;
  Array<double, 1> chp_save,chp0;
  Array<double, 3> MOsave;
  AQsave.resize(config.nbins);
  chp_save.resize(config.nbins);
  MOsave.resize(config.nbins,config.nlayer,config.max_number_of_phases);
  chp0.resize(config.nbins);

  //The method of newton raphson cannot be used due to the high number of variables
  //(nbins*(1+nlayer*nphase)
  //
  //Method to reach equilibrium
  // - initialisation of MOinit and AQinit
  // - initialisation of MO and AQ
  // - compute error_tot = maximum of relative errors
  //   if error_tot < 0.1% the system has converged
  // - if error_tot = error_old2 (error from iteration - 2) the system cannot converge
  //   due to strong variation of the composition.
  //   in that case a method is used so that variations are smaller:
  //     Ap(iter+1)=Ap(iter)*(1-1.0/nh)+1.0/nh*Ap(computed)
  //     1.0/nh weight factor
  //     nh is increased to ensure convergence
  //     error_tot must be inferior to factor * 0.1% to reach equilibrium
  int iiter=0;
  int nh;
  int ntoo_big=0;

  if (config.compute_inorganic and config.compute_organic==false)
    nh=config.nh_inorg_init;
  else if (config.compute_inorganic and config.compute_organic)
    nh=max(config.nh_inorg_init,max(config.nh_aq_init,config.nh_org_init));
  else
    nh=max(config.nh_aq_init,config.nh_org_init);

  //nh=nh-2;

  double factor_max=1.;
  double factor_min=pow(0.5,config.nh_max-1);
  double var_min=0.08;
  double var_max=0.12;
  double var_rej=1.;
  double factor=pow(0.5,nh-1);
  double var;
  double factor_old=factor;
  int m;
  double significant=1.0e-5;
  double chp_max=10.;
  double chp_min=1.e-20;
  //chp_save=chp;
  chp0=chp;
  config.nh_max=7;
  //config.precision=1.e-8;
  for (i=0;i<n;i++)
    {
      surrogate[i].Atot0=surrogate[i].Atot;
      surrogate[i].Ag0=surrogate[i].Ag;
      for(b=0;b<config.nbins;b++)
        {
          if (surrogate[i].hydrophilic)
            surrogate[i].Aaq_bins_init0(b)=surrogate[i].Aaq_bins_init(b);
          if (surrogate[i].hydrophobic)
            for (ilayer=0;ilayer<config.nlayer;++ilayer)
              for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                surrogate[i].Ap_layer_init0(b,ilayer,iphase)=surrogate[i].Ap_layer_init(b,ilayer,iphase);
	  if (surrogate[i].is_solid)
	    surrogate[i].Asol_bins_init0(b)=surrogate[i].Asol_bins_init(b);
        }

         
      surrogate[i].Atot20=surrogate[i].Atot;
      surrogate[i].Ag20=surrogate[i].Ag;
      for(b=0;b<config.nbins;b++)
        {
          if (surrogate[i].hydrophilic)
            surrogate[i].Aaq_bins_init20(b)=surrogate[i].Aaq_bins_init(b);
          if (surrogate[i].hydrophobic)
            for (ilayer=0;ilayer<config.nlayer;++ilayer)
              for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                surrogate[i].Ap_layer_init20(b,ilayer,iphase)=surrogate[i].Ap_layer_init(b,ilayer,iphase);
	  if (surrogate[i].is_solid)
	    surrogate[i].Asol_bins_init20(b)=surrogate[i].Asol_bins_init(b);
        }
    }
  chp20=chp;

  

  for(b=0;b<config.nbins;b++)
    {
      AQsave(b)=AQinit(b);
      for (ilayer=0;ilayer<config.nlayer;++ilayer)
        for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
          MOsave(b,ilayer,iphase)=MOinit(b,ilayer,iphase); 
    }
  //cout << "init0 " << chp << " " << chp_save << endl;
  chp_save=chp;

  bool hygroscopicity_save=config.hygroscopicity;
  bool inorganic_save=config.compute_inorganic;

  if (config.isorropia_ph==false)
    for (b=0;b<config.nbins;b++)
      chp(b)=min(1.,max(chp(b),1.e-6));

  //cout << "index_max: " << index_max << endl;
  while ((error_tot>config.relative_precision and index < index_max) or index<2)
    {
      if (index>2)
        {
          //ensure that the system can converge
          non_convergence=false;              
          if (non_convergence and nh<config.nh_max)
            {
       
              if (abs(vec_error_org(index-1))>100.0 or abs(vec_error_aq(index-1))>100.0 or abs(vec_error_gas(index-1))>100.0)
                {
		  //cout << "rejected " << endl;
                  for (i=0;i<n;i++)
                    {
		      surrogate[i].Atot=surrogate[i].Atot1;
                      surrogate[i].Ag=surrogate[i].Ag1;
                      if (surrogate[i].hydrophobic)
                        for(b=0;b<config.nbins;b++)
                          for (ilayer=0;ilayer<config.nlayer;++ilayer)
                            for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                              surrogate[i].Ap_layer_init(b,ilayer,iphase)=surrogate[i].Ap_layer(b,ilayer,iphase);          
                      if (surrogate[i].hydrophilic)
                        for(b=0;b<config.nbins;b++)
                          surrogate[i].Aaq_bins_init(b)=surrogate[i].Aaq_bins(b);

		      if (surrogate[i].is_solid)
			for(b=0;b<config.nbins;b++)
                          surrogate[i].Asol_bins_init(b)=surrogate[i].Asol_bins(b);
		   
                    }
                  
                  for(b=0;b<config.nbins;b++)
                    {
                      AQinit(b)=AQsave(b);
                      for (ilayer=0;ilayer<config.nlayer;++ilayer)
                        for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                          MOinit(b,ilayer,iphase)=MOsave(b,ilayer,iphase);
                      chp(b)=chp_save(b);
                    }
                }
              //++nh;
	      //factor=pow(0.5,nh-1);
              factor_max=max(factor_max/2,factor_min);
              factor=max(factor/2,factor_min);
              for (i=max(index-50,0);i<index;i++)
                {
                  vec_error_org(i)=0.0;
                  vec_error_aq(i)=0.0;
		  vec_error_gas(i)=0.0;
		  vec_error_chp(i)=0.0;
                  vec_error_compo(i)=0.0;
                }
              iiter=0;

              for(b=0;b<config.nbins;b++)
                {
                  AQinit(b)=AQsave(b);
                  for (ilayer=0;ilayer<config.nlayer;++ilayer)
                    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                      MOinit(b,ilayer,iphase)=MOsave(b,ilayer,iphase); 
                }
            }
        }

      for (b=0;b<config.nbins;++b)
        for (ilayer=0;ilayer<config.nlayer;++ilayer)
          for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
            for (i=0;i<n;++i)
              {
                surrogate[i].Ap_layer(b,ilayer,iphase)=
                  surrogate[i].Ap_layer_init(b,ilayer,iphase);
                surrogate[i].gamma_org_layer0(b,ilayer,iphase)=surrogate[i].gamma_org_layer(b,ilayer,iphase);
              }

      for (b=0;b<config.nbins;++b)
        for (i=0;i<n;i++)
          {
            surrogate[i].Aaq_bins(b)=surrogate[i].Aaq_bins_init(b);
            surrogate[i].gamma_aq_bins_old(b)=surrogate[i].gamma_aq_bins(b);
	    surrogate[i].Asol_bins(b)=surrogate[i].Asol_bins_init(b);
          }

      for (i=0;i<n;i++)
	surrogate[i].Ag1=surrogate[i].Ag;

      for (i=0;i<n;i++)
	surrogate[i].Atot1=surrogate[i].Atot;

      //cout << "NO3: " << sum(surrogate[config.iNO3m].Aaq_bins_init)/surrogate[config.iNO3m].MM*surrogate[config.iHNO3].MM+surrogate[config.iHNO3].Ag << " " << nh << endl;;
      //cout << "NO3i: " << sum(surrogate[config.iNO3m].Aaq_bins_init0)/surrogate[config.iNO3m].MM*surrogate[config.iHNO3].MM << " " << surrogate[config.iHNO3].Ag0 << " " << nh << endl;
   
      if (LWCtot>config.LWClimit)
        density_aqueous_phase_ssh(config, surrogate, LWC, Temperature);
          
      for (b=0;b<config.nbins;b++)
	chp_save(b)=chp(b);
      factor_old=factor;
      
      for (b=0;b<config.nbins;b++)
        {
          AQsave(b)=AQinit(b);
          for (ilayer=0;ilayer<config.nlayer;++ilayer)
            for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
              MOsave(b,ilayer,iphase)=MOinit(b,ilayer,iphase);
        }
      //cout << " AVANT " << sum(surrogate[config.iNO3m].Aaq_bins_init)+surrogate[config.iHNO3].Ag/63*62 << " " << surrogate[config.iHNO3].Ag << " " << surrogate[config.iNO3m].Aaq_bins_init <<  endl;

      if (index<20 and index>15)
	{
	  chp20=chp;
	  for (i=0;i<n;i++)
	    {
	      surrogate[i].Atot20=surrogate[i].Atot;
	      surrogate[i].Ag20=surrogate[i].Ag;
	      for(b=0;b<config.nbins;b++)
		{
		  if (surrogate[i].hydrophilic)
		    surrogate[i].Aaq_bins_init20(b)=surrogate[i].Aaq_bins_init(b);
		  if (surrogate[i].hydrophobic)
		    for (ilayer=0;ilayer<config.nlayer;++ilayer)
		      for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
			surrogate[i].Ap_layer_init20(b,ilayer,iphase)=surrogate[i].Ap_layer_init(b,ilayer,iphase);
		  if (surrogate[i].is_solid)
		    surrogate[i].Asol_bins_init20(b)=surrogate[i].Asol_bins_init(b);
		}
	    }
	}

      if (config.first_evaluation_activity_coefficients==false)
        {	  
          twostep_tot_ssh(config, surrogate, config.tequilibrium, AQinit, AQ, conc_inorganic,
                          ionic, ionic_organic, organion, chp, LWC, MOinit, MO,
                          MOW, Temperature, RH, MMaq, true, factor, t, deltat, index);       
          if (config.compute_saturation and config.compute_organic)
            phase_repartition_ssh(config,surrogate,Temperature,MOinit,MO,MOW,factor);	  
        }
      else
        {
          twostep_tot_ssh(config, surrogate, config.tequilibrium, AQinit, AQ, conc_inorganic,
                          ionic, ionic_organic, organion, chp, LWC, MOinit, MO,
                          MOW, Temperature, RH, MMaq, false, factor, t, deltat, index);
          if (config.compute_saturation and config.compute_organic)
            phase_repartition_ssh(config,surrogate,Temperature,MOinit,MO,MOW,factor);	  
        }
      
      //cout << surrogate[config.iCa].Aaq_bins_init << endl;

      if (iiter>200)
        {
          factor_max=max(factor_max/2,factor_min);
          factor=max(factor/2,factor_min);
          iiter=0;
        }
      
      var=0.;
      m=0;
      double negligeable=config.MOmin;
      for (i=0;i<n;i++)
        if ((surrogate[i].is_ion and config.compute_inorganic) or i==config.iH2O or (surrogate[i].is_organic and config.compute_organic))
          if (surrogate[i].hydrophilic)
            //if (surrogate[i].hydrophilic and surrogate[i].is_inorganic_precursor==false)
            for (b=0;b<config.nbins;b++)
              //if (error_spec(i)/factor_old>relprec)
              {
                if (0.5*surrogate[i].gamma_aq_bins(b)+0.5*surrogate[i].gamma_aq_bins_old(b)>1.0e-6 and (surrogate[i].Aaq_bins_init(b)>negligeable or surrogate[i].Aaq_bins(b)>negligeable or i==config.iHp or i==config.iOHm))
                  {
                    m++;
                    var+=abs((surrogate[i].gamma_aq_bins_old(b)-surrogate[i].gamma_aq_bins(b))/(0.5*surrogate[i].gamma_aq_bins(b)+0.5*surrogate[i].gamma_aq_bins_old(b)));
                  }
              }

      if (config.compute_organic)
        for (i=0;i<n;i++)
          if (surrogate[i].is_organic or i==config.iH2O)
            if (surrogate[i].hydrophobic)
              for (b=0;b<config.nbins;b++)
                for (ilayer=0;ilayer<config.nlayer;++ilayer)
                  for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                    //if (error_spec(i)/factor_old>relprec)
                    {
                      if (0.5*surrogate[i].gamma_org_layer(b,ilayer,iphase)+0.5*surrogate[i].gamma_org_layer0(b,ilayer,iphase)>1.0e-6 and (surrogate[i].Ap_layer_init(b,ilayer,iphase)>negligeable or surrogate[i].Ap_layer(b,ilayer,iphase)>negligeable))
                        {
                          m++;
                          var+=abs((surrogate[i].gamma_org_layer0(b,ilayer,iphase)-surrogate[i].gamma_org_layer(b,ilayer,iphase))/(0.5*surrogate[i].gamma_org_layer(b,ilayer,iphase)+0.5*surrogate[i].gamma_org_layer0(b,ilayer,iphase)));
                        }
                    }

      if (m>0)
        var=var/m;      

      if (var>var_max and index>0 and factor>factor_min)
        {                      
          factor=max(factor_old*max(var_max/var,0.1),factor_min);
          /*
            if (var>var_rej)
            for (i=0;i<n;i++)
            {
            surrogate[i].Atot=surrogate[i].Atot1;
            surrogate[i].Ag=surrogate[i].Ag1;
            if (surrogate[i].hydrophobic)
            for(b=0;b<config.nbins;b++)
            for (ilayer=0;ilayer<config.nlayer;++ilayer)
            for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
            surrogate[i].Ap_layer_init(b,ilayer,iphase)=surrogate[i].Ap_layer(b,ilayer,iphase);          
            if (surrogate[i].hydrophilic)
            for(b=0;b<config.nbins;b++)
            surrogate[i].Aaq_bins_init(b)=surrogate[i].Aaq_bins(b);             
            }*/
        }
      else if (var<var_min and index>0 and factor<1.)
        {
          if (var>0.)
            factor=min(factor_old*min(var_max/var,10.),factor_max);
          else
            factor=min(factor_old*10.,factor_max);
        }

      factor=max(min(factor,factor_max),factor_min);

      //redistribute concentrations to ensure that the volume of layers are constant      
      redistribution_ssh(config, surrogate,MOinit,MO);

      water_concentration_ssh(config, surrogate, Temperature, RH);
      
      tau_dif_ssh(config, surrogate, number, Vsol,Temperature);
      tau_kmt_ssh(config, surrogate, Temperature, number);
      characteristic_time_ssh(config, surrogate, MOinit, AQinit, LWCtot); 
      if (LWCtot>config.LWClimit)
	characteristic_time_aq_ssh(config, surrogate, Temperature, chp, LWC, AQinit, MOinit, MMaq, ionic);

      //compute the new diameters of particle due to the growth of particles by condensation
      compute_diameters_ssh(config, surrogate, Vsol, number, LWC, LWCtot);
      if (config.explicit_representation)
	compute_morphology_ssh(config, Vsol, number);
		  
      //Computation of error_tot
      vec_error_org(index)=0.0;
      vec_error_aq(index)=0.0;
      vec_error_gas(index)=0.0;
      vec_error_chp(index)=0.0;
      vec_error_compo(index)=0.0;
      double errloc;
      for (i=0;i<n;i++)
	if (surrogate[i].is_organic or (surrogate[i].is_inorganic_precursor and surrogate[i].is_solid==false))
          {
            if (surrogate[i].Atot1>1.0e-10)
              {
                errloc=(surrogate[i].Atot-surrogate[i].Atot1)/surrogate[i].Atot1/factor_old;
                if (abs(errloc)>abs(vec_error_gas(index)))
                  vec_error_gas(index)=errloc;
              }
            if (surrogate[i].Ag-surrogate[i].Ag1>config.precision and surrogate[i].Ag1>1.e-3)
              {
                errloc=(surrogate[i].Ag-surrogate[i].Ag1)/surrogate[i].Ag1/factor_old;
                if (abs(errloc)>abs(vec_error_gas(index)))
                  vec_error_gas(index)=errloc;
              }
          }
      
      for (i=0;i<n;i++)
	if (surrogate[i].is_ion)
	  {
	    double btot=sum(surrogate[i].Aaq_bins_init);
	    double btot1=sum(surrogate[i].Aaq_bins);
	    if (i==config.iNO3m)
	      {
		btot+=surrogate[config.iHNO3].Ag/surrogate[config.iHNO3].MM*surrogate[config.iNO3m].MM;
		btot1+=surrogate[config.iHNO3].Ag1/surrogate[config.iHNO3].MM*surrogate[config.iNO3m].MM;
	      }
	    else if (i==config.iClm)
	      {
		btot+=surrogate[config.iHCl].Ag/surrogate[config.iHCl].MM*surrogate[i].MM;
		btot1+=surrogate[config.iHCl].Ag1/surrogate[config.iHCl].MM*surrogate[i].MM;
	      }
	    else if (i==config.iNH4p)
	      {
		btot+=surrogate[config.iNH3].Ag/surrogate[config.iNH3].MM*surrogate[i].MM;
		btot1+=surrogate[config.iNH3].Ag1/surrogate[config.iNH3].MM*surrogate[i].MM;
	      }
	    else if (i==config.iSO4mm)
	      {
		btot+=sum(surrogate[config.iHSO4m].Aaq_bins_init)/surrogate[config.iHSO4m].MM*surrogate[i].MM;
		btot1+=sum(surrogate[config.iHSO4m].Aaq_bins)/surrogate[config.iHSO4m].MM*surrogate[i].MM;		
	      }
	    else
	      btot1=0.;
	      
	    
	    if (btot1>1.0e-10)
              {
		//Add strict convergence criteria on the total mass of inorganics when chemistry is accounted for to avoid numerical creation of loss of matter
		// max(60./deltat,1.0) create a stricter parameter for lower time step
                errloc=(btot-btot1)/(1.0e-5)/factor_old*max(60./deltat,1.0);
                if (abs(errloc)>abs(vec_error_gas(index)))
                  vec_error_gas(index)=errloc;
              }
	  }
	

	
      for (b=0;b<config.nbins;++b)
        {
          AQ(b)=max(AQ(b),config.MOmin);         
          for (ilayer=0;ilayer<config.nlayer;++ilayer)
            for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
              //if (abs(MO(b,ilayer,iphase)-MOinit(b,ilayer,iphase))>config.precision)
              {
                if (MO(b,ilayer,iphase)>1.0e-5*config.Vlayer(ilayer))
                  {
                    errloc=(MO(b,ilayer,iphase)-MOinit(b,ilayer,iphase))/MO(b,ilayer,iphase)/factor_old;
                    if (abs(errloc)>abs(vec_error_org(index)))
                      vec_error_org(index)=errloc;
                  }
		    
                //if (var<=var_rej)
                MOinit(b,ilayer,iphase)=max(MO(b,ilayer,iphase),config.MOmin*config.Vlayer(ilayer));
              }


          if (config.compute_inorganic)
            //if ((surrogate[config.iHp].gamma_aq_bins(b)*chp(b)<1.e-7 or surrogate[config.iHp].gamma_aq_bins_old(b)*chp_save(b)<1.e-7 or
	    //    surrogate[config.iHp].gamma_aq_bins_old(b)*chp(b)>0.1 or surrogate[config.iHp].gamma_aq_bins_old(b)*chp_save(b)>0.1) and AQ(b)>negligeable)
              //if (AQ(b)>1.e-5)
              {
                errloc=(min(max(chp(b),chp_min),chp_max)-min(max(chp_save(b),chp_min),chp_max))/min(max(chp_min,chp(b)),chp_max)/factor_old;
		//if (abs(errloc)>0.3)
		//  cout << "errloc1 " << errloc << " " << AQ(b) << endl;
                if (abs(errloc)>abs(vec_error_chp(index)))
                  vec_error_chp(index)=errloc;
              }
          
          //if (abs(AQinit(b)-AQ(b))>config.precision)
	  if (AQ(b)>1.e-5) //config.MOmin)
	      {
		errloc=(AQ(b)-AQinit(b))/AQ(b)/factor_old;
		if (abs(errloc)>abs(vec_error_aq(index)))
		  vec_error_aq(index)=errloc;
		
		if (chp(b)>0. and config.compute_inorganic)
		  {
		    errloc=(min(max(chp(b),chp_min),chp_max)-min(max(chp_save(b),chp_min),chp_max))/min(max(chp_min,chp(b)),chp_max)/factor_old;
		    //if (abs(errloc)>0.3)
		    //cout << "errloc2 " << errloc << endl;
		    if (abs(errloc)>abs(vec_error_chp(index)))
		      vec_error_chp(index)=errloc;
		  }
		//vec_error_aq(index)=max(vec_error_aq(index),vec_error_chp(index));
	      }
       
	  for (i=0;i<n;i++)
	    if (surrogate[i].hydrophilic and i!=config.iHp and i!=config.iOHm and i!=config.iHSO4m and i!=config.iSO4mm) // and i!=config.iCa and i!=config.iCO3mm and i!=config.iHCO3m)
	      if (surrogate[i].Aaq_bins(b)>0 and (surrogate[i].Aaq_bins(b)>significant or surrogate[i].Aaq_bins_init(b)>significant))// and abs(surrogate[i].Aaq_bins_init(b)-surrogate[i].Aaq_bins(b))>config.precision)
		{
		  /*if (abs(surrogate[i].Aaq_bins_init(b)-surrogate[i].Aaq_bins(b))/surrogate[i].Aaq_bins(b)>vec_error_aq(index))
		    cout << surrogate[i].name << " " << abs(surrogate[i].Aaq_bins_init(b)-surrogate[i].Aaq_bins(b))/surrogate[i].Aaq_bins(b) << " " << surrogate[i].gamma_aq_bins(b) << b << " " << surrogate[i].Aaq_bins_init(b) << " " << surrogate[i].Aaq_bins(b) << " " << b << endl;*/
		  errloc=(surrogate[i].Aaq_bins_init(b)-surrogate[i].Aaq_bins(b))/surrogate[i].Aaq_bins(b)/factor_old;
		  /*
		    if (abs(errloc)>1.e20)
		    {
		    cout << "errloc: " << surrogate[i].name << " " << errloc << " " << surrogate[i].Aaq_bins_init(b) << " " << surrogate[i].Aaq_bins(b) << " " << factor << endl;
		    //exit(1);
		    }*/
		  
		}
	  // }
	      
	  for (i=0;i<n;i++)
	    if (surrogate[i].hydrophobic)
	      for (ilayer=0;ilayer<config.nlayer;++ilayer)
		for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		  if (surrogate[i].Ap_layer(b,ilayer,iphase)>config.Vlayer(ilayer)*significant)
		    {
		      errloc=(surrogate[i].Ap_layer_init(b,ilayer,iphase)-surrogate[i].Ap_layer(b,ilayer,iphase))/surrogate[i].Ap_layer(b,ilayer,iphase)/factor_old;
		      if (abs(errloc)>abs(vec_error_compo(index)))
			vec_error_compo(index)=errloc;
		    }

	  for (i=0;i<n;i++)
	    if (surrogate[i].is_solid)
	      if (surrogate[i].Asol_bins(b)>negligeable) //config.MOmin)
		if (abs(surrogate[i].Asol_bins_init(b)-surrogate[i].Asol_bins(b))>negligeable)
		{
		  errloc=(surrogate[i].Asol_bins_init(b)-surrogate[i].Asol_bins(b))/surrogate[i].Asol_bins(b)/factor_old;
		  if (index>config.max_iter-5)
		    cout << index << " errloc: " << surrogate[i].name << " " << errloc << " " << surrogate[i].Asol_bins_init(b) << " " << surrogate[i].Asol_bins(b) << endl;
		  if (abs(errloc)>abs(vec_error_compo(index)))
		    vec_error_compo(index)=errloc;
		}
          //if (AQ(b)>1.0e-5)
          //  vec_error_chp(index)=max(vec_error_chp(index),abs(chp(b)-chp_save(b))/chp(b));
          //vec_error_aq(index)=max(vec_error_aq(index),vec_error_chp(index));
          //if (var<=var_rej)
          AQinit(b)=max(AQ(b),config.MOmin);
        }

      if (iiter>2)
        {             
          if (vec_error_org(index)*vec_error_org(index-1)<0.0 and abs(vec_error_aq(index))<config.relative_precision and abs(vec_error_org(index)-vec_error_org(index-2))<0.01*abs(vec_error_org(index)))
            ntoo_big=ntoo_big+10;
          else if (vec_error_aq(index)*vec_error_aq(index-1)<0.0 and abs(vec_error_org(index))<config.relative_precision and abs(vec_error_aq(index)-vec_error_aq(index-2))<0.01*abs(vec_error_aq(index)))                        
            ntoo_big=ntoo_big+10;
          else if (vec_error_org(index)*vec_error_org(index-1)<0.0 and abs(vec_error_aq(index))>config.relative_precision and vec_error_aq(index)*vec_error_aq(index-1)<0.0 and abs(vec_error_org(index))>config.relative_precision and abs(vec_error_org(index)-vec_error_org(index-2))<0.01*abs(vec_error_org(index)) and abs(vec_error_aq(index)-vec_error_aq(index-2))<0.01*abs(vec_error_aq(index))) 
            ntoo_big=ntoo_big+10;
          else if (vec_error_chp(index)*vec_error_chp(index-1)<0.0 and abs(vec_error_chp(index)-vec_error_chp(index-2))<0.01*abs(vec_error_chp(index))) // and abs(vec_error_chp(index))>config.relative_precision) 
            ntoo_big=ntoo_big+10;
          
          
          int a=0;
          //cout << error3/factor_old << " " << error3_old << endl;
          if (config.compute_inorganic)
            if (vec_error_chp(index)*vec_error_chp(index-1)<0.0 and abs(vec_error_chp(index))>config.relative_precision) 
              {
                //if (max(abs(error3),abs(error3_old))/(abs(error3)+abs(error3_old))<0.6)
                
                factor=min(factor,factor_old*(max(abs(vec_error_chp(index)),abs(vec_error_chp(index-1)))/(abs(vec_error_chp(index)+vec_error_chp(index))))); 
                ntoo_big=ntoo_big+1;
              }                     
          
          /*
            if (error4/factor_old>relprec)
            for (i=0;i<n;i++)
            if (error_spec(i)*error_spec_old(i)<0.0 and abs(error_spec(i))/factor_old>relprec)     
            {                       
            //factor=min(factor,0.5*factor_old); //
            //factor=min(factor,factor_old*(max(abs(error_spec(i)),abs(error_spec_old(i)))/(abs(error_spec(i))+abs(error_spec_old(i)))));                              
            //ntoo_big++;
            a=a+1;
            }*/

          //factpr=min(factor,factor_old*(ionic

          /*
            if (error4>var_rej)
            {
            //factor=min(factor,factor_old*var_rej/error4);
            //a=a+1;
            ntoo_big++; //(error4/factor_old/var_rej);
            }*/               
          //if (a>0)
          //  ntoo_big++;
          factor=max(factor,factor_min);
        }                  

          
      if (ntoo_big>200 and factor_max>factor_min)
        {
          factor_max=max(factor_max/2,factor_min);
          
          ntoo_big=0;                  
        }

      /*
        if (vec_error_compo(index)>1e5 and iiter>5)
        {
        factor_max=max(factor_max/2,factor_min);
        factor=max(min(factor,factor_old/2),factor_min);
        }*/
                  

      //vec_error_org(index)=max(vec_error_org(index),0.001*config.relative_precision);
      //vec_error_aq(index)=max(vec_error_aq(index),0.001*config.relative_precision);
      //vec_error_gas(index)=max(vec_error_gas(index),0.001*config.relative_precision);
      /*
        if (vec_error_org(index)>0.)
        vec_error_org(index)=max(vec_error_org(index),0.001*config.relative_precision*factor);
        else
        vec_error_org(index)=min(vec_error_org(index),-0.001*config.relative_precision*factor);

        if (vec_error_aq(index)>0.)
        vec_error_aq(index)=max(vec_error_aq(index),0.001*config.relative_precision*factor);
        else
        vec_error_aq(index)=min(vec_error_aq(index),-0.001*config.relative_precision*factor);

        if (vec_error_gas(index)>0.)
        vec_error_gas(index)=max(vec_error_gas(index),0.001*config.relative_precision*factor);
        else
        vec_error_gas(index)=min(vec_error_gas(index),-0.001*config.relative_precision*factor);*/
   
      //if (config.imethod==3)
      //	error_tot=vec_error_gas(index);
      //else
      
      error_tot=max(max(max(abs(vec_error_org(index)),abs(vec_error_aq(index))),abs(vec_error_gas(index))),abs(vec_error_compo(index)));
      //for (b=0;b<config.nbins;++b)
      // if (AQ(b)>1.e-5)

      if (config.compute_inorganic and error_tot>0.)
        error_tot=max(error_tot,abs(vec_error_chp(index)));
      //cout << index << " " << vec_error_org(index) << " " << abs(vec_error_aq(index)) << " " << abs(vec_error_gas(index)) << " " << abs(vec_error_compo(index)) << endl;
      //cout << factor_max << " " << factor << " " << chp << endl;
      /*
      if (int(index/100)*100==index and index>0)
	{
	  //     if (config.compute_inorganic)
	  for (b=0;b<config.nbins;++b)	  
	    if ((surrogate[config.iHp].gamma_aq_bins(b)*chp(b)<1.e-7 or surrogate[config.iHp].gamma_aq_bins_old(b)*chp_save(b)<1.e-7 or
		 surrogate[config.iHp].gamma_aq_bins_old(b)*chp(b)>0.1 or surrogate[config.iHp].gamma_aq_bins_old(b)*chp_save(b)>0.1) and AQ(b)>1.0e-5)
              {
                errloc=(min(chp(b),chp_max)-min(chp_save(b),chp_max))/min(chp(b),chp_max)/factor_old;
                if (abs(errloc)>0.9*abs(error_tot))
		  cout << "pH error " << b << " " << AQinit(b) << " " << chp(b) << " " << chp_save(b) << " " << errloc << endl;
		    //vec_error_chp(index)=errloc;
              }
	  cout << index << " " << error_tot << " " << deltat << " " << vec_error_chp(index) << " " << vec_error_org(index) << " " << vec_error_aq(index) << " " << vec_error_gas(index) << " " << vec_error_compo(index) << endl;
	}*/
      ++index;
      ++iiter;


      //exit(0);
      /*
      for (i=0;i<n;i++)
	if (surrogate[i].is_solid)
	  {
	    cout << surrogate[i].name << " " << surrogate[i].Asol_bins_init << endl;
	    for (b=0;b<config.nbins;b++)
	      cout << "error : " << b << " " << (surrogate[i].Asol_bins_init(b)-surrogate[i].Asol_bins(b))/factor << endl;
	      }*/

    }
  //cout << index << " " << deltat/index <<endl;
  
  //cout << "final " << index << " " << error_tot << " " << deltat << " " << vec_error_chp(index) << " " << vec_error_org(index) << " " << vec_error_aq(index) << " " << vec_error_gas(index) << " " << vec_error_compo(index) << endl;
  if (error_tot>config.relative_precision and index_max<config.max_iter)
    {
      reject_step=true;
      for (i=0;i<n;i++)
	{
	  surrogate[i].Atot=surrogate[i].Atot0;
	  surrogate[i].Ag=surrogate[i].Ag0;
	  if (surrogate[i].hydrophobic)
	    for(b=0;b<config.nbins;b++)
	      for (ilayer=0;ilayer<config.nlayer;++ilayer)
		for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		  surrogate[i].Ap_layer_init(b,ilayer,iphase)=surrogate[i].Ap_layer_init0(b,ilayer,iphase);          
	  if (surrogate[i].hydrophilic)
	    for(b=0;b<config.nbins;b++)
	      surrogate[i].Aaq_bins_init(b)=surrogate[i].Aaq_bins_init0(b);

	  if (surrogate[i].is_solid)
	    for(b=0;b<config.nbins;b++)
	      surrogate[i].Asol_bins_init(b)=surrogate[i].Asol_bins_init0(b);		   
	}
      chp0=chp;
    }
  else
    reject_step=false;
  
  if (error_tot>config.relative_precision and index_max==config.max_iter)
    {
      cout << "The model did not converged... " << RH << " " << Temperature << " " << factor << " " << config.MOmin << endl;
      cout << "error tot" << error_tot << endl;
      if (config.compute_inorganic and config.compute_organic)
        cout << "coupled-both " << vec_error_chp(index-1) << " " << vec_error_aq(index-1) << " " << vec_error_org(index-1) << " " << vec_error_gas(index-1) << " " << vec_error_compo(index-1) << endl;
      else if (config.compute_inorganic)
        cout << "coupled-inorg " << vec_error_chp(index-1) << " " << vec_error_aq(index-1) << " " << vec_error_org(index-1) << " " << vec_error_gas(index-1) << " " << vec_error_compo(index-1) << endl;
      else if (config.compute_organic)
        cout << "coupled-org " << vec_error_chp(index-1) << " " << vec_error_aq(index-1) << " " << vec_error_org(index-1) << " " << vec_error_gas(index-1) << " " << vec_error_compo(index-1) << endl;
      /*
        for (b=0;b<config.nbins;b++)
	cout << b << " " << abs(chp(b)-chp_save(b))/chp(b)/factor << " " << chp(b) << " " << chp_save(b) << endl;
        cout << config.nphase << " " << nh << endl;
        cout << surrogate[config.iH2O].time_aq << endl;
        cout << "error" << error_tot << endl;
        cout << "avant" << MOinit << endl;
        cout << "avant" << AQinit << endl;
        cout << "apres" << MO << endl;
        cout << "apres" << AQ << endl;*/
      cout << "chp: " << chp << " " << chp_save << endl;
      cout << AQinit << endl;
      cout << surrogate[config.iH2O].Aaq_bins_init << endl;
      cout << surrogate[config.iHp].gamma_aq_bins << endl;
      
      for (i=0;i<n;i++)
        {
          if (sum(surrogate[i].Aaq_bins_init)>1.0e-10)
            //if (surrogate[i].Atot>0.) 
            cout << surrogate[i].name << " " << surrogate[i].Aaq_bins_init << " " << surrogate[i].Aaq_bins << endl;
          if (config.compute_organic)
            if (sum(surrogate[i].Ap_layer_init)>1.0e-10)
	      {
		ilayer=0;
		iphase=0;
		for (b=0;b<config.nbins;b++)
		  cout << surrogate[i].name << " " << b << " " << surrogate[i].Ap_layer_init(b,ilayer,iphase) << " " << (surrogate[i].Ap_layer_init(b,ilayer,iphase)-surrogate[i].Ap_layer(b,ilayer,iphase))/surrogate[i].Ap_layer(b,ilayer,iphase)/factor_old << endl;
		cout << "gas " << surrogate[i].Ag << " " << (surrogate[i].Ag-surrogate[i].Ag1)/surrogate[i].Ag/factor_old << endl;
	      }
          if (surrogate[i].is_inorganic_precursor and config.compute_inorganic)
            cout << surrogate[i].name << " " << surrogate[i].Ag0 << endl;
        }
      cout << config.diameters << endl;
      cout << chp << endl;
      cout << t << " " << deltat << " " << config.max_iter << " " << factor << " " << index << endl;
      //throw string("stop");
      //exit(0);
    }
}


void solve_implicit_aqorg_repart_ssh(model_config config, vector<species> &surrogate,
                                     Array<double, 3> &MOinit, Array<double, 3> &MOW, Array<double, 1> &number,
                                     Array<double, 1> &Vsol,
                                     Array<double, 1> &LWC, Array<double, 1> &AQinit, Array<double, 1> &ionic,
                                     Array<double, 1> &chp,
                                     double &Temperature, double &RH,
                                     Array<double, 1> &AQ, Array<double, 3> &MO,
                                     Array<double, 1> &conc_inorganic, Array<double, 1> &ionic_organic, Array<double, 1> &organion, Array<double, 1> &MMaq,
                                     double t, double deltat)
{
  int b,ilayer,iphase,i;
  int n=surrogate.size();
  Array<double, 1> vec_error_org,vec_error_aq,vec_error_compo; 
  vec_error_org.resize(config.max_iter);
  vec_error_aq.resize(config.max_iter);
  vec_error_compo.resize(config.max_iter);
  bool non_convergence;
  double LWCtot=0.0;
  for (b=0;b<config.nbins;++b)
    LWCtot+=LWC(b);
  Array<double, 1> AQsave;
  Array<double, 1> chp_save;
  Array<double, 3> MOsave;
  AQsave.resize(config.nbins);
  chp_save.resize(config.nbins);
  MOsave.resize(config.nbins,config.nlayer,config.max_number_of_phases);

  //The method of newton raphson cannot be used due to the high number of variables
  //(nbins*(1+nlayer*nphase)
  //
  //Method to reach equilibrium
  // - initialisation of MOinit and AQinit
  // - initialisation of MO and AQ
  // - compute error_tot = maximum of relative errors
  //   if error_tot < 0.1% the system has converged
  // - if error_tot = error_old2 (error from iteration - 2) the system cannot converge
  //   due to strong variation of the composition.
  //   in that case a method is used so that variations are smaller:
  //     Ap(iter+1)=Ap(iter)*(1-1.0/nh)+1.0/nh*Ap(computed)
  //     1.0/nh weight factor
  //     nh is increased to ensure convergence
  //     error_tot must be inferior to factor * 0.1% to reach equilibrium
  int nh;

  if (config.compute_inorganic and config.compute_organic==false)
    nh=config.nh_inorg_init;
  else if (config.compute_inorganic and config.compute_organic)
    nh=max(config.nh_inorg_init,max(config.nh_aq_init,config.nh_org_init));
  else
    nh=max(config.nh_aq_init,config.nh_org_init);
  
  config.nh_max=7;
  int m;
  chp_save=chp;
  
  for (i=0;i<n;i++)
    {
      surrogate[i].Atot0=surrogate[i].Atot;
      surrogate[i].Ag0=surrogate[i].Ag;
      for(b=0;b<config.nbins;b++)
        {
          if (surrogate[i].hydrophilic)
            surrogate[i].Aaq_bins_init0(b)=surrogate[i].Aaq_bins_init(b);
          if (surrogate[i].hydrophobic)
            for (ilayer=0;ilayer<config.nlayer;++ilayer)
              for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                surrogate[i].Ap_layer_init0(b,ilayer,iphase)=surrogate[i].Ap_layer_init(b,ilayer,iphase); 
        }
    }

  for(b=0;b<config.nbins;b++)
    {
      AQsave(b)=AQinit(b);
      for (ilayer=0;ilayer<config.nlayer;++ilayer)
        for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
          MOsave(b,ilayer,iphase)=MOinit(b,ilayer,iphase); 
    }

  bool hygroscopicity_save=config.hygroscopicity;
  bool inorganic_save=config.compute_inorganic;
  config.compute_inorganic=false;

  vec_error_org=0.;
  vec_error_aq=0.;
  vec_error_compo=0.;

  for (b=0;b<config.nbins;b++)
    {
      double factor_max=1.;
      double factor_min=pow(0.5,config.nh_max-1);
      double var_min=0.08;
      double var_max=0.12;
      double var_rej=1.;
      double factor=pow(0.5,nh-1);
      double var;
      double factor_old=factor;
      int iiter=0;
      int ntoo_big=0;
      double error_tot=10.0;
      int index=0;
  
      chp(b)=min(1.,max(chp(b),1.e-6));
      while (error_tot>config.relative_precision and index < config.max_iter) 
        {
          if (index>2)
            {
              //ensure that the system can converge
              non_convergence=false;
              /*
                if (iiter>20)
                for (i=max(index-20,2);i<index-1;i++)
                if (((abs(vec_error_org(i)-vec_error_org(index-1))/vec_error_org(index-1)/factor<1.0e-4 and vec_error_org(index-1)/factor>config.relative_precision) and
                (abs(vec_error_org(i-1)-vec_error_org(index-2))/vec_error_org(index-2)/factor<1.0e-4) and (abs(vec_error_org(i-2)-vec_error_org(index-3))/vec_error_org(index-3)/factor<1.0e-4))
                or (abs(vec_error_aq(i)-vec_error_aq(index-1))/vec_error_aq(index-1)/factor<1.0e-4 and vec_error_aq(index-1)/factor>config.relative_precision and 
                (abs(vec_error_aq(i-1)-vec_error_aq(index-2))/vec_error_aq(index-2)/factor<1.0e-4) and (abs(vec_error_aq(i-2)-vec_error_aq(index-3))/vec_error_aq(index-3)/factor<1.0e-4))
                or (abs(vec_error_gas(i)-vec_error_gas(index-1))/vec_error_gas(index-1)/factor<1.0e-4 and vec_error_gas(index-1)/factor>config.relative_precision and 
                (abs(vec_error_gas(i-1)-vec_error_gas(index-2))/vec_error_gas(index-2)/factor<1.0e-4) and (abs(vec_error_gas(i-2)-vec_error_gas(index-3))/vec_error_gas(index-3)/factor<1.0e-4)) 
                or vec_error_org(index-1)/factor>10.0 or vec_error_aq(index-1)/factor>10.0 or vec_error_gas(index-1)/factor>10.0)
                non_convergence=true;*/

              /*
                if (iiter>100)
                if (abs(vec_error_aq(index-1))>config.relative_precision)
                {
                double a=1;
                double b=1;
                for (i=index-100;i<index-90;i++)
                a*=vec_error_aq(i);
              
                
                for (i=index-11;i<index-1;i++)
                b*=vec_error_aq(i);
                  
                if (abs(b)>abs(a)*0.999)
                non_convergence=true;
                }

                if (iiter>100)
                if (abs(vec_error_org(index-1))>config.relative_precision)
                {
                double a=1;
                double b=1;
                for (i=index-100;i<index-90;i++)
                a*=vec_error_org(i);
                
                for (i=index-11;i<index-1;i++)
                b*=vec_error_org(i);
                  
                if (abs(b)>abs(a)*0.999)
                non_convergence=true;
                }
	  
                if (iiter>100)
                if (abs(vec_error_gas(index-1))>config.relative_precision)
                {
                double a=1;
                double b=1;
                for (i=index-100;i<index-90;i++)
                a*=vec_error_gas(i);
                
                for (i=index-11;i<index-1;i++)
                b*=vec_error_gas(i);
                  
                if (abs(b)>abs(a)*0.999)
                non_convergence=true;
                }*/
              
              if (non_convergence and nh<config.nh_max)
                { 
                  if (abs(vec_error_org(index-1))>100.0 or abs(vec_error_aq(index-1))>100.0)
                    {
                      for (i=0;i<n;i++)
                        {
                          surrogate[i].Atot=surrogate[i].Atot1;
                          surrogate[i].Ag=surrogate[i].Ag1;
                          if (surrogate[i].hydrophobic)
                            for (ilayer=0;ilayer<config.nlayer;++ilayer)
                              for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                                surrogate[i].Ap_layer_init(b,ilayer,iphase)=surrogate[i].Ap_layer(b,ilayer,iphase);          
                          if (surrogate[i].hydrophilic)
                            surrogate[i].Aaq_bins_init(b)=surrogate[i].Aaq_bins(b);             
                        }
                      
                      AQinit(b)=AQsave(b);
                      for (ilayer=0;ilayer<config.nlayer;++ilayer)
                        for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                          MOinit(b,ilayer,iphase)=MOsave(b,ilayer,iphase);
                      chp(b)=chp_save(b);
                    }
                  //++nh;
                  //factor=pow(0.5,nh-1);
                  factor_max=max(factor_max/2,factor_min);
                  factor=max(factor/2,factor_min);
                  for (i=max(index-50,0);i<index;i++)
                    {
                      vec_error_org(i)=0.0;
                      vec_error_aq(i)=0.0;
                      vec_error_compo(i)=0.0;
                    }
                  iiter=0;

                }
            }

          for (ilayer=0;ilayer<config.nlayer;++ilayer)
            for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
              for (i=0;i<n;++i)
                {
                  surrogate[i].Ap_layer(b,ilayer,iphase)=
                    surrogate[i].Ap_layer_init(b,ilayer,iphase);
                  surrogate[i].gamma_org_layer0(b,ilayer,iphase)=surrogate[i].gamma_org_layer(b,ilayer,iphase);
                }

          for (i=0;i<n;i++)
            {
              surrogate[i].Aaq_bins(b)=surrogate[i].Aaq_bins_init(b);
              surrogate[i].gamma_aq_bins_old(b)=surrogate[i].gamma_aq_bins(b);
            }

          for (i=0;i<n;i++)
            surrogate[i].Ag1=surrogate[i].Ag;

          for (i=0;i<n;i++)
            surrogate[i].Atot1=surrogate[i].Atot;

          //cout << "NO3: " << sum(surrogate[config.iNO3m].Aaq_bins_init)/surrogate[config.iNO3m].MM*surrogate[config.iHNO3].MM+surrogate[config.iHNO3].Ag << " " << nh << endl;;
          //cout << "NO3i: " << sum(surrogate[config.iNO3m].Aaq_bins_init0)/surrogate[config.iNO3m].MM*surrogate[config.iHNO3].MM << " " << surrogate[config.iHNO3].Ag0 << " " << nh << endl;
   
          if (LWCtot>config.LWClimit)
            density_aqueous_phase_ssh(config, surrogate, LWC, Temperature);
          
          chp_save(b)=chp(b);
          factor_old=factor;
          AQsave(b)=AQinit(b);
          for (ilayer=0;ilayer<config.nlayer;++ilayer)
            for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
              MOsave(b,ilayer,iphase)=MOinit(b,ilayer,iphase);

          if (config.first_evaluation_activity_coefficients==false)
            {	  
              twostep_aqorg_repart_ssh(config, surrogate, config.tequilibrium, AQinit, AQ, conc_inorganic,
                                       ionic, ionic_organic, organion, chp, LWC, MOinit, MO,
                                       MOW, Temperature, RH, MMaq, true, factor, t, deltat, index, b);       
              if (config.compute_saturation and config.compute_organic)
                phase_repartition_ssh(config,surrogate,Temperature,MOinit,MO,MOW,factor);	  
            }
          else
            {
              twostep_aqorg_repart_ssh(config, surrogate, config.tequilibrium, AQinit, AQ, conc_inorganic,
                                       ionic, ionic_organic, organion, chp, LWC, MOinit, MO,
                                       MOW, Temperature, RH, MMaq, false, factor, t, deltat, index, b);
              if (config.compute_saturation and config.compute_organic)
                phase_repartition_ssh(config,surrogate,Temperature,MOinit,MO,MOW,factor);
            }

          if (iiter>200)
            {
              factor_max=max(factor_max/2,factor_min);
              factor=max(factor/2,factor_min);
              iiter=0;
            }
      
          var=0.;
          m=0;
          double negligeable=1.e-5; //config.MOmin;
          
          for (i=0;i<n;i++)
            if ((surrogate[i].is_ion and config.compute_inorganic) or i==config.iH2O or (surrogate[i].is_organic and config.compute_organic))
              if (surrogate[i].hydrophilic)
                //if (surrogate[i].hydrophilic and surrogate[i].is_inorganic_precursor==false)
                //if (error_spec(i)/factor_old>relprec)
                {
                  if (0.5*surrogate[i].gamma_aq_bins(b)+0.5*surrogate[i].gamma_aq_bins_old(b)>1.0e-6 and (surrogate[i].Aaq_bins_init(b)>negligeable or surrogate[i].Aaq_bins(b)>negligeable or i==config.iHp or i==config.iOHm))
                    {
                      m++;
                      var+=abs((surrogate[i].gamma_aq_bins_old(b)-surrogate[i].gamma_aq_bins(b))/(0.5*surrogate[i].gamma_aq_bins(b)+0.5*surrogate[i].gamma_aq_bins_old(b)));
                    }
                }
          
          if (config.compute_organic)
            for (i=0;i<n;i++)
              if (surrogate[i].is_organic or i==config.iH2O)
                if (surrogate[i].hydrophobic)
                  for (ilayer=0;ilayer<config.nlayer;++ilayer)
                    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                      //if (error_spec(i)/factor_old>relprec)
                      {
                        if (0.5*surrogate[i].gamma_org_layer(b,ilayer,iphase)+0.5*surrogate[i].gamma_org_layer0(b,ilayer,iphase)>1.0e-6 and (surrogate[i].Ap_layer_init(b,ilayer,iphase)>negligeable or surrogate[i].Ap_layer(b,ilayer,iphase)>negligeable))
                          {
                            m++;
                            var+=abs((surrogate[i].gamma_org_layer0(b,ilayer,iphase)-surrogate[i].gamma_org_layer(b,ilayer,iphase))/(0.5*surrogate[i].gamma_org_layer(b,ilayer,iphase)+0.5*surrogate[i].gamma_org_layer0(b,ilayer,iphase)));
                          }
                      }

          if (m>0)
            var=var/m;


          if (var>var_max and index>0 and factor>factor_min)
            {                      
              factor=max(factor_old*max(var_max/var,0.1),factor_min);
              /*
                if (var>var_rej)
                for (i=0;i<n;i++)
                {
                surrogate[i].Atot=surrogate[i].Atot1;
                surrogate[i].Ag=surrogate[i].Ag1;
                if (surrogate[i].hydrophobic)
                for(b=0;b<config.nbins;b++)
                for (ilayer=0;ilayer<config.nlayer;++ilayer)
                for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                surrogate[i].Ap_layer_init(b,ilayer,iphase)=surrogate[i].Ap_layer(b,ilayer,iphase);          
                if (surrogate[i].hydrophilic)
                for(b=0;b<config.nbins;b++)
                surrogate[i].Aaq_bins_init(b)=surrogate[i].Aaq_bins(b);             
                }*/
            }
          else if (var<var_min and index>0 and factor<1.)
            {
              if (var>0.)
                factor=min(factor_old*min(var_max/var,10.),factor_max);
              else
                factor=min(factor_old*10.,factor_max);
            }

          factor=max(min(factor,factor_max),factor_min);

          //redistribute concentrations to ensure that the volume of layers are constant      
          redistribution_ssh(config, surrogate,MOinit,MO);

          water_concentration_ssh(config, surrogate, Temperature, RH);
      
          tau_dif_ssh(config, surrogate, number, Vsol,Temperature);
          tau_kmt_ssh(config, surrogate, Temperature, number);
	  characteristic_time_ssh(config, surrogate, MOinit, AQinit, LWCtot); 
	  if (LWCtot>config.LWClimit)
	    characteristic_time_aq_ssh(config, surrogate, Temperature, chp, LWC, AQinit, MOinit, MMaq, ionic);

          //compute the new diameters of particle due to the growth of particles by condensation
          compute_diameters_ssh(config, surrogate, Vsol, number, LWC, LWCtot);
          if (config.explicit_representation)
            compute_morphology_ssh(config, Vsol, number);
		  
          //Computation of error_tot
          vec_error_org(index)=0.0;
          vec_error_aq(index)=0.0;
          vec_error_compo(index)=0.0;
          double errloc;
              
          AQ(b)=max(AQ(b),config.MOmin);         
          for (ilayer=0;ilayer<config.nlayer;++ilayer)
            for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
              //if (abs(MO(b,ilayer,iphase)-MOinit(b,ilayer,iphase))>config.precision)
              {
                if (MO(b,ilayer,iphase)>1.0e-5)
                  {
                    errloc=(MO(b,ilayer,iphase)-MOinit(b,ilayer,iphase))/MO(b,ilayer,iphase)/factor_old;
                    if (abs(errloc)>abs(vec_error_org(index)))
                      vec_error_org(index)=errloc;
                  }
		    
                //if (var<=var_rej)
                MOinit(b,ilayer,iphase)=max(MO(b,ilayer,iphase),config.MOmin*config.Vlayer(ilayer));
              }
          
          if (abs(AQinit(b)-AQ(b))>config.precision)
            {
              if (AQ(b)>1.0e-5)
                {
                  errloc=(AQ(b)-AQinit(b))/AQ(b)/factor_old;
                  if (abs(errloc)>abs(vec_error_aq(index)))
                    vec_error_aq(index)=errloc;
                  
                  //vec_error_aq(index)=max(vec_error_aq(index),vec_error_chp(index));
                }
       
              /*
              for (i=0;i<n;i++)
                if (surrogate[i].hydrophilic and i!=config.iHp and i!=config.iHp)
                  if (surrogate[i].Aaq_bins(b)>0 and (surrogate[i].Aaq_bins(b)>config.MOmin or surrogate[i].Aaq_bins_init(b)>config.MOmin))// and abs(surrogate[i].Aaq_bins_init(b)-surrogate[i].Aaq_bins(b))>config.precision)
                    {
                      errloc=(surrogate[i].Aaq_bins_init(b)-surrogate[i].Aaq_bins(b))/surrogate[i].Aaq_bins(b)/factor_old;
                      if (abs(errloc)>abs(vec_error_aq(index)))
                        vec_error_compo(index)=errloc;                      
                    }*/

              for (i=0;i<n;i++)
                if (surrogate[i].hydrophilic and i!=config.iHp and i!=config.iOHm)
                  if (surrogate[i].Aaq_bins(b)>0 and (surrogate[i].Aaq_bins(b)>negligeable or surrogate[i].Aaq_bins_init(b)>negligeable))// and abs(surrogate[i].Aaq_bins_init(b)-surrogate[i].Aaq_bins(b))>config.precision)
                    {
                      /*if (abs(surrogate[i].Aaq_bins_init(b)-surrogate[i].Aaq_bins(b))/surrogate[i].Aaq_bins(b)>vec_error_aq(index))
                        cout << surrogate[i].name << " " << abs(surrogate[i].Aaq_bins_init(b)-surrogate[i].Aaq_bins(b))/surrogate[i].Aaq_bins(b) << " " << surrogate[i].gamma_aq_bins(b) << b << " " << surrogate[i].Aaq_bins_init(b) << " " << surrogate[i].Aaq_bins(b) << " " << b << endl;*/
                      errloc=(surrogate[i].Aaq_bins_init(b)-surrogate[i].Aaq_bins(b))/surrogate[i].Aaq_bins(b)/factor_old;
                      if (abs(errloc)>abs(vec_error_compo(index)))
                        vec_error_compo(index)=errloc;
                      /*
                        if (abs(errloc)>1.e20)
                        {
                        cout << "errloc: " << surrogate[i].name << " " << errloc << " " << surrogate[i].Aaq_bins_init(b) << " " << surrogate[i].Aaq_bins(b) << " " << factor << endl;
                        //exit(1);
                        }*/
                      
                    }
	    }
	      

          for (i=0;i<n;i++)
	    if (surrogate[i].hydrophobic)
	      for (ilayer=0;ilayer<config.nlayer;++ilayer)
		for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		  if (surrogate[i].Ap_layer(b,ilayer,iphase)>config.Vlayer(ilayer)*negligeable)
		    {
		      errloc=(surrogate[i].Ap_layer_init(b,ilayer,iphase)-surrogate[i].Ap_layer(b,ilayer,iphase))/surrogate[i].Ap_layer(b,ilayer,iphase)/factor_old;
		      if (abs(errloc)>abs(vec_error_compo(index)))
			vec_error_compo(index)=errloc;
		    }
      
     
          AQinit(b)=max(AQ(b),config.MOmin);       

          if (iiter>2)
            {             
              if (vec_error_org(index)*vec_error_org(index-1)<0.0 and abs(vec_error_aq(index))<config.relative_precision and abs(vec_error_org(index)-vec_error_org(index-2))<0.01*abs(vec_error_org(index)))
                ntoo_big=ntoo_big+10;
              else if (vec_error_aq(index)*vec_error_aq(index-1)<0.0 and abs(vec_error_org(index))<config.relative_precision and abs(vec_error_aq(index)-vec_error_aq(index-2))<0.01*abs(vec_error_aq(index)))                        
                ntoo_big=ntoo_big+10;
              else if (vec_error_org(index)*vec_error_org(index-1)<0.0 and abs(vec_error_aq(index))>config.relative_precision and vec_error_aq(index)*vec_error_aq(index-1)<0.0 and abs(vec_error_org(index))>config.relative_precision and abs(vec_error_org(index)-vec_error_org(index-2))<0.01*abs(vec_error_org(index)) and abs(vec_error_aq(index)-vec_error_aq(index-2))<0.01*abs(vec_error_aq(index))) 
                ntoo_big=ntoo_big+10;           
          
              factor=max(factor,factor_min);
            }                  

          
          if (ntoo_big>200 and factor_max>factor_min)
            {
              factor_max=max(factor_max/2,factor_min);
              ntoo_big=0;                  
            }

          error_tot=max(max(abs(vec_error_org(index)),abs(vec_error_aq(index))),abs(vec_error_compo(index)));

          ++index;
          ++iiter;
        }

      //cout << index << " " << error_tot << endl;
      if (error_tot>config.relative_precision)
        {
          cout << "The model did not converged... " << RH << " " << Temperature << endl;
          cout << "error tot" << error_tot << endl;
          cout << "coupled-repart " << vec_error_aq(index-1) << " " << vec_error_org(index-1) << " " << vec_error_compo(index-1) << " " << b << endl;
          
          for (i=0;i<n;i++)
            {
              if (surrogate[i].Aaq_bins_init(b)>1.0e-10)
                //if (surrogate[i].Atot>0.) 
                cout << surrogate[i].name << " " << surrogate[i].Aaq_bins_init(b) << endl;
              if (surrogate[i].Ap_layer_init(b)>1.0e-10)
                cout << surrogate[i].name << " " << surrogate[i].Ap_layer_init(b) << endl;
            }
          cout << chp(b) << endl;
          cout << t << " " << deltat << " " << config.max_iter << " " << factor << " " << index << endl;
          //throw string("stop");
          //exit(0);
        }
    }
  //cout << index << endl;
  config.compute_inorganic=inorganic_save;
}



void solve_implicit_ssh(model_config config, vector<species> &surrogate,
                        Array<double, 3> &MOinit, Array<double, 3> &MOW, Array<double, 1> &number,
                        Array<double, 1> &Vsol,
                        Array<double, 1> &LWC, Array<double, 1> &AQinit, Array<double, 1> &ionic,
                        Array<double, 1> &chp,
                        double &Temperature, double &RH,
                        Array<double, 1> &AQ, Array<double, 3> &MO,
                        Array<double, 1> &conc_inorganic, Array<double, 1> &ionic_organic, Array<double, 1> &organion, Array<double, 1> &MMaq,
                        double t, double deltat, int &index, bool &reject_step, Array<double, 1> &chp20)
{
  int b,ilayer,iphase,i;
  int n=surrogate.size();
  /*
    double error_tot=10.0;
    int index=0;
    Array<double, 1> vec_error_org,vec_error_aq,vec_error_gas,vec_error_chp,vec_error_compo; 
    vec_error_org.resize(config.max_iter);
    vec_error_aq.resize(config.max_iter);
    vec_error_gas.resize(config.max_iter);
    vec_error_chp.resize(config.max_iter);
    vec_error_compo.resize(config.max_iter);
    vec_error_org=0.;
    vec_error_aq=0.;
    vec_error_gas=0.;
    vec_error_chp=0.;
    vec_error_compo=0.;
    bool non_convergence;*/
  double LWCtot=0.0;
  for (b=0;b<config.nbins;++b)
    LWCtot+=LWC(b);
  reject_step=false;
  /*
    Array<double, 1> AQsave;
    Array<double, 1> chp_save;
    Array<double, 3> MOsave;
    AQsave.resize(config.nbins);
    chp_save.resize(config.nbins);
    MOsave.resize(config.nbins,config.nlayer,config.max_number_of_phases);*/

  //The method of newton raphson cannot be used due to the high number of variables
  //(nbins*(1+nlayer*nphase)
  //
  //Method to reach equilibrium
  // - initialisation of MOinit and AQinit
  // - initialisation of MO and AQ
  // - compute error_tot = maximum of relative errors
  //   if error_tot < 0.1% the system has converged
  // - if error_tot = error_old2 (error from iteration - 2) the system cannot converge
  //   due to strong variation of the composition.
  //   in that case a method is used so that variations are smaller:
  //     Ap(iter+1)=Ap(iter)*(1-1.0/nh)+1.0/nh*Ap(computed)
  //     1.0/nh weight factor
  //     nh is increased to ensure convergence
  //     error_tot must be inferior to factor * 0.1% to reach equilibrium
  /*
    int iiter=0;
    int nh;
    int ntoo_big=0;

    if (config.compute_inorganic and config.compute_organic==false)
    nh=config.nh_inorg_init;
    else if (config.compute_inorganic and config.compute_organic)
    nh=max(config.nh_inorg_init,max(config.nh_aq_init,config.nh_org_init));
    else
    nh=max(config.nh_aq_init,config.nh_org_init);
  
    config.nh_max=7;

    double factor_max=1.;
    double factor_min=pow(0.5,config.nh_max-1);
    double var_min=0.08;
    double var_max=0.12;
    double var_rej=1.;
    double factor=pow(0.5,nh-1);
    double var;
    double factor_old=factor;
    int m;
    chp_save=chp;*/

  if (config.aqorg_repart)
    {
      solve_implicit_aqorg_repart_ssh(config, surrogate, MOinit, MOW, number, Vsol, LWC, AQinit, ionic, chp, Temperature, RH, AQ, MO,
                                      conc_inorganic, ionic_organic, organion, MMaq, t, deltat);
    }
    
  if (config.imethod>=3 and config.hygroscopicity)
    {
      solve_implicit_water_coupled_ssh(config, surrogate, MOinit, MOW, number, Vsol, LWC, AQinit, ionic, chp, Temperature, RH, AQ, MO,
                                       conc_inorganic, ionic_organic, organion, MMaq, t, 0);
    }
    
  bool hygroscopicity_save=config.hygroscopicity;
  bool compute_organic_save=config.compute_organic;
  if (config.imethod==3)
    {
      config.hygroscopicity=false;
      if (config.compute_inorganic)
        config.compute_organic=false;
    }

  int index1=index;
  int index2=index;
    
  solve_implicit_coupled_ssh(config, surrogate, MOinit, MOW, number, Vsol, LWC, AQinit, ionic, chp, Temperature, RH, AQ, MO,
                             conc_inorganic, ionic_organic, organion, MMaq, t, deltat, index1, reject_step, chp20);
  index=index1;

  if (config.imethod==3 and config.compute_inorganic and compute_organic_save)
    {
      config.compute_organic=true;
      config.compute_inorganic=false;
      bool reject_step2;
      //config.compute_hygroscopicity=true;
      if (config.compute_organic)
        solve_implicit_coupled_ssh(config, surrogate, MOinit, MOW, number, Vsol, LWC, AQinit, ionic, chp, Temperature, RH, AQ, MO,
                                   conc_inorganic, ionic_organic, organion, MMaq, t, deltat, index2, reject_step2, chp20);
      if (reject_step2)
	reject_step=true;
      config.compute_inorganic=true;
      index=max(index1,index2);
    }
  
  
  
  if (config.imethod==3)
    {
      config.hygroscopicity=hygroscopicity_save;
      config.compute_organic=compute_organic_save;
      if (config.hygroscopicity)
	solve_implicit_water_coupled_ssh(config, surrogate, MOinit, MOW, number, Vsol, LWC, AQinit, ionic, chp, Temperature, RH, AQ, MO,
                                         conc_inorganic, ionic_organic, organion, MMaq, t, 0);
    }
  
  if (config.compute_saturation and config.first_evaluation_of_saturation==false and config.compute_organic)
    {
      number_org_phases_ssh(config,surrogate,Temperature,MOinit,MOW);
      compute_kp_org_ssh(config, surrogate, MOinit, Temperature, MOW);
      tau_dif_ssh(config, surrogate, number, Vsol, Temperature);
      tau_kmt_ssh(config, surrogate, Temperature, number);
      characteristic_time_ssh(config, surrogate, MOinit, AQinit, LWCtot);
      if (LWCtot>config.LWClimit)
        characteristic_time_aq_ssh(config, surrogate, Temperature, chp, LWC, AQinit, MOinit, MMaq, ionic);
    }

}
void initialisation_ssh(model_config &config, vector<species> &surrogate,
                        Array<double,3> &MOinit, Array<double, 3> &MO, Array<double, 3> &MOW,
                        Array<double,1> &AQinit,Array<double,1> &AQ, Array<double,1> &MMaq,
                        double &LWCtot, Array<double,1> &LWC, Array<double, 1> &chp, Array<double, 1> &ionic,
                        Array<double, 1> &ionic_organic, Array<double, 1> &organion, Array<double, 1> &conc_inorganic,
			Array <double, 1> &number, Array<double, 1> &Vsol,
                        double &Temperature, double &RH)
{  
  int b,ilayer,iphase,i;
  int n=surrogate.size();
  Array <double, 1> conc_org;
  conc_org.resize(config.nbins);
  double MOWloc=1.;
  double T0=298.15;
  double deltaH_over_RT0,deltaCp0_over_R;
  double R=8.314; //ideal gas constant (J/K/mol)
  /*  if (config.compute_inorganic)
      for (b=0;b<config.nbins;b++)
      if (chp(b)<=1.e-8)
      chp(b)=1.e-2;*/
  config.wat_min=0.;

  config.viscosity_layer=0.;
  if (config.compute_viscosity)
    compute_pure_viscosity(config, surrogate, Temperature);

  if (config.isorropia_ph)
     for (b=0;b<config.nbins;++b)
       chp(b)=max(chp(b),1.e-7);

  for (i=0;i<n;i++)
    {
      surrogate[i].gamma_aq_bins=1.;          
      if (surrogate[i].hydrophilic==false)
	surrogate[i].Aaq_bins_init=0;
      if (surrogate[i].hydrophobic==false)
	surrogate[i].Ap_layer_init=0;
      surrogate[i].velocity=pow(2.11714271498563e4*Temperature/surrogate[i].MM,0.5);
      surrogate[i].knui=pow(3.0*surrogate[i].KDiffusion_air/surrogate[i].velocity,2);
      for (b=0;b<config.nbins;b++)
        {
          surrogate[i].fac_corr_ph(b)=1.;
          surrogate[i].LR(b)=1.0;
          surrogate[i].SRMR(b)=1.0;
        }
      if (surrogate[i].is_solid==false)
	surrogate[i].Asol_bins_init=0;
    }
  
  for (b=0;b<config.nbins;b++)
    MMaq(b)=18.;
  
  for (i=0;i<n;i++)
    if (surrogate[i].is_organic or i==config.iH2O)
      {
	if (surrogate[i].nonvolatile==false)
	  if(surrogate[i].hydrophobic or i==config.iH2O)
	    {        
	      if (surrogate[i].kp_from_experiment)
		surrogate[i].kpi=surrogate[i].Kp_exp_org_ssh(Temperature);
	      else if (surrogate[i].kp_from_experiment==false)
		surrogate[i].kpi=surrogate[i].Kp_eff_org_ssh(Temperature, MOWloc);
	    }


	if (i==config.iH2O)
	  surrogate[i].kaqi=760.0*8.202e-5*Temperature/(MOWloc*1.0e6*surrogate[i].Psat_ssh(Temperature));
	else
	  if (surrogate[i].hydrophilic)
	    for (b=0;b<config.nbins;b++)
	      {
		double gamma=pow(10,-0.511*pow(298.0/Temperature,1.5)*pow(ionic(b),0.5)/(1.0+pow(ionic(b),0.5)));
		if (config.compute_aqueous_phase_properties or chp(b)==0.)
		  if (surrogate[i].nonvolatile)
		    surrogate[i].veckaqi(b)=1000.;
		  else
		    surrogate[i].veckaqi(b)=surrogate[i].Kpart_aq_ssh(Temperature,MOWloc);              
		else
		  if (surrogate[i].nonvolatile)
		    {
		      surrogate[i].veckaqi(b)=1000.;
		      if (surrogate[i].aqt==2) //diacid
			{
			  surrogate[i].vecfioni1(b)=(surrogate[i].Kacidity1/(gamma*chp(b)))/
			    (1.0+surrogate[i].Kacidity1/(gamma*chp(b))*(1.0+surrogate[i].Kacidity2/(gamma*chp(b))));
			  surrogate[i].vecfioni2(b)=(surrogate[i].Kacidity1/(gamma*chp(b)))*(surrogate[i].Kacidity2/(gamma*chp(b)))/
			    (1.0+surrogate[i].Kacidity1/(gamma*chp(b))*(1.0+surrogate[i].Kacidity2/(gamma*chp(b))));
			}
		      else if (surrogate[i].aqt==1) //monoacid
			surrogate[i].vecfioni1(b)=(surrogate[i].Kacidity1/(gamma*chp(b)))/(1.0+surrogate[i].Kacidity1/(gamma*chp(b)));			  
		    }
		  else
		    {		      
		      if (surrogate[i].aqt==2) //diacid
			{
			  surrogate[i].veckaqi(b)=surrogate[i].Kpart_aq_ssh(Temperature,MOWloc)*
			    (1.0+surrogate[i].Kacidity1/(gamma*chp(b))*
			     (1.0+surrogate[i].Kacidity2/(gamma*chp(b))));
			  surrogate[i].vecfioni1(b)=(surrogate[i].Kacidity1/(gamma*chp(b)))/
			    (1.0+surrogate[i].Kacidity1/(gamma*chp(b))*(1.0+surrogate[i].Kacidity2/(gamma*chp(b))));
			  surrogate[i].vecfioni2(b)=(surrogate[i].Kacidity1/(gamma*chp(b)))*(surrogate[i].Kacidity2/(gamma*chp(b)))/
			    (1.0+surrogate[i].Kacidity1/(gamma*chp(b))*(1.0+surrogate[i].Kacidity2/(gamma*chp(b))));
			}
		      else if (surrogate[i].aqt==1) //monoacid
			{
			  surrogate[i].veckaqi(b)=surrogate[i].Kpart_aq_ssh(Temperature,MOWloc)*(1.0+surrogate[i].Kacidity1/(gamma*chp(b)));
			  surrogate[i].vecfioni1(b)=(surrogate[i].Kacidity1/(gamma*chp(b)))/(1.0+surrogate[i].Kacidity1/(gamma*chp(b)));
			}
		      else if (surrogate[i].aqt==3) //aldehyde
			surrogate[i].veckaqi(b)=surrogate[i].Kpart_aq_ssh(Temperature,MOWloc)*(1.0+surrogate[i].Koligo_aq*pow(chp(b)/pow(10,-surrogate[i].pHref),surrogate[i].beta));
		      else
			surrogate[i].veckaqi(b)=surrogate[i].Kpart_aq_ssh(Temperature,MOWloc);
		    }
	      }
      }
  
    else if (surrogate[i].is_inorganic_precursor and config.compute_inorganic and surrogate[i].is_solid==false)
      {
        if (i==config.iH2SO4)
          {
            surrogate[i].kaqi=1.0e10;
            surrogate[i].keqi=surrogate[i].Kequilibrium_ssh(Temperature);
          }
        else if (i==config.iNH3)
          {
            deltaH_over_RT0=-13.79;
            deltaCp0_over_R=5.39;
            surrogate[i].kaqi=
              surrogate[i].Henry*exp(-deltaH_over_RT0*(T0/Temperature-1.0)-deltaCp0_over_R*(1.+log(T0/Temperature)-T0/Temperature))
              *R*Temperature/(1000.*1.0e6*1.013e5);
            surrogate[i].keqi=surrogate[i].Kequilibrium_ssh(Temperature);
          }
        else if (i==config.iHNO3)
          {
            deltaH_over_RT0=-29.17;
            deltaCp0_over_R=-16.83;
            surrogate[i].kaqi=
              surrogate[i].Henry*exp(-deltaH_over_RT0*(T0/Temperature-1.0)-deltaCp0_over_R*(1.+log(T0/Temperature)-T0/Temperature))
              *R*Temperature/(1000.*1.0e6*1.013e5);
            surrogate[i].keqi=surrogate[i].Kequilibrium_ssh(Temperature);
          }
	else if (i==config.iCO2)
	  {
	    surrogate[i].kaqi=exp(2.50E+02+0.0457081*Temperature-1.59e+04/298-4.05E+01*log(Temperature)+1.54E+06/Temperature/Temperature)*R*Temperature/(1000.*1.0e6*1.013e5);
	    Kequilibrium_co2_ssh(config, surrogate, Temperature);	      
	  }
        else if (i==config.iHCl)
          {
            
            deltaH_over_RT0=-30.20;
            deltaCp0_over_R=-19.91;
            surrogate[i].kaqi=surrogate[i].Henry*exp(-deltaH_over_RT0*(T0/Temperature-1.0)-deltaCp0_over_R*(1.+log(T0/Temperature)-T0/Temperature))
              *R*Temperature/(1000.*1.0e6*1.013e5);
            surrogate[i].keqi=surrogate[i].Kequilibrium_ssh(Temperature);
          }
      }

  for (i=0;i<n;i++)
    if (surrogate[i].is_solid)
      surrogate[i].keq=surrogate[i].Ksol*exp(-surrogate[i].deltaH*(298./Temperature-1.0)-surrogate[i].dCp*(1.+log(298./Temperature)-298/Temperature));

  double Pwater=surrogate[config.iH2O].Psat_ssh(Temperature)*RH;
  surrogate[config.iH2O].Atot=(Pwater/760.0*1.013e5)*surrogate[config.iH2O].MM*1.0e6/
    (8.314*Temperature);

  if (config.activity_model=="unifac")
    {
      double tval1=1.0/298.15-1.0/Temperature;
      double tval2=298.15/Temperature-1.0+log(Temperature/298.15);
      /*
        config.Inter2_aq.resize(config.nfunc_aq,config.nfunc_aq);
        config.Inter2_org.resize(config.nfunc_org,config.nfunc_org);
        config.Inter2_tot.resize(config.nfunc_tot,config.nfunc_tot);*/
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

      if (config.compute_viscosity)
	if (config.temperature_dependancy)
	  for (j=0;j<config.nfunc_tot;j++)
	    for (k=0;k<config.nfunc_tot;k++)          
	      config.logInter2_tot(j,k)=-config.Inter_tot(j,k)/Temperature+config.InterB_tot(j,k)*tval1+config.InterC_tot(j,k)*tval2;         
	else
	  for (j=0;j<config.nfunc_tot;j++)
	    for (k=0;k<config.nfunc_tot;k++)          
	      config.logInter2_tot(j,k)=-config.Inter_tot(j,k)/Temperature;

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

      if (config.compute_viscosity)
	if (config.temperature_dependancy)
	  for (j=0;j<config.nfunc_org;j++)
	    for (k=0;k<config.nfunc_org;k++)
	      config.logInter2_org(j,k)=-config.Inter_org(j,k)/Temperature+config.InterB_org(j,k)*tval1+config.InterC_org(j,k)*tval2;         
	else
	  for (j=0;j<config.nfunc_aq;j++)
	    for (k=0;k<config.nfunc_aq;k++)          
	      config.logInter2_org(j,k)=-config.Inter_org(j,k)/Temperature;	

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
  
  for (b=0;b<config.nbins;++b)		  
    for (ilayer=0;ilayer<config.nlayer;++ilayer)
      config.nphase(b,ilayer)=1;
  config.Ke=1.010e-14*exp(-22.52*(298./Temperature-1.0)+26.92*(1+log(298./Temperature)-298./Temperature));

  if (config.compute_inorganic)
    {      
      config.LWClimit=-1.0;
      if (config.compute_long_and_medium_range_interactions)
	config.compute_aqueous_phase_properties=true;
      else
	config.compute_aqueous_phase_properties=false;

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
  else
    config.compute_aqueous_phase_properties=false;

  for (b=0;b<config.nbins;++b)		  
    for (ilayer=0;ilayer<config.nlayer;++ilayer)
      for (iphase=0;iphase<config.max_number_of_phases;++iphase)
        MO(b,ilayer,iphase)=0.0;

  density_aqueous_phase_ssh(config, surrogate, LWC, Temperature);
  compute_diameters_ssh(config, surrogate, Vsol, number, LWC, LWCtot);

  //initialisation of some parameters
  if (config.hygroscopicity)
    for (b=0;b<config.nbins;++b)
      {
	double AQ2=0.;
	AQinit(b)=LWC(b);
	for (i=0;i<n;++i)
	  if (surrogate[i].hydrophilic)
	    {
	      AQinit(b)+=surrogate[i].Aaq_bins_init(b);
	      if (i!=config.iH2O)
		AQ2+=surrogate[i].Aaq_bins_init(b);
	    }
	AQinit(b)=max(AQinit(b),config.MOmin);
	AQ2=max(AQ2,config.MOmin);
	if (surrogate[config.iH2O].Aaq_bins_init(b)<0.01*AQ2)
	  {
	    AQ2=AQinit(b);
	    surrogate[config.iH2O].Aaq_bins_init(b)=AQ2; //*0.01;
	  }
	
      }

  for (b=0;b<config.nbins;b++)
    {
      conc_org(b)=LWC(b);  
      for (i=0;i<n;++i)
        if((surrogate[i].is_organic or i==config.iH2O) and surrogate[i].hydrophilic)
          conc_org(b)+=surrogate[i].Aaq_bins_init(b);
      conc_org(b)=max(conc_org(b),config.MOmin);
    }
 
  if (LWCtot>config.LWClimit)
    density_aqueous_phase_ssh(config, surrogate, LWC, Temperature);
  compute_diameters_ssh(config, surrogate, Vsol, number, LWC, LWCtot);
 
  for (b=0;b<config.nbins;++b)
    {
      ionic_organic(b)=0.0;
      organion(b)=0.0;
      conc_inorganic(b)=0.0;
      if (config.compute_inorganic)
        {  
          LWC(b)=0.0;
          if (AQinit(b)>0.0)
            {
              double inorg1;
              inorg1=0.0;
              for (i=0;i<n;++i) 
                if (surrogate[i].is_organic==false and surrogate[i].is_inorganic_precursor==false and i!=config.iHp and i!=config.iOHm and i!=config.iH2O)
                  {
                    inorg1-=surrogate[i].Aaq_bins_init(b)/surrogate[i].MM*surrogate[i].charge*config.AQrho(b)/AQinit(b);
                  }
              chp(b)=0.5*(inorg1+pow(pow(inorg1,2)+4*config.Ke,0.5));

	      if (chp(b)<1.e-14)
		chp(b)=1.0e-3;
	    }
	  else
	    chp(b)=1.0e-7;
	}
      else
	{
	  for (i=0;i<n;++i)
	    if (surrogate[i].is_organic==false and i!=config.iH2O)
	      conc_inorganic(b)+=surrogate[i].Aaq_bins_init(b);
	}
    }

  density_aqueous_phase_ssh(config, surrogate, LWC, Temperature);
  compute_diameters_ssh(config, surrogate, Vsol, number, LWC, LWCtot);
  for (b=0;b<config.nbins;++b)		  
    for (ilayer=0;ilayer<config.nlayer;++ilayer)
      for (i=0;i<n;++i)
	if (surrogate[i].hydrophobic)
	  for (iphase=config.nphase(b,ilayer);iphase<config.max_number_of_phases;++iphase)
	    {
	      surrogate[i].Ap_layer_init(b,ilayer,config.nphase(b,ilayer))+=surrogate[i].Ap_layer_init(b,ilayer,iphase);
	      surrogate[i].Ap_layer_init(b,ilayer,iphase)=0.0;
	    }
  
  for (b=0;b<config.nbins;++b)		  
    for (ilayer=0;ilayer<config.nlayer;++ilayer)
      for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
	/*if (MOinit(b,ilayer,iphase)>0.0)
	  {
          double error=100.0;
          int index=0.0;
          double temp;
          while (error>config.relative_precision and index<10)
          {
          if (config.hygroscopicity)
          {   
          activity_coefficients_dyn_sat_ssh(config, surrogate, Temperature, MOW, b, ilayer);
          temp=surrogate[config.iH2O].MM/MOW(b,ilayer,iphase)*MOinit(b,ilayer,iphase)*RH
          /surrogate[config.iH2O].gamma_org_layer(b,ilayer,iphase);
          if (surrogate[config.iH2O].Ap_layer_init(b,ilayer,iphase) > 1.e-20) // YK
          error=(temp-surrogate[config.iH2O].Ap_layer_init(b,ilayer,iphase))/surrogate[config.iH2O].Ap_layer_init(b,ilayer,iphase);
          surrogate[config.iH2O].Ap_layer_init(b,ilayer,iphase)=temp;
          }
          else
          error = 0.0;
          index++;
          MOinit(b,ilayer,iphase)=0.0;
          for (i=0;i<n;++i)
          if (surrogate[i].hydrophobic)
          MOinit(b,ilayer,iphase)+=surrogate[i].Ap_layer_init(b,ilayer,iphase);
          }
	  }
          else*/
        {
          //Remplace concentrations of organic aerosols equal to zero by MOmin
          //(to ensure formation of organic aerosols)
          MOinit(b,ilayer,iphase)=max(MOinit(b,ilayer,iphase),
                                      config.MOmin*config.Vlayer(ilayer));
          for (i=0;i<n;++i)
            surrogate[i].gamma_org_layer(b,ilayer,iphase)=1.0;
        }

  density_aqueous_phase_ssh(config, surrogate, LWC, Temperature);
  compute_diameters_ssh(config, surrogate, Vsol, number, LWC, LWCtot); 
  for (i=0;i<n;++i)
    if (surrogate[i].hydrophilic or surrogate[i].is_inorganic_precursor)
      for (b=0;b<config.nbins;++b)
	surrogate[i].time_aq(b)=2.0*config.tequilibrium;

  //Computation of Ag for water
  water_concentration_ssh(config, surrogate, Temperature, RH);
  compute_diameters_ssh(config, surrogate, Vsol, number, LWC, LWCtot);
  tau_kmt_ssh(config, surrogate, Temperature, number);

  //Initialize Kaq to avoid initialization issues
  for (i=0;i<n;i++)
    surrogate[i].Kaq=config.kp_low_volatility*2;
  activity_coefficients_dyn_aq_ssh(config, surrogate, Temperature,AQinit,MOinit,
                                   conc_inorganic, ionic, ionic_organic,
                                   organion,chp,LWC,MMaq,0.0,0.,0);	   
  if (config.compute_inorganic)
    for (b=0;b<config.nbins;++b)
      {
	double error=100.0;
	int index=0;
	double temp;
	double factor=1.0;
	double XH2O=0.0;
	while (error>config.relative_precision and index<1)
	  {
	    density_aqueous_phase_ssh(config, surrogate, LWC, Temperature);
	    config.rho_aqueous=config.AQrho(b);

            for (i=0;i<n;++i)
              surrogate[i].Aaq=surrogate[i].Aaq_bins_init(b);
            
            conc_org(b)=LWC(b);	  
            for (i=0;i<n;i++)
              if (surrogate[i].is_organic or i==config.iH2O)
                conc_org(b)+=surrogate[i].Aaq_bins_init(b);
            conc_org(b)=max(conc_org(b),config.MOmin);

	    compute_ionic_strenght2_ssh(config, surrogate, Temperature, AQinit(b), conc_inorganic(b), ionic(b), chp(b),
                                        organion(b), ionic_organic(b), conc_org(b), factor, 1.e-14);

	    /*
		  
	      activity_coefficients_aq_ssh(config,surrogate,Temperature,0.0,MMaq(b),XH2O,conc_org(b));
	      if (config.compute_long_and_medium_range_interactions)
	      activity_coefficients_LR_MR_ssh(config, surrogate, Temperature, 0.0, ionic(b));
	      if (index==0)
              surrogate[config.iH2O].gamma_aq_bins(b)=1.0;
	      else
	      surrogate[config.iH2O].gamma_aq_bins(b)=surrogate[config.iH2O].gamma_aq;*/	   	   
            
	    compute_kp_aq_ssh(config, surrogate, Temperature, ionic, chp, AQinit, MMaq, b, b+1);
	    /*
	      temp=surrogate[config.iH2O].Aaq_bins_init(b);
	      surrogate[config.iH2O].Aaq_bins_init(b)=surrogate[config.iH2O].Atot*surrogate[config.iH2O].Kaq(b)*AQinit(b)/(1.0+surrogate[config.iH2O].Kaq(b)*AQinit(b));
            
	      if (surrogate[config.iH2O].Aaq_bins_init(b)>0.)
              error=abs(temp-surrogate[config.iH2O].Aaq_bins_init(b))/surrogate[config.iH2O].Aaq_bins_init(b);
	      else
	      error=1000.;*/
	  	    
	    if (AQinit(b)>0.0)
	      {
		double inorg1;
		inorg1=0.0;
		for (i=0;i<n;++i) 
		  if (surrogate[i].is_organic==false and surrogate[i].is_inorganic_precursor==false and i!=config.iHp and i!=config.iOHm and i!=config.iH2O)
		    inorg1-=surrogate[i].Aaq_bins_init(b)/surrogate[i].MM*surrogate[i].charge*config.AQrho(b)/AQinit(b);
		   
		chp(b)=0.5*(inorg1+pow(pow(inorg1,2)+4*config.Ke,0.5));
		if (chp(b)<1.e-14)
		  chp(b)=1.0e-2;
		
		//chp(b)=min(chp(b),1.0);
		//surrogate[config.iHp].Aaq_bins_init(b)=0.;
	      }
	    else
	      chp(b)=1.0e-7;

            //surrogate[config.iH2O].Aaq_bins_init(b)=temp;           
	    //surrogate[config.iHp].Aaq_bins_init(b)=chp(b)*AQinit(b)/config.AQrho(b);

	    AQinit(b)=0.0;
	    for (i=0;i<n;++i)
	      if (surrogate[i].hydrophilic)
                AQinit(b)+=surrogate[i].Aaq_bins_init(b);
            AQinit(b)=max(AQinit(b),config.MOmin);

	    index++;
	  }

	for (i=0;i<n;++i)
	  if (surrogate[i].hydrophilic)
	    surrogate[i].gamma_aq_bins(b)=surrogate[i].gamma_aq;
      }
  else
    if (LWCtot>config.LWClimit)
      activity_coefficients_dyn_aq_ssh(config, surrogate, Temperature,AQinit,MOinit,
                                       conc_inorganic, ionic, ionic_organic,
                                       organion,chp,LWC,MMaq,0.0,0.,0);
  for (i=0;i<n;i++)
    {
      surrogate[i].Ap_layer=surrogate[i].Ap_layer_init;
      surrogate[i].Ap_layer_init0=surrogate[i].Ap_layer_init;
      surrogate[i].Aaq_bins=surrogate[i].Aaq_bins_init;
      surrogate[i].Aaq_bins_init0=surrogate[i].Aaq_bins_init;
      surrogate[i].Jdn=0.;
      surrogate[i].k1=0.;   
      surrogate[i].Jdn_aq=0.;
      surrogate[i].k1_aq=0.;      
    }
  MO=MOinit;
  AQ=AQinit;
  
  for (b=0;b<config.nbins;++b)		  
    for (ilayer=0;ilayer<config.nlayer;++ilayer)
      for (iphase=0;iphase<config.max_number_of_phases;++iphase)
        if(MOW(b,ilayer,iphase)<1)
          MOW(b,ilayer,iphase) = 200.;

  config.gamma_MR_ions_bins=-1.0;

  for (b=0;b<config.nbins;++b)
    if (chp(b)>0)
      surrogate[config.iOHm].Aaq_bins_init(b)=config.Ke/chp(b)*LWC(b)/1000.;
    else
      surrogate[config.iOHm].Aaq_bins_init(b)=0.;

  if (config.compute_viscosity)
    activity_coefficients_dyn_org_ssh(config, surrogate, Temperature, MOW);

  characteristic_time_ssh(config, surrogate, MOinit, AQinit, LWCtot); 
  if (LWCtot>config.LWClimit)
    characteristic_time_aq_ssh(config, surrogate, Temperature, chp, LWC, AQinit, MOinit, MMaq, ionic);

  /*if (config.isorropia_ph==false)
    for (b=0;b<config.nbins;b++)
      chp(b)=min(1.,max(chp(b),1.e-6));
  */

}

void dynamic_system_ssh(model_config &config, vector<species> &surrogate,
                        Array<double, 3> &MOinit, Array<double, 3> &MOW, Array<double, 1> &number,
                        Array<double, 1> &Vsol,
                        Array<double, 1> &LWC, Array<double, 1> &AQinit, Array<double, 1> &ionic,
                        Array<double, 1> &chp,
                        double &Temperature, double &RH, double &deltatmax)
{
  //Method use to solve the system dynamically:
  //for each compounds:
  // for a layer if characteristic time to reach equilibrium < tequilibrium concentration of the
  //   compound in the layer is at equilibrium with the gas phase
  //   if it is not the case, concentrations are computed dynamically
  //
  // - initialisation of concentrations at equilibrium
  // - for each time step:
  //    - computation of dynamic evolution
  //    - computation of new equilibrium
  int b,ilayer,iphase,i;
  int n=surrogate.size();
  double deltat1=config.deltatmin;
  double deltat2=deltat1;
  double t=0;
  Array<double,3> MO;
  //double error_aq=10.0;
  //double error_org=10.0;
  //double error_tot=10.0;
  MO.resize(config.nbins,config.nlayer,config.max_number_of_phases);  
  Array<double,1> AQ,conc_inorganic,ionic_organic, organion,MMaq,chp1,chp0,chp20;
  AQ.resize(config.nbins);
  ionic_organic.resize(config.nbins);
  conc_inorganic.resize(config.nbins);
  organion.resize(config.nbins);
  MMaq.resize(config.nbins);
  chp1.resize(config.nbins);
  chp0.resize(config.nbins);
  chp20.resize(config.nbins);
  int icycle;
  //config.MOmin=1.e-30;

  //for (b=0;b<config.nbins;b++)
  //  cout << number(b)*4./3*3.14159*pow(config.diameters(b)*0.5e-6,3)*1000*1.e9*1.e-5 << endl;
  /*
  if (config.imethod>=1)    
    config.tequilibrium=0;*/

  for (b=0;b<config.nbins;++b)		  
    for (ilayer=0;ilayer<config.nlayer;++ilayer)
      config.nphase(b,ilayer)=config.max_number_of_phases;

  if (config.compute_inorganic)
    {
      config.LWClimit=-1.0;
      if (config.compute_long_and_medium_range_interactions)
	config.compute_aqueous_phase_properties=true;
      else
	config.compute_aqueous_phase_properties=false;
    }
  else
    config.compute_aqueous_phase_properties=false;
  
  //double error_old,error_old2;
  double LWCtot=0.0;
  for (b=0;b<config.nbins;++b)
    LWCtot+=LWC(b);

  initialisation_ssh(config, surrogate, MOinit, MO, MOW, AQinit, AQ, MMaq, LWCtot, LWC, chp, ionic, ionic_organic,
                     organion, conc_inorganic, number, Vsol, Temperature, RH);
  chp0=chp;
  chp1=chp;
  
  compute_diameters_ssh(config, surrogate, Vsol, number, LWC, LWCtot);
  if (config.explicit_representation)
    compute_morphology_ssh(config, Vsol, number);

  //Computation of characteristic times to reach equilibrium
  tau_dif_ssh(config, surrogate, number, Vsol, Temperature);
  tau_kmt_ssh(config, surrogate, Temperature, number);
  compute_kp_org_ssh(config, surrogate, MOinit, Temperature, MOW);
  if (LWCtot>config.LWClimit)
    compute_kp_aq_ssh(config, surrogate, Temperature, ionic, chp, AQinit, MMaq, 0, config.nbins);

  characteristic_time_ssh(config, surrogate, MOinit, AQinit, LWCtot); 
  if (LWCtot>config.LWClimit)
    characteristic_time_aq_ssh(config, surrogate, Temperature, chp, LWC, AQinit, MOinit, MMaq, ionic);    
  
  if (config.imethod>=1)    
    {          
      if (config.compute_saturation and config.compute_organic)
	{
	  number_org_phases_ssh(config,surrogate,Temperature,MOinit,MOW);
	  tau_dif_ssh(config, surrogate, number, Vsol, Temperature);
	  tau_kmt_ssh(config, surrogate, Temperature, number);
	  compute_kp_org_ssh(config, surrogate, MOinit, Temperature, MOW);
	  characteristic_time_ssh(config, surrogate, MOinit, AQinit, LWCtot); 
	  if (LWCtot>config.LWClimit)
	    characteristic_time_aq_ssh(config, surrogate, Temperature, chp, LWC, AQinit, MOinit, MMaq, ionic);
	}

      /*
	if (config.imethod>=1)
	{
	solve_implicit_ssh(config, surrogate, MOinit, MOW, number, Vsol, LWC, AQinit, ionic, chp, Temperature, RH, AQ, MO,
	conc_inorganic, ionic_organic, organion, MMaq, t, 0);
	}*/
      
      deltat1=config.deltatmin;
      if (config.imethod>=2)
        deltat1=deltatmax;
      
      double atol=1.e-6*config.EPSER;
      double rtol=config.EPSER;
      double wk;          
      int index=config.max_iter;
      bool reject_step;
      //Array<double,1> diameters0;
      //diameters0.resize(config.nbins);      
      while (t<deltatmax)
        {
	  deltat1=min(deltatmax-t,deltat1);
	  deltat2=deltat1;
	  int index_save=index;
	  if (deltat1>config.deltatmin)
	    {
	      index=config.max_iter; //min(index*10,config.max_iter);
	    }
	  else
	    index=config.max_iter;

	  //diameters0=config.diameters;
	  chp0=chp;
	  solve_implicit_ssh(config, surrogate, MOinit, MOW, number, Vsol, LWC, AQinit, ionic, chp, Temperature, RH, AQ, MO,
                             conc_inorganic, ionic_organic, organion, MMaq, t, deltat1, index, reject_step, chp20);
	  
	  //cout << config.diameters << endl;
	  if (reject_step)
	    {
	      // cout << "rejected!!!" << endl;
	      index=index_save;
	      deltat1=max(deltat1/2,config.deltatmin);
	      //exit(0);
	      chp=chp0;
	    }
	  else
	    {	      
	      t+=deltat1;
	      
	      for (i=0;i<n;i++)
		if (surrogate[i].is_inorganic_precursor and surrogate[i].is_solid==false)
		  {
		    surrogate[i].Atot=surrogate[i].Ag;
		    if (i==config.iNH3)
		      {
			surrogate[i].Atot+=sum(surrogate[config.iNH4p].Aaq_bins_init)/surrogate[config.iNH4p].MM*surrogate[i].MM;
			surrogate[config.iNH4p].Atot=surrogate[i].Atot;
		      }
		    else if (i==config.iHNO3)
		      {
			surrogate[i].Atot+=sum(surrogate[config.iNO3m].Aaq_bins_init)/surrogate[config.iNO3m].MM*surrogate[i].MM;
			surrogate[config.iNO3m].Atot=surrogate[i].Atot;
		      }
		    else if (i==config.iHCl)
		      {
			surrogate[i].Atot+=sum(surrogate[config.iClm].Aaq_bins_init)/surrogate[config.iClm].MM*surrogate[i].MM;
			surrogate[config.iClm].Atot=surrogate[i].Atot;
		      }
		    else if (i==config.iH2SO4)
		      {
			surrogate[i].Atot+=sum(surrogate[config.iSO4mm].Aaq_bins_init)/surrogate[config.iSO4mm].MM*surrogate[i].MM+
			  sum(surrogate[config.iHSO4m].Aaq_bins_init)/surrogate[config.iHSO4m].MM*surrogate[i].MM;
			surrogate[config.iSO4mm].Atot=surrogate[i].Atot;
			surrogate[config.iHSO4m].Atot=surrogate[i].Atot;
		      }

		    surrogate[i].Atot0=surrogate[i].Ag0;
		    if (i==config.iNH3)
		      {
			surrogate[i].Atot0+=sum(surrogate[config.iNH4p].Aaq_bins_init0)/surrogate[config.iNH4p].MM*surrogate[i].MM;
			surrogate[config.iNH4p].Atot0=surrogate[i].Atot0;
		      }
		    else if (i==config.iHNO3)
		      {
			surrogate[i].Atot0+=sum(surrogate[config.iNO3m].Aaq_bins_init0)/surrogate[config.iNO3m].MM*surrogate[i].MM;
			surrogate[config.iNO3m].Atot0=surrogate[i].Atot0;
		      }
		    else if (i==config.iHCl)
		      {
			surrogate[i].Atot0+=sum(surrogate[config.iClm].Aaq_bins_init0)/surrogate[config.iClm].MM*surrogate[i].MM;
			surrogate[config.iClm].Atot0=surrogate[i].Atot0;
		      }
		    else if (i==config.iH2SO4)
		      {
			surrogate[i].Atot0+=sum(surrogate[config.iSO4mm].Aaq_bins_init0)/surrogate[config.iSO4mm].MM*surrogate[i].MM+
			  sum(surrogate[config.iHSO4m].Aaq_bins_init0)/surrogate[config.iHSO4m].MM*surrogate[i].MM;
			surrogate[config.iSO4mm].Atot0=surrogate[i].Atot0;
			surrogate[config.iHSO4m].Atot0=surrogate[i].Atot0;
		      }
		  }

	  
	      
	      double error_max=0.;
	      /*for (b=0;b<config.nbins;++b)
		{
		  error_max=max(error_max,abs(diameters0(b)-config.diameters(b))/diameters0(b)/0.1);
		}*/
	      if (config.EPSER_hp>0.)
		for (b=0;b<config.nbins;++b)
		  error_max=max(error_max,min(abs(chp0(b)-chp(b)),abs(chp20(b)-chp(b)))/chp(b)/config.EPSER_hp);
	      
	      for (i=0;i<n;i++)
		if (i!=config.iHp and i!=config.iHSO4m and i!=config.iCa and i!=config.iCO3mm and i!=config.iHCO3m and i!=config.iOHm)
		  {
		    if (surrogate[i].is_organic or surrogate[i].is_inorganic_precursor)
		      {
			wk=atol+rtol*surrogate[i].Ag;
			wk=max(wk,config.EPSER*1.e-5*surrogate[i].Atot);
			error_max=max(error_max,min(abs(surrogate[i].Ag-surrogate[i].Ag20),abs(surrogate[i].Ag-surrogate[i].Ag0))/wk);
			//error_max=max(error_max,abs(surrogate[i].Ag-surrogate[i].Ag0)/wk);
			//if ((surrogate[i].Ag-surrogate[i].Ag0)/wk>0.63)
			//   cout << surrogate[i].name << " gas " << error_max << endl;
		      }

		    if (surrogate[i].is_organic and surrogate[i].hydrophobic)
		      for (ilayer=0;ilayer<config.nlayer;++ilayer)
			for (b=0;b<config.nbins;++b)
			  for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
			    {
			      wk=atol+rtol*surrogate[i].Ap_layer_init(b,ilayer,iphase);
			      wk=max(wk,config.EPSER*1.e-5*surrogate[i].Atot);
			      error_max=max(error_max,min(abs(surrogate[i].Ap_layer_init(b,ilayer,iphase)-surrogate[i].Ap_layer_init0(b,ilayer,iphase)),
							  abs(surrogate[i].Ap_layer_init(b,ilayer,iphase)-surrogate[i].Ap_layer_init20(b,ilayer,iphase)))/wk);
			      //error_max=max(error_max,abs(surrogate[i].Ap_layer_init(b,ilayer,iphase)-surrogate[i].Ap_layer_init0(b,ilayer,iphase))/wk);
			      //if ((surrogate[i].Ap_layer_init(b,ilayer,iphase)-surrogate[i].Ap_layer_init0(b,ilayer,iphase))/wk>0.63)
			      //if (surrogate[i].name=="GLYOHOH")
			      //  cout << surrogate[i].name << " org " << error_max << " " << surrogate[i].Ap_layer_init(b,ilayer,iphase) << " " << surrogate[i].Ap_layer_init0(b,ilayer,iphase) << endl;
			    }

		    /*
		      if ((surrogate[i].is_organic and surrogate[i].hydrophilic) or (surrogate[i].is_ion or config.iHp!=i))
		      {
		      wk=atol+rtol*sum(surrogate[i].Aaq_bins_init);
		      error_max=max(error_max,(sum(surrogate[i].Aaq_bins_init)-sum(surrogate[i].Aaq_bins_init0))/wk);
		      }*/
	      
		    for (b=0;b<config.nbins;++b)
		      {
			double a1=surrogate[i].Aaq_bins_init(b);
			double a0=surrogate[i].Aaq_bins_init0(b);
			double a20=surrogate[i].Aaq_bins_init20(b);
			if (i==config.iSO4mm)
			  {
			    a1+=surrogate[config.iHSO4m].Aaq_bins_init(b)/surrogate[config.iHSO4m].MM*surrogate[config.iSO4mm].MM;
			    a0+=surrogate[config.iHSO4m].Aaq_bins_init20(b)/surrogate[config.iHSO4m].MM*surrogate[config.iSO4mm].MM;
			    a20+=surrogate[config.iHSO4m].Aaq_bins_init20(b)/surrogate[config.iHSO4m].MM*surrogate[config.iSO4mm].MM;
			  }
		  
			wk=atol+rtol*a1;
			wk=max(wk,config.EPSER*1.e-5*surrogate[i].Atot);
			error_max=max(error_max,min(abs(a1-a20),abs(a1-a0))/wk);
			//if ((surrogate[i].Aaq_bins_init(b)-surrogate[i].Aaq_bins_init0(b))/wk>0.63)
			//  cout << surrogate[i].name << " aq  " << error_max << endl;
		      }	
		  }

	      /*	      double error_max2;
	      for (b=0;b<config.nbins;++b)
		{
		  error_max2=min(abs(chp0(b)-chp(b)),abs(chp20(b)-chp(b)))/chp(b)/0.01;
		  //if (error_max2>0.9*error_max)
		  //  cout << b << " " << error_max2 << " " << min(abs(chp0(b)-chp(b)),abs(chp20(b)-chp(b)))/config.EPSER << " " << min(abs(chp0(b)-chp(b)),abs(chp20(b)-chp(b))) << " " << chp(b) << " " << chp0(b) << " " << chp20(b) << endl;
		  }*/
	      //cout << t << " " << error_max << " " << deltat1 << endl;
	      /*double error_max2=0.0;
	      deltat2=0.8/pow(error_max,0.5)*deltat1;
	      if (deltat2<0.5*deltat1 or deltat2<1.e-3)
		for (i=0;i<n;i++)
		  if (i!=config.iHp and i!=config.iHSO4m and i!=config.iCa and i!=config.iCO3mm and i!=config.iHCO3m and i!=config.iOHm)
		    {
		      if (surrogate[i].is_organic or surrogate[i].is_inorganic_precursor)
			{
			  wk=atol+rtol*surrogate[i].Ag;
			  wk=max(wk,config.EPSER*1.e-5*surrogate[i].Atot);
			  error_max2=abs(surrogate[i].Ag-surrogate[i].Ag0)/wk;
			  if (error_max2>0.9*error_max)
			    cout << surrogate[i].name << " gas " << error_max2 << " " << surrogate[i].Ag << " " << surrogate[i].Ag0 << endl;
			}

		      if (surrogate[i].is_organic and surrogate[i].hydrophobic)
			for (ilayer=0;ilayer<config.nlayer;++ilayer)
			  for (b=0;b<config.nbins;++b)
			    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
			      {
				wk=atol+rtol*surrogate[i].Ap_layer_init(b,ilayer,iphase);
				wk=max(wk,config.EPSER*1.e-5*surrogate[i].Atot);
				error_max2=abs(surrogate[i].Ap_layer_init(b,ilayer,iphase)-surrogate[i].Ap_layer_init0(b,ilayer,iphase))/wk;
				if (error_max2>0.9*error_max)
				  cout << surrogate[i].name << " org " << error_max2 << " " << surrogate[i].Ap_layer_init(b,ilayer,iphase) << " " << surrogate[i].Ap_layer_init0(b,ilayer,iphase) << endl;
			      }
	      
		      for (b=0;b<config.nbins;++b)
			{
			  double a1=surrogate[i].Aaq_bins_init(b);
			  double a0=surrogate[i].Aaq_bins_init0(b);
			  if (i==config.iSO4mm)
			    {
			      a1+=surrogate[config.iHSO4m].Aaq_bins_init(b)/surrogate[config.iHSO4m].MM*surrogate[config.iSO4mm].MM;
			      a0+=surrogate[config.iHSO4m].Aaq_bins_init0(b)/surrogate[config.iHSO4m].MM*surrogate[config.iSO4mm].MM;
			    }
		  
			  wk=rtol*a1;
			  wk=max(wk,config.EPSER*1.e-5*surrogate[i].Atot);
			  error_max2=abs(a1-a0)/wk;
			  if (error_max2>0.9*error_max)
			    cout << surrogate[i].name << " aq " << b << " " << error_max2 << " " << surrogate[i].Aaq_bins_init(b) << " " << surrogate[i].Aaq_bins_init0(b) << " " << config.EPSER*1.e-5*surrogate[i].Atot << endl;
		  
			}	
			}*/
	      /*
		for (i=0;i<n;i++)
		if ((surrogate[i].is_organic and surrogate[i].hydrophilic) or (surrogate[i].is_ion and i!=config.iHp))
		for (b=0;b<config.nbins;++b)
		{
		wk=atol+rtol*surrogate[i].Aaq_bins_init(b);
		if ((surrogate[i].Aaq_bins_init(b)-surrogate[i].Aaq_bins_init0(b))/wk>=0.9*error_max)
		cout << surrogate[i].name << " " << b << " " << surrogate[i].Aaq_bins_init(b) << " " << surrogate[i].Aaq_bins_init0(b) << endl;
		}*/
	  
	      //cout << error_max << endl;
	      if (error_max>0)
		{
		  deltat2=0.8/pow(error_max,0.5)*deltat1;
		  deltat1=max(min(deltat2,deltat1*10),deltat1*0.1);
		}
	      else
		deltat1=10*deltat1;

	      //cout << error_max << " " << deltat1 << endl;
	      deltat1=max(deltat1,config.deltatmin);
	  
	    }
	  for (i=0;i<n;i++)
	    {
	      surrogate[i].Aaq_bins=surrogate[i].Aaq_bins_init;
	      surrogate[i].Ap_layer=surrogate[i].Ap_layer_init;
	      /*
		if (surrogate[i].name=="SOAlP")
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
		cout << atot << " " << surrogate[i].Ag << " " << surrogate[i].Ap_layer_init << " " << surrogate[i].Atot << endl;
		}*/
	    }
	}
      /*
	for (i=0;i<n;i++)
	if (sum(surrogate[i].Aaq_bins_init)>0)
	cout << surrogate[i].name << " out " << surrogate[i].Aaq_bins_init << endl;
	else if (surrogate[i].Ag>0. and surrogate[i].is_inorganic_precursor)
	cout << surrogate[i].name << " out " << surrogate[i].Ag << endl;*/
    }
  else
    if (config.coupled_phases or 
	(RH>=config.RHcoupling and surrogate[config.iH2O].hydrophilic and surrogate[config.iH2O].hydrophobic))
      {
	// If the syst is coupled particulate hydrophilic compounds concentrations
	//must be computed simultaneously with the particulate hydrophobic compounds concentrations
	//cout << "The system is coupled." << endl;

	//Initialisation of equibrium      
	if (config.compute_saturation and config.compute_organic)
	  {
	    number_org_phases_ssh(config,surrogate,Temperature,MOinit,MOW);
	    compute_kp_org_ssh(config, surrogate, MOinit, Temperature, MOW);
	    tau_dif_ssh(config, surrogate, number, Vsol, Temperature);
	    tau_kmt_ssh(config, surrogate, Temperature, number);
	    
	    characteristic_time_ssh(config, surrogate, MOinit, AQinit, LWCtot); 
	    if (LWCtot>config.LWClimit)
	      characteristic_time_aq_ssh(config, surrogate, Temperature, chp, LWC, AQinit, MOinit, MMaq, ionic);
	  }
   
	if (config.coupling_organic_inorganic or config.compute_organic==false 
	    or config.compute_inorganic==false)
	  solve_local_equilibriums_coupled_ssh(config, surrogate, MOinit, MOW, number, Vsol, LWC, AQinit, ionic, chp, Temperature, RH, AQ, MO,
					       conc_inorganic, ionic_organic, organion, MMaq);
	else
	  for (icycle=0;icycle<config.number_of_org_inorg_cycles;icycle++)
	    {
	      config.compute_organic=false;
	      solve_local_equilibriums_coupled_ssh(config, surrogate, MOinit, MOW, number, Vsol, LWC, AQinit, ionic, chp, Temperature, RH, AQ, MO,
						   conc_inorganic, ionic_organic, organion, MMaq);
	      config.compute_inorganic=false; 
	      config.compute_organic=true;
	      solve_local_equilibriums_coupled_ssh(config, surrogate, MOinit, MOW, number, Vsol, LWC, AQinit, ionic, chp, Temperature, RH, AQ, MO,
						   conc_inorganic, ionic_organic, organion, MMaq);
	      config.compute_inorganic=true;
	    }

	characteristic_time_ssh(config, surrogate, MOinit, AQinit, LWCtot); 
	if (LWCtot>config.LWClimit)
	  characteristic_time_aq_ssh(config, surrogate, Temperature, chp, LWC, AQinit, MOinit, MMaq, ionic);

	//Dynamic evolution
	while (t<deltatmax)
	  {	  
	    deltat1=min(deltatmax-t,deltat1);
		  
	    //save the old time step in deltat2          
	    deltat2=deltat1;
	    //compute the dynamic evolution for dt=deltat1
	    double tequilibrium;
	    if (config.explicit_representation)
	      tequilibrium=0.0;
	    else
	      tequilibrium=config.tequilibrium;

	    config.to_be_rejected=false;
	    if (config.first_evaluation_activity_coefficients==false)
	      dynamic_tot_ssh(config,surrogate,MOinit,MO,MOW,AQinit,AQ,conc_inorganic,
			      ionic,ionic_organic,organion,chp,chp1,chp0,LWC,MMaq,Temperature,RH,
			      deltat1,tequilibrium,true);
	    else
	      dynamic_tot_ssh(config,surrogate,MOinit,MO,MOW,AQinit,AQ,conc_inorganic,
			      ionic,ionic_organic,organion,chp,chp1,chp0,LWC,MMaq,Temperature,RH,
			      deltat1,tequilibrium,false); 

	    //compute the new time step so that changes are small
	    adapstep_ssh(config,surrogate,Temperature,config.tequilibrium,deltat1,t,deltatmax,config.deltatmin,
			 MOinit,MO,LWCtot,AQinit,AQ,LWC,conc_inorganic,chp,chp1,chp0,number);	  
		  
	    if (deltat1<0.999*deltat2 or config.to_be_rejected) //if the new time step is inferior to the old one
	      {
		//cout << "rejected " << t << " " << deltat2 << " " << config.to_be_rejected << endl;
		//the old time step is rejected
		for (i=0;i<n;++i)
		  {
		    surrogate[i].Ag=surrogate[i].Ag0;

		    if (surrogate[i].hydrophobic)
		      for (ilayer=0;ilayer<config.nlayer;++ilayer)
			for (b=0;b<config.nbins;++b)
			  for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
			    {
			      surrogate[i].Ap_layer_init(b,ilayer,iphase)=
				surrogate[i].Ap_layer_init0(b,ilayer,iphase);
			      surrogate[i].Ap_layer(b,ilayer,iphase)=
				surrogate[i].Ap_layer_init(b,ilayer,iphase);
			    }

		    if (LWCtot>config.LWClimit)
		      if (surrogate[i].hydrophilic)
			for (b=0;b<config.nbins;++b)
			  {
			    surrogate[i].Aaq_bins_init(b)=surrogate[i].Aaq_bins_init0(b);		  
			    surrogate[i].Aaq_bins(b)=surrogate[i].Aaq_bins(b);	
			  }
		    surrogate[i].Atot=surrogate[i].Atot0;
		  }

		for (b=0;b<config.nbins;++b)
		  {
		    AQinit(b)=LWC(b);
		    for (i=0;i<n;++i)
		      if(surrogate[i].hydrophilic)
			AQinit(b)+=surrogate[i].Aaq_bins_init(b);
		    AQinit(b)=max(AQinit(b),config.MOmin);		      
		    chp(b)=chp0(b);
		  }

		for (b=0;b<config.nbins;++b)
		  for (ilayer=0;ilayer<config.nlayer;++ilayer)
		    {
		      for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
			{			   
			  MOinit(b,ilayer,iphase)=0.0;
			  for (i=0;i<n;++i)	
			    if(surrogate[i].hydrophobic)
			      MOinit(b,ilayer,iphase)+=surrogate[i].Ap_layer_init(b,ilayer,iphase);
			}
                    
		      double sum=0.0;
		      for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
			sum+=MOinit(b,ilayer,iphase);
                    
		      if (sum>0.0)
			for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
			  MOinit(b,ilayer,iphase)=max(MOinit(b,ilayer,iphase),
						      config.MOmin*config.Vlayer(ilayer)*MOinit(b,ilayer,iphase)/sum);
		      else
			for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
			  if(iphase==0)
			    MOinit(b,ilayer,iphase)=config.MOmin*config.Vlayer(ilayer);
			  else
			    MOinit(b,ilayer,iphase)=0.0;
		    }               
	      }
	    else                      //the time step is accepted 
	      {	
		t+=deltat2;	      	     
		for (b=0;b<config.nbins;++b)
		  {
		    for (ilayer=0;ilayer<config.nlayer;++ilayer)
		      for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
			for (i=0;i<n;++i)
			  if (surrogate[i].hydrophobic and
			      (surrogate[i].is_organic or i==config.iH2O))
			    {
			      surrogate[i].Ap_layer_init0(b,ilayer,iphase)=
				surrogate[i].Ap_layer_init(b,ilayer,iphase);
			      surrogate[i].Ap_layer_init(b,ilayer,iphase)=
				surrogate[i].Ap_layer(b,ilayer,iphase);
			    }

		    if (config.compute_inorganic)
		      {
			for (i=0;i<n;++i)
			  if (surrogate[i].hydrophilic)
			    surrogate[i].Aaq_bins_init(b)=surrogate[i].Aaq_bins(b);
		      }
		    else
		      for (i=0;i<n;++i)
			if (surrogate[i].hydrophilic and surrogate[i].is_organic)
			  surrogate[i].Aaq_bins_init(b)=surrogate[i].Aaq_bins(b);   
		  }

		//redistribute concentrations to ensure that the volume of layers are constant
		redistribution_ssh(config, surrogate,MOinit,MO);

		for (b=0;b<config.nbins;++b)
		  {
		    for (ilayer=0;ilayer<config.nlayer;++ilayer)
		      for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
			MOinit(b,ilayer,iphase)=MO(b,ilayer,iphase);
		    AQinit(b)=AQ(b);
		  }

		if (LWCtot>config.LWClimit)
		  density_aqueous_phase_ssh(config, surrogate, LWC, Temperature);

		//compute the new diameters of particles
		compute_diameters_ssh(config, surrogate, Vsol, number, LWC, LWCtot);
		if (config.explicit_representation)
		  compute_morphology_ssh(config, Vsol,number);

		//Computation of characteristic times to reach equilibrium
		tau_dif_ssh(config, surrogate, number, Vsol,Temperature);
		tau_kmt_ssh(config, surrogate, Temperature, number);
		characteristic_time_ssh(config, surrogate, MOinit, AQinit, LWCtot);
		if (LWCtot>config.LWClimit)
		  characteristic_time_aq_ssh(config, surrogate, Temperature, chp, LWC, AQinit, MOinit, MMaq, ionic);
			  
		for (b=0;b<config.nbins;++b)
		  {
		    for (ilayer=0;ilayer<config.nlayer;++ilayer)
		      for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
			MOinit(b,ilayer,iphase)=MO(b,ilayer,iphase);
		    AQinit(b)=AQ(b);
		  }

		//computation of concentrations at equilibrium	      
	      
		if (config.coupling_organic_inorganic or config.compute_organic==false 
		    or config.compute_inorganic==false)
		  solve_local_equilibriums_coupled_ssh(config, surrogate, MOinit, MOW, number, Vsol, LWC, AQinit, ionic, chp, Temperature, RH, AQ, MO,
						       conc_inorganic, ionic_organic, organion, MMaq);
		else
		  for (icycle=0;icycle<config.number_of_org_inorg_cycles;icycle++)
		    {
		      config.compute_organic=false;
		      solve_local_equilibriums_coupled_ssh(config, surrogate, MOinit, MOW, number, Vsol, LWC, AQinit, ionic, chp, Temperature, RH, AQ, MO,
							   conc_inorganic, ionic_organic, organion, MMaq);
		      config.compute_inorganic=false;
		      config.compute_organic=true;
		      solve_local_equilibriums_coupled_ssh(config, surrogate, MOinit, MOW, number, Vsol, LWC, AQinit, ionic, chp, Temperature, RH, AQ, MO,
							   conc_inorganic, ionic_organic, organion, MMaq);
		      config.compute_inorganic=true;
		    }
	      }
	  }
      }
    else
      {
      
	// If the syst is not coupled particulate hydrophilic compounds concentrations
	//can be computed independently of the particulate hydrophobic compounds concentrations
	//cout << "The system is not coupled." << endl;

	//Initialisation of equibrium
	//compute phase separation
	if (config.compute_saturation and config.compute_organic)
	  {
	    number_org_phases_ssh(config,surrogate,Temperature,MOinit,MOW);
	    compute_kp_org_ssh(config, surrogate, MOinit, Temperature, MOW);
	    tau_dif_ssh(config, surrogate, number, Vsol, Temperature);
	    tau_kmt_ssh(config, surrogate, Temperature, number);
	    
	    characteristic_time_ssh(config, surrogate, MOinit, AQinit, LWCtot); 
	    if (LWCtot>config.LWClimit)
	      characteristic_time_aq_ssh(config, surrogate, Temperature, chp, LWC, AQinit, MOinit, MMaq, ionic); 
	  }
      
	if (config.coupling_organic_inorganic or config.compute_organic==false 
	    or config.compute_inorganic==false)
	  solve_local_equilibriums_uncoupled_ssh(config, surrogate, MOinit, MOW, number, Vsol, LWC, AQinit, ionic, chp, Temperature, RH, AQ, MO,
						 conc_inorganic, ionic_organic, organion, MMaq);
	else
	  for (icycle=0;icycle<config.number_of_org_inorg_cycles;icycle++)
	    {
	      config.compute_organic=false;
	      solve_local_equilibriums_uncoupled_ssh(config, surrogate, MOinit, MOW, number, Vsol, LWC, AQinit, ionic, chp, Temperature, RH, AQ, MO,
						     conc_inorganic, ionic_organic, organion, MMaq);
	      config.compute_inorganic=false;
	      config.compute_organic=true;
	      solve_local_equilibriums_uncoupled_ssh(config, surrogate, MOinit, MOW, number, Vsol, LWC, AQinit, ionic, chp, Temperature, RH, AQ, MO,
						     conc_inorganic, ionic_organic, organion, MMaq);
	      config.compute_inorganic=true;
	    }
	  
	//Dynamic evolution
	while (t<deltatmax)
	  {
	    deltat1=min(deltatmax-t,deltat1);
	    deltat2=deltat1;
		  
	    //compute the dynamic evolution for dt=deltat1
	    /*
	    if (config.first_evaluation_activity_coefficients==false)
	      {
		dynamic_org_ssh(config,surrogate,MOinit,MO,AQinit,
				MOW,Temperature,deltat1,config.tequilibrium, true);             
		if (LWCtot>config.LWClimit)
		  dynamic_aq_ssh(config,surrogate,AQinit,AQ,MOinit,conc_inorganic,ionic,ionic_organic,
				 organion,chp,chp1,chp0,LWC,MMaq,MOW,Temperature,deltat1,config.tequilibrium, true);
	      }
	    else
	      {
		dynamic_org_ssh(config,surrogate,MOinit,MO,AQinit,
				MOW,Temperature,deltat1,config.tequilibrium, false);
		if (LWCtot>config.LWClimit)
		  dynamic_aq_ssh(config,surrogate,AQinit,AQ,MOinit,conc_inorganic,ionic,ionic_organic,
				 organion,chp,chp1,chp0,LWC,MMaq,MOW,Temperature,deltat1,config.tequilibrium, false);
	      }*/
	    if (config.first_evaluation_activity_coefficients==false)
	      dynamic_tot_ssh(config,surrogate,MOinit,MO,MOW,AQinit,AQ,conc_inorganic,
			      ionic,ionic_organic,organion,chp,chp1,chp0,LWC,MMaq,Temperature,RH,
			      deltat1,config.tequilibrium,true);
	    else
	      dynamic_tot_ssh(config,surrogate,MOinit,MO,MOW,AQinit,AQ,conc_inorganic,
			      ionic,ionic_organic,organion,chp,chp1,chp0,LWC,MMaq,Temperature,RH,
			      deltat1,config.tequilibrium,false); 
          
		  
	    //compute the new time step so that changes are small
	    adapstep_ssh(config,surrogate,Temperature,config.tequilibrium,deltat1,t,deltatmax,config.deltatmin,
			 MOinit,MO,LWCtot,AQinit,AQ,LWC,conc_inorganic,chp,chp1,chp0,number);
		  
	    if (deltat1<0.95*deltat2) //if the new time step is inferior to the old one
	      {             
		//the old time step is rejected
		for (i=0;i<n;++i)
		  {
		    surrogate[i].Ag=surrogate[i].Ag0;

		    if (surrogate[i].hydrophobic)
		      for (ilayer=0;ilayer<config.nlayer;++ilayer)
			for (b=0;b<config.nbins;++b)
			  for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
			    {
			      surrogate[i].Ap_layer_init(b,ilayer,iphase)=
				surrogate[i].Ap_layer_init0(b,ilayer,iphase);
			      surrogate[i].Ap_layer(b,ilayer,iphase)=
				surrogate[i].Ap_layer_init(b,ilayer,iphase);
			    }

		    if (LWCtot>config.LWClimit)
		      if (surrogate[i].hydrophilic)
			for (b=0;b<config.nbins;++b)
			  {
			    surrogate[i].Aaq_bins_init(b)=surrogate[i].Aaq_bins_init0(b);		  
			    surrogate[i].Aaq_bins(b)=surrogate[i].Aaq_bins(b);	
			  }
		  }

		for (b=0;b<config.nbins;++b)
		  {
		    AQinit(b)=LWC(b);
		    for (i=0;i<n;++i)
		      if(surrogate[i].hydrophilic)
			AQinit(b)+=surrogate[i].Aaq_bins_init(b);
		    AQinit(b)=max(AQinit(b),config.MOmin);		      
		    chp(b)=chp0(b);
		  }

		for (b=0;b<config.nbins;++b)
		  for (ilayer=0;ilayer<config.nlayer;++ilayer)
		    {
		      for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
			{			   
			  MOinit(b,ilayer,iphase)=0.0;
			  for (i=0;i<n;++i)	
			    if(surrogate[i].hydrophobic)
			      MOinit(b,ilayer,iphase)+=surrogate[i].Ap_layer_init(b,ilayer,iphase);
			}
                    
		      double sum=0.0;
		      for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
			sum+=MOinit(b,ilayer,iphase);
                    
		      if (sum>0.0)
			for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
			  MOinit(b,ilayer,iphase)=max(MOinit(b,ilayer,iphase),
						      config.MOmin*config.Vlayer(ilayer)*MOinit(b,ilayer,iphase)/sum);
		      else
			for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
			  if(iphase==0)
			    MOinit(b,ilayer,iphase)=config.MOmin*config.Vlayer(ilayer);
			  else
			    MOinit(b,ilayer,iphase)=0.0;
		    }                  		                  
              
	      }
	    else                     //the time step is acepted 
	      {
		t+=deltat2;
		//cout << "dt: " << deltat2 << endl;
		for (b=0;b<config.nbins;++b)
		  {
		    for (i=0;i<n;++i)
		      if (surrogate[i].hydrophobic)
			for (ilayer=0;ilayer<config.nlayer;++ilayer)
			  for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
			    surrogate[i].Ap_layer_init(b,ilayer,iphase)=
			      surrogate[i].Ap_layer(b,ilayer,iphase);

		    if (config.compute_inorganic)
		      {
			for (i=0;i<n;++i)
			  if (surrogate[i].hydrophilic)
			    surrogate[i].Aaq_bins_init(b)=surrogate[i].Aaq_bins(b);
		      }
		    else
		      if (LWCtot>config.LWClimit)
			for (i=0;i<n;++i)
			  if (surrogate[i].hydrophilic and surrogate[i].is_organic)
			    surrogate[i].Aaq_bins_init(b)=surrogate[i].Aaq_bins(b);
		  }

		//redistribute concentrations to ensure that the volume of layers are constant
		redistribution_ssh(config, surrogate,MOinit,MO);
		for (b=0;b<config.nbins;++b)
		  {
		    for (ilayer=0;ilayer<config.nlayer;++ilayer)
		      for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
			MOinit(b,ilayer,iphase)=MO(b,ilayer,iphase);
		    if (LWCtot>config.LWClimit)
		      AQinit(b)=AQ(b);
		  }

		//compute the new diameters of particles
		if (LWCtot>config.LWClimit)
		  density_aqueous_phase_ssh(config, surrogate, LWC, Temperature);

		compute_diameters_ssh(config, surrogate, Vsol, number, LWC, LWCtot);
		if (config.explicit_representation)
		  compute_morphology_ssh(config, Vsol, number);

		//Computation of characteristic times to reach equilibrium
		tau_dif_ssh(config, surrogate, number, Vsol, Temperature);
		tau_kmt_ssh(config, surrogate, Temperature, number);
		characteristic_time_ssh(config, surrogate, MOinit, AQinit, LWCtot);
		if (LWCtot>config.LWClimit)
		  characteristic_time_aq_ssh(config, surrogate, Temperature, chp, LWC, AQinit, MOinit, MMaq, ionic);
              
              
		//computation of concentrations at equilibrium	      
		if (config.coupling_organic_inorganic or config.compute_organic==false 
		    or config.compute_inorganic==false)
		  solve_local_equilibriums_uncoupled_ssh(config, surrogate, MOinit, MOW, number, Vsol, LWC, AQinit, ionic, chp, Temperature, RH, AQ, MO,
							 conc_inorganic, ionic_organic, organion, MMaq);
		else
		  for (icycle=0;icycle<config.number_of_org_inorg_cycles;icycle++)
		    {
		      config.compute_organic=false;
		      solve_local_equilibriums_uncoupled_ssh(config, surrogate, MOinit, MOW, number, Vsol, LWC, AQinit, ionic, chp, Temperature, RH, AQ, MO,
							     conc_inorganic, ionic_organic, organion, MMaq);
		      config.compute_inorganic=false;
		      config.compute_organic=true;
		      solve_local_equilibriums_uncoupled_ssh(config, surrogate, MOinit, MOW, number, Vsol, LWC, AQinit, ionic, chp, Temperature, RH, AQ, MO,
							     conc_inorganic, ionic_organic, organion, MMaq);
		      config.compute_inorganic=true;
		    }              
	      } 
	  }
      }
}

