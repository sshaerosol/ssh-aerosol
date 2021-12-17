

//!!-----------------------------------------------------------------------
//!!     Copyright (C) 2019 CEREA (ENPC) - INERIS
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
  double MMaq;
  int i; //,it;
  int n=surrogate.size();
  for (i=0;i<n;++i)
    if (surrogate[i].is_organic==false and i!=config.iH2O)
      conc_inorganic+=surrogate[i].Aaq;
  bool compute_activity_coefficients=true;
    
  //if (config.compute_inorganic)
  //  config.coupled_phases=true;

  //cout << RH << " " << config.RHcoupling << " " << surrogate[config.iH2O].hydrophilic << " " <<  surrogate[config.iH2O].hydrophobic << endl;

  
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
              config.first_evaluation_activity_coefficients=true;                            
              //config.compute_organic=false; 

              double RHsave=RH;
              //RH=max(RH,0.2);
              
              while ((index_iter < config.max_iter) and (abs(error1)/factor_old > config.precision or abs(error2)/factor_old > config.precision or
                                                         ((abs(error3)/factor_old>relprec  or abs(error4)/factor_old>relprec) and AQ>100*config.MOmin) or
                                                         config.first_evaluation_activity_coefficients==true or RH>RHsave))
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

		  error4=0.0;
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
			    surrogate[i].Aaq=surrogate[i].Aaq_old;                           
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
	      while ((index_iter < max_iter2) and ((abs(error1)/factor > config.precision)
						   or (abs(error2)/factor > config.precision)))
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

		      error4=0.0;
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
          while ((index_iter < config.max_iter) and (abs(error3) > config.precision))
            {
              if (config.first_evaluation_activity_coefficients==false or index_iter==0)
                compute_activity_coefficients=true;
              else
                compute_activity_coefficients=false;
              
              error_org_ssh(config,surrogate,MO,MOW,Temperature,error3,derivative,RH,1.0,
                            all_hydrophobic, compute_activity_coefficients);
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
	  while ((index_iter < max_iter2) and (abs(error)/factor > config.precision))
	    {
	      if (config.first_evaluation_activity_coefficients==false or index_iter==0)
		compute_activity_coefficients=true;
	      else
		compute_activity_coefficients=false;
	      
	      error_org_ssh(config,surrogate,MO,MOW,Temperature,error,derivative,RH,factor,
			    all_hydrophobic, compute_activity_coefficients);
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

		  error4=0.0;
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

		  error4=0.0;
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
	      //int nmin=0;
	      //double AQmin=AQ;
	      while ((index_iter < max_iter2) and (abs(error2)/factor > config.precision))
		{
		  if (config.first_evaluation_activity_coefficients==false or index_iter==0)
		    compute_activity_coefficients=true;
		  else
		    compute_activity_coefficients=false;

		  error_aq_ssh(config,surrogate,AQ,LWC,conc_inorganic,ionic,chp,MMaq,
			       Temperature,error2,derivative,RH,
			       organion,ionic_organic,factor,compute_activity_coefficients);

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

		      error4=0.0;
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
	  // minimize the function error with a method of newton raphson
	  // error3 = MOinit - sum of concentrations of organic compounds in the organic phase
	  //           - hygroscopicity of the organic phase  
	  while ((index_iter < config.max_iter) and (abs(error3) > config.precision))
	    {
	      if (config.first_evaluation_activity_coefficients==false or index_iter==0)
		compute_activity_coefficients=true;
	      else
		compute_activity_coefficients=false;
              
	      error_org_ssh(config,surrogate,MO,MOW,Temperature,error3,derivative,RH,1.0,
			    all_hydrophobic, compute_activity_coefficients);
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
    initialisation_eq_ssh(config,surrogate,Temperature,RH,ionic, chp, AQinit,false);
  else
    initialisation_eq_ssh(config,surrogate,Temperature,RH,ionic,chp,AQinit,true);   
  
  if (config.chemistry)
    for (it=0;it<config.nt;it++)
      { 
        solve_equilibrium_ssh(config, surrogate, MOinit,MOW, LWC, AQinit, ionic, chp,
                              Temperature, RH);
        solve_chemistry_ssh(config, surrogate, MOinit,MOW, LWC, AQinit, ionic, chp,
                            Temperature, RH, deltat/config.nt, compute_activity_coefficients);
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
                    i!=config.iHp)
                  error_aq=max(error_aq,abs(surrogate[i].Aaq_bins_init(b)
                                            -surrogate[i].Aaq_bins(b))
                               /surrogate[i].Aaq_bins(b));

            }
        }

      vec_error_org(index)=error_org;
      vec_error_aq(index)=error_aq;
		  
      //Computation of characteristic times to reach equilibrium
      tau_dif_ssh(config, surrogate, number, Vsol);
      tau_kmt_ssh(config, surrogate, Temperature, number);
      /*      if (config.explicit_representation)
              compute_morphology_ssh(config, Vsol);*/

      characteristic_time_ssh(config, surrogate, MOinit, AQinit, LWCtot); 
      if (LWCtot>config.LWClimit)
        characteristic_time_aq_ssh(config, surrogate, Temperature, chp, LWC, AQinit, MOinit);  

      ++index;
    }


  if (config.compute_saturation and config.first_evaluation_of_saturation==false and config.compute_organic)
    {
      number_org_phases_ssh(config,surrogate,Temperature,MOinit,MOW);
      tau_dif_ssh(config, surrogate, number, Vsol);
      tau_kmt_ssh(config, surrogate, Temperature, number);
      compute_kp_org_ssh(config, surrogate, MOinit, Temperature, MOW);     
      characteristic_time_ssh(config, surrogate, MOinit, AQinit, LWCtot);
      if (LWCtot>config.LWClimit)
        characteristic_time_aq_ssh(config, surrogate, Temperature, chp, LWC, AQinit, MOinit);
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
            if (surrogate[i].hydrophilic and i!=config.iHp)
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
      tau_dif_ssh(config, surrogate, number, Vsol);
      tau_kmt_ssh(config, surrogate, Temperature, number);
      characteristic_time_ssh(config, surrogate, MOinit, AQinit, LWCtot); 
      if (LWCtot>config.LWClimit)
        characteristic_time_aq_ssh(config, surrogate, Temperature, chp, LWC, AQinit, MOinit);      

      ++index;  
    }

  //cout << "niter: " << index << " " << chp(0) << " " << chp_save(0) << endl; 
  
  if (config.compute_saturation and config.first_evaluation_of_saturation==false and config.compute_organic)
    {
      number_org_phases_ssh(config,surrogate,Temperature,MOinit,MOW);
      tau_dif_ssh(config, surrogate, number, Vsol);
      tau_kmt_ssh(config, surrogate, Temperature, number);
      compute_kp_org_ssh(config, surrogate, MOinit, Temperature, MOW);
      characteristic_time_ssh(config, surrogate, MOinit, AQinit, LWCtot);
      if (LWCtot>config.LWClimit)
        characteristic_time_aq_ssh(config, surrogate, Temperature, chp, LWC, AQinit, MOinit);
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



void solve_semi_implicit_coupled_ssh(model_config config, vector<species> &surrogate,
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
  Array<double, 1> vec_error_org,vec_error_aq,vec_error_gas; 
  vec_error_org.resize(config.max_iter);
  vec_error_aq.resize(config.max_iter);
  vec_error_gas.resize(config.max_iter);
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
  config.nh_max=10;
  
  double factor=pow(0.5,nh-1);

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

  while (error_tot>config.relative_precision*factor and index < config.max_iter) 
    {
      if (index>2)
        {
          //ensure that the system can converge
          non_convergence=false;   
          if (iiter>20)
            for (i=max(index-20,2);i<index-1;i++)
              if (((abs(vec_error_org(i)-vec_error_org(index-1))/vec_error_org(index-1)/factor<1.0e-4 and vec_error_org(index-1)/factor>config.precision) and
                   (abs(vec_error_org(i-1)-vec_error_org(index-2))/vec_error_org(index-2)/factor<1.0e-4) and (abs(vec_error_org(i-2)-vec_error_org(index-3))/vec_error_org(index-3)/factor<1.0e-4))
                  or (abs(vec_error_aq(i)-vec_error_aq(index-1))/vec_error_aq(index-1)/factor<1.0e-4 and vec_error_aq(index-1)/factor>config.precision and 
                      (abs(vec_error_aq(i-1)-vec_error_aq(index-2))/vec_error_aq(index-2)/factor<1.0e-4) and (abs(vec_error_aq(i-2)-vec_error_aq(index-3))/vec_error_aq(index-3)/factor<1.0e-4))
		  or (abs(vec_error_gas(i)-vec_error_gas(index-1))/vec_error_gas(index-1)/factor<1.0e-4 and vec_error_gas(index-1)/factor>config.precision and 
                      (abs(vec_error_gas(i-1)-vec_error_gas(index-2))/vec_error_gas(index-2)/factor<1.0e-4) and (abs(vec_error_gas(i-2)-vec_error_gas(index-3))/vec_error_gas(index-3)/factor<1.0e-4)) 
                  or vec_error_org(index-1)/factor>10.0 or vec_error_aq(index-1)/factor>10.0 or vec_error_gas(index-1)/factor>10.0)
                non_convergence=true;

          if (iiter>100)
            if (vec_error_aq(index-1)/factor>config.precision)
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
            if (vec_error_org(index-1)/factor>config.precision)
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
	  
          if (iiter>100)
            if (vec_error_gas(index-1)/factor>config.precision)
              {
                double a=1;
                double b=1;
                for (i=index-100;i<index-90;i++)
                  a*=vec_error_gas(i);
                
                for (i=index-11;i<index-1;i++)
                  b*=vec_error_gas(i);
                  
                if (b>a*0.999)
                  non_convergence=true;
              }
              
          if (non_convergence and nh<config.nh_max)
            { 
              if (vec_error_org(index-1)/factor>100.0 or vec_error_aq(index-1)/factor>100.0 or vec_error_gas(index-1)/factor>100.0)
                {
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
                    }
                }
              ++nh;
	      factor=pow(0.5,nh-1);
              for (i=max(index-50,0);i<index;i++)
                {
                  vec_error_org(i)=-1.0;
                  vec_error_aq(i)=-1.0;
		  vec_error_gas(i)=-1.0;
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

      for (i=0;i<n;i++)
	surrogate[i].Ag1=surrogate[i].Ag;

      for (i=0;i<n;i++)
	surrogate[i].Atot1=surrogate[i].Atot;

      //cout << "NO3: " << sum(surrogate[config.iNO3m].Aaq_bins_init)/surrogate[config.iNO3m].MM*surrogate[config.iHNO3].MM+surrogate[config.iHNO3].Ag << " " << nh << endl;;
      //cout << "NO3i: " << sum(surrogate[config.iNO3m].Aaq_bins_init0)/surrogate[config.iNO3m].MM*surrogate[config.iHNO3].MM << " " << surrogate[config.iHNO3].Ag0 << " " << nh << endl;

      iiter++;   
      if (LWCtot>config.LWClimit)
        density_aqueous_phase_ssh(config, surrogate, LWC, Temperature);
          
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
      
      tau_dif_ssh(config, surrogate, number, Vsol);
      tau_kmt_ssh(config, surrogate, Temperature, number);
      

      //compute the new diameters of particle due to the growth of particles by condensation
      compute_diameters_ssh(config, surrogate, Vsol, number, LWC, LWCtot);
      if (config.explicit_representation)
	compute_morphology_ssh(config, Vsol, number);
		  
      //Computation of error_tot
      vec_error_org(index)=0.0;
      vec_error_aq(index)=0.0;
      vec_error_gas(index)=0.0;

      for (i=0;i<n;i++)
	if (surrogate[i].is_organic or (surrogate[i].is_inorganic_precursor and surrogate[i].is_solid==false))
	  if (surrogate[i].Atot1>1.0e-10)
	    vec_error_gas(index)=max(vec_error_gas(index),abs(surrogate[i].Atot-surrogate[i].Atot1)
				     /surrogate[i].Atot1);	
      
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
		    {
		      /*
			if (abs(surrogate[i].Ap_layer_init(b,ilayer,iphase)
			-surrogate[i].Ap_layer(b,ilayer,iphase))
			/surrogate[i].Ap_layer(b,ilayer,iphase)>vec_error_org(index))
			cout << surrogate[i].name << " " << abs(surrogate[i].Ap_layer_init(b,ilayer,iphase)
			-surrogate[i].Ap_layer(b,ilayer,iphase))
			/surrogate[i].Ap_layer(b,ilayer,iphase) << endl;*/
		      vec_error_org(index)=max(vec_error_org(index),
					       abs(surrogate[i].Ap_layer_init(b,ilayer,iphase)
						   -surrogate[i].Ap_layer(b,ilayer,iphase))
					       /surrogate[i].Ap_layer(b,ilayer,iphase));
		    }

                MOinit(b,ilayer,iphase)=max(MO(b,ilayer,iphase),config.MOmin*config.Vlayer(ilayer));
              }

          if (AQ(b)>1.0e-5)
            vec_error_aq(index)=max(vec_error_aq(index),abs(AQ(b)-AQinit(b))/AQ(b));

	  if (chp(b)>0.)
	    {
	      /*if (abs(chp(b)-chp_save(b))/chp(b)>0.2)
		cout << "error ph " << abs(chp(b)-chp_save(b))/chp(b) << endl;*/
	      vec_error_aq(index)=max(vec_error_aq(index),abs(chp(b)-chp_save(b))/chp(b));
	    }
       
          for (i=0;i<n;i++)
            if (surrogate[i].hydrophilic and i!=config.iHp)
	      if (surrogate[i].Aaq_bins(b)>1.0e-5)
		{
		  vec_error_aq(index)=max(vec_error_aq(index),abs(surrogate[i].Aaq_bins_init(b)
								  -surrogate[i].Aaq_bins(b))
					  /surrogate[i].Aaq_bins(b));
		  /*
		    if (abs(surrogate[i].Aaq_bins_init(b)
		    -surrogate[i].Aaq_bins(b))
		    /surrogate[i].Aaq_bins(b)>1.5)
		    cout << surrogate[i].name << " " << abs(surrogate[i].Aaq_bins_init(b)
		    -surrogate[i].Aaq_bins(b))
		    /surrogate[i].Aaq_bins(b) <<endl;*/
		    
		}
			  
          AQinit(b)=AQ(b); //max(AQ(b),config.MOmin);

	  
        }

      vec_error_org(index)=max(vec_error_org(index),0.001*config.precision*factor);
      vec_error_aq(index)=max(vec_error_aq(index),0.001*config.precision*factor);
      vec_error_gas(index)=max(vec_error_gas(index),0.001*config.precision*factor);
      
      error_tot=max(max(vec_error_org(index),vec_error_aq(index)),vec_error_gas(index));

      /*
	cout << error_tot << " " << sum(surrogate[config.iNO3m].Aaq_bins_init) << " " << vec_error_aq(index) << " " << factor << endl;
	cout << chp << endl;
	cout << sum(surrogate[config.iNO3m].Aaq_bins_init) << " "<< sum(surrogate[config.iNH4p].Aaq_bins_init) << " " << sum(surrogate[config.iHSO4m].Aaq_bins_init) << " " << sum(surrogate[config.iH2O].Aaq_bins_init) << endl;
      
	for (i=0;i<n;i++)
	if (surrogate[i].hydrophilic and i!=config.iHp)
	for (b=0;b<config.nbins;b++)
	if (surrogate[i].Aaq_bins(b)>1.0e-5)
	if (abs(surrogate[i].Aaq_bins_init(b)-surrogate[i].Aaq_bins(b))/surrogate[i].Aaq_bins(b)>=error_tot)
	cout << "error : " << surrogate[i].name << " " << b << " " << surrogate[i].Aaq_bins_init(b) << " " << surrogate[i].Aaq_bins(b) << endl;*/

      /*
	for (b=0;b<config.nbins;++b)
	if (abs(chp(b)-chp_save(b))/chp(b)>=error_tot)
	cout << "error_ph " << b << endl;*/
      
      //cout << config.gamma_MR_ions << endl;
      //Computation of characteristic times to reach equilibrium

      /*
	characteristic_time_ssh(config, surrogate, MOinit, AQinit, LWCtot); 
	if (LWCtot>config.LWClimit)
	characteristic_time_aq_ssh(config, surrogate, Temperature, chp, LWC, AQinit, MOinit);      */
      ++index;  
    }

  //cout << "niter: " << index << " " << chp(0) << " " << chp_save(0) << endl; 
  
  if (config.compute_saturation and config.first_evaluation_of_saturation==false and config.compute_organic)
    {
      number_org_phases_ssh(config,surrogate,Temperature,MOinit,MOW);
      tau_dif_ssh(config, surrogate, number, Vsol);
      tau_kmt_ssh(config, surrogate, Temperature, number);
      compute_kp_org_ssh(config, surrogate, MOinit, Temperature, MOW);
      characteristic_time_ssh(config, surrogate, MOinit, AQinit, LWCtot);
      if (LWCtot>config.LWClimit)
        characteristic_time_aq_ssh(config, surrogate, Temperature, chp, LWC, AQinit, MOinit);
    }

  //cout << index << " " << error_tot << endl;
  if (error_tot>config.relative_precision*factor)
    {
      cout << "The model did not converged..." << endl;
      cout << "coupled" << endl;
      for (b=0;b<config.nbins;b++)
	cout << b << " " << abs(chp(b)-chp_save(b))/chp(b)/factor << endl;
      cout << config.nphase << " " << nh << endl;
      cout << surrogate[config.iH2O].time_aq << endl;
      cout << "error" << error_tot << endl;
      cout << "avant" << MOinit << endl;
      cout << "avant" << AQinit << endl;
      cout << "apres" << MO << endl;
      cout << "apres" << AQ << endl;
      for (i=0;i<n;i++)
	if (surrogate[i].Atot>0.) 
	  cout << surrogate[i].name << " " << sum(surrogate[i].Ap_layer_init) << " " << endl;
      cout << chp << endl;
      cout << t << " " << deltat << " " << config.max_iter << " " << factor << " " << index << endl;
      //exit(0);
    }
  
}
void initialisation_ssh(model_config &config, vector<species> &surrogate,
                        Array<double,3> &MOinit, Array<double, 3> &MO, Array<double, 3> &MOW,
                        Array<double,1> &AQinit,Array<double,1> &AQ, Array<double,1> &MMaq,
                        double &LWCtot, Array<double,1> &LWC, Array<double, 1> &chp, Array<double, 1> &ionic,
                        Array<double, 1> &ionic_organic, Array<double, 1> &organion, Array<double, 1> &conc_inorganic,
                        double &Temperature, double &RH)
{  
  int b,ilayer,iphase,i;
  int n=surrogate.size();
  Array <double, 1> conc_org;
  conc_org.resize(config.nbins);
  double MOWloc=1.;

  /*  if (config.compute_inorganic)
      for (b=0;b<config.nbins;b++)
      if (chp(b)<=1.e-8)
      chp(b)=1.e-2;*/

  for (i=0;i<n;i++)
    {
      if (surrogate[i].hydrophilic==false)
	surrogate[i].Aaq_bins_init=0;
      if (surrogate[i].hydrophobic==false)
	surrogate[i].Ap_layer_init=0;
      surrogate[i].velocity=pow(2.11714271498563e4*Temperature/surrogate[i].MM,0.5);
      surrogate[i].knui=pow(3.0*surrogate[i].KDiffusion_air/surrogate[i].velocity,2);
      for (b=0;b<config.nbins;b++)
	surrogate[i].fac_corr_ph(b)=1.;
    }

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
		  surrogate[i].veckaqi(b)=surrogate[i].Kpart_aq_ssh(Temperature,MOWloc);              
		else
		  if (surrogate[i].nonvolatile)
		    surrogate[i].veckaqi(b)=1000.;
		  else
		    {
		      if (surrogate[i].aqt==2) //diacid
			{
			  surrogate[i].veckaqi(b)=surrogate[i].Kpart_aq_ssh(Temperature,MOWloc)*
			    (1.0+surrogate[i].Kacidity1/(pow(gamma,2)*chp(b))*
			     (1.0+surrogate[i].Kacidity2/(pow(gamma,2)*chp(b))));
			  surrogate[i].vecfioni1(b)=(surrogate[i].Kacidity1/(pow(gamma,2)*chp(b)))/
			    (1.0+surrogate[i].Kacidity1/(pow(gamma,2)*chp(b))*(1.0+surrogate[i].Kacidity2/(pow(gamma,2)*chp(b))));
			  surrogate[i].vecfioni2(b)=(surrogate[i].Kacidity1/(pow(gamma,2)*chp(b)))*(surrogate[i].Kacidity2/(pow(gamma,2)*chp(b)))/
			    (1.0+surrogate[i].Kacidity1/(pow(gamma,2)*chp(b))*(1.0+surrogate[i].Kacidity2/(pow(gamma,2)*chp(b))));
			}
		      else if (surrogate[i].aqt==1) //monoacid
			{
			  surrogate[i].veckaqi(b)=surrogate[i].Kpart_aq_ssh(Temperature,MOWloc)*(1.0+surrogate[i].Kacidity1/(pow(gamma,2)*chp(b)));
			  surrogate[i].vecfioni1(b)=(surrogate[i].Kacidity1/(pow(gamma,2)*chp(b)))/(1.0+surrogate[i].Kacidity1/(pow(gamma,2)*chp(b)));
			}
		      else if (surrogate[i].aqt==3) //aldehyde
			surrogate[i].veckaqi(b)=surrogate[i].Kpart_aq_ssh(Temperature,MOWloc)*(1.0+surrogate[i].Koligo_aq*pow(gamma*chp(b)/pow(10,-surrogate[i].pHref),surrogate[i].beta));
		      else
			surrogate[i].veckaqi(b)=surrogate[i].Kpart_aq_ssh(Temperature,MOWloc);
		    }
	      }
      }

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
  else
    config.compute_aqueous_phase_properties=false;

  for (b=0;b<config.nbins;++b)		  
    for (ilayer=0;ilayer<config.nlayer;++ilayer)
      for (iphase=0;iphase<config.max_number_of_phases;++iphase)
        MO(b,ilayer,iphase)=0.0;

  //initialisation of some parameters
  for (b=0;b<config.nbins;++b)
    {
      AQinit(b)=LWC(b);
      for (i=0;i<n;++i)
        if (surrogate[i].hydrophilic)
          AQinit(b)+=surrogate[i].Aaq_bins_init(b);
	 
      AQinit(b)=max(AQinit(b),config.MOmin);
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
                if (surrogate[i].is_organic==false and surrogate[i].is_inorganic_precursor==false and i!=config.iHp and i!=config.iH2O)
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
	if (MOinit(b,ilayer,iphase)>0.0)
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
	else
	  {
	    //Remplace concentrations of organic aerosols equal to zero by MOmin
	    //(to ensure formation of organic aerosols)
	    MOinit(b,ilayer,iphase)=max(MOinit(b,ilayer,iphase),
					config.MOmin*config.Vlayer(ilayer));
	    for (i=0;i<n;++i)
	      surrogate[i].gamma_org_layer(b,ilayer,iphase)=1.0;
	  }
  
  for (i=0;i<n;++i)
    if (surrogate[i].hydrophilic or surrogate[i].is_inorganic_precursor)
      for (b=0;b<config.nbins;++b)
	surrogate[i].time_aq(b)=2.0*config.tequilibrium;

  //Computation of Ag for water
  water_concentration_ssh(config, surrogate, Temperature, RH);

  if (config.compute_inorganic)
    for (b=0;b<config.nbins;++b)
      {
	double error=100.0;
	int index=0;
	double temp;
	double factor=1.0;
	double XH2O=0.0;
	while (error>config.relative_precision and index<10)
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
                                        organion(b), ionic_organic(b), conc_org(b), factor);

	    /*
		  
	      activity_coefficients_aq_ssh(config,surrogate,Temperature,0.0,MMaq(b),XH2O,conc_org(b));
	      if (config.compute_long_and_medium_range_interactions)
	      activity_coefficients_LR_MR_ssh(config, surrogate, Temperature, 0.0, ionic(b));
	      if (index==0)
              surrogate[config.iH2O].gamma_aq_bins(b)=1.0;
	      else
	      surrogate[config.iH2O].gamma_aq_bins(b)=surrogate[config.iH2O].gamma_aq;*/
	    
	    activity_coefficients_dyn_aq_ssh(config, surrogate, Temperature,AQinit,MOinit,
					     conc_inorganic, ionic, ionic_organic,
					     organion,chp,LWC,MMaq,0.0,0.,0);	    
            
	    compute_kp_aq_ssh(config, surrogate, Temperature, ionic, chp, MMaq);
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
		  if (surrogate[i].is_organic==false and surrogate[i].is_inorganic_precursor==false and i!=config.iHp and i!=config.iH2O)
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

  /*
    for (i=0;i<n;i++)
    if (surrogate[i].Ag+sum(surrogate[i].Aaq_bins_init)>0.)
    cout << surrogate[i].name << " init2 " << surrogate[i].Ag << " " << sum(surrogate[i].Aaq_bins_init) << endl;  
    cout << chp << endl;
    cout << LWC << endl;*/
  
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
  Array<double,1> AQ,conc_inorganic,ionic_organic, organion,MMaq,chp1,chp0;
  AQ.resize(config.nbins);
  ionic_organic.resize(config.nbins);
  conc_inorganic.resize(config.nbins);
  organion.resize(config.nbins);
  MMaq.resize(config.nbins);
  chp1.resize(config.nbins);
  chp0.resize(config.nbins);
  int icycle;
  /*
    cout << "SO4: " << surrogate[config.iH2SO4].Ag+(sum(surrogate[config.iSO4mm].Aaq_bins_init)/surrogate[config.iSO4mm].MM+sum(surrogate[config.iHSO4m].Aaq_bins_init)/surrogate[config.iHSO4m].MM)*surrogate[config.iH2SO4].MM << endl;
    cout << "NO3: " << surrogate[config.iHNO3].Ag+sum(surrogate[config.iNO3m].Aaq_bins_init)/surrogate[config.iNO3m].MM*surrogate[config.iHNO3].MM << endl;
    cout << "NH4: " << surrogate[config.iNH3].Ag+sum(surrogate[config.iNH4p].Aaq_bins_init)/surrogate[config.iNH4p].MM*surrogate[config.iNH3].MM << endl;
    cout << "PNO3: " << sum(surrogate[config.iNO3m].Aaq_bins_init)/surrogate[config.iNO3m].MM*surrogate[config.iHNO3].MM << endl;
    cout << "PNH4: " << sum(surrogate[config.iNH4p].Aaq_bins_init)/surrogate[config.iNH4p].MM*surrogate[config.iNH3].MM << endl;*/
  if (config.imethod>=1)
    config.tequilibrium=0;

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
                     organion, conc_inorganic, Temperature, RH);

  compute_diameters_ssh(config, surrogate, Vsol, number, LWC, LWCtot);
  if (config.explicit_representation)
    compute_morphology_ssh(config, Vsol, number);

  //Computation of characteristic times to reach equilibrium
  tau_dif_ssh(config, surrogate, number, Vsol);
  tau_kmt_ssh(config, surrogate, Temperature, number);
  compute_kp_org_ssh(config, surrogate, MOinit, Temperature, MOW);
  if (LWCtot>config.LWClimit)
    compute_kp_aq_ssh(config, surrogate, Temperature, ionic, chp, MMaq);

  characteristic_time_ssh(config, surrogate, MOinit, AQinit, LWCtot); 
  if (LWCtot>config.LWClimit)
    characteristic_time_aq_ssh(config, surrogate, Temperature, chp, LWC, AQinit, MOinit);

  if (config.imethod>=1)    
    {     
      if (config.imethod==1)
	{
	  if (config.compute_saturation and config.compute_organic)
	    {
	      number_org_phases_ssh(config,surrogate,Temperature,MOinit,MOW);
	      tau_dif_ssh(config, surrogate, number, Vsol);
	      tau_kmt_ssh(config, surrogate, Temperature, number);
	      compute_kp_org_ssh(config, surrogate, MOinit, Temperature, MOW);
	      characteristic_time_ssh(config, surrogate, MOinit, AQinit, LWCtot); 
	      if (LWCtot>config.LWClimit)
		characteristic_time_aq_ssh(config, surrogate, Temperature, chp, LWC, AQinit, MOinit);
	    }
      
	  solve_semi_implicit_coupled_ssh(config, surrogate, MOinit, MOW, number, Vsol, LWC, AQinit, ionic, chp, Temperature, RH, AQ, MO,
					  conc_inorganic, ionic_organic, organion, MMaq, t, 0);
	}
      /*
	for (i=0;i<n;i++)
	if (sum(surrogate[i].Aaq_bins_init)>0)
	cout << surrogate[i].name << " in " << surrogate[i].Aaq_bins_init << endl;
	else if (surrogate[i].Ag>0. and surrogate[i].is_inorganic_precursor)
	cout << surrogate[i].name << " in " << surrogate[i].Ag << endl;*/
      
      deltat1=config.deltatmin;
      if (config.imethod==2)
	deltat1=deltatmax;
      
      double atol=1.e-6*config.EPSER;
      double rtol=config.EPSER;
      double wk;          
      
      while (t<deltatmax)
        {
	  //cout << "evol: " << t << " " << deltat1 << endl;
	  deltat1=min(deltatmax-t,deltat1);
	  deltat2=deltat1;
	  solve_semi_implicit_coupled_ssh(config, surrogate, MOinit, MOW, number, Vsol, LWC, AQinit, ionic, chp, Temperature, RH, AQ, MO,
					  conc_inorganic, ionic_organic, organion, MMaq, t, deltat1);
	  t+=deltat1;

	  double error_max=0.;
	  for (i=0;i<n;i++)
	    {
	      if (surrogate[i].is_organic or surrogate[i].is_inorganic_precursor)
		{
		  wk=atol+rtol*surrogate[i].Ag;
		  error_max=max(error_max,(surrogate[i].Ag-surrogate[i].Ag0)/wk);
		}

	      if (surrogate[i].is_organic and surrogate[i].hydrophobic)
		for (ilayer=0;ilayer<config.nlayer;++ilayer)
		  for (b=0;b<config.nbins;++b)
		    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
		      {
			wk=atol+rtol*surrogate[i].Ap_layer_init(b,ilayer,iphase);
			error_max=max(error_max,(surrogate[i].Ap_layer_init(b,ilayer,iphase)-surrogate[i].Ap_layer_init0(b,ilayer,iphase))/wk);
		      }

	      /*
		if ((surrogate[i].is_organic and surrogate[i].hydrophilic) or (surrogate[i].is_ion or config.iHp!=i))
		{
		wk=atol+rtol*sum(surrogate[i].Aaq_bins_init);
		error_max=max(error_max,(sum(surrogate[i].Aaq_bins_init)-sum(surrogate[i].Aaq_bins_init0))/wk);
		}*/
	      
	      for (b=0;b<config.nbins;++b)
		{
		  wk=atol+rtol*surrogate[i].Aaq_bins_init(b);
		  error_max=max(error_max,(surrogate[i].Aaq_bins_init(b)-surrogate[i].Aaq_bins_init0(b))/wk);
		}
	    }

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
	  deltat1=max(deltat1,config.deltatmin);
	}
      for (i=0;i<n;i++)
	{
	  surrogate[i].Aaq_bins=surrogate[i].Aaq_bins_init;
	  surrogate[i].Ap_layer=surrogate[i].Ap_layer_init;
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
	      //cout << atot << " " << surrogate[i].Ag << " " << surrogate[i].Ap_layer_init << " " << surrogate[i].Atot << endl;
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
	    tau_dif_ssh(config, surrogate, number, Vsol);
	    tau_kmt_ssh(config, surrogate, Temperature, number);
	    compute_kp_org_ssh(config, surrogate, MOinit, Temperature, MOW);
	    characteristic_time_ssh(config, surrogate, MOinit, AQinit, LWCtot); 
	    if (LWCtot>config.LWClimit)
	      characteristic_time_aq_ssh(config, surrogate, Temperature, chp, LWC, AQinit, MOinit);
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
	  characteristic_time_aq_ssh(config, surrogate, Temperature, chp, LWC, AQinit, MOinit);

	//Dynamic evolution
	while (t<deltatmax)
	  {	  
	    deltat1=min(deltatmax-t,deltat1);
	    //cout << "evol: " << t << " " << deltat1 << endl; //" " << chp(0) << " " << surrogate[config.iNO3m].Aaq_bins_init(0) << " " << surrogate[config.iH2O].Aaq_bins_init(0) << endl;
	    /*for (i=0;i<n;i++)
	      if (surrogate[i].Ag+sum(surrogate[i].Aaq_bins_init)+sum(surrogate[i].Ap_layer_init)>0.)
	      cout << surrogate[i].name << " " << surrogate[i].Ag << " " << sum(surrogate[i].Aaq_bins_init) << " " << sum(surrogate[i].Ap_layer_init) << endl;*/
		  
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
			      ionic,ionic_organic,organion,chp,chp1,chp0,LWC,MMaq,Temperature,
			      deltat1,tequilibrium,true);
	    else
	      dynamic_tot_ssh(config,surrogate,MOinit,MO,MOW,AQinit,AQ,conc_inorganic,
			      ionic,ionic_organic,organion,chp,chp1,chp0,LWC,MMaq,Temperature,
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
		tau_dif_ssh(config, surrogate, number, Vsol);
		tau_kmt_ssh(config, surrogate, Temperature, number);
		characteristic_time_ssh(config, surrogate, MOinit, AQinit, LWCtot);
		if (LWCtot>config.LWClimit)
		  characteristic_time_aq_ssh(config, surrogate, Temperature, chp, LWC, AQinit, MOinit);
			  
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
	    tau_dif_ssh(config, surrogate, number, Vsol);
	    tau_kmt_ssh(config, surrogate, Temperature, number);
	    compute_kp_org_ssh(config, surrogate, MOinit, Temperature, MOW);
	    characteristic_time_ssh(config, surrogate, MOinit, AQinit, LWCtot); 
	    if (LWCtot>config.LWClimit)
	      characteristic_time_aq_ssh(config, surrogate, Temperature, chp, LWC, AQinit, MOinit); 
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
	      }
          
		  
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
		tau_dif_ssh(config, surrogate, number, Vsol);
		tau_kmt_ssh(config, surrogate, Temperature, number);
		characteristic_time_ssh(config, surrogate, MOinit, AQinit, LWCtot);
		if (LWCtot>config.LWClimit)
		  characteristic_time_aq_ssh(config, surrogate, Temperature, chp, LWC, AQinit, MOinit);
              
              
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

  /*cout << "out: " << chp << " " << LWC << " " << AQinit << endl;  
    for (i=0;i<n;i++)
    cout << surrogate[i].name << " " << sum(surrogate[i].Aaq_bins_init) << " " << surrogate[i].Ag << endl;*/
  /*
    cout << chp << endl;
    cout << LWC << endl;
    for (i=0;i<n;i++)
    if (surrogate[i].Ag+sum(surrogate[i].Aaq_bins_init)>0.)
    cout << surrogate[i].name << " out " << surrogate[i].Ag << " " << sum(surrogate[i].Aaq_bins_init) << " " << sum(surrogate[i].gamma_aq_bins) << endl;
    cout << config.gamma_MR_ions << endl;*/
  /*
    cout << "out: " << endl;
    cout << "SO4: " << surrogate[config.iH2SO4].Ag+(sum(surrogate[config.iSO4mm].Aaq_bins_init)/surrogate[config.iSO4mm].MM+sum(surrogate[config.iHSO4m].Aaq_bins_init)/surrogate[config.iHSO4m].MM)*surrogate[config.iH2SO4].MM << endl;
    cout << "NO3: " << surrogate[config.iHNO3].Ag+sum(surrogate[config.iNO3m].Aaq_bins_init)/surrogate[config.iNO3m].MM*surrogate[config.iHNO3].MM << endl;
    cout << "NH4: " << surrogate[config.iNH3].Ag+sum(surrogate[config.iNH4p].Aaq_bins_init)/surrogate[config.iNH4p].MM*surrogate[config.iNH3].MM << endl;
    cout << "PNO3: " << sum(surrogate[config.iNO3m].Aaq_bins_init)/surrogate[config.iNO3m].MM*surrogate[config.iHNO3].MM << endl;
    cout << "PNH4: " << sum(surrogate[config.iNH4p].Aaq_bins_init)/surrogate[config.iNH4p].MM*surrogate[config.iNH3].MM << endl;*/
  
}

