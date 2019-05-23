

#include "properties.cxx"
#include "equilibrium.cxx"
#include "dynamic.cxx"
using namespace soap;

void solve_system(model_config &config, vector<species>& surrogate,
                  double &MOinit,double &MOW,
                  double &LWC, double &AQinit, double &ionic, double &chp,
                  double &Temperature, double &RH)
{
  double organion=0.0;
  double conc_inorganic=0.0;
  double ionic_organic=0.0;
  double MMaq;
  int i;
  int n=surrogate.size();
  for (i=0;i<n;++i)
    if (surrogate[i].is_organic==false and i!=config.iH2O)
      conc_inorganic+=surrogate[i].Aaq;

  bool compute_activity_coefficients=true;

  if (LWC>config.LWClimit)
    initialisation_eq(config,surrogate,Temperature,ionic, chp, false);
  else
    initialisation_eq(config,surrogate,Temperature,ionic,chp,true);
  
  if (config.coupled_phases or 
      (RH>=config.RHcoupling and surrogate[config.iH2O].hydrophilic and surrogate[config.iH2O].hydrophobic))
    {
      // If the syst is coupled particulate hydrophilic compounds concentrations
      //must be computed simultaneously with the particulate hydrophobic compounds concentrations
      if (LWC>config.LWClimit)
        {
          bool all_hydrophobic=false;
          double error1=1000.0;
          double error2=1000.0;
          double deriv_error1_MO,deriv_error1_AQ,deriv_error2_MO,deriv_error2_AQ;
          double derivative=0.0;
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
              Array<double, 1> vec_error_org,vec_error_aq,vec_error_chp,vec_error_compaq;
              vec_error_aq.resize(config.max_iter);	      
              vec_error_org.resize(config.max_iter);
	      vec_error_compaq.resize(config.max_iter);
	      vec_error_chp.resize(config.max_iter);
              bool non_convergence;
              int iiter=0;
              double chp_old=chp;
	      double error3=1000.0;
	      double error4=1000.0;
	      double relprec=1.0e-3;

              while ((index_iter < config.max_iter) and (abs(error2)*nh > config.precision*MO or abs(error1)*nh > config.precision*AQ or abs(error3)*nh>relprec or abs(error4)*nh>relprec))
                {

                  if (config.first_evaluation_activity_coefficients==false or index_iter==0)
                    compute_activity_coefficients=true;
                  else
                    compute_activity_coefficients=false;

		  for (i=0;i<n;++i)
		    {
		      if (surrogate[i].hydrophilic)
			surrogate[i].Aaq_old=surrogate[i].Aaq;
		      if (surrogate[i].hydrophobic)
			surrogate[i].Ap_old=surrogate[i].Ap;
		    }

                  if (iiter>20)
                    {
                      non_convergence=false;
                      for (i=max(index_iter-20,0);i<index_iter-1;i++)
                        if ((abs((vec_error_org(i)-abs(error1))/error1)*nh<1.0e-3 and abs(error1)*nh>config.precision) or
                            (abs((vec_error_aq(i)-abs(error2))/error2)*nh<1.0e-3 and abs(error2)*nh>config.precision) or
			    (abs((vec_error_chp(i)-abs(error3))/error3)*nh<1.0e-3 and abs(error3)*nh>relprec) or
			    (abs((vec_error_compaq(i)-abs(error4))/error4)*nh<1.0e-3 and abs(error4)*nh>relprec) or
			    error1>1.0 or error2>1.0 or error3>1.0 or error4>1.0)
                          non_convergence=true;

                      if (non_convergence and nh<config.nh_max)
                        {
                          iiter=0;
                          ++nh;
                          for (i=max(index_iter-10,0);i<index_iter;i++)
                            {
                              vec_error_org(i)=-1.0;
                              vec_error_aq(i)=-1.0;
			      vec_error_chp(i)=-1.0;
			      vec_error_compaq(i)=-1.0;
                            }
                        }
                    }
                  iiter++;
		  error3=abs(chp-chp_old)/chp_old;
                  chp_old=chp;
		  for (i=0;i<n;i++)
		    if (surrogate[i].hydrophilic)
		      surrogate[i].Aaq_old=surrogate[i].Aaq;

                  error_coupled(config,surrogate,MO,MOW,MMaq,AQ,LWC,conc_inorganic,ionic,ionic_organic,
                                chp, organion,
                                Temperature,RH,
                                error1,deriv_error1_MO,deriv_error1_AQ,
                                error2,deriv_error2_MO,deriv_error2_AQ,1.0/nh, compute_activity_coefficients);
                          
                  AQ=AQ-error2;
                  MO=MO-error1;
		  error4=0.0;
		  for (i=0;i<n;i++)
		    if (surrogate[i].hydrophilic)
		      if (surrogate[i].Aaq_old>0.0)
			error4=max(error4,abs(surrogate[i].Aaq-surrogate[i].Aaq_old)/surrogate[i].Aaq_old);

                  vec_error_org(index_iter)=abs(error1);
                  vec_error_aq(index_iter)=abs(error2);
		  vec_error_chp(index_iter)=abs(error3);
		  vec_error_compaq(index_iter)=abs(error4);
                  ++index_iter;
                }
            }
          else
            while ((index_iter < config.max_iter) and ((abs(error1) > config.precision)
                                                       or (abs(error2) > config.precision)))
              {
                if (config.first_evaluation_activity_coefficients==false or index_iter==0)
                  compute_activity_coefficients=true;
                else
                  compute_activity_coefficients=false;
                
                error_coupled(config,surrogate,MO,MOW,MMaq,AQ,LWC,conc_inorganic,ionic,ionic_organic,
                              chp, organion,
                              Temperature,RH,
                              error1,deriv_error1_MO,deriv_error1_AQ,
                              error2,deriv_error2_MO,deriv_error2_AQ,1.0, compute_activity_coefficients);
                //solve the system with a method of newton raphson
                newton_raphson_coupled(MO,AQ,error1,deriv_error1_MO,deriv_error1_AQ,
                                       error2,deriv_error2_MO,deriv_error2_AQ); 
                ++index_iter;
              }

          if (config.compute_saturation and MO > 0.0 and config.compute_organic)
            saturation(config,surrogate,all_hydrophobic,LWC,ionic,conc_inorganic,ionic_organic,
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
              
              error_org(config,surrogate,MO,MOW,Temperature,error3,derivative,RH,
                        all_hydrophobic, compute_activity_coefficients);
              //solve the system with a method of newton raphson
              if (derivative != 0.0 and MO-error3/derivative >= 0.0)
                MO=MO-error3/derivative;
              else
                MO=MO-error3;

              ++index_iter;
            }
       		  
          for (i=0;i<n;++i)
	    if (surrogate[i].is_organic)
	      surrogate[i].Aaq=0.0;

          if (config.compute_saturation and MO > 0.0 and config.compute_organic)
            saturation(config,surrogate,all_hydrophobic,LWC,ionic,conc_inorganic,ionic_organic,
                       chp,organion,Temperature,RH);
        }
    }
  else
    {
      // If the syst is not coupled particulate hydrophilic compounds concentrations
      //can be computed independently of the particulate hydrophobic compounds concentrations
      bool all_hydrophobic;
      double MO=MOinit;
      if (LWC>config.LWClimit)
        {
          all_hydrophobic=false;
          double error=1000.0;
          double derivative=0.0;
          int index_iter=0;
		  
          // minimize the function error with a method of newton raphson
          // error = MOinit - sum of concentrations of organic compounds in the organic phase
          //           - hygroscopicity of the organic phase 
          while ((index_iter < config.max_iter) and (abs(error) > config.precision))
            {
              if (config.first_evaluation_activity_coefficients==false or index_iter==0)
                compute_activity_coefficients=true;
              else
                compute_activity_coefficients=false;

              error_org(config,surrogate,MO,MOW,Temperature,error,derivative,RH,
                        all_hydrophobic, compute_activity_coefficients);
              //solve the system with a method of newton raphson
              if (derivative != 0.0 and MO-error/derivative >= 0.0)
                MO=MO-error/derivative;
              else
                MO=MO-error;

              ++index_iter;
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
              Array<double, 1> vec_error_aq,vec_error_chp,vec_error_compaq;
              vec_error_aq.resize(config.max_iter);
	      vec_error_chp.resize(config.max_iter);
	      vec_error_compaq.resize(config.max_iter);
              bool non_convergence;
              int iiter=0;
	      double chp_old=chp;
              double error3=1000.0;
	      double error4=1000.0;
              while ((index_iter < config.max_iter) and (abs(error2)*nh > config.precision or abs(error3)*nh>1.0e-4 or abs(error4)*nh>1.0e-4))
                {
                  if (config.first_evaluation_activity_coefficients==false or index_iter==0)
                    compute_activity_coefficients=true;
                  else
                    compute_activity_coefficients=false;

                  if (iiter>20)
                    {
                      non_convergence=false;
                      for (i=max(index_iter-10,0);i<index_iter-1;i++)
                        if ((abs((vec_error_aq(i)-abs(error2))/error2)*nh<1.0e-3 and abs(error2)*nh>config.precision) or 
			    (abs((vec_error_chp(i)-abs(error3))/error3)*nh<1.0e-3 and abs(error3)*nh>1.0e-4) or
			    (abs((vec_error_compaq(i)-abs(error4))/error4)*nh<1.0e-3 and abs(error4)*nh>1.0e-4) or
			    error2>1.0 or error3>1.0)
                          non_convergence=true;

                      if (non_convergence and nh<config.nh_max)
                        {
                          ++nh;
                          iiter=0;
                          for (i=max(index_iter-10,0);i<index_iter;i++)
			    {
			      vec_error_aq(i)=-1.0;
			      vec_error_chp(i)=-1.0;
			      vec_error_compaq(i)=-1.0;
			    }
                        }
                    }
                  iiter++;

		  error3=abs(chp-chp_old)/chp_old;
		  for (i=0;i<n;i++)
		    if (surrogate[i].hydrophilic)
		      surrogate[i].Aaq_old=surrogate[i].Aaq;

                  chp_old=chp;		  
		  double deriv_error1_MO,deriv_error1_AQ,deriv_error2_MO,deriv_error2_AQ,error1;
		  
                  error_aq(config,surrogate,AQ,LWC,conc_inorganic,ionic,chp,MMaq,
                           Temperature,error2,derivative,RH,
                           organion,ionic_organic,1.0/nh, compute_activity_coefficients);
                          
                  AQ=AQ-error2;

                  vec_error_aq(index_iter)=abs(error2);
		  vec_error_chp(index_iter)=abs(error3);
		  error4=0.0;
		  for (i=0;i<n;i++)
		    if (surrogate[i].hydrophilic)
		      if (surrogate[i].Aaq_old>0.0)
			error4=max(error4,abs(surrogate[i].Aaq-surrogate[i].Aaq_old)/surrogate[i].Aaq_old);

		  vec_error_compaq(index_iter)=abs(error4);
                  ++index_iter;
                }	      
            }
          else
            while ((index_iter < config.max_iter) and (abs(error2) > config.precision))
              {
                if (config.first_evaluation_activity_coefficients==false or index_iter==0)
                  compute_activity_coefficients=true;
                else
                  compute_activity_coefficients=false;

                error_aq(config,surrogate,AQ,LWC,conc_inorganic,ionic,chp,MMaq,
                         Temperature,error2,derivative,RH,
                         organion,ionic_organic,1.0,compute_activity_coefficients);

                //solve the system with a method of newton raphson
                if (derivative != 0.0 and MO-error2/derivative >= 0.0)
                  AQ=AQ-error2/derivative;
                else
                  AQ=AQ-error2;

                ++index_iter;
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
              
              error_org(config,surrogate,MO,MOW,Temperature,error3,derivative,RH,
                        all_hydrophobic, compute_activity_coefficients);
              if (derivative != 0.0 and MO-error3/derivative >= 0.0)
                MO=MO-error3/derivative;
              else
                MO=MO-error3;

              ++index_iter;
            }
          
          for (i=0;i<n;++i)
	    if (surrogate[i].is_organic)
	      surrogate[i].Aaq=0.0;
        }
	  
      if (config.compute_saturation and MO > 0.0 and config.compute_organic)
        saturation(config,surrogate,all_hydrophobic,LWC,ionic,conc_inorganic,ionic_organic,
                   chp,organion,Temperature,RH);
    }

  for(i=0;i<n;i++)
    if (surrogate[i].is_organic)
      surrogate[i].Ag=max(surrogate[i].Atot-surrogate[i].Aaq
                          -surrogate[i].Ap, 0.0);
       
      
}

void global_equilibrium(model_config &config, vector<species>& surrogate,
                        double &MOinit,double &MOW,
                        double &LWC, double &AQinit, double &ionic, double &chp,
                        double &Temperature, double &RH)
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
    solve_system(config, surrogate, MOinit, MOW, LWC, AQinit, ionic, 
                 chp, Temperature, RH);
  else
    for (icycle=0;icycle<config.number_of_org_inorg_cycles;icycle++)
      {
        config.compute_organic=false;
        solve_system(config, surrogate, MOinit, MOW, LWC, AQinit, ionic, 
                     chp, Temperature, RH);
        config.compute_inorganic=false;
	

        config.compute_organic=true;
        solve_system(config, surrogate, MOinit, MOW, LWC, AQinit, ionic, 
                     chp, Temperature, RH);
        config.compute_inorganic=true;
      }
}

void solve_local_equilibriums_uncoupled(model_config config, vector<species> &surrogate,
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
              if ((vec_error_org(index-1) > 0.0) and
                  (vec_error_org(index-2) > 0.0) and
                  (vec_error_org(index-3) > 0.0))
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
              if (vec_error_aq(index-1) > 0.0 and
                  vec_error_aq(index-2) > 0.0 and 
                  vec_error_aq(index-3) > 0.0)
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
          if (config.first_evaluation_activity_coefficients==false)
            {
              equilibrium_org(config, surrogate, config.tequilibrium, MOinit,MO,
                              Temperature, MOW, true,1.0/nh_org);
              if (config.compute_saturation and config.compute_organic)
                phase_repartition(config,surrogate,Temperature,MOinit,MO,MOW,1.0/nh_org);
            }
          else
            {
              equilibrium_org(config, surrogate, config.tequilibrium, MOinit,MO,
                              Temperature, MOW, false,1.0/nh_org);
              if (config.compute_saturation and config.compute_organic)
                phase_repartition(config,surrogate,Temperature,MOinit,MO,MOW,1.0/nh_org);
            }

          //redistribute concentrations to ensure that the volume of layers are constant
          redistribution(config, surrogate,MOinit,MO);
        }


      if (error_aq>1.0/nh_aq*config.relative_precision and LWCtot>config.LWClimit)
        density_aqueous_phase(config, surrogate, LWC, Temperature);
		  
      //if the system has not converged for aqueous concentrations
      if (error_aq>1.0/nh_aq*config.relative_precision and LWCtot>config.LWClimit)
        if (config.first_evaluation_activity_coefficients==false)
          equilibrium_aq(config, surrogate, config.tequilibrium, AQinit, AQ,
                         MOinit,conc_inorganic, ionic, ionic_organic, organion,
                         chp, LWC, Temperature, RH, MMaq, true,1.0/nh_aq);
        else
          equilibrium_aq(config, surrogate, config.tequilibrium, AQinit, AQ,
                         MOinit,conc_inorganic, ionic, ionic_organic, organion,
                         chp, LWC, Temperature, RH, MMaq, false,1.0/nh_aq);

      water_concentration(config, surrogate, Temperature, RH);

      //compute the new diameters of particles
      compute_diameters(config, surrogate, Vsol, number, LWC, LWCtot);		  		  
      if (config.explicit_representation)
	compute_morphology(config, Vsol, number);
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
                      surrogate[i].Ap_layer(b,ilayer,iphase)>1.0e-5)
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
                if (surrogate[i].hydrophilic and surrogate[i].Aaq_bins(b)>1.0e-5 and i!=config.iHp)
                  error_aq=max(error_aq,abs(surrogate[i].Aaq_bins_init(b)
                                            -surrogate[i].Aaq_bins(b))
                               /surrogate[i].Aaq_bins(b));

            }
        }

      vec_error_org(index)=error_org;
      vec_error_aq(index)=error_aq;
		  
      //Computation of characteristic times to reach equilibrium
      tau_dif(config, surrogate, number, Vsol);
      tau_kmt(config, surrogate, Temperature, number);

      characteristic_time(config, surrogate, MOinit, AQinit, LWCtot); 
      if (LWCtot>config.LWClimit)
        characteristic_time_aq(config, surrogate, LWC, AQinit, MOinit);

      ++index;

    }

  if (config.compute_saturation and config.first_evaluation_of_saturation==false and config.compute_organic)
    {
      number_org_phases(config,surrogate,Temperature,MOinit,MOW);
      tau_dif(config, surrogate, number, Vsol);
      tau_kmt(config, surrogate, Temperature, number);
      compute_kp_org(config, surrogate, MOinit, Temperature, MOW);     
      characteristic_time(config, surrogate, MOinit, AQinit, LWCtot);
      if (LWCtot>config.LWClimit)
        characteristic_time_aq(config, surrogate, LWC, AQinit, MOinit);
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

void solve_local_equilibriums_coupled(model_config config, vector<species> &surrogate,
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
              if (vec_error_org(index-1) > 0.0 and
                  vec_error_org(index-2) > 0.0 and
                  vec_error_org(index-3) > 0.0 and
                  vec_error_aq(index-1) > 0.0 and
                  vec_error_aq(index-2) > 0.0 and
                  vec_error_aq(index-3) > 0.0)
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
        density_aqueous_phase(config, surrogate, LWC, Temperature);
          
      for (b=0;b<config.nbins;b++)
	chp_save(b)=chp(b);

      if (config.first_evaluation_activity_coefficients==false)
        {	  
          equilibrium_tot(config, surrogate, config.tequilibrium, AQinit, AQ, conc_inorganic,
                          ionic, ionic_organic, organion, chp, LWC, MOinit, MO,
                          MOW, Temperature, RH, MMaq, true, 1.0/nh);       
          if (config.compute_saturation and config.compute_organic)
            phase_repartition(config,surrogate,Temperature,MOinit,MO,MOW,1.0/nh);	  
        }
      else
        {
          equilibrium_tot(config, surrogate, config.tequilibrium, AQinit, AQ, conc_inorganic,
                          ionic, ionic_organic, organion, chp, LWC, MOinit, MO,
                          MOW, Temperature, RH, MMaq, false, 1.0/nh);
          if (config.compute_saturation and config.compute_organic)
            phase_repartition(config,surrogate,Temperature,MOinit,MO,MOW,1.0/nh);
	}

      //redistribute concentrations to ensure that the volume of layers are constant      
      redistribution(config, surrogate,MOinit,MO);

      water_concentration(config, surrogate, Temperature, RH);


      //compute the new diameters of particle due to the growth of particles by condensation
      compute_diameters(config, surrogate, Vsol, number, LWC, LWCtot);
      if (config.explicit_representation)
	compute_morphology(config, Vsol, number);
		  
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

          // Error in the aqueous phase
          if (AQ(b)>1.0e-5)
            vec_error_aq(index)=max(vec_error_aq(index),abs(AQ(b)-AQinit(b))/AQ(b));
	  vec_error_aq(index)=max(vec_error_aq(index),abs(chp(b)-chp_save(b))/chp(b));

          for (i=0;i<n;i++)
            if (surrogate[i].hydrophilic and surrogate[i].Aaq_bins(b)>1.0e-5 and i!=config.iHp)
              vec_error_aq(index)=max(vec_error_aq(index),abs(surrogate[i].Aaq_bins_init(b)
                                                              -surrogate[i].Aaq_bins(b))
                                      /surrogate[i].Aaq_bins(b));
          AQinit(b)=max(AQ(b),config.MOmin);
        }
      error_tot=max(vec_error_org(index),vec_error_aq(index));

      //Computation of characteristic times to reach equilibrium
      tau_dif(config, surrogate, number, Vsol);
      tau_kmt(config, surrogate, Temperature, number);
      characteristic_time(config, surrogate, MOinit, AQinit, LWCtot); 
      if (LWCtot>config.LWClimit)
        characteristic_time_aq(config, surrogate, LWC, AQinit, MOinit);

      ++index;  
    }
  
  if (config.compute_saturation and config.first_evaluation_of_saturation==false and config.compute_organic)
    {
      number_org_phases(config,surrogate,Temperature,MOinit,MOW);
      tau_dif(config, surrogate, number, Vsol);
      tau_kmt(config, surrogate, Temperature, number);
      compute_kp_org(config, surrogate, MOinit, Temperature, MOW);
      characteristic_time(config, surrogate, MOinit, AQinit, LWCtot);
      if (LWCtot>config.LWClimit)
        characteristic_time_aq(config, surrogate, LWC, AQinit, MOinit);
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
    }
}

void initialisation(model_config &config, vector<species> &surrogate,
		    Array<double,3> &MOinit, Array<double, 3> &MO, Array<double, 3> &MOW,
		    Array<double,1> &AQinit, Array<double,1> &MMaq,
		    double &LWCtot, Array<double,1> &LWC, Array<double, 1> &chp, Array<double, 1> &ionic,
		    Array<double, 1> &ionic_organic, Array<double, 1> &organion, Array<double, 1> &conc_inorganic,
		    double &Temperature, double &RH)
{  
  int b,ilayer,iphase,i;
  int n=surrogate.size();
  
  for (b=0;b<config.nbins;++b)		  
    for (ilayer=0;ilayer<config.nlayer;++ilayer)
       config.nphase(b,ilayer)=1;

  if (config.compute_inorganic)
    {      
      config.LWClimit=-1.0;
      if (config.compute_long_and_medium_range_interactions)
	config.compute_aqueous_phase_properties=true;
      else
	config.compute_aqueous_phase_properties=false;

      for (i=0;i<n;i++)
	if (surrogate[i].is_inorganic_precursor)
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
 
 if (LWCtot>config.LWClimit)
   density_aqueous_phase(config, surrogate, LWC, Temperature);
 
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
	     double Ke=1.0e-14;
	     inorg1=0.0;
              for (i=0;i<n;++i) 
                if (surrogate[i].is_organic==false and surrogate[i].is_inorganic_precursor==false and i!=config.iHp and i!=config.iH2O)
                  {
                    inorg1-=surrogate[i].Aaq_bins_init(b)/surrogate[i].MM*surrogate[i].charge*config.AQrho(b)/AQinit(b);
                  }
              chp(b)=0.5*(inorg1+pow(pow(inorg1,2)+4*Ke,0.5));
	      if (chp(b)==0.0)
		chp(b)=1.0e-7;
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
                    activity_coefficients_dyn_sat(config, surrogate, Temperature, MOW, b, ilayer);
                    temp=surrogate[config.iH2O].MM/MOW(b,ilayer,iphase)*MOinit(b,ilayer,iphase)*RH
                      /surrogate[config.iH2O].gamma_org_layer(b,ilayer,iphase);
                    if(surrogate[config.iH2O].Ap_layer_init(b,ilayer,iphase) > 0.0)
                      {
                       error=(temp-surrogate[config.iH2O].Ap_layer_init(b,ilayer,iphase))/surrogate[config.iH2O].Ap_layer_init(b,ilayer,iphase);
                       surrogate[config.iH2O].Ap_layer_init(b,ilayer,iphase)=temp;
                      }
		    else
		      surrogate[config.iH2O].Ap_layer_init(b,ilayer,iphase)=0.0;
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
	    density_aqueous_phase(config, surrogate, LWC, Temperature);
	    config.rho_aqueous=config.AQrho(b);

            for (i=0;i<n;++i)
                surrogate[i].Aaq=surrogate[i].Aaq_bins_init(b);


	    compute_ionic_strenght2(config,surrogate,
                                    AQinit(b), conc_inorganic(b), ionic(b), chp(b),
				    organion(b), ionic_organic(b), factor);
		  
	    activity_coefficients_aq(config,surrogate,Temperature,0.0,MMaq(b),XH2O);
	    if (config.compute_long_and_medium_range_interactions)
	      activity_coefficients_LR_MR(config, surrogate, Temperature, 0.0, ionic(b));
	    
	    temp=surrogate[config.iH2O].Aaq_bins_init(b);
	    error=(temp-surrogate[config.iH2O].Aaq_bins_init(b))/surrogate[config.iH2O].Aaq_bins_init(b);
	  	    
	    if (AQinit(b)>0.0)
	      {
		double inorg1;
		double Ke=1.0e-14;
		inorg1=0.0;
		for (i=0;i<n;++i) 
		  if (surrogate[i].is_organic==false and surrogate[i].is_inorganic_precursor==false and i!=config.iHp and i!=config.iH2O)
		    inorg1-=surrogate[i].Aaq_bins_init(b)/surrogate[i].MM*surrogate[i].charge*config.AQrho(b)/AQinit(b);
		   
		chp(b)=0.5*(inorg1+pow(pow(inorg1,2)+4*Ke,0.5));
		if (chp(b)==0.0)
		  chp(b)=1.0e-7;
	      }
	    else
	      chp(b)=1.0e-7;

	    surrogate[config.iH2O].Aaq_bins_init(b)=temp;
	    surrogate[config.iHp].Aaq_bins_init(b)=chp(b)*AQinit(b)/config.AQrho(b);

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
      activity_coefficients_dyn_aq(config, surrogate, Temperature,AQinit,MOinit,
				   conc_inorganic, ionic, ionic_organic,
				   organion,chp,LWC,MMaq,1.0);
  
  //Computation of Ag for water
  water_concentration(config, surrogate, Temperature, RH);
}

void dynamic_system(model_config &config, vector<species> &surrogate,
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
  double error_aq=10.0;
  double error_org=10.0;
  double error_tot=10.0;
  MO.resize(config.nbins,config.nlayer,config.max_number_of_phases);
  int index=0;
  Array<double,1> AQ,conc_inorganic,ionic_organic, organion,MMaq,chp1,chp0;
  AQ.resize(config.nbins);
  ionic_organic.resize(config.nbins);
  conc_inorganic.resize(config.nbins);
  organion.resize(config.nbins);
  MMaq.resize(config.nbins);
  chp1.resize(config.nbins);
  chp0.resize(config.nbins);
  int icycle;

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
  
  double error_old,error_old2;
  double LWCtot=0.0;
  for (b=0;b<config.nbins;++b)
    LWCtot+=LWC(b);

  initialisation(config, surrogate, MOinit, MO, MOW, AQinit, MMaq, LWCtot, LWC, chp, ionic, ionic_organic,
		 organion, conc_inorganic, Temperature, RH);

  compute_diameters(config, surrogate, Vsol, number, LWC, LWCtot);

  if (config.explicit_representation)
    compute_morphology(config, Vsol, number);

  //Computation of characteristic times to reach equilibrium
  tau_dif(config, surrogate, number, Vsol);
  tau_kmt(config, surrogate, Temperature, number);
  compute_kp_org(config, surrogate, MOinit, Temperature, MOW);
  if (LWCtot>config.LWClimit)
    compute_kp_aq(config, surrogate, Temperature, ionic, chp, MMaq);

  characteristic_time(config, surrogate, MOinit, AQinit, LWCtot); 
  if (LWCtot>config.LWClimit)
    characteristic_time_aq(config, surrogate, LWC, AQinit, MOinit);

  if (config.coupled_phases or 
      (RH>=config.RHcoupling and surrogate[config.iH2O].hydrophilic and surrogate[config.iH2O].hydrophobic))
    {
      // If the syst is coupled particulate hydrophilic compounds concentrations
      //must be computed simultaneously with the particulate hydrophobic compounds concentrations

      //Initialisation of equibrium
      if (config.compute_saturation and config.compute_organic)
        {
          number_org_phases(config,surrogate,Temperature,MOinit,MOW);
          tau_dif(config, surrogate, number, Vsol);
          tau_kmt(config, surrogate, Temperature, number);
          compute_kp_org(config, surrogate, MOinit, Temperature, MOW);
          characteristic_time(config, surrogate, MOinit, AQinit, LWCtot); 
          if (LWCtot>config.LWClimit)
            characteristic_time_aq(config, surrogate, LWC, AQinit, MOinit);
        }

      if (config.coupling_organic_inorganic or config.compute_organic==false 
          or config.compute_inorganic==false)
        solve_local_equilibriums_coupled(config, surrogate, MOinit, MOW, number, Vsol, LWC, AQinit, ionic, chp, Temperature, RH, AQ, MO,
                                         conc_inorganic, ionic_organic, organion, MMaq);
      else
        for (icycle=0;icycle<config.number_of_org_inorg_cycles;icycle++)
          {
            config.compute_organic=false;
            solve_local_equilibriums_coupled(config, surrogate, MOinit, MOW, number, Vsol, LWC, AQinit, ionic, chp, Temperature, RH, AQ, MO,
                                             conc_inorganic, ionic_organic, organion, MMaq);
            config.compute_inorganic=false; 
            config.compute_organic=true;
            solve_local_equilibriums_coupled(config, surrogate, MOinit, MOW, number, Vsol, LWC, AQinit, ionic, chp, Temperature, RH, AQ, MO,
                                             conc_inorganic, ionic_organic, organion, MMaq);
            config.compute_inorganic=true;
          }

      characteristic_time(config, surrogate, MOinit, AQinit, LWCtot); 
      if (LWCtot>config.LWClimit)
        characteristic_time_aq(config, surrogate, LWC, AQinit, MOinit);
   
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
       
	  if (config.first_evaluation_activity_coefficients==false)
	    dynamic_tot(config,surrogate,MOinit,MO,MOW,AQinit,AQ,conc_inorganic,
			ionic,ionic_organic,organion,chp,chp1,chp0,LWC,MMaq,Temperature,
			deltat1,tequilibrium,true);
	  else
	    dynamic_tot(config,surrogate,MOinit,MO,MOW,AQinit,AQ,conc_inorganic,
			ionic,ionic_organic,organion,chp,chp1,chp0,LWC,MMaq,Temperature,
			deltat1,tequilibrium,false);	    	  
		  
          //compute the new time step so that changes are small
          adapstep(config,surrogate,config.tequilibrium,deltat1,t,deltatmax,config.deltatmin,
                   MOinit,MO,LWCtot,AQinit,AQ,LWC,conc_inorganic,chp,chp1,chp0);
		  
          if (deltat1<0.999*deltat2) //if the new time step is inferior to the old one
            {                       //the old time step is rejected	     
              for (i=0;i<n;++i)
                {
                  surrogate[i].Ag=surrogate[i].Ag0;
                  for (ilayer=0;ilayer<config.nlayer;++ilayer)
                    for (b=0;b<config.nbins;++b)
                      for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                        surrogate[i].Ap_layer_init(b,ilayer,iphase)=
                          surrogate[i].Ap_layer_init0(b,ilayer,iphase);

                  for (b=0;b<config.nbins;++b)
                    surrogate[i].Aaq_bins_init(b)=surrogate[i].Aaq_bins_init0(b);				  
		}
		  
	      for (b=0;b<config.nbins;++b)
		{
		  AQinit(b)=LWC(b);
		  for (i=0;i<n;++i)
		    if(surrogate[i].hydrophilic)
		      AQinit(b)+=surrogate[i].Aaq_bins(b);
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
              redistribution(config, surrogate,MOinit,MO);

              for (b=0;b<config.nbins;++b)
                {
                  for (ilayer=0;ilayer<config.nlayer;++ilayer)
                    for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
                      MOinit(b,ilayer,iphase)=MO(b,ilayer,iphase);
                  AQinit(b)=AQ(b);
                }

              if (LWCtot>config.LWClimit)
                density_aqueous_phase(config, surrogate, LWC, Temperature);

              //compute the new diameters of particles
              compute_diameters(config, surrogate, Vsol, number, LWC, LWCtot);
	      if (config.explicit_representation)
		compute_morphology(config, Vsol,number);	     

              //Computation of characteristic times to reach equilibrium
              tau_dif(config, surrogate, number, Vsol);
              tau_kmt(config, surrogate, Temperature, number);
              characteristic_time(config, surrogate, MOinit, AQinit, LWCtot);
              if (LWCtot>config.LWClimit)
                characteristic_time_aq(config, surrogate, LWC, AQinit, MOinit);
			  
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
                solve_local_equilibriums_coupled(config, surrogate, MOinit, MOW, number, Vsol, LWC, AQinit, ionic, chp, Temperature, RH, AQ, MO,
                                                 conc_inorganic, ionic_organic, organion, MMaq);
              else
                for (icycle=0;icycle<config.number_of_org_inorg_cycles;icycle++)
                  {
                    config.compute_organic=false;
                    solve_local_equilibriums_coupled(config, surrogate, MOinit, MOW, number, Vsol, LWC, AQinit, ionic, chp, Temperature, RH, AQ, MO,
                                                     conc_inorganic, ionic_organic, organion, MMaq);
                    config.compute_inorganic=false;
                    config.compute_organic=true;
                    solve_local_equilibriums_coupled(config, surrogate, MOinit, MOW, number, Vsol, LWC, AQinit, ionic, chp, Temperature, RH, AQ, MO,
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
  
      //Initialisation of equibrium
      //compute phase separation
      if (config.compute_saturation and config.compute_organic)
        {
          number_org_phases(config,surrogate,Temperature,MOinit,MOW);
          tau_dif(config, surrogate, number, Vsol);

          tau_kmt(config, surrogate, Temperature, number);
          compute_kp_org(config, surrogate, MOinit, Temperature, MOW);
          characteristic_time(config, surrogate, MOinit, AQinit, LWCtot); 
          if (LWCtot>config.LWClimit)
            characteristic_time_aq(config, surrogate, LWC, AQinit, MOinit); 
        }

      if (config.coupling_organic_inorganic or config.compute_organic==false 
          or config.compute_inorganic==false)
        solve_local_equilibriums_uncoupled(config, surrogate, MOinit, MOW, number, Vsol, LWC, AQinit, ionic, chp, Temperature, RH, AQ, MO,
                                         conc_inorganic, ionic_organic, organion, MMaq);
      else
        for (icycle=0;icycle<config.number_of_org_inorg_cycles;icycle++)
          {
            config.compute_organic=false;
            solve_local_equilibriums_uncoupled(config, surrogate, MOinit, MOW, number, Vsol, LWC, AQinit, ionic, chp, Temperature, RH, AQ, MO,
                                               conc_inorganic, ionic_organic, organion, MMaq);
            config.compute_inorganic=false;
            config.compute_organic=true;
            solve_local_equilibriums_uncoupled(config, surrogate, MOinit, MOW, number, Vsol, LWC, AQinit, ionic, chp, Temperature, RH, AQ, MO,
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
              dynamic_org(config,surrogate,MOinit,MO,AQinit,
                          MOW,Temperature,deltat1,config.tequilibrium, true);
              if (LWCtot>config.LWClimit)
                dynamic_aq(config,surrogate,AQinit,AQ,MOinit,conc_inorganic,ionic,ionic_organic,
                           organion,chp,chp1,chp0,LWC,MMaq,Temperature,deltat1,config.tequilibrium, true);
            }
          else
            {
              dynamic_org(config,surrogate,MOinit,MO,AQinit,
                          MOW,Temperature,deltat1,config.tequilibrium, false);
              if (LWCtot>config.LWClimit)
                dynamic_aq(config,surrogate,AQinit,AQ,MOinit,conc_inorganic,ionic,ionic_organic,
                           organion,chp,chp1,chp0,LWC,MMaq,Temperature,deltat1,config.tequilibrium, false);
            }
	  
          //compute the new time step so that changes are small
          adapstep(config,surrogate,config.tequilibrium,deltat1,t,deltatmax,config.deltatmin,
                   MOinit,MO,LWCtot,AQinit,AQ,LWC,conc_inorganic,chp,chp1,chp0);

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
                          surrogate[i].Ap_layer_init(b,ilayer,iphase)=
                            surrogate[i].Ap_layer_init0(b,ilayer,iphase);

		  for (b=0;b<config.nbins;++b)
		    {
		      AQinit(b)=LWC(b);
		      for (i=0;i<n;++i)
			if(surrogate[i].hydrophilic)
			  AQinit(b)+=surrogate[i].Aaq_bins(b);
		      AQinit(b)=max(AQinit(b),config.MOmin);		      
		      chp(b)=chp0(b);
		    }

                  if (LWCtot>config.LWClimit)
                    if (surrogate[i].hydrophilic)
                      for (b=0;b<config.nbins;++b)
                        surrogate[i].Aaq_bins_init(b)=surrogate[i].Aaq_bins_init0(b);
		  
		  for (b=0;b<config.nbins;++b)
		    for (ilayer=0;ilayer<config.nlayer;++ilayer)
		      {
			for (iphase=0;iphase<config.nphase(b,ilayer);++iphase)
			  {			   
			    MOinit(b,ilayer,iphase)=0.0;
			    for (i=0;i<n;++i)			
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
              
            }
          else                     //the time step is accepted 
            {
              t+=deltat2;
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
              redistribution(config, surrogate,MOinit,MO);

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
                density_aqueous_phase(config, surrogate, LWC, Temperature);

              compute_diameters(config, surrogate, Vsol, number, LWC, LWCtot);
	      if (config.explicit_representation)
		compute_morphology(config, Vsol, number);

              //Computation of characteristic times to reach equilibrium
              tau_dif(config, surrogate, number, Vsol);
              tau_kmt(config, surrogate, Temperature, number);
              characteristic_time(config, surrogate, MOinit, AQinit, LWCtot);
              if (LWCtot>config.LWClimit)
                characteristic_time_aq(config, surrogate, LWC, AQinit, MOinit);

              //computation of concentrations at equilibrium
              if (config.coupling_organic_inorganic or config.compute_organic==false 
                  or config.compute_inorganic==false)
                solve_local_equilibriums_uncoupled(config, surrogate, MOinit, MOW, number, Vsol, LWC, AQinit, ionic, chp, Temperature, RH, AQ, MO,
                                                   conc_inorganic, ionic_organic, organion, MMaq);
              else
                for (icycle=0;icycle<config.number_of_org_inorg_cycles;icycle++)
                  {
                    config.compute_organic=false;
                    solve_local_equilibriums_uncoupled(config, surrogate, MOinit, MOW, number, Vsol, LWC, AQinit, ionic, chp, Temperature, RH, AQ, MO,
                                                       conc_inorganic, ionic_organic, organion, MMaq);
                    config.compute_inorganic=false;
                    config.compute_organic=true;
                    solve_local_equilibriums_uncoupled(config, surrogate, MOinit, MOW, number, Vsol, LWC, AQinit, ionic, chp, Temperature, RH, AQ, MO,
                                                       conc_inorganic, ionic_organic, organion, MMaq);
                    config.compute_inorganic=true;
                  }
            }
        }
    }
}

