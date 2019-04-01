
void aiomfac(Array<double, 1> &X_solvents, Array<double, 1> &gamma_LR_solvents,
             Array<double, 1> &gamma_MR_solvents, Array<double,1> &molar_mass_solvents,
             Array<double, 1> &molar_mass_groups,
             Array<double, 1> &molality_ions, Array<double, 1> &gamma_LR_ions,
             Array<double, 1> &gamma_MR_ions, Array<double, 1> &charges_ions, 
             double &Temperature, double &ionic, int &ngroups,
             Array<double, 2> &groups_solvents,
             Array<double, 2> &b1ki, Array<double, 2> &b2ki, Array<double, 2> &b1ca,
             Array<double, 2> &b2ca, Array<double, 2> &b3ca, Array<double, 2> &c1ca,
             Array<double, 2> &c2ca, Array<double, 2> &Rcc, Array<double, 3> &Qcca)
{
  int nsolvents=X_solvents.size();
  int nions=molality_ions.size();
 
  double rho=997.0;  //volumic mass in kg.m-3 
  double permittivity=78.54;  //relative static permittivity of water"

  int i,j,k,l,isol;    
  double logmax=20.;

  if (ionic>0.0)
    {      
      ////////////////////////////////////////
      // LONG - RANGE computation           //
      ////////////////////////////////////////

      double sqrt_ionic=pow(ionic,0.5); 
      double A=1.327757e5*pow(rho,0.5)/(pow(permittivity*Temperature,3.0/2));
      double b=6.369696*pow(rho,0.5)/(pow(permittivity*Temperature,0.5));      
      double A2=A*(1+b*sqrt_ionic-1.0/(1+b*sqrt_ionic)
                  -2.0*log(1.0+b*sqrt_ionic));
      double A3=A*sqrt_ionic/(1+b*sqrt_ionic);

      for (k=0;k<nsolvents;++k)
        gamma_LR_solvents(k)=exp(2*molar_mass_solvents(k)/pow(b,3)*A2); 
	  
      for (i=0;i<nions;++i)
        gamma_LR_ions(i)=exp(-pow(charges_ions(i),2)*A3);

      //////////////////////////////////////////
      // MEDIUM - RANGE computation           //
      //////////////////////////////////////////
	  
      Array<double, 2> Bki,IBki_deriv;
      Array<double, 2> Bca,IBca_deriv;
      Array<double, 2> Cca,ICca_deriv;
      Bki.resize(ngroups,nions);
      IBki_deriv.resize(ngroups,nions);
      Bca.resize(nions,nions);
      IBca_deriv.resize(nions,nions);
      Cca.resize(nions,nions);
      ICca_deriv.resize(nions,nions);
      Bki=0.;
      Bca=0.;      
      double I1=exp(-1.2*sqrt_ionic);
      double I2=I1/sqrt_ionic;

      for (k=0;k<ngroups;++k)
        for (i=0;i<nions;++i)
          {
            Bki(k,i)=b1ki(k,i)+b2ki(k,i)*I1;
            IBki_deriv(k,i)=-0.6*b2ki(k,i)*I2;
          }

      
      
      for (i=0;i<nions;++i)
        for (j=0;j<nions;++j)
          if(charges_ions(i)*charges_ions(j)<0.0)
            {
              Bca(i,j)=b1ca(i,j)+b2ca(i,j)*exp(-b3ca(i,j)*sqrt_ionic);
              IBca_deriv(i,j)=-0.5*b3ca(i,j)*b2ca(i,j)*
                exp(-b3ca(i,j)*sqrt_ionic)/sqrt_ionic;
              Cca(i,j)=c1ca(i,j)*exp(-c2ca(i,j)*sqrt_ionic);
              ICca_deriv(i,j)=-0.5*c2ca(i,j)*Cca(i,j)/sqrt_ionic;
            }
          else
            {
              Bca(i,j)=0.0;
              IBca_deriv(i,j)=0.0;
              Cca(i,j)=0.0;
	      ICca_deriv(i,j)=0.0;
            }

      Array<double, 1> Xk;
      Array<double, 1> log_gammak;
      log_gammak.resize(ngroups);
      Xk.resize(ngroups);
	  
      double sumXk=0.0;
      for (k=0;k<ngroups;++k)
        {
          Xk(k)=0.0;
          for (i=0;i<nsolvents;++i)
            Xk(k)+=X_solvents(i)*groups_solvents(i,k);      
         
          sumXk+=Xk(k);
        }     

      if (sumXk>0.0)
        for (k=0;k<ngroups;++k)
          Xk(k)/=sumXk;    
      else
        cout << "sumXk is zero" << endl;

      double molar_mass_average=0.0;
      for (k=0;k<ngroups;++k)
        molar_mass_average+=Xk(k)*molar_mass_groups(k); //CHANGED Based on a corrigeum published in ACP by Zuend. BEFORE: molar_mass_average+=X_solvents(k)*molar_mass_solvents(k); 

      double mi_zi=0.0;
      for (i=0;i<nions;++i)
        mi_zi+=molality_ions(i)*abs(charges_ions(i));
	  
      for (k=0;k<ngroups;++k)
        {
          log_gammak(k)=0.0;	 
          for (i=0;i<nions;++i)
            {
              log_gammak(k)+=Bki(k,i)*molality_ions(i);	   
	   
              double a=0.;   
              for (l=0;l<ngroups;++l)            
                a+=(Bki(l,i)+IBki_deriv(l,i)*ionic)*Xk(l)/molar_mass_average;
              	  
              if(charges_ions(i)>0.0)
                {
                                 
                  for (j=0;j<nions;++j)
                    if(charges_ions(j)<0.0)
                      a+=(Bca(i,j)+IBca_deriv(i,j)*ionic+mi_zi*(2*Cca(i,j)+ICca_deriv(i,j)*ionic))*molality_ions(j);

                 for (j=i;j<nions;++j)
                    if(charges_ions(j)>0.0)
                      {
                        a+=Rcc(i,j)*molality_ions(j);
                        for (l=0;l<nions;++l)
                          if(charges_ions(l)<0.0)
                            a+=2.0*Qcca(i,j,l)*molality_ions(j)*molality_ions(l);
                      }
                }	
              log_gammak(k)-=a*molality_ions(i)*molar_mass_groups(k);
            }
        }

      for (l=0;l<nsolvents;++l)
        {
          gamma_MR_solvents(l)=0.0;
          for (k=0;k<ngroups;++k)
            gamma_MR_solvents(l)+=groups_solvents(l,k)*log_gammak(k);

	  gamma_MR_solvents(l)=min(gamma_MR_solvents(l),logmax);
	  gamma_MR_solvents(l)=max(gamma_MR_solvents(l),-logmax);	  
          gamma_MR_solvents(l)=exp(gamma_MR_solvents(l));
        }

      double loggamma;
      for (i=0;i<nions;++i)
        if(charges_ions(i)>0.0)
          {
            loggamma=0;
            for (k=0;k<ngroups;++k)
              loggamma+=1.0/molar_mass_average*Bki(k,i)*Xk(k);           
			
            for (j=0;j<nions;++j)
              {
                double a=0.0;
                for (k=0;k<ngroups;++k)                
                  a+=IBki_deriv(k,j)*Xk(k);

                a*=pow(charges_ions(i),2)/(2*molar_mass_average);
			          
                if(charges_ions(j)<0.0)
                  a+=Bca(i,j)+Cca(i,j)*mi_zi;

                if(charges_ions(j)>0.0)
                  {
                    a+=Rcc(i,j);
                    for (l=0;l<nions;++l)
                      if(charges_ions(l)<0.0)
                        {
                          a+=molality_ions(l)*
                            (0.5*pow(charges_ions(i),2)*IBca_deriv(j,l)+
                             +(Cca(j,l)*abs(charges_ions(i))+ICca_deriv(j,l)*0.5*pow(charges_ions(i),2)*mi_zi)+			                                               			                                 
                             Qcca(i,j,l));
                        }
                  }
                loggamma+=a*molality_ions(j);
              }

	    loggamma=min(loggamma,logmax);
	    loggamma=max(loggamma,-logmax);
            gamma_MR_ions(i)=exp(loggamma);
          }
        else
          {
            loggamma=0;
            for (k=0;k<ngroups;++k)
              loggamma+=1.0/molar_mass_average*Bki(k,i)*Xk(k);
            
            for (j=0;j<nions;++j)
              {
                double a=0.;
                for (k=0;k<ngroups;++k)              
                  a+=IBki_deriv(k,j)*Xk(k);
                a*=pow(charges_ions(i),2)/(2*molar_mass_average);
			   
                if(charges_ions(j)>0.0)
                  {
                    a+=Bca(j,i)+Cca(j,i)*mi_zi;

                    for (l=0;l<nions;++l)
                      {
                        if(charges_ions(l)<0.0)
                          {
                            a+=molality_ions(l)*
                              (0.5*pow(charges_ions(i),2)*IBca_deriv(j,l)+
                               (Cca(j,l)*abs(charges_ions(i))+ICca_deriv(j,l)*0.5*pow(charges_ions(i),2)*mi_zi));
                          }
                        else if(charges_ions(l)>0.0)
                          a+=Qcca(j,l,i)*molality_ions(l);
                      }
                  }
                loggamma+=a*molality_ions(j);
              }
			
	    loggamma=min(loggamma,logmax);
	    loggamma=max(loggamma,-logmax);
            gamma_MR_ions(i)=exp(loggamma);
          }
    }
  else
    {
      for (i=0;i<nions;++i)
        {
          gamma_MR_ions(i)=1.0;
          gamma_LR_ions(i)=1.0;
        }

      for (i=0;i<nsolvents;++i)
        {
          gamma_MR_solvents(i)=1.0;
          gamma_LR_solvents(i)=1.0;
        } 
    }
}
