//!!-----------------------------------------------------------------------
//!!     Copyright (C) 2019 CEREA (ENPC) - INERIS
//!!     SSH-aerosol is distributed under the GNU General Public License v3
//!!-----------------------------------------------------------------------

#include <math.h>

void aiomfac_ssh(Array<double, 1> &X_solvents, Array<double, 1> &gamma_LR_solvents,
		 Array<double, 1> &gamma_MR_solvents, Array<double,1> &molar_mass_solvents,
		 Array<double, 1> &molar_mass_groups,
		 Array<double, 1> &molality_ions, Array<double, 1> &gamma_LR_ions,
		 Array<double, 1> &gamma_MR_ions, Array<double, 1> &charges_ions, 
		 double &Temperature, double &ionic, int &ngroups,
		 Array<double, 2> &groups_solvents,
		 Array<double, 2> &b1ki, Array<double, 2> &b2ki, Array<double, 2> &b1ca,
		 Array<double, 2> &b2ca, Array<double, 2> &b3ca, Array<double, 2> &c1ca,
		 Array<double, 2> &c2ca, Array<double, 2> &Rcc, Array<double, 3> &Qcca,
		 int &jH2O, int &iH, bool compute_organic)
{
  int nsolvents=X_solvents.size();
  int nions=molality_ions.size();
 
  //double rho=config.rho_aqueous;
  double rho=997.0; //1000.0;  //volumic mass in kg.m-3 
  double permittivity=78.54;  //relative static permittivity of water"

  int i,j,k,l;    
  double logmax=12.;
   
  /*cout << "A" << groups_solvents << endl;
    cout << "B" << b1ki << endl;
    cout << "C" << b2ki << endl;
    cout << "D" << b1ca << endl;*/
  //b1ca=0.0;
  //cout << "E" << b2ca << endl;
  //b2ca=0.0;
  //cout << "F" << b3ca << endl;
  //b3ca=0.0;
  //cout << "G" << c1ca << endl;
  //c1ca=0.0;
  //cout << "H" << c2ca << endl;
  //c2ca=0.0;
  //Rcc=0.0;
  //cout << "I" << Rcc << endl;
  //Qcca=0.0;
  //cout << "J" << Qcca << endl;
  if (ionic>0.0)
    {      
      //cout << "ionic:" << ionic << endl;
      ////////////////////////////////////////
      // LONG - RANGE computation           //
      ////////////////////////////////////////

      double sqrt_ionic=pow(ionic,0.5);
      double sqrt_rho=pow(rho,0.5);
      double A=1.327757e5*sqrt_rho/(pow(permittivity*Temperature,3.0/2));
      double b=6.369696*sqrt_rho/(pow(permittivity*Temperature,0.5));      
      double A2=A*(1+b*sqrt_ionic-1.0/(1+b*sqrt_ionic)
                   -2.0*log(1.0+b*sqrt_ionic))/pow(b,3);
      double A3=A*sqrt_ionic/(1+b*sqrt_ionic);

      static Array<double, 2> Bki,IBki_deriv;
      static Array<double, 2> Bca,IBca_deriv;
      static Array<double, 2> Cca,ICca_deriv;
      static Array<double, 1> Afaci,sumIBkXk;
      static Array<double, 1> Xk;
      static Array<double, 1> log_gammak;
      static int iRcc,jRcc,jQcc;
      //cout << Afaci << " " << Afaci.size() << endl;
      if (int(Afaci.size())!=nions)
        {
          Afaci.resize(nions);
          sumIBkXk.resize(nions);
          Bki.resize(ngroups,nions);
          IBki_deriv.resize(ngroups,nions);
          Bca.resize(nions,nions);
          IBca_deriv.resize(nions,nions);
          Cca.resize(nions,nions);
          ICca_deriv.resize(nions,nions);
          log_gammak.resize(ngroups);
          Xk.resize(ngroups);
          for (i=0;i<nions;++i)
            if(charges_ions(i)>0.0)
              {   
                for (j=i;j<nions;++j)
                  if(charges_ions(j)>0.0)
                    {
                      if (Rcc(i,j)!=0)
                        {
                          //cout << Rcc(i,j) << " " << i << " " << j << endl;
                          iRcc=i;
                          jRcc=j;
                          for (l=0;l<nions;++l)
                            if(charges_ions(l)<0.0)
                              if (Qcca(i,j,l)>0.0)
                                jQcc=l;
                        }
                    }           
              }

          
        }

      //cout << Rcc(iRcc,jRcc) << " " << Qcca(iRcc,jRcc,jQcc) << endl;
      
      for (k=0;k<nsolvents;++k)
        gamma_LR_solvents(k)=exp(2.*molar_mass_solvents(k)*A2);

      gamma_MR_ions=1.;
      gamma_LR_ions=1.;
      if (iH<0)
        for (i=0;i<nions;++i)
          gamma_LR_ions(i)=exp(-pow(charges_ions(i),2)*A3);
      else
        gamma_LR_ions(iH)=exp(-pow(charges_ions(iH),2)*A3);

      //////////////////////////////////////////
      // MEDIUM - RANGE computation           //
      //////////////////////////////////////////
	  
      
      Bki=0.;
      Bca=0.;      
      double I1=exp(-1.2*sqrt_ionic);
      double I2=I1/sqrt_ionic;

      for (k=0;k<ngroups;++k)
        for (i=0;i<nions;++i)
          {
            Bki(k,i)=b1ki(k,i)+b2ki(k,i)*I1;
            IBki_deriv(k,i)=-0.6*b2ki(k,i)*I2;
	    //if (k==0 and i==1)
	    //  cout << "Bbbbb " << b1ki(k,i) << " " << b2ki(k,i) << " " << b1ki(k,i)+b2ki(k,i)*exp(-1.2*sqrt_ionic) << endl;
          }

      
      //Cca=c1ca*exp(-c2ca*sqrt_ionic);
      for (i=0;i<nions;++i)
        for (j=0;j<nions;++j)
          if(charges_ions(i)*charges_ions(j)<0.0)
            {
              double expb3ca=exp(-b3ca(i,j)*sqrt_ionic);
              Bca(i,j)=b1ca(i,j)+b2ca(i,j)*expb3ca;
              IBca_deriv(i,j)=-0.5*b3ca(i,j)*b2ca(i,j)*
                expb3ca/sqrt_ionic;
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

      //cout << Xk << endl;
      double molar_mass_average=0.0;
      for (k=0;k<ngroups;++k)
        molar_mass_average+=Xk(k)*molar_mass_groups(k); //CHANGED Based on a corrigeum published in ACP by Zuend. BEFORE: molar_mass_average+=X_solvents(k)*molar_mass_solvents(k); 

      double mi_zi=0.0;
      for (i=0;i<nions;++i)
        mi_zi+=molality_ions(i)*abs(charges_ions(i));
      
      int kbeg,kend;
      
      if (compute_organic==false)
        {
          kbeg=jH2O;
          kend=jH2O+1;
        }
      else
        {
          //cout << "compute organic is true" << endl;
          kbeg=0;
          kend=ngroups;
        }
      //cout << molar_mass_groups(jH2O) << endl;
      //cout << kbeg << " " << kend << endl;
      log_gammak=0.0;
      
      Afaci=0.;
      for (i=0;i<nions;++i)
        {
          if(charges_ions(i)>0.0)
            {
                                 
              for (j=0;j<nions;++j)
                if(charges_ions(j)<0.0)
                  {
                    //cout << Bca(i,j) << " " << IBca_deriv(i,j) << " " << 2*Cca(i,j)+ICca_deriv(i,j)*ionic << endl;
                    Afaci(i)+=(Bca(i,j)+IBca_deriv(i,j)*ionic+mi_zi*(2*Cca(i,j)+ICca_deriv(i,j)*ionic))*molality_ions(j);
                  }

              if (i==iRcc)
                {
                  j=jRcc;
                  l=jQcc;
                  Afaci(i)+=Rcc(i,j)*molality_ions(j)+2.0*Qcca(i,j,l)*molality_ions(j)*molality_ions(l);  
                }           
            }
        }


      for (k=kbeg;k<kend;++k)
        {
          //log_gammak(k)=0.0;	 
          for (i=0;i<nions;++i)
            {
              
              log_gammak(k)+=Bki(k,i)*molality_ions(i);	                 
              
              double a=0.;                 
              for (l=0;l<ngroups;++l)            
                a+=(Bki(l,i)+IBki_deriv(l,i)*ionic)*Xk(l)/molar_mass_average;
              
              if(charges_ions(i)>0.0)
                a+=Afaci(i);                              
             
              log_gammak(k)-=a*molality_ions(i)*molar_mass_groups(k);
            }
          //cout << log_gammak(k) << endl;
        }
      //cout << molality_ions << endl;
      //cout << ionic << endl;

      /*
        cout << "Gr " << groups_solvents << endl;
        cout << "lng " << log_gammak << endl;*/
      //cout << Bki << endl;

      for (l=0;l<nsolvents;++l)
        {
          gamma_MR_solvents(l)=0.0;
          for (k=0;k<ngroups;++k)
            gamma_MR_solvents(l)+=groups_solvents(l,k)*log_gammak(k);

          if (gamma_MR_solvents(l)>logmax) gamma_MR_solvents(l)=logmax;
          if (gamma_MR_solvents(l)<-logmax) gamma_MR_solvents(l)=-logmax;
          gamma_MR_solvents(l)=exp(gamma_MR_solvents(l));
          //cout << "MR " << l << " " << gamma_MR_solvents(l) << endl;
        }

      sumIBkXk=0.;
      for (j=0;j<nions;++j)
        {
          //double a=0.;
          for (k=0;k<ngroups;++k)              
            sumIBkXk(j)+=IBki_deriv(k,j)*Xk(k);
        }

      double loggamma;
      for (i=0;i<nions;++i)
        if(charges_ions(i)>0.0 and (iH<0 or i==iH))
          {
            loggamma=0;
            for (k=0;k<ngroups;++k)
              loggamma+=1.0/molar_mass_average*Bki(k,i)*Xk(k);           
			
            for (j=0;j<nions;++j)
              {
                double a=sumIBkXk(j);

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

            if (loggamma>logmax) loggamma=logmax;
            if (loggamma<-logmax) loggamma=-logmax; 
            gamma_MR_ions(i)=exp(loggamma);
          }
        else if (charges_ions(i)<0.0 and iH<0)
          {
            //cout << i << " " << charges_ions(i) << endl;
            loggamma=0;
            for (k=0;k<ngroups;++k)
              loggamma+=1.0/molar_mass_average*Bki(k,i)*Xk(k);
            
            //if (i==1)
            //  cout << "A: " << loggamma << " " << molar_mass_average << endl;            
            for (j=0;j<nions;++j)
              {
                double a=sumIBkXk(j);
                a*=pow(charges_ions(i),2)/(2*molar_mass_average);
                /*
                if (i==6)
                cout << "A1: " << a*molality_ions(j) << endl;*/
			   
                if(charges_ions(j)>0.0)
                  {
                    a+=Bca(j,i)+Cca(j,i)*mi_zi;
                    /*
                    if (i==6 and molality_ions(j)>0.)
                    cout << "A2: " << j << " " << a*molality_ions(j) << endl;*/

                    for (l=0;l<nions;++l)
                      {
                        if(charges_ions(l)<0.0)
                          {
                            a+=molality_ions(l)*
                              (0.5*pow(charges_ions(i),2)*IBca_deriv(j,l)+
                               (Cca(j,l)*abs(charges_ions(i))+ICca_deriv(j,l)*0.5*pow(charges_ions(i),2)*mi_zi));
                            /*
                            if (i==6)                              
                              if (molality_ions(l)>0. and molality_ions(j)>0.)
                                cout << "A3 :" << j << " " << l << " " << molality_ions(l)*
                                  (0.5*pow(charges_ions(i),2)*IBca_deriv(j,l)+
                                  (Cca(j,l)*abs(charges_ions(i))+ICca_deriv(j,l)*0.5*pow(charges_ions(i),2)*mi_zi))*molality_ions(j) << endl;*/
                          }
                        else if(charges_ions(l)>0.0)
                          {
                            a+=Qcca(j,l,i)*molality_ions(l);
                            /*
                            if (i==6)
                              if (molality_ions(l)>0. and molality_ions(j)>0.)
                              cout << "A4: " << j << " " << l << " " << Qcca(j,l,i)*molality_ions(l)*molality_ions(j) << endl;*/
                          }
                      }
                  }
                loggamma+=a*molality_ions(j);                
                /*if (i==1)
                  if(charges_ions(j)>0.0 and molality_ions(j)>0.)                  
                  cout << j << " " << a*molality_ions(j) << " " << molality_ions(j) << endl;*/
              }

            /*
            if (i==6)
            cout << charges_ions(i) << " " << i << " " << loggamma << endl;*/
            if (loggamma>logmax) loggamma=logmax;
            if (loggamma<-logmax) loggamma=-logmax;
            gamma_MR_ions(i)=exp(loggamma);
          }
      //cout << molality_ions << endl;
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
