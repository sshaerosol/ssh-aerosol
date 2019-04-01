#ifndef UNIFAC_CXX
#define UNIFAC_CXX

void unifac(int &nmols, int &nfunc, Array<double, 2> &groups, Array<double, 1> &Xmol,
	    Array<double, 2> &Inter, Array<double, 2> &InterB, Array<double, 2> &InterC, 
	    Array<double,1> &RG, Array <double, 1> &QG,
	    Array<double, 1> &Rparam, Array <double, 1> &Qparam, Array <double,1> &Lparam,
	    double &Z, double &Temperature, Array<double,1> &gamma, 
	    bool temperature_dependancy)
{  
  double xmin=1.0e-11;
  double xtot=0.0;
  int i,j,k;
  for (i=0;i<nmols;i++)
    {      
      Xmol(i)=max(Xmol(i),xmin);
      xtot+=Xmol(i);
      gamma(i)=1.0;
    }
  
  if (xtot>0.0)
    {
      for (i=0;i<nmols;i++)
        Xmol(i)/=xtot;
  
      Array <double, 1> theta,phi;
      theta.resize(nmols);
      phi.resize(nmols);
      double Lmean=0.0;
      double sumtheta=0.0;
      double sumphi=0.0;
      for (i=0;i<nmols;i++)
        {
          theta(i)=Xmol(i)*Qparam(i);
          sumtheta+=Xmol(i)*Qparam(i);
          phi(i)=Xmol(i)*Rparam(i);
          sumphi+=Xmol(i)*Rparam(i);
          Lmean+=Xmol(i)*Lparam(i);
        }

      
      //Computation of residual gamma
      for (i=0;i<nmols;i++)
        {
          theta(i)/=sumtheta;
          phi(i)/=sumphi; 
          gamma(i)=exp(log(phi(i)/Xmol(i))+Z/2*Qparam(i)*log(theta(i)/phi(i))+Lparam(i)-phi(i)/Xmol(i)*Lmean);                    
        }      
      
      //Compute surface fraction
      Array <double, 1> surface_fraction;
      Array <double, 2> surface_fraction_mol;
      surface_fraction.resize(nfunc);
      surface_fraction_mol.resize(nfunc,nmols);
      surface_fraction=0.0;
      surface_fraction_mol=0.0;
      double sum_surf=0.0;
      for (i=0;i<nmols;i++)
        {
          double sum_surf_mol=0.0;
          for (j=0;j<nfunc;j++)
            if (groups(j,i)>0.0)
              {
                surface_fraction(j)+=QG(j)*Xmol(i)*groups(j,i);            
                surface_fraction_mol(j,i)+=QG(j)*groups(j,i);
                sum_surf+=QG(j)*Xmol(i)*groups(j,i);  
                sum_surf_mol+=QG(j)*groups(j,i);  
              }

          for (j=0;j<nfunc;j++)      
            surface_fraction_mol(j,i)/=sum_surf_mol;
        }
    
      for (j=0;j<nfunc;j++)      
        surface_fraction(j)/=sum_surf;
      
      Array<double,2> Inter2;
      Inter2.resize(nfunc,nfunc);
      if (temperature_dependancy)
	for (j=0;j<nfunc;j++)
	  for (k=0;k<nfunc;k++)          
	    Inter2(j,k)=exp(-Inter(j,k)/Temperature+InterB(j,k)*(1.0/298.15-1.0/Temperature)+InterC(j,k)*(298.15/Temperature-1.0+log(Temperature/298.15)));         
      else
	for (j=0;j<nfunc;j++)
	  for (k=0;k<nfunc;k++)          
	    Inter2(j,k)=exp(-Inter(j,k)/Temperature);   

      Array <double, 1> group_activity;
      group_activity.resize(nfunc);
      Array <double, 2> group_activity_mol,sum2mol;
      group_activity_mol.resize(nfunc,nmols);
      sum2mol.resize(nfunc,nmols);
      
      Array <double, 1> sum2;
      sum2.resize(nfunc);      
      for (j=0;j<nfunc;j++)
        {
          sum2(j)=0.0;
          for (k=0;k<nfunc;k++)
            sum2(j)+=surface_fraction(k)*Inter2(k,j);
        }

      for (j=0;j<nfunc;j++)
        {
          group_activity(j)=1.0-log(sum2(j));                             
          for (k=0;k<nfunc;k++)            
            group_activity(j)-=surface_fraction(k)*Inter2(j,k)/sum2(k);           
        }

      for (i=0;i<nmols;i++)
        for (j=0;j<nfunc;j++)
          if (groups(j,i)>0.0)
            {
              sum2mol(j,i)=0.0;
              for (k=0;k<nfunc;k++)
                if (groups(k,i)>0.0)
                  sum2mol(j,i)+=surface_fraction_mol(k,i)*Inter2(k,j);              
            }

      for (i=0;i<nmols;i++)
        for (j=0;j<nfunc;j++)
          if (groups(j,i)>0.0)
            {
              group_activity_mol(j,i)=1.0-log(sum2mol(j,i));                     
              for (k=0;k<nfunc;k++)      
                if (groups(k,i)>0.0)
                  group_activity_mol(j,i)-=surface_fraction_mol(k,i)*Inter2(j,k)/sum2mol(k,i);           
            }
      
      for (i=0;i<nmols;i++)
        {
          double sum1=0.0;         
          for (j=0;j<nfunc;j++)
            if (groups(j,i)>0.0)
              sum1+=groups(j,i)*QG(j)*(group_activity(j)-group_activity_mol(j,i));                

          gamma(i)*=exp(sum1);	  
        }
    }
}

void tab_unifac(int &nmols, int &nfunc, Array<double, 2> &groups, Array<double, 1> &Xmol,
                Array<double, 2> &Inter, Array<double, 2> &InterB, Array<double, 2> &InterC, 
		Array<double,1> &RG, Array <double, 1> &QG,
                Array<double, 1> &Rparam, Array <double, 1> &Qparam, Array <double,1> &Lparam,
                Array<double, 1> & spec_activity, double &Z, double &Temperature, Array<double,1> &gamma,
		bool temperature_dependancy)
{
  
  double xmin=1.0e-11;
  double xtot=0.0;
  int i,j,k;
  for (i=0;i<nmols;i++)
    {      
      Xmol(i)=max(Xmol(i),xmin);
      xtot+=Xmol(i);
      gamma(i)=1.0;
    }
  
  if (xtot>0.0)
    {
      for (i=0;i<nmols;i++)
        Xmol(i)/=xtot;
  
      Array <double, 1> theta,phi;
      theta.resize(nmols);
      phi.resize(nmols);
      double Lmean=0.0;
      double sumtheta=0.0;
      double sumphi=0.0;
      for (i=0;i<nmols;i++)
        {
          theta(i)=Xmol(i)*Qparam(i);
          sumtheta+=Xmol(i)*Qparam(i);
          phi(i)=Xmol(i)*Rparam(i);
          sumphi+=Xmol(i)*Rparam(i);
          Lmean+=Xmol(i)*Lparam(i);
        }

      
      //Computation of residual gamma
      for (i=0;i<nmols;i++)
        {
          theta(i)/=sumtheta;
          phi(i)/=sumphi;
          gamma(i)=exp(log(phi(i)/Xmol(i))+Z/2*Qparam(i)*log(theta(i)/phi(i))+Lparam(i)-phi(i)/Xmol(i)*Lmean);          
        }      
      
      //Compute surface fraction
      Array <double, 1> surface_fraction;      
      surface_fraction.resize(nfunc);      
      surface_fraction=0.0;      
      double sum_surf=0.0;
      for (i=0;i<nmols;i++)
        {         
          for (j=0;j<nfunc;j++)
            if (groups(j,i)>0.0)
              {
                surface_fraction(j)+=QG(j)*Xmol(i)*groups(j,i);                            
                sum_surf+=QG(j)*Xmol(i)*groups(j,i);                  
              }          
        }
    
      for (j=0;j<nfunc;j++)      
        surface_fraction(j)/=sum_surf;
      
      Array<double,2> Inter2;
      Inter2.resize(nfunc,nfunc);        
      if (temperature_dependancy)
	for (j=0;j<nfunc;j++)
	  for (k=0;k<nfunc;k++)          
	    Inter2(j,k)=exp(-Inter(j,k)/Temperature+InterB(j,k)*(1.0/298.15-1.0/Temperature)+InterC(j,k)*(298.15/Temperature-1.0+log(Temperature/298.15)));         
      else
	for (j=0;j<nfunc;j++)
	  for (k=0;k<nfunc;k++)          
	    Inter2(j,k)=exp(-Inter(j,k)/Temperature); 
      
      Array <double, 1> group_activity;
      group_activity.resize(nfunc);      
      
      Array <double, 1> sum2;
      sum2.resize(nfunc);      
      for (j=0;j<nfunc;j++)
        {
          sum2(j)=0.0;
          for (k=0;k<nfunc;k++)
            sum2(j)+=surface_fraction(k)*Inter2(k,j);
        }

      for (j=0;j<nfunc;j++)
        {         
          group_activity(j)=1.0-log(sum2(j));                             
          for (k=0;k<nfunc;k++)            
            group_activity(j)-=surface_fraction(k)*Inter2(j,k)/sum2(k);                              
        }
      
      for (i=0;i<nmols;i++)
        {
          double sum1=0.0;
          for (j=0;j<nfunc;j++)
            if (groups(j,i)>0.0)
              {                          
                sum1+=groups(j,i)*QG(j)*group_activity(j);                               
              }

          gamma(i)*=exp(sum1-spec_activity(i));         	  
        }
    }
}

void unifac_aq(int &nmols, int &nions, int &nfunc, Array<double, 2> &groups, Array<double, 1> &Xmol,
               Array<double, 1> &Xions,
               Array<double, 2> &Inter, Array<double, 2> &InterB, Array<double, 2> &InterC, 
	       Array<double,1> &RG, Array <double, 1> &QG,
               Array<double, 1> &Rparam, Array <double, 1> &Qparam, Array <double,1> &Lparam,
               Array <double, 1> & RGions, Array <double, 1> & QGions, Array <double, 1> & Lions,
               double &Z, double &Temperature, Array<double,1> &gamma, Array<double,1> &gamma_ions,
	       bool temperature_dependancy)
{ 
  double xmin=1.0e-11;
  double xtot=0.0;
  int i,j,k;
  for (i=0;i<nmols;i++)
    {      
      Xmol(i)=max(Xmol(i),xmin);
      xtot+=Xmol(i);
      gamma(i)=1.0;
    }

  for (i=0;i<nions;i++)
    {      
      Xions(i)=max(Xions(i),xmin);
      xtot+=Xions(i);
      gamma_ions(i)=1.0;
    }
  
  if (xtot>0.0)
    {
      for (i=0;i<nmols;i++)
        Xmol(i)/=xtot;
      
      for (i=0;i<nions;i++)
        Xions(i)/=xtot;
  
      Array <double, 1> theta,theta_ions,phi,phi_ions;
      theta.resize(nmols);
      phi.resize(nmols);
      theta_ions.resize(nions);
      phi_ions.resize(nions);
      double Lmean=0.0;
      double sumtheta=0.0;
      double sumphi=0.0;
      for (i=0;i<nmols;i++)
        {
          theta(i)=Xmol(i)*Qparam(i);
          sumtheta+=Xmol(i)*Qparam(i);
          phi(i)=Xmol(i)*Rparam(i);
          sumphi+=Xmol(i)*Rparam(i);
          Lmean+=Xmol(i)*Lparam(i);
        }

      for (i=0;i<nions;i++)
        {
          theta_ions(i)=Xions(i)*QGions(i);
          sumtheta+=Xions(i)*QGions(i);
          phi_ions(i)=Xions(i)*RGions(i);
          sumphi+=Xions(i)*RGions(i);
          Lmean+=Xions(i)*Lions(i);
        }
      
      //Computation of residual gamma
      // Assumes that SR interactions with involves ions are zero as in UNIFAC and LIFAC
      for (i=0;i<nmols;i++)
        {
          theta(i)/=sumtheta;
          phi(i)/=sumphi;
          gamma(i)=exp(log(phi(i)/Xmol(i))+Z/2*Qparam(i)*log(theta(i)/phi(i))+Lparam(i)-phi(i)/Xmol(i)*Lmean);         
        }

      for (i=0;i<nions;i++)
        {
          theta_ions(i)/=sumtheta;
          phi_ions(i)/=sumphi;
          gamma_ions(i)=exp(log(phi_ions(i)/Xions(i))+Z/2*QGions(i)*log(theta_ions(i)/phi_ions(i))+Lions(i)-phi_ions(i)/Xions(i)*Lmean)/
                exp(log(RGions(i)/Rparam(nmols-1))+1.0-RGions(i)/Rparam(nmols-1)+Z/2*QGions(i)*(log(Rparam(nmols-1)*QGions(i)/RGions(i)/Qparam(nmols-1))
                    -1.0+RGions(i)*Qparam(nmols-1)/Rparam(nmols-1)/QGions(i)));
        }
      
      //Compute surface fraction
      Array <double, 1> surface_fraction,surface_fraction_ions;
      Array <double, 2> surface_fraction_mol;
      surface_fraction.resize(nfunc);
      surface_fraction_ions.resize(nions);
      surface_fraction_mol.resize(nfunc,nmols);      
      surface_fraction=0.0;
      surface_fraction_mol=0.0;
      double sum_surf=0.0;
      for (i=0;i<nmols;i++)
        {
          double sum_surf_mol=0.0;
          for (j=0;j<nfunc;j++)
            if (groups(j,i)>0.0)
              {
                surface_fraction(j)+=QG(j)*Xmol(i)*groups(j,i);            
                surface_fraction_mol(j,i)+=QG(j)*groups(j,i);
                sum_surf+=QG(j)*Xmol(i)*groups(j,i);  
                sum_surf_mol+=QG(j)*groups(j,i);  
              }
          
          for (j=0;j<nfunc;j++)      
            surface_fraction_mol(j,i)/=sum_surf_mol;
        }

      for (i=0;i<nions;i++) 
        {
          surface_fraction_ions(i)=QGions(i)*Xions(i);
          sum_surf+=QGions(i)*Xions(i);         
        }
    
      for (j=0;j<nfunc;j++)      
        surface_fraction(j)/=sum_surf;

      for (i=0;i<nions;i++)         
        surface_fraction_ions(i)/=sum_surf;
      
      Array<double,2> Inter2;
      Inter2.resize(nfunc,nfunc);
      if (temperature_dependancy)
	for (j=0;j<nfunc;j++)
	  for (k=0;k<nfunc;k++)          
	    Inter2(j,k)=exp(-Inter(j,k)/Temperature+InterB(j,k)*(1.0/298.15-1.0/Temperature)+InterC(j,k)*(298.15/Temperature-1.0+log(Temperature/298.15)));         
      else
	for (j=0;j<nfunc;j++)
	  for (k=0;k<nfunc;k++)          
	    Inter2(j,k)=exp(-Inter(j,k)/Temperature);       
      
      Array <double, 1> group_activity;
      group_activity.resize(nfunc);
      Array <double, 2> group_activity_mol,sum2mol;
      group_activity_mol.resize(nfunc,nmols);
      sum2mol.resize(nfunc,nmols);
      
      Array <double, 1> sum2;
      sum2.resize(nfunc);            
      for (j=0;j<nfunc;j++)
        {
          sum2(j)=0.0;
          for (k=0;k<nfunc;k++)
            sum2(j)+=surface_fraction(k)*Inter2(k,j);

          for (k=0;k<nions;k++)
            sum2(j)+=surface_fraction_ions(k);
        }

      for (j=0;j<nfunc;j++)
        {
          group_activity(j)=1.0-log(sum2(j));                             
          for (k=0;k<nfunc;k++)            
            group_activity(j)-=surface_fraction(k)*Inter2(j,k)/sum2(k);           
          
          for (k=0;k<nions;k++)
            group_activity(j)-=surface_fraction_ions(k);    
        }

      for (i=0;i<nmols;i++)
        for (j=0;j<nfunc;j++)
          if (groups(j,i)>0.0)
            {
              sum2mol(j,i)=0.0;
              for (k=0;k<nfunc;k++)
                if (groups(k,i)>0.0)
                  sum2mol(j,i)+=surface_fraction_mol(k,i)*Inter2(k,j);              
            }

      for (i=0;i<nmols;i++)
        for (j=0;j<nfunc;j++)
          if (groups(j,i)>0.0)
            {
              group_activity_mol(j,i)=1.0-log(sum2mol(j,i));                     
              for (k=0;k<nfunc;k++)      
                if (groups(k,i)>0.0)
                  group_activity_mol(j,i)-=surface_fraction_mol(k,i)*Inter2(j,k)/sum2mol(k,i);           
            }
      
      for (i=0;i<nmols;i++)
        {
          double sum1=0.0;
          for (j=0;j<nfunc;j++)
            if (groups(j,i)>0.0)           
              sum1+=groups(j,i)*QG(j)*(group_activity(j)-group_activity_mol(j,i));                

          gamma(i)*=exp(sum1);
        }
    }
}

void tab_unifac_aq(int &nmols, int &nions, int &nfunc, Array<double, 2> &groups, Array<double, 1> &Xmol,
                   Array<double, 1> &Xions,
                   Array<double, 2> &Inter, Array<double, 2> &InterB, Array<double, 2> &InterC, 
		   Array<double,1> &RG, Array <double, 1> &QG,
                   Array<double, 1> &Rparam, Array <double, 1> &Qparam, Array <double,1> &Lparam,
                   Array <double, 1> & RGions, Array <double, 1> & QGions, Array <double, 1> & Lions,
                   Array<double, 1> & spec_activity, 
                   double &Z, double &Temperature, Array<double,1> &gamma, Array<double,1> &gamma_ions,
		   bool temperature_dependancy)
{
  
  double xmin=1.0e-11;
  double xtot=0.0;
  int i,j,k;
  for (i=0;i<nmols;i++)
    {      
      Xmol(i)=max(Xmol(i),xmin);
      xtot+=Xmol(i);
      gamma(i)=1.0;
    }

  for (i=0;i<nions;i++)
    {      
      Xions(i)=max(Xions(i),xmin);
      xtot+=Xions(i);
      gamma_ions(i)=1.0;
    }
  
  if (xtot>0.0)
    {
      for (i=0;i<nmols;i++)
        Xmol(i)/=xtot;
      
      for (i=0;i<nions;i++)
        Xions(i)/=xtot;
  
      Array <double, 1> theta,theta_ions,phi,phi_ions;
      theta.resize(nmols);
      phi.resize(nmols);
      theta_ions.resize(nions);
      phi_ions.resize(nions);
      double Lmean=0.0;
      double sumtheta=0.0;
      double sumphi=0.0;
      for (i=0;i<nmols;i++)
        {
          theta(i)=Xmol(i)*Qparam(i);
          sumtheta+=Xmol(i)*Qparam(i);
          phi(i)=Xmol(i)*Rparam(i);
          sumphi+=Xmol(i)*Rparam(i);
          Lmean+=Xmol(i)*Lparam(i);
        }

      for (i=0;i<nions;i++)
        {
          theta_ions(i)=Xions(i)*QGions(i);
          sumtheta+=Xions(i)*QGions(i);
          phi_ions(i)=Xions(i)*RGions(i);
          sumphi+=Xions(i)*RGions(i);
          Lmean+=Xions(i)*Lions(i);
        }
      
      //Computation of residual gamma
      // Assumes that SR interactions with involves ions are zero as in UNIFAC and LIFAC
      for (i=0;i<nmols;i++)
        {
          theta(i)/=sumtheta;
          phi(i)/=sumphi;
          gamma(i)=exp(log(phi(i)/Xmol(i))+Z/2*Qparam(i)*log(theta(i)/phi(i))+Lparam(i)-phi(i)/Xmol(i)*Lmean);          	  
        }

      for (i=0;i<nions;i++)
        {
          theta_ions(i)/=sumtheta;
          phi_ions(i)/=sumphi;
          gamma_ions(i)=exp(log(phi_ions(i)/Xions(i))+Z/2*QGions(i)*log(theta_ions(i)/phi_ions(i))+Lions(i)-phi_ions(i)/Xions(i)*Lmean)/
                exp(log(RGions(i)/Rparam(nmols-1))+1.0-RGions(i)/Rparam(nmols-1)+Z/2*QGions(i)*(log(Rparam(nmols-1)*QGions(i)/RGions(i)/Qparam(nmols-1))
                    -1.0+RGions(i)*Qparam(nmols-1)/Rparam(nmols-1)/QGions(i)));
        }                            

      //Compute surface fraction
      Array <double, 1> surface_fraction,surface_fraction_ions;      
      surface_fraction.resize(nfunc);    
      surface_fraction_ions.resize(nions);
      surface_fraction=0.0;      
      double sum_surf=0.0;
      for (i=0;i<nmols;i++)
        {         
          for (j=0;j<nfunc;j++)
            if (groups(j,i)>0.0)
              {
                surface_fraction(j)+=QG(j)*Xmol(i)*groups(j,i);                            
                sum_surf+=QG(j)*Xmol(i)*groups(j,i);                  
              }          
        }

      for (i=0;i<nions;i++) 
        {
          surface_fraction_ions(i)=QGions(i)*Xions(i);
          sum_surf+=QGions(i)*Xions(i);         
        }

      for (j=0;j<nfunc;j++)      
        surface_fraction(j)/=sum_surf;

      for (i=0;i<nions;i++)         
        surface_fraction_ions(i)/=sum_surf;
      
      Array<double,2> Inter2;
      Inter2.resize(nfunc,nfunc);    
      if (temperature_dependancy)
	for (j=0;j<nfunc;j++)
	  for (k=0;k<nfunc;k++)          
	    Inter2(j,k)=exp(-Inter(j,k)/Temperature+InterB(j,k)*(1.0/298.15-1.0/Temperature)+InterC(j,k)*(298.15/Temperature-1.0+log(Temperature/298.15)));         
      else
	for (j=0;j<nfunc;j++)
	  for (k=0;k<nfunc;k++)          
	    Inter2(j,k)=exp(-Inter(j,k)/Temperature); 
      
      Array <double, 1> group_activity;
      group_activity.resize(nfunc);      
      
      Array <double, 1> sum2;
      sum2.resize(nfunc);      
      for (j=0;j<nfunc;j++)
        {
          sum2(j)=0.0;
          for (k=0;k<nfunc;k++)
            sum2(j)+=surface_fraction(k)*Inter2(k,j);

          for (k=0;k<nions;k++)
            sum2(j)+=surface_fraction_ions(k);
        }

      for (j=0;j<nfunc;j++)
        {         
          group_activity(j)=1.0-log(sum2(j));                             
          for (k=0;k<nfunc;k++)            
            group_activity(j)-=surface_fraction(k)*Inter2(j,k)/sum2(k);      

          for (k=0;k<nions;k++)
            group_activity(j)-=surface_fraction_ions(k);                            
        }
      
      for (i=0;i<nmols;i++)
        {
          double sum1=0.0;
          for (j=0;j<nfunc;j++)
            if (groups(j,i)>0.0)            
              sum1+=groups(j,i)*QG(j)*group_activity(j);                                             
          gamma(i)*=exp(sum1-spec_activity(i));           
        }
    }
}

#endif
