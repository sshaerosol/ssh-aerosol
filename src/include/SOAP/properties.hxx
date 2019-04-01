#ifndef PROPERTIES_HXX

#define _unifac unifac_
#define _isoropia isoropia_

  extern "C"
  {
    void _unifac(int*, int*, double*, double*, double*, double*, double*,
				 double*, double*, double*);
 
    void _isoropia(double*,double*,double*,double*,double*,double*,
                   double*,double*, double*,double*, double*);     

  }

#define PROPERTIES_HXX
#endif
