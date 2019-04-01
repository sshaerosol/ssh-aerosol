//////////////
// INCLUDES //

#include <vector>
#include "species.h"
#include "solving.cxx"
using namespace soap;

void soap_main(double LWC, double RH, double Temperature, 
               double ionic, double chp, double& LWCorg,
               double& deltat,
               double DSD[], double csol[], double liquid[],
               int ns_aer, int neq, double q[], double qaero[], double qgas[],
               double lwc_Nsize[], double ionic_Nsize[], double chp_Nsize[],               double liquid_Nsize[], int nbin);

