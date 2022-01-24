#ifndef MC_SPHERE_HEADER
#define MC_SPHERE_HEADER

#include"global.h"
#include"math_utilities.h"
#include"mc_pyrochlore.h"

void mc_sphere(string lattice, double spin, int R, int64_t nsamples, int64_t nburn, string start_config, 
		   string mcmove, double temp, int ntemps, double hx, double hy, double hz, 
		   double J, double Dip, double D,
		   double & eavg, 
		   double &mxavg, double &myavg, double &mzavg, double &e2avg, 
		   double &mx2avg, double &my2avg, double &mz2avg, bool &measure_corrs);
#endif
