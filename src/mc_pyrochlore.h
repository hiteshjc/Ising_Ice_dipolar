#ifndef MC_PYRO_HEADER
#define MC_PYRO_HEADER

#include"global.h"
#include"math_utilities.h"
#include<omp.h>

using namespace std; 

class QMC_Info
{
	public:
	int 			nsites;
	std::vector<double> 	configx, configy, configz;
        RMatrix 		sxsxcorrs, sysycorrs, szszcorrs;
      	RMatrix 		sxsycorrs, sxszcorrs, syszcorrs;
	RMatrix 		sxsxcorrstot, sysycorrstot, szszcorrstot;
	RMatrix 		sxsycorrstot, sxszcorrstot, syszcorrstot;
	double 			etot,e2tot,e4tot;
	double 			mxtot,mx2tot,mx4tot;
	double 			mytot,my2tot,my4tot;
	double 			mztot,mz2tot,mz4tot;
	double 			eavg,e2avg,e4avg;
	double 			mxavg,mx2avg,mx4avg;
	double 			myavg,my2avg,my4avg;
	double 			mzavg,mz2avg,mz4avg;
	double 			nmeas;
	double 			spheatpersitenodim; 
	double 			spheat;
	double			spheatpersite;
	double			temp,tempKelvin,beta;
	double 			energy,mx,my,mz;

	void init(int nsites, double tempKelvin)
	{
		this->nsites=nsites;
		this->tempKelvin=tempKelvin;
        	this->temp=tempKelvin*0.08621738;
		this->beta=1.0/temp;
		/*this->sxsxcorrstot.resize(nsites,nsites);
		this->sysycorrstot.resize(nsites,nsites);
		this->szszcorrstot.resize(nsites,nsites);
	
		this->sxsycorrstot.resize(nsites,nsites);
		this->sxszcorrstot.resize(nsites,nsites);
		this->syszcorrstot.resize(nsites,nsites);

		#pragma omp parallel for	
		for (int s=0;s<nsites*nsites;s++) 
		{
			this->sxsxcorrstot[s]=0.0;this->sysycorrstot[s]=0.0;this->szszcorrstot[s]=0.0;
			this->sxsycorrstot[s]=0.0;this->sxszcorrstot[s]=0.0;this->syszcorrstot[s]=0.0;
		}*/
		this->etot=0.0;this->e2tot=0.0;this->e4tot=0.0;
		this->mxtot=0.0;this->mx2tot=0.0;this->mx4tot=0.0;
		this->mytot=0.0;this->my2tot=0.0;this->my4tot=0.0;
		this->mztot=0.0;this->mz2tot=0.0;this->mz4tot=0.0;
		this->nmeas=0;
	}
	
	void update_totals()
	{
		this->etot+=this->energy;this->e2tot+=pow(this->energy,2.0);this->e4tot+=pow(this->energy,4.0);
		this->mxtot+=this->mx;this->mytot+=this->my;this->mztot+=this->mz;
		this->mx2tot+=pow(this->mx,2.0);this->my2tot+=pow(this->my,2.0);this->mz2tot+=pow(this->mz,2.0);
		this->mx4tot+=pow(this->mx,4.0);this->my4tot+=pow(this->my,4.0);this->mz4tot+=pow(this->mz,4.0);
		this->nmeas+=1;
	}

	void update_total_correlations()
	{
		#pragma omp parallel for
		for (int s=0;s<nsites*nsites;s++) this->sxsxcorrstot[s]+=this->sxsxcorrs[s];
		#pragma omp parallel for
		for (int s=0;s<nsites*nsites;s++) this->sysycorrstot[s]+=this->sysycorrs[s];
		#pragma omp parallel for
		for (int s=0;s<nsites*nsites;s++) this->szszcorrstot[s]+=this->szszcorrs[s];
		#pragma omp parallel for
		for (int s=0;s<nsites*nsites;s++) this->sxsycorrstot[s]+=this->sxsycorrs[s];
		#pragma omp parallel for
		for (int s=0;s<nsites*nsites;s++) this->sxszcorrstot[s]+=this->sxszcorrs[s];
		#pragma omp parallel for
		for (int s=0;s<nsites*nsites;s++) this->syszcorrstot[s]+=this->syszcorrs[s];
	}

	void average()
	{
		//cout<<"Nmeas = "<<this->nmeas<<endl;
		this->eavg=this->etot/this->nmeas;this->e2avg=this->e2tot/this->nmeas;this->e4avg=this->e4tot/this->nmeas;
		this->mxavg=this->mxtot/this->nmeas;this->myavg=this->mytot/this->nmeas;this->mzavg=this->mztot/this->nmeas;
		this->mx2avg=this->mx2tot/this->nmeas;this->my2avg=this->my2tot/this->nmeas;this->mz2avg=this->mz2tot/this->nmeas;
		this->mx4avg=this->mx4tot/this->nmeas;this->my4avg=this->my4tot/this->nmeas;this->mz4avg=this->mz4tot/this->nmeas;
		this->spheatpersitenodim=(this->e2avg-(this->eavg*this->eavg))/(this->temp*this->temp*double(this->nsites));
		this->spheat=(1119.67107046)*(this->e2avg-(this->eavg*this->eavg))/(this->tempKelvin*this->tempKelvin);
		this->spheatpersite=this->spheat/double(this->nsites);
		/*for (int s=0;s<nsites*nsites;s++) this->sxsxcorrstot[s]=this->sxsxcorrstot[s]/double(this->nmeas);
		for (int s=0;s<nsites*nsites;s++) this->sysycorrstot[s]=this->sysycorrstot[s]/double(this->nmeas);
		for (int s=0;s<nsites*nsites;s++) this->szszcorrstot[s]=this->szszcorrstot[s]/double(this->nmeas);
		for (int s=0;s<nsites*nsites;s++) this->sxsycorrstot[s]=this->sxsycorrstot[s]/double(this->nmeas);
		for (int s=0;s<nsites*nsites;s++) this->sxszcorrstot[s]=this->sxszcorrstot[s]/double(this->nmeas);
		for (int s=0;s<nsites*nsites;s++) this->syszcorrstot[s]=this->syszcorrstot[s]/double(this->nmeas);*/
	}
	
};

void iterative_pyrochlore(double spin, int L,int nsamples, int nburn, string start_config, 
		   	  double hx, double hy, double hz, 
		   	  double J1, double J2, double J3, double J4, double Jnnn,
		   	  double disorder_strength,
		   	  double gxy, double gz, 
		   	  double & eavg, 
		   	  double &mxavg, double &myavg, double &mzavg, double &e2avg, 
		   	  double &mx2avg, double &my2avg, double &mz2avg, bool &measure_corrs);

void mc_pyrochlore(double spin, int L,int nsamples, int nburn, string start_config, 
		   string mcmove, double temp, double hx, double hy, double hz, 
		   double J1, double J2, double J3, double J4, double Jnnn, double disorder,
		   double gxy, double gz, 
		   double & eavg, 
		   double &mxavg, double &myavg, double &mzavg, double &e2avg, 
		   double &mx2avg, double &my2avg, double &mz2avg, bool &measure_corrs);


void mc_pyrochlore_ground_state(double spin, int L, int nburn, string start_config, 
		   string mcmove, double temp, double hx, double hy, double hz, 
		   double J1, double J2, double J3, double J4, double Jnnn,
		   double gxy, double gz, 
		   double & eavg, 
		   double &mxavg, double &myavg, double &mzavg, double &e2avg, 
		   double &mx2avg, double &my2avg, double &mz2avg, bool &measure_corrs);

void mc_pyrochlore_pt(double spin, int L, int64_t nsamples, int64_t nburn, 
		   double temp, int ntemps, double hx, double hy, double hz, 
		   double J1, double J2, double D, double g, double alpha, int kcut, double loop_prob, double muprime);
#endif
