#include"mc_sphere.h"
#include"mc_finite_D.h"
#include"printing_functions.h"
#include"matrix_functions.h"
#include"number_functions.h"

using namespace std;

////////////////////////////////////////////////////////////////////////
void make_sphere(double R,
		 std::vector< std::vector<double> > &fullcoords,
		 std::vector< std::vector<int> > &neighbors)
{
	fullcoords.clear();
	neighbors.clear();
        int site=0;
	for (int n1=-R;n1<=R;n1++)
	{
		for (int n2=-R;n2<=R;n2++)
		{
			for (int n3=R;n3<=R;n3++)
			{
				double x = double(n1);
				double y = double(n2);
				double z = double(n3);
				if (x*x + y*y + z*z <= (R*R))
				{
					std::vector<double> fullcoordsentry;
					fullcoordsentry.push_back(x);fullcoordsentry.push_back(y);fullcoordsentry.push_back(z);
					fullcoords.push_back(fullcoordsentry);
					neighbors.push_back(std::vector<int>());
					site+=1;
				}
			}
		}
	}
	
	int nsites=site;
	cout<<"Number of sites recorded is = "<<nsites<<endl;
	for (int i=0;i<nsites;i++)
	{
		double xi=fullcoords[i][0];
		double yi=fullcoords[i][1];
		double zi=fullcoords[i][2];
		for (int j=0;j<nsites;j++)
		{
			double xj=fullcoords[j][0];
			double yj=fullcoords[j][1];
			double zj=fullcoords[j][2];
			double dx=xi-xj;
			double dy=yi-yj;
			double dz=zi-zj;
			double dist=sqrt(dx*dx + dy*dy + dz*dz);  // distance between points 
			if (dist>0.99 and dist<1.01)  neighbors[i].push_back(j);  //      nearest neighbor on cube
		}
	}
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double total_jdipolarDh_energy(double &spin,
		      double J,
		      double Dip,
		      double D,
		      double hx,
		      double hy,
		      double hz,
		      std::vector<double> &configx, 
		      std::vector<double> &configy, 
		      std::vector<double> &configz, 
		      std::vector< std::vector<int> > &neighbors, 
		      std::vector< std::vector<double> > &fullcoords) 
{
	double e1=0.0;
	int nsites=neighbors.size();
	for (int i=0;i<nsites;i++)
	{
		for (int j=0;j<neighbors[i].size();j++)
		{
			int k=neighbors[i][j];
			e1+=-J*((configx[i]*configx[k])+(configy[i]*configy[k])+(configz[i]*configz[k]));			
		}
	}
	for (int i=0;i<nsites;i++)
	{
		for (int j=0;j<nsites;j++)
		{
			if (i!=j)
			{
				double rijx=fullcoords[j][0]-fullcoords[i][0];
				double rijy=fullcoords[j][1]-fullcoords[i][1];
				double rijz=fullcoords[j][2]-fullcoords[i][2];
				double r=sqrt(rijx*rijx + rijy*rijy + rijz*rijz);
				double dot1=(configx[i]*rijx + configy[i]*rijy + configz[i]*rijz);
				double dot2=(configx[j]*rijx + configy[j]*rijy + configz[j]*rijz);
				e1+=(Dip*((configx[i]*configx[j])+(configy[i]*configy[j])+(configz[i]*configz[j]))/(r*r*r));			
				e1-=(Dip*(3.0*dot1*dot2)/(r*r*r*r*r));
			}			
		}
	}
	e1=(e1/2.0)*spin*spin;  // Pairs counted twice so divide by 2
	
	double e2=0.0;
	for (int i=0;i<nsites;i++) e2-=((hx*configx[i])+(hy*configy[i])+(hz*configz[i]));
	e2=e2*spin;
	
	double e3=0.0;
	for (int i=0;i<nsites;i++) e3-=(pow(configx[i],4.0)+pow(configy[i],4.0)+pow(configz[i],4.0));
	e3=e3*D*spin*spin*spin*spin;

	return e1+e2+e3;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<double> total_magnetization_no_g(double 				&spin,
					std::vector<double> 			&configx, 
					std::vector<double> 			&configy, 
					std::vector<double> 			&configz)
{
	double mx=0;
	double my=0;
	double mz=0;
	std::vector<double> m;
	int nsites=int(configx.size());
	for (int n=0;n<nsites;n++)
	{
		double sx=configx[n];
		double sy=configy[n];
		double sz=configz[n];
		mx+=(sx);
		my+=(sy);
		mz+=(sz);
	}
	m.push_back(mx*spin);
	m.push_back(my*spin);
	m.push_back(mz*spin);
	return m;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double local_jdipolarDh_energy(double &spin,
		      int site,
		      double sx,
		      double sy,
		      double sz,
		      double J,
		      double Dip,
		      double D,
		      double hx,
		      double hy,
		      double hz,
		      std::vector<double> &configx, 
		      std::vector<double> &configy, 
		      std::vector<double> &configz, 
		      std::vector< std::vector<int> > &neighbors, 
		      std::vector< std::vector<double> > &fullcoords) 
{
	double e1=0.0;
	int nsites=neighbors.size();
	for (int j=0;j<neighbors[site].size();j++)
	{
			int k=neighbors[site][j];
			e1+=-J*((sx*configx[k])+(sy*configy[k])+(sz*configz[k]));			
	}
		
	for (int j=0;j<nsites;j++)
	{
			if (site!=j)
			{
				double rijx=fullcoords[j][0]-fullcoords[site][0];
				double rijy=fullcoords[j][1]-fullcoords[site][1];
				double rijz=fullcoords[j][2]-fullcoords[site][2];
				double r=sqrt(rijx*rijx + rijy*rijy + rijz*rijz);
				double dot1=(sx*rijx + sy*rijy + sz*rijz);
				double dot2=(configx[j]*rijx + configy[j]*rijy + configz[j]*rijz);
				e1+=(Dip*((sx*configx[j])+(sy*configy[j])+(sz*configz[j]))/(r*r*r));			
				e1-=(Dip*(3.0*dot1*dot2)/(r*r*r*r*r));
			}			
	}
	e1=(e1)*spin*spin;  // Pairs counted twice so divide by 2
	
	double e2=0.0;
	e2-=((hx*sx)+(hy*sy)+(hz*sz));
	e2=e2*spin;
	
	double e3=0.0;
	e3-=(pow(sx,4.0)+pow(sy,4.0)+pow(sz,4.0));
	e3=e3*D*spin*spin*spin*spin;

	return e1+e2+e3;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void mc_sphere(string lattice, double spin, int R, int64_t nsamples, int64_t nburn, string start_config, 
		   string mcmove, double temp, int ntemps, double hx, double hy, double hz, 
		   double J, double Dip, double D,
		   double & eavg, 
		   double &mxavg, double &myavg, double &mzavg, double &e2avg, 
		   double &mx2avg, double &my2avg, double &mz2avg, bool &measure_corrs)
{
	/////////////////////////////////////////////////////////////////////////
	// MC related quantities
	// Units - Assume J in meV (millielectron volts)
        //         Assume h in T   (tesla)
        //         Convert h to meV units in total_h_energy and local_h_energy
        // Convert temperature in Kelvin to temperature in meV
	// Set 16xLxLxL pyrochlore lattice
	int nsites;
       	double kB=1.38064852*1e-23;
        double NA=6.02214179*1e23;
        double JpermeV=1.60218*1e-22; 
	double tempKelvin=temp;
        temp=tempKelvin*0.08621738;
	/////////////////////////////////////////////////////////////////////////
	// MC related quantities
	// Set nburn
	nburn=nsamples;
	cout<<"Nsamples = "<<nsamples<<endl;
	cout<<"Nburn    = "<<nburn<<endl;
	double e4avg,mx4avg,my4avg,mz4avg;
	
	/////////////////////////////////////////////////////////////////////////
	std::vector< std::vector<int> >    neighbors;	
	std::vector< std::vector<double> > fullcoords;	
        make_sphere(R, fullcoords, neighbors);
	
	/////////////////////////////////////////////////////////////////////////
	// Make random configuration of spins or selected type 
	std::vector<QMC_Info> infos;
	cout<<"Ntemps   = "<<ntemps<<endl;
	double exponent=pow(100.0,1.0/double(ntemps)); // Tmax/Tmin=30
	cout<<"exponent = "<<exponent<<endl;
	for (int b=0;b<ntemps;b++)
	{
		QMC_Info qmc;
		qmc.init(nsites,tempKelvin*pow(exponent,double(b)));
		infos.push_back(qmc);
		cout<<"Beta = "<<infos[b].beta<<endl;
		make_random_config(nsites,infos[b].configx,infos[b].configy,infos[b].configz); 
	}

	
	for (int b=0;b<ntemps;b++)
	{
		double 		    energy=total_jdipolarDh_energy(spin, J, Dip, D, hx, hy, hz, 
								    infos[b].configx,infos[b].configy,infos[b].configz, 
								    neighbors,fullcoords);
		std::vector<double> magnetization=total_magnetization_no_g(spin, infos[b].configx, infos[b].configy, infos[b].configz);
		infos[b].energy=energy;
		infos[b].mx=magnetization[0];infos[b].my=magnetization[1];infos[b].mz=magnetization[2];
		cout<<"Total J energy + total h +D energy ="<<infos[b].energy<<endl;
		cout<<"Total magnetization is             ="<<endl; print_vec_acc(magnetization,true);
	}

	cout<<endl;	
	cout<<endl;	
	cout<<"========================================================"<<endl;
	cout<<"Energy history "<<endl;	
	/////////////////////////////////////////////////////////////////////////
	// Accept reject Metropolis
	double accept=0.0;
	double reject=0.0;
	double nreplicatries=0.0;
	double nswaps=0.0;
        for (int64_t n=0; n<(nsamples+(nburn));n++)
	{
	   if (n%nsites!=0 or ntemps==1) // Do usual Metropolis MC
	   {
			///////////////////////////////////////////////////////////////////////
			// Usual Moves of a serial Metropolis Monte Carlo
			///////////////////////////////////////////////////////////////////////
			std::vector<double> rnd1,rnd2,rnd3,rnd4;
			std::vector<int> rndints;
		        // random numbers generated in advance
			for (int b=0;b<ntemps;b++)
			{
				rnd1.push_back(uniform_rnd());
				rnd3.push_back(uniform_rnd());
				rnd4.push_back(uniform_rnd());
				rndints.push_back(uniform_rand_int(0,nsites));
			}
			# pragma omp parallel for
			for (int b=0;b<ntemps;b++)
			{
					// Very Small moves needed at low temperatures to increase acceptance rates
					double move_size=min(0.3,0.1*(infos[0].beta)/(infos[b].beta));  
					int site=rndints[b];
					double sxnew,synew,sznew;
					double sx=infos[b].configx[site];double sy=infos[b].configy[site];double sz=infos[b].configz[site];						if (mcmove=="conical") conical_move_continuous_spin_rnds_provided(move_size,rnd1[b],rnd2[b],sx,sy,sz,sxnew,synew,sznew);			
					// Normalize new direction
					double norm=sqrt(sxnew*sxnew + synew*synew + sznew*sznew); 
					sxnew=sxnew/norm; synew=synew/norm; sznew=sznew/norm;

					double local_energy1=local_jdipolarDh_energy(spin,site, sx,sy,sz, J, Dip, D, hx,hy, hz,
						       infos[b].configx,infos[b].configy,infos[b].configz,neighbors,fullcoords); 
					double local_energy2=local_jdipolarDh_energy(spin,site, sxnew,synew,sznew, J, Dip, D, hx,hy, hz,
						       infos[b].configx,infos[b].configy,infos[b].configz,neighbors,fullcoords);
					double mxdiff=(sxnew-sx)*spin;	
					double mydiff=(synew-sy)*spin;
					double mzdiff=(sznew-sz)*spin;
					double ediff=(local_energy2-local_energy1);
					double beta; 
					beta=infos[b].beta;
					double prob=exp(-beta*ediff);
					double rand=rnd3[b]; // random number previously generated
					if (rand<prob) // Metropolis Accept-reject for a given temperature
					{
						if (b==0 and n>(nburn)) accept=accept+1.0;
						infos[b].configx[site]=sxnew;infos[b].configy[site]=synew;infos[b].configz[site]=sznew;
						infos[b].energy+=ediff;
						infos[b].mx+=mxdiff;infos[b].my+=mydiff;infos[b].mz+=mzdiff;
					}
					else
					{
					       	if (b==0 and n>(nburn)) reject=reject+1.0;
					}
			}
	}	
	else if (ntemps>1)  //50 % chance of doing replica exchange
	{
		nreplicatries+=1.0;	
		///////////////////////////////////////////////////////////////////////
		////////// Attempted exchange moves of parallel tempering
		///////////////////////////////////////////////////////////////////////
		// Metropolis move done, try swapping every 2 sweeps
		for (int which1=0;which1<ntemps;which1++)
		{
			int which2=(which1+1)%(ntemps);
			// Slight bias, if i try to swap serially ?
			double rand=uniform_rnd();	
			double beta_i=infos[which1].beta;
			double beta_j=infos[which2].beta;
			double energy_i=infos[which1].energy;
			double energy_j=infos[which2].energy;
			double power=(beta_j-beta_i)*(energy_i - energy_j);
			double ratio=exp(-power);
			if (rand<ratio)
			{
				// SWAP quantities which are being saved
				std::vector<double> temp;
				double temp_energy, temp_mx, temp_my, temp_mz;
				temp=infos[which1].configx;	
				infos[which1].configx=infos[which2].configx;	
				infos[which2].configx=temp;	
				
				temp=infos[which1].configy;	
				infos[which1].configy=infos[which2].configy;	
				infos[which2].configy=temp;	

				temp=infos[which1].configz;	
				infos[which1].configz=infos[which2].configz;	
				infos[which2].configz=temp;	
				
				temp_energy=infos[which1].energy;	
				infos[which1].energy=infos[which2].energy;	
				infos[which2].energy=temp_energy;	
				
				temp_mx=infos[which1].mx;	
				infos[which1].mx=infos[which2].mx;	
				infos[which2].mx=temp_mx;	
				
				temp_my=infos[which1].my;	
				infos[which1].my=infos[which2].my;	
				infos[which2].my=temp_my;	
				
				temp_mz=infos[which1].mz;	
				infos[which1].mz=infos[which2].mz;	
				infos[which2].mz=temp_mz;
				if (which1==0 or which2==0) nswaps+=1.0;	
				//nswaps+=1.0;	
			}
		}
	}
	if (n%nsites==0 or n==nburn+nsamples-1) cout<<boost::format("%+ .5f") %infos[0].energy<<endl;
	}
	cout<<"========================================================"<<endl;
	cout<<endl;
	accept=accept/(accept+reject);
	nswaps=nswaps/nreplicatries;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

