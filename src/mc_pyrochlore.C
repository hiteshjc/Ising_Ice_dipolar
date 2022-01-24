#include"mc_pyrochlore.h"
#include"mc_finite_D.h"
#include"printing_functions.h"
#include"matrix_functions.h"
#include"number_functions.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////			
void check_surroundings(std::vector< std::vector<int> > &neighbors,int site, std::vector<double> &configz, 
			std::vector<int> &touched_sites, bool &aborted)
{
	aborted=false;
	int type1=0;
	int type2=0;
	// Look for a neighbor with opposite spin and on a different up tetrahedron		
	for (int i=0;i<neighbors[site].size();i++)
	{
		int nbr=neighbors[site][i];
		if ( (nbr/4 != site/4) and abs(configz[nbr]-configz[site])<0.001)  // neighbors on different up tetrahedra and have specified spin
		{
			type1+=1;
		}
		if ( (nbr/4 == site/4) and abs(configz[nbr]-configz[site])<0.001)  // neighbors on same up tetrahedra and have specified spin
		{
			type2+=1;
		}
	}
	if (type1!=1 and type2!=1) 
	{
		//touched_sites.clear();
		//touched_sites.push_back(site);
		aborted = true;
	}
}
///////////////////////////////////////////////////////////////////////////			
void spin_nbrs_on_different_up_tet(std::vector< std::vector<int> > &neighbors,int currentsite,double inputsz,std::vector<double> &configz,std::vector<int> &tmp_nbrs, 
				   std::vector<int> &tempsites, std::vector<int> &touched_sites, bool &aborted)
{
	tmp_nbrs.clear();
	aborted=false;
	// Look for a neighbor with opposite spin and on a different up tetrahedron		
	for (int i=0;i<neighbors[currentsite].size();i++)
	{
		int nbr=neighbors[currentsite][i];
		if ( (nbr/4 != currentsite/4) and abs(configz[nbr]-inputsz)<0.001)  // neighbors on different up tetrahedra and have specified spin
		{
			tmp_nbrs.push_back(nbr);
		}
	}
	// If not 2 neighbors end the loop search
	if (tmp_nbrs.size()!=2) 
	{
		//touched_sites.clear();
		//touched_sites.push_back(tempsites[0]);
		aborted = true;
	}
}
///////////////////////////////////////////////////////////////////////////			
void spin_nbrs_on_same_up_tet(std::vector< std::vector<int> > &neighbors,int currentsite,double inputsz,std::vector<double> &configz,std::vector<int> &tmp_nbrs, 
				   std::vector<int> &tempsites, std::vector<int> &touched_sites, bool &aborted)
{
	tmp_nbrs.clear();
	aborted=false;
	// Look for a neighbor with opposite spin and same up tetrahedron		
	for (int i=0;i<neighbors[currentsite].size();i++)
	{
		int nbr=neighbors[currentsite][i];
		if ( (nbr/4==currentsite/4) and abs(configz[nbr]-inputsz)<0.001)  // neighbors on same up tetrahedra and have specified spin
		{
			tmp_nbrs.push_back(nbr);
		}
	}
	// If not 2 neighbors end the loop search
	if (tmp_nbrs.size()!=2) 
	{
		//touched_sites.clear();
		//touched_sites.push_back(tempsites[0]);
		aborted = true;
	}
}
///////////////////////////////////////////////////////////////////////////			

void check_loop_closed_on_adding_last_site(RMatrix &dmat, std::vector<int> &tempsites, bool &found_loop, int &closing_point, std::vector<int> &touched_sites)
{
	closing_point=-1;
	found_loop=false;
	int chosen=tempsites[tempsites.size()-1];
	// Check if a loop has been closed
	for (int i=0;i<tempsites.size()-2;i++) // Do not check last and second last site
	{
		if (dmat(chosen,tempsites[i])>0.353 and dmat(chosen,tempsites[i])<0.354) // NN condition 
		{
			found_loop=true; // If the new site is on the same up tetrahedron as some other site previously in the string
			closing_point=i;
			//cout<<"LOOP COMPLETE"<<endl;
			//i=tempsites.size(); // DO NOT quit for-loop, rather WAIT for the site which is part of the loop, and not the string
		}
	}
	
	if (found_loop==true)
	{
		touched_sites.clear();
		for (int i=closing_point;i<tempsites.size();i++) //"Erase" rest of the loop, and save only from the closing point to the final site 
		{
			touched_sites.push_back(tempsites[i]);
		}

	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_spinxyz_given_t(int t,double sz, double &spinx, double &spiny, double &spinz)
{
	double factor=1.0/sqrt(3.0);
	if (t==0)
	{
		spinx=sz*factor;
		spiny=sz*factor;
		spinz=sz*factor;
	}	
	if (t==1)
	{
		spinx=sz*factor;
		spiny=-sz*factor;
		spinz=-sz*factor;
	}	
	if (t==2)
	{
		spinx=-sz*factor;
		spiny=sz*factor;
		spinz=-sz*factor;
	}	
	if (t==3)
	{
		spinx=-sz*factor;
		spiny=-sz*factor;
		spinz=sz*factor;
	}	

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void total_dipolar_energy_pyrochlore(	    double 				coupling,
					    double 				alpha, 
					    std::vector<double> 		&configz, 
					    RMatrix				&Drmatrix,
					    RMatrix				&Dkmatrix,
					    double 				&total_e)
{
	double invrootpi=0.56418958354;
	int nsites=configz.size();
	double total_e_temp=0.0;
	//#pragma omp parallel for default(shared) reduction(+:total_e_temp)
	for (int i=0;i<nsites;i++)
	{
		for (int j=0;j<nsites;j++)          
		{
			total_e_temp+=(0.5*configz[i]*configz[j]*Drmatrix(i,j));
		}
	}
	for (int i=0;i<nsites;i++)
	{
		for (int j=0;j<nsites;j++)          // In this way of doing things the factor of 1/2 is not needed, the j=i term will be treated separately which comes from k space
		{
			total_e_temp+=(configz[i]*configz[j]*Dkmatrix(i,j));
		}
	}
	for (int i=0;i<nsites;i++) total_e_temp+=((-2.0/3.0)*invrootpi*alpha*alpha*alpha); //- a constant piece that is alpha dependent, total energy is alpha independent
	total_e=total_e_temp*coupling;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ediff_dipolar_pyrochlore(         double 				coupling, 
				       std::vector<int> 		&touched_sites,
				       std::vector<double> 		&configz, 
				       std::vector<double> 		&newconfigz, 
				       RMatrix				&Drmatrix,
				       RMatrix				&Dkmatrix,
				       double 				&ediff)
{
	ediff=0;
	int nsites=configz.size();
	double local_en=0.0;

	#pragma omp parallel for default(shared) reduction(+:local_en)
	for (int p=0;p<touched_sites.size();p++)
	{
		int site=touched_sites[p];	
		// i=j terms do not contribute
		//#pragma omp parallel for default(shared) reduction(+:local_en)
		for (int j=0;j<nsites;j++)
		{
			double factor=1.0;
			if (count(touched_sites.begin(),touched_sites.end(),j)==1) factor=0.5; // to avoid double counting of site pairs in touched_sites list
			double diff=(newconfigz[j]*newconfigz[site])-(configz[j]*configz[site]);
			local_en+=(diff*factor*(Drmatrix(site,j)+2.0*Dkmatrix(site,j)));       // Do not account for onsite term since it does not contribute to energy differences
											   // It does contribute to the total energy though (Depends), so BE careful 
		}
	}
	ediff=(local_en)*coupling;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ediff_nn_pyrochlore(              double 				coupling1, 
					double 				coupling2, 
				       std::vector<std::vector<int> >	&neighbors,
				       std::vector<std::vector<int> >	&nneighbors,
				       std::vector<int> 		&touched_sites,
				       std::vector<double> 		&configz, 
				       std::vector<double> 		&newconfigz, 
				       double 				&ediff)
{
	ediff=0;
	int nsites=configz.size();
	double local_en1=0.0;
	double local_en2=0.0;

	for (int p=0;p<touched_sites.size();p++)
	{
		int site=touched_sites[p];	
		//#pragma omp parallel for default(shared) reduction(+:local_en)
		for (int j=0;j<neighbors[site].size();j++)
		{
			double factor=1.0;
			int k=neighbors[site][j];
			if (count(touched_sites.begin(),touched_sites.end(),k)==1) factor=0.5; // to avoid double counting of site pairs in touched_sites list
			double diff=(newconfigz[k]*newconfigz[site])-(configz[k]*configz[site]);  // Double counting if touched sites are nearest neighbors, need to FIX!!
			local_en1+=(diff*factor);
		}
	}
	for (int p=0;p<touched_sites.size();p++)
	{
		int site=touched_sites[p];	
		//#pragma omp parallel for default(shared) reduction(+:local_en)
		for (int j=0;j<nneighbors[site].size();j++)
		{
			double factor=1.0;
			int k=nneighbors[site][j];
			if (count(touched_sites.begin(),touched_sites.end(),k)==1) factor=0.5; // to avoid double counting of site pairs in touched_sites list
			double diff=(newconfigz[k]*newconfigz[site])-(configz[k]*configz[site]);  // Double counting if touched sites are nearest neighbors, need to FIX!!
			local_en2+=(diff*factor);
		}
	}
	ediff=(local_en1)*coupling1*(1.0/3.0);
	ediff+=((local_en2)*coupling2*(1.0/3.0));
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ediff_and_mdiff_field_pyrochlore( std::vector< std::vector<int> >  &ijkt,  
				       double 				g, 
				       double 				hx, 
				       double 				hy, 
				       double 				hz, 
				       std::vector<int> 		&touched_sites,
				       std::vector<double> 		&configz, 
				       std::vector<double> 		&newconfigz, 
				       double 				&ediff, 
				       double 				&mxdiff,
				       double 				&mydiff,
				       double 				&mzdiff)
{
        double mu_b=5.7883818012*0.01; // meV/Tesla
	double f=1.0/sqrt(3.0);
	g=g/sqrt(3.0);
	ediff=0.0;
	mxdiff=0.0;
	mydiff=0.0;
	mzdiff=0.0;
	
	int nsites=configz.size();
	double local_en=0.0;

	for (int p=0;p<touched_sites.size();p++)
	{
		int site=touched_sites[p];	
		double diff=(newconfigz[site])-(configz[site]);
		int t=ijkt[site][3];
		if(t==0) 
		{
			local_en+=((hx+hy+hz)*diff);
			mxdiff+=(diff*f);
			mydiff+=(diff*f);
			mzdiff+=(diff*f);
		}
		if(t==1) 
		{
			local_en+=((hx-hy-hz)*diff);
			mxdiff+=(diff*f);
			mydiff+=(-diff*f);
			mzdiff+=(-diff*f);
		}
		if(t==2) 
		{
			local_en+=((-hx+hy-hz)*diff);
			mxdiff+=(-diff*f);
			mydiff+=(diff*f);
			mzdiff+=(-diff*f);
		}
		if(t==3) 
		{
			local_en+=((-hx-hy+hz)*diff);
			mxdiff+=(-diff*f);
			mydiff+=(-diff*f);
			mzdiff+=(diff*f);
		}
	}
	ediff=-(local_en)*g*mu_b;
}


/////////////////////////////////////////////////////////////////////////////
void single_flip_move(std::vector<std::vector<int> > &neighbors,
                      std::vector<std::vector<int> > &nneighbors,
		      std::vector<std::vector<int> > &ijkt,
		      double J1, double J2, double g, double D,double rnn, 
		      RMatrix &Drmatrix, RMatrix &Dkmatrix,
		      double hx, double hy, double hz,
		      std::vector<double> &configz, std::vector<double> &newconfigz,
		      std::vector<int>    &touched_sites,
		      double &ediff, double &mxdiff,double &mydiff, double &mzdiff, bool &aborted)

{
	aborted=false;
	int nsites=configz.size();

	ediff=0.0;
	mxdiff=0.0;
	mydiff=0.0;
	mzdiff=0.0;

	touched_sites.clear();
	int site=uniform_rand_int(0,nsites);
	touched_sites.push_back(site);  // Single site picked

	#pragma omp parallel for
	for (int i=0;i<nsites;i++) newconfigz[i]=configz[i];

	newconfigz[site]=-configz[site];  // Single spin flipped

	// NN part + NNN part
	double ediffnn=0.0;
	ediff_nn_pyrochlore(J1, J2, neighbors, nneighbors, touched_sites,configz,newconfigz,ediffnn);
	
	// Dipolar part
	double ediffdipolar=0.0;
	if (abs(D)>1.0e-8) ediff_dipolar_pyrochlore(D*rnn*rnn*rnn,touched_sites,configz,newconfigz,Drmatrix,Dkmatrix,ediffdipolar);
	
	// Field part
	double edifffield=0.0;
	ediff_and_mdiff_field_pyrochlore(ijkt, g, hx,hy,hz, touched_sites,configz,newconfigz,edifffield,mxdiff,mydiff,mzdiff);

	ediff=ediffnn+ediffdipolar+edifffield;
}

/////////////////////////////////////////////////////////////////////////////
void loop_move(RMatrix &dmat, std::vector<std::vector<int> > &neighbors,
                      std::vector<std::vector<int> > &nneighbors,
		      std::vector<std::vector<int> > &ijkt,
		      double J1, double J2, double g, double D,double rnn, 
		      RMatrix &Drmatrix, RMatrix &Dkmatrix,
		      double hx, double hy, double hz,
		      std::vector<double> &configz, std::vector<double> &newconfigz,
		      std::vector<int>    &touched_sites,
		      double &ediff, double &mxdiff,double &mydiff, double &mzdiff,bool &aborted)

{
	int nsites=configz.size();

	ediff=0.0;
	mxdiff=0.0;
	mydiff=0.0;
	mzdiff=0.0;

	// Make "short loop" - algorithm adapted from Berkema and Newman
	touched_sites.clear();
	
	int startsite=uniform_rand_int(0,nsites); // Start with a random site
	double sz=configz[startsite];

	aborted=false;
	check_surroundings(neighbors,startsite,configz,touched_sites,aborted);

	std::vector<int> tempsites;
	tempsites.push_back(startsite);
	bool found_loop=false; 
	int closing_point=startsite;

	int algotype=uniform_rand_int(0,2); // two types of sequences - this will give 0 or 1 
					    // same up tet, other tet, same tet, other tet - algorithm 0
					    // OR
					    // other up tet (same down tet), same tet, other tet,... algorithm 1
	//int algotype=0;
	//cout<<"algo type = "<<algotype<<endl;
	// Algorithm condition zero 		    
	if (algotype==0)
	{
		while (found_loop==false and aborted==false)
		{
			int chosen;
			std::vector<int> tmp_nbrs;
			int currentsite1=tempsites[tempsites.size()-1]; // The very last site that was added
			spin_nbrs_on_different_up_tet(neighbors,currentsite1,-sz,configz,tmp_nbrs,tempsites,touched_sites,aborted);
			// Now choose a neighbor at random from this list, if the size of this list > 0 , else abort
			if (aborted==false)
			{
				chosen=tmp_nbrs[int(uniform_rnd()*tmp_nbrs.size() - 1.0e-12)];
				tempsites.push_back(chosen);	
				check_loop_closed_on_adding_last_site(dmat,tempsites,found_loop,closing_point,touched_sites);
				if (found_loop==false and aborted==false)
				{
					// Now add a site on the same tetrahedron
					int currentsite2=tempsites[tempsites.size()-1];
					tmp_nbrs.clear();
					spin_nbrs_on_same_up_tet(neighbors,currentsite2,sz,configz,tmp_nbrs,tempsites,touched_sites,aborted);
					// Look for a neighbor with opposite spin and on a different up tetrahedron		
					if (aborted==false)
					{	
						chosen=tmp_nbrs[int(uniform_rnd()*tmp_nbrs.size() - 1.0e-12)];
						tempsites.push_back(chosen);
						check_loop_closed_on_adding_last_site(dmat,tempsites,found_loop,closing_point,touched_sites);
					}
				}
			}	
		}
	}
	else // Algorithm type 1
	{
		while (found_loop==false and aborted==false)
		{
			int chosen;
			std::vector<int> tmp_nbrs;
			int currentsite1=tempsites[tempsites.size()-1]; // The very last site that was added
			spin_nbrs_on_same_up_tet(neighbors,currentsite1,-sz,configz,tmp_nbrs,tempsites,touched_sites,aborted);
			// Now choose a neighbor at random from this list, if the size of this list > 0 , else abort
			if (aborted==false)
			{
				chosen=tmp_nbrs[int(uniform_rnd()*tmp_nbrs.size() - 1.0e-12)];
				tempsites.push_back(chosen);	
				check_loop_closed_on_adding_last_site(dmat,tempsites,found_loop,closing_point,touched_sites);
				if (found_loop==false and aborted==false)
				{
					// Now add a site on the same tetrahedron
					int currentsite2=tempsites[tempsites.size()-1];
					tmp_nbrs.clear();
					spin_nbrs_on_different_up_tet(neighbors,currentsite2,sz,configz,tmp_nbrs,tempsites,touched_sites,aborted);
					// Look for a neighbor with opposite spin and on a different up tetrahedron		
					if (aborted==false)
					{	
						chosen=tmp_nbrs[int(uniform_rnd()*tmp_nbrs.size() - 1.0e-12)];
						tempsites.push_back(chosen);
						check_loop_closed_on_adding_last_site(dmat,tempsites,found_loop,closing_point,touched_sites);
					}
				}
			}	
		}
	}

	if (aborted) return;
	// If a loop move is not possible do not do anything else, just return
	// The next 4 lines are THUS WRONG!
	/*if (aborted==true)  // This means a loop move was not possible, then do a single spin flip 
	{
		//cout<<"Could not find a loop, so doing single spin flip"<<endl;
		touched_sites.clear();
		touched_sites.push_back(startsite);

	}*/
	
	//if (aborted==true)  // Loop move was not possible because a loop could not be found from the chosen startsite
	//{
	//	return;
	//}
	
	#pragma omp parallel for // Setting New config to old config
	for (int i=0;i<nsites;i++) newconfigz[i]=configz[i];

	#pragma omp parallel for  // Flipping the sites that are part of the loop
	for (int i=0;i<touched_sites.size();i++) 
	{
		int site=touched_sites[i];
		newconfigz[site]=-configz[site];
	}

	// NN + NNN part
	double ediffnn=0.0;
	ediff_nn_pyrochlore(J1, J2, neighbors, nneighbors, touched_sites,configz,newconfigz,ediffnn);
	
	// Dipolar part
	double ediffdipolar=0.0;
	if (abs(D)>1.0e-8) ediff_dipolar_pyrochlore(D*rnn*rnn*rnn,touched_sites,configz,newconfigz,Drmatrix,Dkmatrix,ediffdipolar);
	
	// Field part
	double edifffield=0.0;
	ediff_and_mdiff_field_pyrochlore(ijkt, g,hx,hy,hz,touched_sites,configz,newconfigz,edifffield,mxdiff,mydiff,mzdiff);

	ediff=ediffnn+ediffdipolar+edifffield;
}
////////////////////////////////////////////////////////////////////////
void make_pyrochlore(int L,
		     std::vector< std::vector<double> > &fullcoords,
		     std::vector< std::vector<int> > &ijkt,
		     std::vector< std::vector<int> > &neighbors,
		     std::vector< std::vector<int> > &nneighbors, 
		     RMatrix &dmat, RMatrix &rijxmat, RMatrix &rijymat, RMatrix &rijzmat)
{
	cout<<"Making pyrochlore..."<<endl;
	fullcoords.clear();
	ijkt.clear();
	//////////////////////////////////////////////
	//             Pyrochlore lattice
	//////////////////////////////////////////////
	// Make sites on cube of dimensions L,L,L
	// Associate 16 sites with every n1,n2,n3
	//
	// Tetrahedron coords (used in Lucile/Ross paper)
	// r0 = 1/8 (+1,+1,+1)
	// r1 = 1/8 (+1,-1,-1)
	// r2 = 1/8 (-1,+1,-1)
	// r3 = 1/8 (-1,-1,+1)

	double u=1.0/8.0;
	double r0x=+u;double r0y=+u;double r0z=+u;
	double r1x=+u;double r1y=-u;double r1z=-u;
	double r2x=-u;double r2y=+u;double r2z=-u;
	double r3x=-u;double r3y=-u;double r3z=+u;

        // FCC coords....
	// (0.0,0.0,0.0)
	// (0.5,0.5,0.0)
	// (0.5,0.0,0.5)
	// (0.0,0.5,0.5)
        	
        std::vector< std::vector<double> > fcc_coords;
        std::vector<double> coord;
       
	coord.push_back(0.0);
        coord.push_back(0.0);
        coord.push_back(0.0);
	fcc_coords.push_back(coord);
	coord.clear();
        
	coord.push_back(0.5);
        coord.push_back(0.5);
        coord.push_back(0.0);
	fcc_coords.push_back(coord);
	coord.clear();

	coord.push_back(0.5);
        coord.push_back(0.0);
        coord.push_back(0.5);
	fcc_coords.push_back(coord);
	coord.clear();
	
	coord.push_back(0.0);
        coord.push_back(0.5);
        coord.push_back(0.5);
	fcc_coords.push_back(coord);
	
        int site=0;
	for (int nx=0;nx<L;nx++)
	{
		double x=double(nx);
		for (int ny=0;ny<L;ny++)
		{
			double y=double(ny);
			for (int nz=0;nz<L;nz++)
			{
				double z=double(nz);
				for (int p=0; p<fcc_coords.size();p++)
				{
					double xb=fcc_coords[p][0];	
					double yb=fcc_coords[p][1];	
					double zb=fcc_coords[p][2];
	
					std::vector<int>    ijktentry;
					std::vector<double> fullcoordsentry;
					
					///////////////////////////////////////////////////////////////////////////////////////////////
					ijktentry.push_back(nx);ijktentry.push_back(ny);ijktentry.push_back(nz); ijktentry.push_back(0);	
					fullcoordsentry.push_back(x+xb+r0x);fullcoordsentry.push_back(y+yb+r0y);fullcoordsentry.push_back(z+zb+r0z);
					ijkt.push_back(ijktentry);
					fullcoords.push_back(fullcoordsentry);
					site+=1;
					ijktentry.clear();
					fullcoordsentry.clear();
					///////////////////////////////////////////////////////////////////////////////////////////////
	
					ijktentry.push_back(nx);ijktentry.push_back(ny);ijktentry.push_back(nz); ijktentry.push_back(1);	
					fullcoordsentry.push_back(x+xb+r1x);fullcoordsentry.push_back(y+yb+r1y);fullcoordsentry.push_back(z+zb+r1z);
					ijkt.push_back(ijktentry);
					fullcoords.push_back(fullcoordsentry);
					site+=1;
					ijktentry.clear();
					fullcoordsentry.clear();
					///////////////////////////////////////////////////////////////////////////////////////////////
					
					ijktentry.push_back(nx);ijktentry.push_back(ny);ijktentry.push_back(nz); ijktentry.push_back(2);	
					fullcoordsentry.push_back(x+xb+r2x);fullcoordsentry.push_back(y+yb+r2y);fullcoordsentry.push_back(z+zb+r2z);
					ijkt.push_back(ijktentry);
					fullcoords.push_back(fullcoordsentry);
					site+=1;
					ijktentry.clear();
					fullcoordsentry.clear();
					///////////////////////////////////////////////////////////////////////////////////////////////
					
					ijktentry.push_back(nx);ijktentry.push_back(ny);ijktentry.push_back(nz); ijktentry.push_back(3);	
					fullcoordsentry.push_back(x+xb+r3x);fullcoordsentry.push_back(y+yb+r3y);fullcoordsentry.push_back(z+zb+r3z);
					ijkt.push_back(ijktentry);
					fullcoords.push_back(fullcoordsentry);
					site+=1;
					ijktentry.clear();
					fullcoordsentry.clear();
					///////////////////////////////////////////////////////////////////////////////////////////////
				}
			}
		}
	}
	
	cout<<"Number of sites recorded is = "<<site<<endl;
	int nsites=site;
	dmat.resize(nsites,nsites);
	rijxmat.resize(nsites,nsites);
	rijymat.resize(nsites,nsites);
	rijzmat.resize(nsites,nsites);
	for (int i=0;i<nsites;i++)
	{
		double x1=fullcoords[i][0];
		double y1=fullcoords[i][1];
		double z1=fullcoords[i][2];
		for (int j=0;j<nsites;j++)
		{
			double x2=fullcoords[j][0];
			double y2=fullcoords[j][1];
			double z2=fullcoords[j][2];
			
			double nijx=ijkt[j][0]-ijkt[i][0];	
			double nijy=ijkt[j][1]-ijkt[i][1];	
			double nijz=ijkt[j][2]-ijkt[i][2];
			
                        double rijx=x2-x1;	
			double rijy=y2-y1;	
			double rijz=z2-z1;
	
			if (L==1)
			{
				rijx=min(min(abs((x2-x1)-L),abs((x2-x1)+L)),abs(x2-x1));		
				rijy=min(min(abs((y2-y1)-L),abs((y2-y1)+L)),abs(y2-y1));		
				rijz=min(min(abs((z2-z1)-L),abs((z2-z1)+L)),abs(z2-z1));		
				//if ((x2-x1)>=L) rijx=(x2-x1)-L; // Translate x	
				//if ((x1-x2)>=L) rijx=(x2-x1)+L; // Translate x 
				//if ((y2-y1)>=L) rijy=(y2-y1)-L; // Translate y	
				//if ((y1-y2)>=L) rijy=(y2-y1)+L; // Translate y	
				//if ((z2-z1)>=L) rijz=(z2-z1)-L; // Translate z	
				//if ((z1-z2)>=L) rijz=(z2-z1)+L; // Translate z
			}
			else
			{	
				rijx=min(min(abs((x2-x1)-L),abs((x2-x1)+L)),abs(x2-x1));		
				rijy=min(min(abs((y2-y1)-L),abs((y2-y1)+L)),abs(y2-y1));		
				rijz=min(min(abs((z2-z1)-L),abs((z2-z1)+L)),abs(z2-z1));		
				/*if (abs(nijx)>=L/2.0 and nijx>0) rijx=(x2-x1)-L; // Translate x	
				if (abs(nijx)>=L/2.0 and nijx<0) rijx=(x2-x1)+L; // Translate x 
				if (abs(nijy)>=L/2.0 and nijy>0) rijy=(y2-y1)-L; // Translate y	
				if (abs(nijy)>=L/2.0 and nijy<0) rijy=(y2-y1)+L; // Translate y
				if (abs(nijz)>=L/2.0 and nijz>0) rijz=(z2-z1)-L; // Translate z	
				if (abs(nijz)>=L/2.0 and nijz<0) rijz=(z2-z1)+L; // Translate z*/
			}
			rijxmat(i,j)=x1-x2;
			rijymat(i,j)=y1-y2;
			rijzmat(i,j)=z1-z2;
			dmat(i,j)=sqrt((rijx*rijx)+(rijy*rijy)+(rijz*rijz));  // min distance between points in pbc
			//cout<<"dmat(i,j) = "<<dmat(i,j)<<endl;
			if (dmat(i,j)>0.353 and dmat(i,j)<0.354)  neighbors[i].push_back(j);  //      nearest neighbor on pyrochlore
			if (dmat(i,j)>0.612 and dmat(i,j)<0.613)  nneighbors[i].push_back(j); // next nearest neighbor on pyrochlore
		}
	}

     	//cout<<"Here"<<endl;   
	//if (measure_corrs)
	{	
		int six1=0;	
		int twelve1=0;	
		for (int i=0;i<nsites;i++)
		{
			//cout<<" i = "<<i<<", nsize,nnsize ="<<neighbors[i].size()<<"  "<<nneighbors[i].size()<<endl;
			if (neighbors[i].size()==6)   six1+=1;
			if (nneighbors[i].size()==12) twelve1+=1;
		}

		cout<<"Number of sites with  6      nearest neighbors = "<<six1<<endl;
		cout<<"Number of sites with 12 next nearest neighbors = "<<twelve1<<endl;

		/*cout<<"========================================================================================================================="<<endl;
		cout<<" i       j      dij(min dist between i and j in pbc)                                                                     "<<endl;
		cout<<"========================================================================================================================="<<endl;
		for (int i=0;i<nsites;i++)
		{
			for (int j=0;j<nsites;j++) cout<<boost::format("%3d    %3d    %+ .5f") %i %j %dmat(i,j)<<endl;
		}*/
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void make_Dmatrix(int L, int kcut, double alpha, double muprime,std::vector< std::vector<int> > &ijkt,
		  std::vector< std::vector<double> > &fullcoords,
		  RMatrix &rijxmat, RMatrix &rijymat, RMatrix &rijzmat, RMatrix &Drmatrix, RMatrix &Dkmatrix)
{
	int nsites=16*L*L*L;
	Drmatrix.resize(nsites,nsites);
	Dkmatrix.resize(nsites,nsites);
	/*
	//////////////////////////////////////////////////////////////////////////////
	//   Make B and C matrices - ONE time calculation
	//////////////////////////////////////////////////////////////////////////////
	RMatrix bmat(nsites,nsites),cmat(nsites,nsites);
	# pragma omp parallel for
	for (int i=0;i<nsites;i++)
	{
		for (int j=0;j<nsites;j++)
		{
			if (i!=j)
			{
				double rijx=rijxmat(i,j);
				double rijy=rijymat(i,j);
				double rijz=rijzmat(i,j);
				double r=sqrt(rijx*rijx + rijy*rijy + rijz*rijz);
				double b,c;
				get_b_c(r,alpha,b,c);
				bmat(i,j)=b;
				cmat(i,j)=c;
			}
		}
	}*/
	//////////////////////////////////////////////////////////////////////////////
	//   Make k vectors - ONE time calculation
	//////////////////////////////////////////////////////////////////////////////
	std::vector< std::vector<double> > kvecs;
	std::vector<double> kfactors;
	// Make k-sites on box of dimensions L,L,L
	double twopibL=(2.0*3.14159265359)/double(L);
	for (int n1=-kcut;n1<=kcut;n1++)
	{
		double kx=double(n1)*twopibL;
		for (int n2=-kcut;n2<=kcut;n2++)
		{
			double ky=double(n2)*twopibL;
			for (int n3=-kcut;n3<=kcut;n3++)
			{
				double kz=double(n3)*twopibL;
				//if ((n1==0 and n2==0 and n3==0) or ( (n1*n1) + (n2*n2) + (n3*n3)> double(kcut*kcut)+1.0e-6) )
				if ((n1==0 and n2==0 and n3==0))
				{
					//cout<<"Disregarding k vec = 0,0,0 or k outside sphere"<<endl;
				}
				else
				{	
					std::vector<double> kvec;
					kvec.push_back(kx);
					kvec.push_back(ky);
					kvec.push_back(kz);
					kvecs.push_back(kvec);
					double k2=((kx*kx) + (ky*ky) + (kz*kz));
					double kfactor=exp(-k2/(4.0*alpha*alpha))/k2;
					kfactors.push_back(kfactor);
				}
			}
		}	
	}
        cout<<"Made k vectors " <<endl;	
	int numkvecs=kvecs.size();
	double invrootpi=0.56418958354;
	cout<<"numkvecs = "<<numkvecs<<endl;
	//////////////////////////////////////////////////////////////////////////////
	//   Use full coordinates to construct DMatrix(i,j) - ONE TIME calculation
	//////////////////////////////////////////////////////////////////////////////
	RMatrix vecs(4,3);	
	double f=1.0/sqrt(3.0);

	vecs(0,0)=f;vecs(0,1)=f;vecs(0,2)=f;
	vecs(1,0)=f;vecs(1,1)=-f;vecs(1,2)=-f;
	vecs(2,0)=-f;vecs(2,1)=f;vecs(2,2)=-f;
	vecs(3,0)=-f;vecs(3,1)=-f;vecs(3,2)=f;
	
	# pragma omp parallel for
	for (int i=0;i<nsites*nsites;i++) Drmatrix[i]=0.0;
	# pragma omp parallel for
	for (int i=0;i<nsites*nsites;i++) Dkmatrix[i]=0.0;

	// ONE time computation
	// real space part - does not assume minimum image convention
	cout<<"Real space part of Dij"<<endl;	
	for (int n1=-kcut;n1<=kcut;n1++)
	{
		for (int n2=-kcut;n2<=kcut;n2++)
		{
			for (int n3=-kcut;n3<=kcut;n3++)
			{	
				# pragma omp parallel for 
				for (int i=0;i<nsites;i++)
				{
					int ti=ijkt[i][3];
					for (int j=0;j<nsites;j++)
					{
						int tj=ijkt[j][3];
						double rijx=rijxmat(i,j)+n1*L;
						double rijy=rijymat(i,j)+n2*L;
						double rijz=rijzmat(i,j)+n3*L;
						double r=sqrt(rijx*rijx + rijy*rijy + rijz*rijz);
						double b=0,c=0;
						if (abs(r)>1.0e-8) get_b_c(r,alpha,b,c); // avoids the singularity
						double f2=1.0;
						if (ti!=tj) f2=-1.0/3.0;
						double f3=vecs(ti,0)*rijx+vecs(ti,1)*rijy+vecs(ti,2)*rijz;
						double f4=vecs(tj,0)*rijx+vecs(tj,1)*rijy+vecs(tj,2)*rijz;
						if (n1==0 and n2==0 and n3==0) // when doing central cell avoid i=j
						{
							if (i!=j) Drmatrix(i,j)+= (f2*b)-(f3*f4*c); 
						}
						else
						{
							Drmatrix(i,j)+= (f2*b)-(f3*f4*c); 
						}
					}
				 }
			 }
		}
	}
	
	cout<<"K space part of Dij"<<endl;	
	// ONE time computation
	// K space part, K= 0 part removed	
	double twopibv=2.*3.1415926539/double(L*L*L);
	# pragma omp parallel for 
	for (int i=0;i<nsites;i++)
	{
		int ti=ijkt[i][3];
		for (int j=0;j<nsites;j++)
		{
			int tj=ijkt[j][3];
			double rijx=rijxmat(i,j); // which unit cell doesnt matter because of the cosine function
			double rijy=rijymat(i,j);
			double rijz=rijzmat(i,j);
			for (int kint=0;kint<numkvecs;kint++) // k = 0 not part of this list!
			{
				double kx=kvecs[kint][0];
				double ky=kvecs[kint][1];
				double kz=kvecs[kint][2];
				double factori=kx*vecs(ti,0)+ky*vecs(ti,1)+kz*vecs(ti,2);
				double factorj=kx*vecs(tj,0)+ky*vecs(tj,1)+kz*vecs(tj,2);
				double kfactor=kfactors[kint];
				double kdotr=(kx*rijx + ky*rijy + kz*rijz);
				Dkmatrix(i,j)+=(kfactor*cos(kdotr)*factori*factorj*twopibv);  
			}
		}
	}
	// Surface term in Dkmatrix
	# pragma omp parallel for 
	for (int i=0;i<nsites;i++)
	{
		int ti=ijkt[i][3];
		for (int j=0;j<nsites;j++)
		{
			int tj=ijkt[j][3];
			double factor=1.0;
			if (ti!=tj) factor=-1.0/3.0;	
			Dkmatrix(i,j)+=(factor*twopibv*(1.0/(2.0*muprime+1.0)));  
		}
	}
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void local_jh_field(double &spin, int &site, int &t, double &hx, double &hy, double &hz,
		      std::vector<double> &configz, 
		      std::vector< std::vector<int> > &neighbors,
		      std::vector< std::vector<int> > &nneighbors,
		      double J1, double J2,
		      double g,
		      double &eff_field_z)
{
	//cout<<"Calculating local j energy"<<endl;
	//local_energy=0.0;
        double mu_b=5.7883818012*0.01; // meV/Tesla
	eff_field_z=0.0;
	// Disorder free nearest neighbor terms
	for (int j=0;j<neighbors[site].size();j++)
	{
		int k=neighbors[site][j];
		double szn=configz[k];
		eff_field_z+=((J1/3.0)*szn*spin);
	}
	// Disorder free next nearest neighbor terms
	for (int j=0;j<nneighbors[site].size();j++)
	{
		int k=nneighbors[site][j];
		double szn=configz[k];
		eff_field_z+=((J2/3.0)*szn*spin);
	}
	/*
	// Bond disorder nearest neighbor Heisenberg
	for (int j=0;j<neighbors[site].size();j++)
	{
		int k=neighbors[site][j];
		double J=bond_disorder_matrix(site,k);
		eff_field_x+=(J*configx[k]);
		eff_field_y+=(J*configy[k]);
		eff_field_z+=(J*configz[k]);
	}

	if (abs(Jnnn)>1.0e-6)
	{	
		// Next nn Heisenberg
		for (int j=0;j<nneighbors[site].size();j++)
		{
			int k=nneighbors[site][j];
			eff_field_z+=(Jnnn*configz[k]);
		}
	}*/

	g=g/sqrt(3.0);
	// Now accounting for external magnetic field
	if(t==0) eff_field_z=(eff_field_z)-(mu_b*g*(hx+hy+hz));
	if(t==1) eff_field_z=(eff_field_z)-(mu_b*g*(hx-hy-hz));
	if(t==2) eff_field_z=(eff_field_z)-(mu_b*g*(-hx+hy-hz));
	if(t==3) eff_field_z=(eff_field_z)-(mu_b*g*(-hx-hy+hz));
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double local_h_energy(double &spin, int &site, int &t, double &sz,
		      double &hx, double &hy, double &hz, 
		      double g)
{
	g=g/sqrt(3.0);
        double mu_b=5.7883818012*0.01; // meV/Tesla
	if(t==0) return -(mu_b*g*(hx+hy+hz)*sz*spin);
	if(t==1) return -(mu_b*g*(hx-hy-hz)*sz*spin);
	if(t==2) return -(mu_b*g*(-hx+hy-hz)*sz*spin);
	if(t==3) return -(mu_b*g*(-hx-hy+hz)*sz*spin);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double total_h_energy(double &spin, 
		      std::vector<double>  &configz,
		      double &hx, double &hy, double &hz, 
		      double g,
		      std::vector< std::vector<int> > & ijkt)
{
	double e=0;
	g=g/sqrt(3.0);
        double mu_b=5.7883818012*0.01; // meV/Tesla
	int nsites=configz.size();
	for (int n=0; n < nsites; n++)
	{
		int t=ijkt[n][3];
		double sz=configz[n]; 
		if(t==0) e+=((hx+hy+hz)*sz);
		if(t==1) e+=((hx-hy-hz)*sz);
		if(t==2) e+=((-hx+hy-hz)*sz);
		if(t==3) e+=((-hx-hy+hz)*sz);
	}
	return -e*mu_b*g*spin;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double total_j_energy(double &spin, 
		      std::vector<double> &configz, 
		      std::vector< std::vector<int> > &neighbors, 
		      std::vector< std::vector<int> > &nneighbors, 
		      double J1, double J2)
{
	double e=0;
	// Nearest neighbor disorder-free contributions
	for (int i=0;i<neighbors.size();i++)
	{
		double szi=configz[i];
		for (int j=0;j<neighbors[i].size();j++)
		{
			int k=neighbors[i][j];
			double szk=configz[k];
			e=e+((J1/3.0)*szi*szk);
		}
	}

	// Next Nearest neighbor disorder-free contributions
	for (int i=0;i<nneighbors.size();i++)
	{
		double szi=configz[i];
		for (int j=0;j<nneighbors[i].size();j++)
		{
			int k=nneighbors[i][j];
			double szk=configz[k];
			e=e+((J2/3.0)*szi*szk);
		}
	}
	/*
	// Nearest neighbor bond disorder Heisenberg contributions
	for (int i=0;i<neighbors.size();i++)
	{
		double sx1=configx[i];double sy1=configy[i];double sz1=configz[i];
		for (int j=0;j<neighbors[i].size();j++)
		{
			int k=neighbors[i][j];
			double sx2=configx[k];double sy2=configy[k];double sz2=configz[k];
			e=e+(bond_disorder_matrix(i,k)*((sx1*sx2) + (sy1*sy2) + (sz1*sz2)));
		}
	}*/
	// Next nn Heisenberg 
	/*for (int i=0;i<nneighbors.size();i++)
	{
		double sx1=configx[i];double sy1=configy[i];double sz1=configz[i];
		for (int j=0;j<nneighbors[i].size();j++)
		{
			int k=nneighbors[i][j];
			double sx2=configx[k];double sy2=configy[k];double sz2=configz[k];
			e=e+(Jnnn*((sx1*sx2) + (sy1*sy2) + (sz1*sz2)));
		}
	}*/
	return (e/2.0)*spin*spin;  // Pairs counted twice so divide by 2
}
/////////////////////////////////////////////////////////////////////////////////
std::vector<double> local_magnetization(double 				&spin,
					int 				&t,
					double 				&sz)
{
	double  f=1.0/sqrt(3.0);
	double  mx=0;
	double  my=0;
	double  mz=0;
	std::vector<double> m;
	if (t==0) {mx=f;my=f;mz=f;}
	if (t==1) {mx=f;my=-f;mz=-f;}
	if (t==2) {mx=-f;my=f;mz=-f;}
	if (t==3) {mx=-f;my=-f;mz=f;}
	
	m.push_back(mx*sz*spin);
	m.push_back(my*sz*spin);
	m.push_back(mz*sz*spin);
	return m;
}
/////////////////////////////////////////////////////////////////////////////////
std::vector<double> total_magnetization(double 				&spin,
					std::vector<double> 		&configz,
					std::vector< std::vector<int> > &ijkt)
{
	double  f=1.0/sqrt(3.0);
	double mx=0;
	double my=0;
	double mz=0;
	std::vector<double> m;
	int nsites=int(configz.size());
	for (int n=0;n<nsites;n++)
	{
		int t=ijkt[n][3];
		double sz=configz[n];
		if (t==0) {mx+=(f*sz);my+=(f*sz);mz+=(f*sz);}
		if (t==1) {mx+=(f*sz);my+=(-f*sz);mz+=(-f*sz);}
		if (t==2) {mx+=(-f*sz);my+=(f*sz);mz+=(-f*sz);}
		if (t==3) {mx+=(-f*sz);my+=(-f*sz);mz+=(f*sz);}
	}
	m.push_back(mx*spin);
	m.push_back(my*spin);
	m.push_back(mz*spin);
	return m;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void mc_pyrochlore_pt(double spin, int L, int64_t nsamples, int64_t nburn, double temp, 
		   int ntemps, double hx, double hy, double hz, 
		   double J1, double J2, double D, double g, double alphaL, int kcut, double loop_prob, double muprime)
{
	/////////////////////////////////////////////////////////////////////////
	// MC related quantities
	// Units - Assume J in meV (millielectron volts)
        //         Assume h in T   (tesla)
        //         Convert h to meV units in total_h_energy and local_h_energy
        // Convert temperature in Kelvin to temperature in meV
	// Set 16xLxLxL pyrochlore lattice
	double alpha=alphaL/double(L);
	int nsites=16*L*L*L;
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
	
	/////////////////////////////////////////////////////////////////////////
	// Set 16xLxLxL pyrochlore lattice
	std::vector< std::vector<int> > neighbors(nsites);	
	std::vector< std::vector<int> > bondtype(nsites);	
	std::vector< std::vector<int> > nneighbors(nsites);	
	std::vector< std::vector<double> > fullcoords;	
	std::vector< std::vector<int> > ijkt;

	RMatrix dmat,bmat,cmat,Drmatrix,Dkmatrix,rijxmat,rijymat,rijzmat;
        cout<<"Making pyrochlore.." <<endl;	
        make_pyrochlore(L,fullcoords,ijkt,neighbors,nneighbors,dmat,rijxmat,rijymat,rijzmat);
	if (abs(D)>1.0e-8)
	{
        	cout<<"Making DMatrix.." <<endl;	
		make_Dmatrix(L,kcut,alpha,muprime, ijkt,fullcoords,rijxmat,rijymat,rijzmat,Drmatrix,Dkmatrix);
        	cout<<"DMatrix done..." <<endl;	
	}
	/////////////////////////////////////////////////////////////////////////////
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
		infos[b].configz.resize(nsites);
		for (int n=0;n<nsites;n++)
		{
			double rand=uniform_rnd();
			if (rand<0.5)  infos[b].configz[n]=-1.0;
			if (rand>=0.5) infos[b].configz[n]=+1.0;
		}
	}

	double rnn=sqrt(1.0/8.0);
	cout<<"Calculating energies..."<<endl;	
	for (int b=0;b<ntemps;b++)
	{
		double 		    energyj=total_j_energy(spin, infos[b].configz, neighbors, nneighbors,J1, J2);
		double 		    energyh=total_h_energy(spin, infos[b].configz, hx,hy,hz, g, ijkt);
		double 		    energydipolar=0.0;
		double              udr=0.0,udk=0.0,usurf=0.0,uself=0.0;
		if (abs(D)>1.0e-8)
		{
			total_dipolar_energy_pyrochlore(D*rnn*rnn*rnn,alpha,infos[b].configz,Drmatrix,Dkmatrix,energydipolar);
		}
		std::vector<double> magnetization=total_magnetization(spin, infos[b].configz, ijkt);
		infos[b].energy=energyj+energyh+energydipolar;
		infos[b].mx=magnetization[0];infos[b].my=magnetization[1];infos[b].mz=magnetization[2];
		cout<<"Total J       energy                            ="<<energyj<<endl;
		cout<<"Total h       energy                            ="<<energyh<<endl;
		cout<<"Total dipolar energy                            ="<<energydipolar<<endl;
		cout<<"Total J energy + total h + dipolar energy ="<<infos[b].energy<<endl;
		cout<<"Total magnetization is                    ="<<endl; print_vec_acc(magnetization,true);
	}

	cout<<endl;	
	cout<<endl;	
	cout<<"========================================================"<<endl;
	cout<<"Energy history "<<endl;	
	/////////////////////////////////////////////////////////////////////////
	// Accept reject Metropolis
	double ssftries=0.0;
	double looptries=0.0;
	double ssfsuccess=0.0;
	double loopsuccess=0.0;
	double accept=0.0;
	double reject=0.0;
	double nsuccess=0.0;
	//int64_t n=-1;
	while (nsuccess<double(nsamples+nburn))
	{
		//n=n+1;
		int b=0;

		// Perform move
		std::vector<int> touched_sites;
		double ediff,mxdiff,mydiff,mzdiff; 
		double move_rand=uniform_rnd();
		std::vector<double> newconfigz(nsites);
		bool aborted=false;
		bool ssf=false;
		bool loopmove=false;
		if (move_rand<1.0-loop_prob) 
		{
			ssftries+=1.0;
			ssf=true;
			single_flip_move(neighbors,nneighbors,ijkt,J1,J2,g,D,rnn,Drmatrix,Dkmatrix,hx,hy,hz,infos[b].configz,newconfigz,touched_sites,ediff,mxdiff,mydiff,mzdiff,aborted);
		}
		else               
		{
			looptries+=1.0;
			loopmove=true;
			loop_move(dmat,neighbors,nneighbors,ijkt,J1,J2,g,D,rnn,Drmatrix,Dkmatrix,hx,hy,hz,infos[b].configz,newconfigz,touched_sites,ediff,mxdiff,mydiff,mzdiff,aborted);
		}

		if (aborted==false)
		{	
			nsuccess+=1.0;
			if (ssf) ssfsuccess+=1.0;
			if (loopmove) loopsuccess+=1.0;
			// Decide temperature
			double beta; 
			//if (n<nburn/2) beta=0.0 + infos[b].beta*double(n)/double(nburn/2); // Take temperature to infinity and then slowly reduce
			//if (n>=nburn/2) beta=infos[b].beta;
			beta=infos[b].beta;
			
			// Metropolis
			double prob=exp(-beta*ediff);
			double rand=uniform_rnd();
			if (rand<prob) // Metropolis Accept-reject for a given temperature
			{
				if (b==0 and int64_t(nsuccess+1.0e-6)>(nburn)) accept=accept+1.0;
				// Reset configz to new configz, because accepted
				#pragma omp parallel for
				for (int site =0;site<nsites;site++)
				{
					infos[b].configz[site]=newconfigz[site];
				}
				infos[b].energy+=ediff;
				//if (infos[b].energy<elowest) elowest=infos[b].energy;
				infos[b].mx+=mxdiff;infos[b].my+=mydiff;infos[b].mz+=mzdiff;
				//if (n>nburn and n%nsites==0) infos[0].update_totals(); // Update averages - after 1 sweep and after equilibration done 
			}
			else // if rejected - keep earlier config, but update total
			{
				if (b==0 and int64_t(nsuccess+1.0e-6)>(nburn)) reject=reject+1.0;
				//if (n>nburn and n%nsites==0) infos[0].update_totals(); // Update averages - after 1 sweep and after equilibration done 
			}
			if (int64_t(nsuccess+1.0e-6)>nburn and int64_t(nsuccess+1.0e-6)%nsites==0) infos[0].update_totals(); // Update averages - after 1 sweep and after equilibration done 
			if (b==0 and int(nsuccess+1.0e-6)%nsites==0) cout<<infos[0].energy<<","<<infos[0].mx<<","<<infos[0].my<<","<<infos[0].mz<<endl;
		}
	}	
	
	cout<<"========================================================"<<endl;
	cout<<endl;
	accept=accept/(accept+reject);
	// Only the lowest temperature is relevant for averages we are interested in 
	cout<<"Nmeas            = "<<boost::format("%+ .15f") %infos[0].nmeas<<endl;
	cout<<"SSF tries        = "<<boost::format("%+ .15f") %ssftries<<endl;
	cout<<"SSF success      = "<<boost::format("%+ .15f") %ssfsuccess<<endl;
	cout<<"Loop tries       = "<<boost::format("%+ .15f") %looptries<<endl;
	cout<<"Loop success     = "<<boost::format("%+ .15f") %loopsuccess<<endl;
	cout<<"accept           = "<<boost::format("%+ .15f") %accept<<endl;
	cout<<"Mx (last)        = "<<boost::format("%+ .15f") %infos[0].mx<<endl;
	cout<<"My (last)        = "<<boost::format("%+ .15f") %infos[0].my<<endl;
	cout<<"Mz (last)        = "<<boost::format("%+ .15f") %infos[0].mz<<endl;
	cout<<endl;
	cout<<"Taux(last)       = "<<boost::format("%+ .15f") %(infos[0].my*hz - infos[0].mz*hy)<<endl;
	cout<<"Tauy (last)      = "<<boost::format("%+ .15f") %(infos[0].mz*hx - infos[0].mx*hz)<<endl;
	cout<<"Tauz (last)      = "<<boost::format("%+ .15f") %(infos[0].mx*hy - infos[0].my*hx)<<endl;

	cout<<endl;
	infos[0].average();
	cout<<"Mx (avg)         = "<<boost::format("%+ .15f") %infos[0].mxavg<<endl;
	cout<<"My (avg)         = "<<boost::format("%+ .15f") %infos[0].myavg<<endl;
	cout<<"Mz (avg)         = "<<boost::format("%+ .15f") %infos[0].mzavg<<endl;

	cout<<endl;
	cout<<"Taux (avg)       = "<<boost::format("%+ .15f") %(infos[0].myavg*hz - infos[0].mzavg*hy)<<endl;
	cout<<"Tauy (avg)       = "<<boost::format("%+ .15f") %(infos[0].mzavg*hx - infos[0].mxavg*hz)<<endl;
	cout<<"Tauz (avg)       = "<<boost::format("%+ .15f") %(infos[0].mxavg*hy - infos[0].myavg*hx)<<endl;
        
	cout<<endl;
	cout<<"E   (avg)        = "<<boost::format("%+ .15f") %infos[0].eavg<<endl;
	cout<<"E^2 (avg)        = "<<boost::format("%+ .15f") %infos[0].e2avg<<endl;
	cout<<"E^4 (avg)        = "<<boost::format("%+ .15f") %infos[0].e4avg<<endl;
	cout<<"C   (avg)/site   = "<<boost::format("%+ .15f") %infos[0].spheatpersite<<endl;
	cout<<endl;
	cout<<"x                  y                    z             sz tilde"<<endl;
	for (int i=0;i<nsites;i++)
	{
		double x=fullcoords[i][0]; double y=fullcoords[i][1]; double z=fullcoords[i][2];
		double sz=infos[0].configz[i];
		cout<<boost::format("%+ .10f  %+ .10f  %+ .10f  %+ .10f") %x %y %z %sz <<endl;
	}
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

