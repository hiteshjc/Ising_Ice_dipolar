#include"global.h"
#include"mtrand.h"

// Search
#include"search_for.h"

// Hamiltonian
#include"hamiltonian.h"
#include"oleg.h"
// MC
#include"mc.h"
#include"mc_finite_D.h"
#include"mc_pyrochlore.h"
#include"number_functions.h"
using namespace std;


int main(int argc, char *argv[])
{
    time_t start,end;
    double dif;
    bool found;
    int seed;
    MTRand irand;
    string str_ret,neigs_str,nkrylov_str,filename;
    if (argc <= 1)
    {
        cout << "Usage: " << argv[0] << " <Filename>" << endl;
        exit(1);
    }
 
    filename=argv[1];
    search_for(string("seed"),filename,str_ret,found);
    if (found) 
    {
		if (found) {seed=str_to_int(str_ret);}
		else {seed=1;}
    }	
    irand.seed(seed);
    double e,mx,my,mz,e2,mx2,my2,mz2;
    bool Lfound,nstatesfound,measurecorrsfound,disorderfound, nreplicasfound; 
    bool spinfound, tempfound,Dfound, Dpfound, mcmovefound,configfound, nburnfound, nsamplesfound, latticefound, Jnnnfound; 
    bool hxfound,hyfound,hzfound, J1found, J2found, J3found, J4found, Jvalfound, gfound, gxyfound, gzfound, alphaLfound, kcutfound;
    int L,nstates=2;
    double spin, Jval, J1, J2, J3, J4, Jnnn, g, gxy, gz, hx,hy,hz,temp,D,Dp,alphaL;
    string mcmove,start_config,lattice;
    int nreplicas=1; 
    double nsamplesd, nburnd;
    int kcut;
    bool measure_corrs=true;
    double disorder_strength=0.0;
    int Jchoice;
    bool Jchoicefound;

    search_for(string("spin"),filename,str_ret,spinfound);
    if (spinfound){spin=str_to_d(str_ret);} else{spin=0.5;}
    
     search_for(string("L"),filename,str_ret,Lfound);
    if (Lfound){L=str_to_int(str_ret);} else{L=4;}
    
    search_for(string("nsamples"),filename,str_ret,nsamplesfound);
    if (nsamplesfound){nsamplesd=str_to_d(str_ret);} else{nsamplesd=1000000000;}
    
    search_for(string("nburn"),filename,str_ret,nburnfound);
    if (nburnfound){nburnd=str_to_d(str_ret);} else{nburnd=100000000;}

    int64_t nsamples=int64_t(nsamplesd);
    int64_t nburn=int64_t(nburnd);
 
    search_for(string("J1"),filename,str_ret,J1found);
    if (J1found){J1=str_to_d(str_ret);} else{J1=0;}
    
    search_for(string("J2"),filename,str_ret,J2found);
    if (J2found){J2=str_to_d(str_ret);} else{J2=0;}
    
    search_for(string("D"),filename,str_ret,Dfound);
    if (Dfound){D=str_to_d(str_ret);} else{D=0.0;}
    
    search_for(string("alphaL"),filename,str_ret,alphaLfound);
    if (alphaLfound){alphaL=str_to_d(str_ret);} else{alphaL=1;}
    
    search_for(string("kcut"),filename,str_ret,kcutfound);
    if (kcutfound){kcut=str_to_int(str_ret);} else{kcut=10;}
    
    search_for(string("g"),filename,str_ret,gfound);
    if (gfound){g=str_to_d(str_ret);} else{g=0;}
    
    search_for(string("disorder"),filename,str_ret,disorderfound);
    if (disorderfound){disorder_strength=str_to_d(str_ret);} else{disorder_strength=0.0;}
    
    search_for(string("hx"),filename,str_ret,hxfound);
    if (hxfound){hx=str_to_d(str_ret);} else{hx=0;}
    
    search_for(string("hy"),filename,str_ret,hyfound);
    if (hyfound){hy=str_to_d(str_ret);} else{hy=0;}
    
    search_for(string("hz"),filename,str_ret,hzfound);
    if (hzfound){hz=str_to_d(str_ret);} else{hz=0;}
    
    search_for(string("temp"),filename,str_ret,tempfound);
    if (tempfound){temp=str_to_d(str_ret);} else{temp=1.0;}
    
    bool loop_prob_found;
    double loop_prob;
    search_for(string("loop_prob"),filename,str_ret,loop_prob_found);
    if (loop_prob_found){loop_prob=str_to_d(str_ret);} else{loop_prob=0.5;}
    
    bool muprime_found;
    double muprime;
    search_for(string("muprime"),filename,str_ret,muprime_found);
    if (muprime_found){muprime=str_to_d(str_ret);} else{muprime=100000.0;} // Essentially infinite,surface term does not contribute
    									   // muprime=1 is for vacuum
    
    search_for(string("start_config"),filename,str_ret,configfound);
    if (configfound){start_config=str_ret;} else{start_config="random";}
    
    search_for(string("lattice"),filename,str_ret,latticefound);
    if (latticefound){lattice=str_ret;} else{lattice="pyrochlore16";}
    
    search_for(string("measure_corrs"),filename,str_ret,measurecorrsfound);
    if (measurecorrsfound){measure_corrs=str_to_bool(str_ret);} else{measure_corrs=true;}

    cout<<"spin                 = "<<spin<<endl;    
    cout<<"L                    = "<<L<<endl;    
    cout<<"temp (K)             = "<<temp<<endl;    
    cout<<"g                    = "<<g<<endl;    
    cout<<"J1   (meV)           = "<<J1<<endl;    
    cout<<"J2   (meV)           = "<<J2<<endl;    
    cout<<"D   (meV)            = "<<D<<endl;    
    cout<<"hx   (T)             = "<<hx<<endl;    
    cout<<"hy   (T)             = "<<hy<<endl;    
    cout<<"hz   (T)             = "<<hz<<endl;    
    cout<<"nreplicas            = "<<nreplicas<<endl;    
    cout<<"alphaL               = "<<alphaL<<endl;    
    cout<<"kcut                 = "<<kcut<<endl;    
    cout<<"loop_prob            = "<<loop_prob<<endl;    
    cout<<"muprime   		= "<<muprime<<endl;    
    int numprocs=1;
    for (int i=0;i<numprocs;i++)
    {
    	double e,mx,my,mz,e2,mx2,my2,mz2;
    	mc_pyrochlore_pt(spin,L,nsamples,nburn,temp,nreplicas,hx,hy,hz,J1,J2,D,g,alphaL,kcut,loop_prob,muprime);
    	//mc_pyrochlore_pt(lattice, spin,L,Jchoice,nsamples,nburn,start_config,mcmove,temp,nreplicas, hx,hy,hz,J1,J2,J3,J4,Jnnn,disorder_strength,gxy,gz,e,mx,my,mz,e2,mx2,my2,mz2,measure_corrs); //iterative_pyrochlore(spin,L,nsamples,nburn,start_config,hx,hy,hz,J1,J2,J3,J4,Jnnn,disorder_strength,gxy,gz,e,mx,my,mz,e2,mx2,my2,mz2,measure_corrs);
    	//mc_pyrochlore_ground_state(spin,L,nburn,start_config,mcmove,temp,hx,hy,hz,J1,J2,J3,J4,Jnnn,gxy,gz,e,mx,my,mz,e2,mx2,my2,mz2,measure_corrs);
	/*mc_sphere(lattice,spin, R, nsamples, nburn, start_config, 
		   mcmove, temp, ntemps, hx, hy, hz, 
		   J, Dip, D,
		   e, mx, my, mz, e2, mx2, my2, mz2, measure_corrs);*/
    }
}

/////////////////////////////////////////////////////////////////////////////
