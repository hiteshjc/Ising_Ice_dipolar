#include"global.h"
#include"oleg.h"
#include"hamiltonian_spin_functions.h"
#include"number_functions.h"
#include"search_for.h"
#include"printing_functions.h"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////////////
//
//                                 Oleg's MODEL
//
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Spin_Half_Oleg::operator()
                    (int spin_det,
                     std::vector<int> &new_spin_dets, 
                     std::vector< complex<double> > &hints_list)
{
    std::vector< std::vector<int> >:: iterator p;
    
    for (p=this->hexagons.begin();p!=this->hexagons.end();p++)
    {
	int first=(*p)[0]; int second=(*p)[1];
	int third=(*p)[2]; int fourth=(*p)[3];
	int fifth=(*p)[4]; int sixth=(*p)[5];
	calc_hints_xyzxyz(0.0, first, second, third, fourth, fifth, sixth, 
			     spin_det, new_spin_dets, hints_list);

	// Note H = - Sum W
    }
    
    for (p=this->up_triangles.begin();p!=this->up_triangles.end();p++)
    {
	int first=(*p)[0]; int second=(*p)[1];int third=(*p)[2];
	calc_hints_xyz(-1.0, first, second, third, spin_det, new_spin_dets, hints_list);
	// Note H = - Sum W
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void oleg_setup(string filename, 
               Spin_Half_Oleg &oleg)

{    
    bool found=true;
    string str,str_ret;
    std::vector< std::vector<int> > hexagons, up_triangles;
    str="hexagons";
    search_for(str,filename,str_ret,found);
    if (found)
    {
          if (str_ret.substr(0,1)==string("[")) hexagons=convert_string_to_vec_of_vec(str_ret);
    }    
    str="up_triangles";
    search_for(str,filename,str_ret,found);
    if (found)
    {
          if (str_ret.substr(0,1)==string("[")) up_triangles=convert_string_to_vec_of_vec(str_ret);
    }    
    oleg.init(hexagons,up_triangles);
}


