#include<iomanip>
#include<vector>
#include<cmath>
#include<ctime>
#include<cstdio>
#include<cstdlib>
#include<cassert>
#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include <time.h>
#include<deque>
#include<vector>
#include <complex>
#include <valarray>
#include "MMSP.hpp"


using namespace std;

typedef float phi_type; 											
typedef MMSP::sparse<phi_type> store_type; 	

#include "Array.h"
#include "/home/arun/Documents/mmsp-arun/fftw++-2.05/fftw++.h"

using namespace utils;
using namespace Array;
using namespace fftwpp;

#include "global.h"
#include "grid_ind_funcs.h"
#include "grid_dep_vars.hpp"


int main()
{	
    initialize_epsi();
	initialize();
	initialize_alpha_eigen() ;
	define_c_sigma() ;


	//-------Setting the boundary conditions for various grids-------//
	MMSP::b0(grid,0) = MMSP::Dirichlet;
	MMSP::b1(grid,0) = MMSP::Dirichlet;

	MMSP::b0(grid,1) = MMSP::Dirichlet;
	MMSP::b1(grid,1) = MMSP::Dirichlet;

	MMSP::b0(grid,2) = MMSP::Dirichlet;
	MMSP::b1(grid,2) = MMSP::Dirichlet;

	
    srand(time(NULL));
    
    init_fourier();

    for(int t = 0 ; t<steps ; t++)
	{
	
        cout<<t<<endl;

		if(cr > 0.0)    //Implementing the variation of Temperature based on the given cooling rate
		{ 
            update_T_dependent_vars(T)   ;
        }
    
	
		MMSP::b0(update,0) = MMSP::Dirichlet;
		MMSP::b1(update,0) = MMSP::Dirichlet;

		MMSP::b0(update,1) = MMSP::Dirichlet;
		MMSP::b1(update,1) = MMSP::Dirichlet;

		MMSP::b0(update,2) = MMSP::Dirichlet;
		MMSP::b1(update,2) = MMSP::Dirichlet;

        calculate_strains();  
        calc_gradsq_c();
    
        if(t%500==0)   //Generating output files. Refer to outputgen.hpp
		{
            write_output(t);
        }
        
        if(t%5000==0)   //Generating output files. Refer to outputgen.hpp
		{
            string file_name = "variants_dist_" + to_string(t) + ".txt" ;
            ofstream myfile;
            myfile.open(file_name.c_str());
            for(int i = 0 ; i < variant_ids.size() ; i++)
            {
                myfile<<variant_ids[i]<<endl;
            }
            myfile.close();
        }
        
        if(t%500==0  & variant_ids.size() <= 45)
        {
            nucleation(); 
        }
        
		update_node_n();

		swap(grid, update);	
    }

return 0 ;

}







