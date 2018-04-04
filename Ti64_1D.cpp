// Symmetric Toth XMPF written with MMSP classes and no optimizations

#ifndef GRAINGROWTH_UPDATE
#define GRAINGROWTH_UPDATE

#ifdef BGQ
#include<mpi.h>
#endif
#include<iomanip>
#include<vector>
#include<cmath>
#include<ctime>
#include<math.h>
#include<cstdio>
#include<cstdlib>
#include<cassert>
#include "/home/arun/mmsp/include/MMSP.hpp"
#include "Tensor.hpp"
#include "Eigenstrains.hpp"

typedef float phi_type; 											
typedef MMSP::sparse<phi_type> store_type; 									

template <int dim> void makeSplit(MMSP::grid<dim, store_type >& grid);   					
														
#include"Ti64_1D_dec.hpp" 
#include"Ti64_1D_init.hpp"

														
std::vector<std::string> split(const std::string &text, char sep); 					        
void print_progress(const int step, const int steps, const int iterations);

														

#include "tessellate.hpp"

namespace MMSP {

void print_progress(const int step, const int steps, const int iterations) {					
	char* timestring;
	static unsigned long tstart;
	struct tm* timeinfo;

	if (step==0) {
		tstart = time(NULL);
		std::time_t rawtime;
		std::time( &rawtime );
		timeinfo = std::localtime( &rawtime );
		timestring = std::asctime(timeinfo);
		timestring[std::strlen(timestring)-1] = '\0';
		std::cout<<"Pass "<<std::setw(3)<<std::right<<iterations<<": "<<timestring<<" ["<<std::flush;
	} else if (step==steps-1) {
		unsigned long deltat = time(NULL)-tstart;
		std::cout<<"•] "<<std::setw(2)<<std::right<<deltat/3600<<"h:"
										<<std::setw(2)<<std::right<<(deltat%3600)/60<<"m:"
										<<std::setw(2)<<std::right<<deltat%60<<"s"
										<<" (File "<<std::setw(5)<<std::right<<iterations*steps<<")."<<std::endl;
	} else if ((20 * step) % steps == 0) std::cout<<"• "<<std::flush;					 															 
}

std::vector<std::string> split(const std::string &text, char sep) {						
														
  std::vector<std::string> tokens;
  std::size_t start = 0, end = 0;
  while ((end = text.find(sep, start)) != std::string::npos) {
    tokens.push_back(text.substr(start, end - start));								
    start = end + 1;
  }
  tokens.push_back(text.substr(start));		
  return tokens;
}


														//template <int dim> void update(MMSP::grid<dim, store_type >& grid, int steps) ;
template <int dim> void update(MMSP::grid< dim , store_type>& grid, int steps) 
{    														// computation of the model 
	float dx = MMSP::dx(grid, 0);   									// functor
	int id=0;
	#ifdef MPI_VERSION
 	id=MPI::COMM_WORLD.Get_rank();     
	#endif
	
														//thresh is the minimum phi value imposed to minimize mathematical noise in fields
	phi_type thresh = 1.0e-3;
	
	phi_type width = 1.0;											//width of stable interface
	phi_type M = 1.0e-3;											//mobility of field motion
	phi_type dt = dx / (20.0 * M);									 	//timestep calculated as function of resolution, 10x less than Courant maximum
	phi_type gamma = 1.0/3.0;										//energy contribution from gradient component
	phi_type w = 3.0 * (gamma/width);									//energy contribution from double-well component
	phi_type eps_sq = 3.0 * width * gamma;
	phi_type Mc = 1.0 ;
	phi_type fsA = 1.0 ;
	phi_type fsB = 1.0 ;
	phi_type WA = 1.0 ;
	phi_type WB = 1.0 ;
	phi_type R = 1.0 ;
	phi_type T = 1.0 ;
	
														// Spatially-varying simulation parameters calculated
														// Initialize vairiables
	const int num_grains =  fields(grid) - 2 ; 					// Subtracting the composition fields.				// functor
	
	static int iterations = 1;
	
	for (int step = 0; step < steps; step++) {
		if(id == 0)	print_progress(step, steps, iterations);
		ghostswap(grid);    										
		
		MMSP::grid<dim, store_type > update(grid);							

		for (int n = 0; n < nodes(grid); n++){
			vector<int> x = position(grid, n);  										
			
			
			//Updating the composition fields
			
			sparse<phi_type> c;
		    for (int h = 0; h < 2 ; h++) {   				// To extend for better generalization, h < 2 can be replaced by a variable like comp_fields...
			 int index_c = MMSP::index(grid(x), h) ; 								    															  
				    set(c, index_c) = grid(n)[index_c];	
			}
			
			double gphiprime, pphiprime ;
			
			//Calculating gphiprime	 																						   
			for(int h = 2 ; h < length(grid(x)) ; h++){
					for(int j = 2 ; j < length(grid(x)) && j!=h ; j++){
					int index1 = MMSP::index(grid(x), h) ;
					int index2 = MMSP::index(grid(x), j) ;
					gphiprime += 2*grid(n)[index1]*grid(n)[index2]*grid(n)[index2] ;						
			}	
			}
			
			//Calculating pphiprime	 																						   
			for(int h = 2 ; h < length(grid(x)) ; h++){
					int index1 = MMSP::index(grid(x), h) ;
					pphiprime += pow(grid(n)[index1], 3.0) - pow(grid(n)[index1], 2.0) ;	// Need to declare pphiprime and gphiprime					
			}
			
			//Calculating gradients
			//vector<float> gradients_holder(dim) ;
			//vector<phi_type> gradients_holder ;
			MMSP::vector<MMSP:: sparse<float> > gradients_holder = gradient(grid, x);
			
			double gradient_phi = 0.0 ;
			double gradient_C = 0.0	;
				
			//Calculating McdFdC
			double McdFdC[2] ;
			for(int h = 0 ; h < 2 ; h++){
				gradient_phi = 0.0 ;
				gradient_C = gradients_holder[0][h] ;	
				for(int j = 2 ; j < length(grid(x)) ; j++){   //If only the first two fields represent composition, then shouldnt the phase begin from 2?
					gradient_phi += gradients_holder[0][j] ;
				}
				McdFdC[h] =  Mc*gradient_phi*((fsB-fsA)*pphiprime + (WB - WA)*gphiprime) + (Mc*gradient_C*R*T)/(grid(n)[h]*(1-grid(n)[h])) ;
			}
			
			//Updating the phase fields										
			sparse<int> s;
			for (int j = 0; j < dim; j++){
				for (int k = -1; k <= 1; k++) {
				  x[j] += k; //In the following, why are we starting from 2?
				  for (int h = 2; h < length(grid(x)); h++) {   				//grid(x) is similar to grid(n), and returns a sparse variable associated with x (operator() 	
					  														  //is defined in grid for both int and vector<int>). length(sparse s) -> s.length()-> returns 	
					  														  //														  the size of s.
				    int index = MMSP::index(grid(x), h) ; 					//functor - computes grid(x).index(h).  Remember, grid(x) is a sparse variable here. Returns 
				    															  //data[h].index associated with the sparse returned by grid(x). Remember, h is not the id of 
				    															   //the phase field variable. h is simply the position in the array. 
				    set(s, index) = 1;	
				  }
				  x[j] -= k;									// To get x[j] back to the original value. Could be made more clear by defining a temp 															   variable
				}
			}
			phi_type S = phi_type(length(s));  							// ~ float(s.length())

														// if only one field is nonzero,
														// then copy this node to update
			if (S < 2.0){
				update(n) = grid(n);								//This node is copied to the new grid, 'update'.
			} else {
														// compute laplacian of each field
				sparse<phi_type> lap = laplacian(grid, n); 
				
				sparse<phi_type> dFdp;

				for (int h = 0; h < length(s); h++) {	
														// compute variational derivatives
					int hindex = MMSP::index(s, h);   
					phi_type N = num_grains;     						// Why is the number of grains initialized as a float, when int can solve the purpose?
					
					phi_type phi_sq = 0.0; 							// phi_sq is the sum of the squares of the phi field
					phi_type dFall = 0.0;
					
														// Compute phi_sq value by taking the sum of all (defined) phi values
					for (int j = 0; j < length(s); j++) {
						int jindex = MMSP::index(s, j);	
						phi_sq += grid(n)[jindex]*grid(n)[jindex];
					}
					
														// Compute the dFall value by calculating all (defined) dFdp values and summing them
					for (int j = 0; j < length(s); j++) {
						int jindex = MMSP::index(s, j);  
						set(dFdp, jindex) = w * grid(n)[jindex] * (phi_sq - grid(n)[jindex]) - eps_sq * (lap[jindex]);
						dFall += dFdp[jindex];						// If dFdp is a sparse, shouldnt the statement be dFdp[jindex].value? It is a functor. Check 															   MMSP.sparse.hpp
					}
					
					store_type dpdt;
					set(dpdt, hindex) = - (M/N) * (N*dFdp[hindex] - dFall);
					phi_type value = grid(n)[hindex] + dt * dpdt[hindex];			// grid(n)[hindex] refers to the hindex of the sparse variable returned by grid(n)
					if (value > 1.0) value = 1.0;
					if (value < -0.001) value = 0.0;
					if (value > thresh) set(update(n), hindex) = value;	
					else if (grid(n)[hindex] != 0.0) set(update(n), hindex) = 0.0;
				}										//calculate interactions between interacting fields
			}											//perform calculations on non-zero fields
		} 												//loop over nodes		
		swap(grid, update);										//update is copied to grid.
	} 													//loop over steps
	ghostswap(grid);
	++iterations;
}

template <class T> std::ostream& operator<<(std::ostream& o, sparse<T>& s) {
	o<<"    Index  Value\n";
	for (int i=0; i<length(s); ++i) {
		int index = MMSP::index(s, i);
		o<<"    "<<std::setw(5)<<std::right<<index<<"  "<<s[index]<<'\n';
	}
	return o;
}

void generate(int dim, char* filename) {

	int id=0;
	#ifdef MPI_VERSION
 	id=MPI::COMM_WORLD.Get_rank();
	#endif

	std::string search_name(filename);
		std::cout << "Read file name"<<std::endl;
														
	bool planar = (search_name.find("planar") != std::string::npos);
	bool circle = (search_name.find("circle") != std::string::npos);
	bool generated = (search_name.find("file") != std::string::npos);
	bool grad = (search_name.find("grad") != std::string::npos);
			//std::cout << "Booleaned everything"<<std::endl;
	int nx = 100;    											
	int ny = 100;
	int nz = 100;
	int num_grains = 2;

	std::vector<int> domain_scale;										//Each component of domain_scale stores the size of the domain in that particular dimension.

	int bias = 0;
	int num_bins = 0;
	if (generated or grad){
		std::cout << "In loop for grad"<<std::endl;
		std::vector<std::string> splits = split(search_name, '.');
		std::string name_root = splits[0]; 
		std::vector<std::string> metadata = split(name_root, 'x'); 					
		for (int i = 0; i < dim; i++){
			domain_scale.push_back(std::atoi(metadata[i+1].c_str()));     																			  
		}
		num_grains = std::atoi(metadata[dim+1].c_str());  
		if (dim > 0) nx = domain_scale[0];
		if (dim > 1) ny = domain_scale[1];
		if (dim > 2) nz = domain_scale[2];
		if (grad){
			num_bins = std::atoi(metadata[dim+2].c_str());
			bias = std::atoi(metadata[dim+3].c_str());
		}
	}
	else{
		for (int i = 0; i < dim; i++){
			domain_scale.push_back(0);   			
		}
	}

	if(planar){

		if (dim == 1){
			MMSP::grid<1,store_type > grid (4,0,nx); 						//Define the MMSP grid. The input paramters are number of fields, min node in dim 1, max node in dim 1,....
	
			for (int d = 0; d < dim; d++) MMSP::dx(grid, d) = 1.0; 					//Lx/nx ; Thus, Lx = Ly = 100.0 ;
			if (id == 0) std::cout << "Creating 2-grain planar interface."<<std::endl;

		
			makeSplit<1>(grid);
			if (id == 0) std::cout << "Saving..." << std::endl;
			output(grid, filename); 								//write out initialized grid
			if (id == 0) std::cout << "Grid saved as " << filename << std::endl;
		} else if (dim == 2){
			MMSP::grid<2,store_type > grid (2,0,nx,0,ny);
	
			for (int d = 0; d < dim; d++) MMSP::dx(grid, d) = 1.0;
			
			if (id == 0) std::cout << "Creating 2-grain planar interface."<<std::endl;
			makeSplit<2>(grid);
			if (id == 0) std::cout << "Saving..." << std::endl;
			output(grid, filename); 								//write out initialized grid
			if (id == 0) std::cout << "Grid saved as " << filename << std::endl;
		} else if (dim == 3){
			MMSP::grid<3,store_type > grid (2,0,nx,0,ny,0,nz);
	
			for (int d = 0; d < dim; d++) MMSP::dx(grid, d) = 1.0;
			if (id == 0) std::cout << "Creating 2-grain planar interface."<<std::endl;
		
			makeSplit<3>(grid);
			if (id == 0) std::cout << "Saving..." << std::endl;
			output(grid, filename); 								//write out initialized grid
			if (id == 0) std::cout << "Grid saved as " << filename << std::endl;
		}
		

														//std::cout << "Saving..." << std::endl; 
														//output(grid, filename); //write out initialized grid 
														//std::cout << "Grid saved as " << filename << std::endl;
	} 

	        
	else if(circle){
		if(dim == 3 or dim == 1) {
			std::cerr<<"Error: Dimensionality not supported."<<std::endl;
			exit(1);
		}
				
		MMSP::grid<2,store_type > grid (2,0,nx,0,ny);

	
		MMSP::dx(grid, 0) = 1.0;
		MMSP::dx(grid, 1) = 1.0;
		if (id == 0) std::cout << "Creating circular grain-in-grain test grid."<<std::endl;
		
		float radius = 20;   										//Define the radius by default.
		makeCenterCircle<2>(grid, radius);
		
		if (id == 0) std::cout << "Saving..." << std::endl;
		output(grid, filename); 									//write out initialized grid
		if (id == 0) std::cout << "Grid saved as " << filename << std::endl;
	} 
	
	else if(generated)
	{
		std::cout << "In loop for generated"<<std::endl;
		MMSP::grid<2,store_type > grid (num_grains,0,nx,0,ny);  					//Define the MMSP grid. The first input parameter refers to the number of grains.

		if(dim == 3 or dim == 1) {
			std::cerr<<"Error: Dimensionality not supported."<<std::endl;
			exit(1);
		}
		
														//MMSP::grid<2,store_type > grid (num_grains,0,nx,0,ny);

	
		MMSP::dx(grid, 0) = 1.0;									//Lx/nx;
		MMSP::dx(grid, 1) = 1.0;									//Ly/ny;

		if (id == 0) std::cout << "Reading in generated microstructure."<<std::endl;

		treadin<2,store_type>(grid);
		
		if (id == 0) std::cout << "Saving..." << std::endl;
		output(grid, filename); 									//write out initialized grid
		if (id == 0) std::cout << "Grid saved as " << filename << std::endl;
	} 
	
	else if(grad){
		std::cout << "In loop for grad"<<std::endl;
		if (id == 0) std::cout << "Generating biased Voronoi tessellation."<<std::endl;	
														//re-cast nx so that it is definitely integrally divisible by num_bins
		nx = int(nx/num_bins)*num_bins;
		if (dim == 1) {
			MMSP::grid<1,store_type > grid (num_grains,0,nx);
			for (int d = 0; d < dim; d++) MMSP::dx(grid, d) = 1.0;
			tessellate<1,phi_type >(grid, num_grains, num_bins, bias);
			if (id == 0) std::cout << "Saving..." << std::endl;
			output(grid, filename); 								//write out initialized grid
			if (id == 0) std::cout << "Grid saved as " << filename << std::endl;
		} else if (dim == 2) {
			MMSP::grid<2,store_type > grid (num_grains,0,nx,0,ny);
			for (int d = 0; d < dim; d++) MMSP::dx(grid, d) = 1.0;
			tessellate<2,phi_type >(grid, num_grains, num_bins, bias);
			if (id == 0) std::cout << "Saving..." << std::endl;
			output(grid, filename); 								//write out initialized grid
			if (id == 0) std::cout << "Grid saved as " << filename << std::endl;
		} else if (dim == 3) {
			MMSP::grid<3,store_type > grid (num_grains,0,nx,0,ny,0,nz);
			for (int d = 0; d < dim; d++) MMSP::dx(grid, d) = 1.0;
			tessellate<3,phi_type >(grid, num_grains, num_bins, bias);
			if (id == 0) std::cout << "Saving..." << std::endl;
			output(grid, filename); 								//write out initialized grid
			if (id == 0) std::cout << "Grid saved as " << filename << std::endl;
		}
		
	}
	
	else {
		std::cerr<<"Error: No initialization condition selected."<<std::endl;
		exit(1);
	}
	




} 														// namespace MMSP

}


#endif

#include "/home/arun/mmsp/include/MMSP.main.hpp"


