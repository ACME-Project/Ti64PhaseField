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
#include<cstdio>
#include<cstdlib>
#include<cassert>
#include"/home/arun/mmsp/include/MMSP.hpp"
#include "Tensor.hpp"
#include "Eigenstrains.hpp"


typedef float phi_type; 											// Note : This is used extensively in sp-initialize.cpp
typedef MMSP::sparse<phi_type> store_type; 									// Note : This is used extensively in sp-initialize.cpp

template <int dim> void makeSplit(MMSP::grid<dim, store_type >& grid);   					//Always make sure that the function declaration is before the definition. Here, since 															  makeSplit is defined in sp-initialize.cpp, the function is declared before that file is 															  called/referenced. 

														//included after typedefs for proper variable types in the functions of the included files
#include"sp-graingrowth.hpp" 
#include"sp-initialize.hpp"

														//Function declarations 
std::vector<std::string> split(const std::string &text, char sep); 					        // Returns a vector of strings?
void print_progress(const int step, const int steps, const int iterations);

														//included after typedefs for proper variable types in the functions of the included files

#include "tessellate.hpp"





namespace MMSP {

void print_progress(const int step, const int steps, const int iterations) {					//This is the program that prints progress over time
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
	} else if ((20 * step) % steps == 0) std::cout<<"• "<<std::flush;					// std::flush is equivalent to std::endl or the flush that we used in parallel class last 															   semester		
}

std::vector<std::string> split(const std::string &text, char sep) {						//Simply make sure that all the usages of this function appear before after the function 															  definition, or else you might get an undefined reference error. 
														// This functions splits the string into different substrings, separated by the char sep. 
  std::vector<std::string> tokens;
  std::size_t start = 0, end = 0;
  while ((end = text.find(sep, start)) != std::string::npos) {
    tokens.push_back(text.substr(start, end - start));								//text.substr(start, end-start) gives a substring that starts from start character and 															  ranges till end-start.
    start = end + 1;
  }
  tokens.push_back(text.substr(start));		
  return tokens;
}


														//template <int dim> void update(MMSP::grid<dim, store_type >& grid, int steps) ;
template <int dim> void update(MMSP::grid< dim , store_type>& grid, int steps) 
{    														// computation of the model 
	initialize_alpha_eigen();
	float dx = MMSP::dx(grid, 0);   									// functor
	int id=0;
	#ifdef MPI_VERSION
 	id=MPI::COMM_WORLD.Get_rank();     
	#endif
	
														//thresh is the minimum phi value imposed to minimize mathematical noise in fields
	phi_type thresh = 1.0e-3;
	
	phi_type width = 1.0;											//width of stable interface
	phi_type M = 1.0;											//mobility of field motion
	phi_type dt = dx / (20.0 * M);									 	//timestep calculated as function of resolution, 10x less than Courant maximum
	phi_type gamma = 1.0/3.0;										//energy contribution from gradient component
	phi_type w = 3.0 * (gamma/width);									//energy contribution from double-well component
	phi_type eps_sq = 3.0 * width * gamma;
	
	tensor eigen_strain_field[nodes(grid)] ;
	tensor homogeneous_strain ;
	tensor total_strain_field[nodes(grid)] ; 
	tensor elastic_strain_field[nodes(grid)] ; 
	phi_type strain_energy[nodes(grid)];
	
	double Cijkl[6][6] ;

	for(int i = 0 ; i < 6 ; i++)
	{
		for(int j = 0 ; j < 6 ; j++)
		{
			Cijkl[i][j] = 0 ;
		}
	}

	Cijkl[0][0] = 1.0 ;
	Cijkl[0][1] = 1.0 ;
	Cijkl[0][2] = 1.0 ;
	Cijkl[2][2] = 1.0 ;
	Cijkl[3][3] = 1.0 ;
	Cijkl[1][2] = 1.0 ;
	
														// Spatially-varying simulation parameters calculated
														// Initialize vairiables
	const int num_grains =  fields(grid); 									//functor
	
	static int iterations = 1;
	
	for (int step = 0; step < steps; step++) {
		if(id == 0)	print_progress(step, steps, iterations);
		ghostswap(grid);    										//exchange of information 
		
		MMSP::grid<dim, store_type > update(grid);							//A new MMSP grid is defined by the name of update, to which the original grid is passed as 															  initialization parameter. This is equivalent to creating a copy of the original grid. 														  Potential for improvement : lots of duplications which deter easy understanding of MMSP

		for (int n = 0; n < nodes(grid); n++){
			/*
			 * Pseudocode for getting coordinates of interfaces - Write a separate function for this.
			 * 1.) Define array of vectors length equal to number of fields. Is array of vectors possible? 
			 * 2.) Define a dynamic array to store triple junctions.
			 * 3.) At each node, check the number of active field variables.
			 * 4.) If the number of field variables is > 3, add the node id to the vector of all field variables. For example, if node n has 2, 3, and 5 pf variables, add n to array ids 1,2,and 4. 
			 * 5.) Add n to the array of triple junctions.  
			 * 6.) For each phase field variable's node, calculate the prime product of the fields. Segregate the triple points according to their prime products and then calculate the average 
			 * 	   coordinate for each product. 
			 * 7.) Create a curated array in which the average coordinates of TJ associated with a grain are stored.
			 * 
			 * */
			vector<int> x = position(grid, n);  							//position() returns a vector of coordinates at the nth node. position, fields, nodes are all examples of functors.
			
			//for(int n = 0 ; n < nodes(grid) ; n++)
			//{	
				for(int variant = 0 ; variant < 12 ; variant++)  //MAKE SURE TO UNDERSTAND HOW DIMENSIONS AFFECT THIS
				{
					for(int i = 0 ; i <=2 ; i++)
					{
						for(int j = 0 ; j <=2 ; j++)
						{
				
							eigen_strain_field[n].matrix[i][j]+=1.0*eigen_alpha[variant].matrix[i][j] ; //Need to convert the 1.0 into a h(phi)
				
							homogeneous_strain.matrix[i][j]+=1.0*eigen_alpha[variant].matrix[i][j] ;
							
							elastic_strain_field[n].matrix[i][j] = homogeneous_strain.matrix[i][j] - eigen_strain_field[n].matrix[i][j] ;
						}
					}
				}
			//}						
			
			//for(int n = 0 ; n < nodes(grid) ; n++)
			//{	
				strain_energy[n] = elastic_strain_field[n].matrix[0][0]*Cijkl[0][0]*elastic_strain_field[n].matrix[0][0] + 
				elastic_strain_field[n].matrix[0][0]*Cijkl[0][1]*elastic_strain_field[n].matrix[1][1] + elastic_strain_field[n].matrix[0][0]*Cijkl[0][2]*elastic_strain_field[n].matrix[2][2] +  
				elastic_strain_field[n].matrix[2][2]*Cijkl[2][2]*elastic_strain_field[n].matrix[2][2] + elastic_strain_field[n].matrix[1][2]*Cijkl[3][3]*elastic_strain_field[n].matrix[1][2] + 
				Cijkl[0][1]*elastic_strain_field[n].matrix[1][1]*Cijkl[1][2]*elastic_strain_field[n].matrix[2][2] ;
			//}									  
			
														// determine nonzero fields within
														// the neighborhood of this node
													   	// (2 adjacent voxels along each cardinal direction)
			sparse<int> s;
			for (int j = 0; j < dim; j++){
				for (int k = -1; k <= 1; k++) {
				  x[j] += k;
				  for (int h = 0; h < length(grid(x)); h++) {   				//grid(x) is similar to grid(n), and returns a sparse variable associated with x (operator() 															  is defined in grid for both int and vector<int>). length(sparse s) -> s.length()-> returns 															  the size of s.
				    int index = MMSP::index(grid(x), h) ; 					//functor - computes grid(x).index(h).  Remember, grid(x) is a sparse variable here. Returns 															  data[h].index associated with the sparse returned by grid(x). Remember, h is not the id of 															  the phase field variable. h is simply the position in the array. 
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
					value+=(1-2*grid(n)[hindex]*grid(n)[hindex] + pow(grid(n)[hindex], 4))*1.875*strain_energy[n] ;
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
//	Use case for this code is as follows: generate() will create MMSP grid files which can be used to store field data.
//	The type of grid data to be loaded is selected by the name of the filename that you declare, for example: 
//		 planar.000000.dat to generate a two-grain, planar interface
//		 circle.000000.dat to generate a radially-symmetric grain-in-grain structure
//		 filexXxYxZxN.000000.dat to read an appropriatelly generated sharp-interface input file 
//			(dimensions X, Y, Z and with N grains eg filex100x100x15.000.dat for a grid of 100x100 with 15 grains)
//		 gradxXxYxZxNxBxb.000000.dat to generate a Voronoi tessellated domain of dimensions X, Y, Z and N grain seeds dispersed in B bins with a bias of b more grains in each consecutive bin

	int id=0;
	#ifdef MPI_VERSION
 	id=MPI::COMM_WORLD.Get_rank();
	#endif

	std::string search_name(filename);
		std::cout << "Read file name"<<std::endl;
														//Check if the filename has any of the following phrases. The string search function 															  returns npos if no matches are found. npos, in this context, refers to the end of file. 
														// search the filename for keywords to direct the initial condition generation routine
	bool planar = (search_name.find("planar") != std::string::npos);
	bool circle = (search_name.find("circle") != std::string::npos);
	bool generated = (search_name.find("file") != std::string::npos);
	bool grad = (search_name.find("grad") != std::string::npos);
			std::cout << "Booleaned everything"<<std::endl;
	int nx = 100;    											//Default values of nx,ny, and nz.
	int ny = 100;
	int nz = 100;
	int num_grains = 2;

	std::vector<int> dimensions;										//Each component of dimensions stores the number of grid points in that dimension.
	//if (generated){// determining the dimensions based off filename data

	int bias = 0;
	int num_bins = 0;
	//std::vector<int> dimensions;
	if (generated or grad){// determining the dimensions based off filename data
		std::cout <<"In loop for grad"<<std::endl;
		std::vector<std::string> splits = split(search_name, '.');
		std::string name_root = splits[0]; //"filex#####x#####
		std::vector<std::string> metadata = split(name_root, 'x'); 					//This splits creates individual strings in the title separated by x. metadata contains nx, 															  ny, nz, number of grains and the bias. 
		for (int i = 0; i < dim; i++){
			dimensions.push_back(std::atoi(metadata[i+1].c_str()));     				//atoi interprets a string of numbers as integer values. c_str converts the std::string of 															  C++ into a C string.
		}
		num_grains = std::atoi(metadata[dim+1].c_str());  
		if (dim > 0) nx = dimensions[0];
		if (dim > 1) ny = dimensions[1];
		if (dim > 2) nz = dimensions[2];
		if (grad){
			num_bins = std::atoi(metadata[dim+2].c_str());
			bias = std::atoi(metadata[dim+3].c_str());
		}
	}
	else{
		for (int i = 0; i < dim; i++){
			dimensions.push_back(0);   			
		}
	}

	if(planar){

		if (dim == 1){
			MMSP::grid<1,store_type > grid (2,0,nx); 						//Define the MMSP grid. The first input parameter refers to the number of grains.
	
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
	} else if(generated){
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
	} else if(grad){
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

#include"MMSP.main.hpp"


