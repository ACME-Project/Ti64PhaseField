// Initialize.cpp

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include "MersenneTwister.h"

template<int dim, typename T>
void treadin(MMSP::grid<dim, store_type >& grid);   

/*


There are two classes of MMSP used by this program, the sparse class and the grid class. The sparse class defines the data structure that
will store the list of the phase field variables and their respective values. The grid class represents the simulation domain's grid (or the portion of the grid computed
by the processor). Both the classes have an attribute named data. While the sparse class's data is an array of struct 'item' objects, each having an index and a value, the grid
class's data is an array of class 'sparse' objects. The former array's size is equal to the number of phase field variables, while the latter array's is equal to the number of cells. 

*/

//makeSplit - Initializes the microstructure into two grains and the initial interface is at 1/2 the distance of the 1st dimension.
template <int dim> void makeSplit(MMSP::grid<dim, store_type >& grid) {
	for(int n = 0; n < nodes(grid); n++){    //nodes(grid) returns the number of nodes in that grid. 
		MMSP::vector<int> x = position(grid, n);   //vector x has number of components equal to the number of dimensions, position(grid, n) returns the n coordinates 
												   //corresponding to the n-dimensions , at the node n of the grid 'grid'. x[0] thus corresponds to the 												     coordinate of the dimension 1, 
												   //and x1(grid) would correspond to the length of the grid in dimension no.1 
		store_type newGrain;
		if (x[0] < (x1(grid)/2))
		{ //if on the left side of the grid
			set(newGrain, 1) = 1.0;
			/*set() is a member function of the sparse class.  There is a template defined as follows : 
template <int ind, typename T> T& set(const target<0, ind, sparse<T> >& s, int index) {	return s.set(index);} . Thus, set(newGrain, index) is equivalent to newGrain.set(index). This
function updates the data variable of the sparse object newGrain and appends the index to the list of phase field variables. It then returns the address for the value instance of that
data array member, to which we write the requisite value. */
		}
		else{ //else on the right side of the grid
			set(newGrain, 0) = 1.0;
		}
		grid(n) = newGrain;
		/* grid(n) , returns the address for the sparse data structure associated with the node n of the grid, to which we write our computed sparse variable. grid(n) is a typical 
		example of a functor. */
	}
}

/*makeCenterCircle initializes a circular nucleus, with the center being the center of the grid. I believe that g1 and g0 represent the coordinates of the two extremes of the grid, for 
each dimension. */

template <int dim> void makeCenterCircle(MMSP::grid<dim, store_type >& grid, float radius) {
	for(int n = 0; n < nodes(grid); n++){
		MMSP::vector<int> x = position(grid, n);
		float sum = 0;
		for (int i = 0; i < dim; i++)
			sum += pow( ((x[i]-(g1(grid,i)-g0(grid,i)))/2.0), 2);
		float distance = sqrt(sum);
		
		bool inCircle = distance <= radius;
		store_type newGrain;
		if (inCircle){
			set(newGrain, 1) = 1.0;
		}
		else{
			set(newGrain, 0) = 1.0;
		}
		grid(n) = newGrain;
	}
}

/*

NOTES ABOUT THE TREADIN FUNCTION

The nodes of the MMSP grid and those of the input file should be identical. At the least, the number of nodes need to be the same.

Since grain_list contain the grain ids in the order of the aapearance in the input dat file, the id map of the reconstructed microstructure is different from the id_map of the input 
microstructure.

*/

// Sparse fielding function to import a microstructure from a grain-ID indexed matrix
// currently only supported for 2D
template<int dim, typename T>
void treadin(MMSP::grid<dim, store_type >& grid) {	
	// This loop reads in the input files and creates the arrays all_list and grain_list necessary for the grid
	// all_list contains the direct data from the Taylor output.
	// grain_list is a processed, sorted list of grain ID numbers without repeats.
	// As written, the data file input is not automatic and certain things (like filename and array sizes) are hard coded.
	int x_dim = g1(grid,0) - g0(grid,0); // I think g1 and g0 return the distance in that direction. 
	int y_dim = g1(grid,1) - g0(grid,1);
	std::ostringstream fileNameStream;	//String buffer?
	fileNameStream << x_dim << "x" << y_dim << ".txt";
	std::string filename = fileNameStream.str();
	std::ifstream taylor_dat(filename.c_str());

	std::vector<int> all_list;
	std::vector<int> grain_list;
	
	std::cout<<"Begin loading dataset."<<std::endl;
	
	std::string line;
	if (taylor_dat.is_open()){
		while (getline(taylor_dat, line)){
			int number;
			std::stringstream buffer(line); 
			buffer >> number;
			all_list.push_back(int(number));
			bool in_list = false;
			for(int i = 0; i < grain_list.size(); i++) {
				if (number == grain_list[i])
					in_list = true;
			}
			if (! in_list)
				grain_list.push_back(number);
		}
		taylor_dat.close();
		std::cout<<"Taylor data imported successfully."<<std::endl;
		
	}
	else{
		std::cerr<<"ERROR: "<<filename<<" not found!"<<std::endl;
  		std::exit(-1);
	}

	// Loop through the grid and assign new, shorter grain id values
	for (int n = 0; n < nodes(grid); n++) {
		const MMSP::vector<int> x=position(grid,n);
		int list_pos = x[1] * x_dim + x[0];
		int grain_ID = all_list[list_pos];   //We can do this because grain ids in the input file are assumed to be arranged in a sequential manner.
		int new_ID = grain_ID;
		for (int i = 0; i < grain_list.size(); i++) {
			if (grain_list[i] == grain_ID) {
				new_ID = i;
			}
		}

		store_type newGrain;
		set(newGrain, new_ID) = 1.;
		grid(n) = newGrain;
	}

	std::cout<<"Completed Microstructure Recreation."<<std::endl;
} // sparse treadin 

