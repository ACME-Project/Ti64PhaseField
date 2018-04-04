// Initialize.cpp

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
//#include "MersenneTwister.h"

template<int dim, typename T>
void treadin(MMSP::grid<dim, store_type >& grid);   

template <int dim> void makeSplit(MMSP::grid<dim, store_type >& grid) {
	for(int n = 0; n < nodes(grid); n++){    
		MMSP::vector<int> x = position(grid, n);   
		store_type newGrain;
		if (x[0] < (x1(grid)/2))
		{
			set(newGrain, 0) = 0.075;
			set(newGrain, 1) = 0.025 ;
			set(newGrain, 3) = 1.0;
		}
		else if( x[0] >= (x1(grid)/2))
		{			  
			set(newGrain, 0) = 0.035 ;   		// Does 0 here mean 0th element or element with the index 0?
			set(newGrain, 1) = 0.050 ;
			set(newGrain, 2) = 1.0;
		}
		grid(n) = newGrain;
	}
}

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


template<int dim, typename T>
void treadin(MMSP::grid<dim, store_type >& grid) {	
	int x_dim = g1(grid,0) - g0(grid,0);  
	int y_dim = g1(grid,1) - g0(grid,1);
	std::ostringstream fileNameStream;	
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

	for (int n = 0; n < nodes(grid); n++) {
		const MMSP::vector<int> x=position(grid,n);
		int list_pos = x[1] * x_dim + x[0];
		int grain_ID = all_list[list_pos];   
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

