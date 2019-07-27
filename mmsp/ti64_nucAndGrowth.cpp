/* Simultaneous nucleation and growth of alpha variants. Anisotropy is introduced through anisotropic gradient energy coefficients.
 * 
 * Mode of computation: Serial
 * 
 * Code developed by: Arun Baskaran
 */

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
#include"MMSP.hpp"

typedef float phi_type; 											
typedef MMSP::sparse<phi_type> store_type; 	
using namespace std;


//----------- Data structure to store a variant's six eigenstrains (three, in 2D) under the Voigt notation---------//
class sfts
{
	
	public:
	
	double e1, e2, e3, e4, e5, e6 ;
	
	sfts()
	{
		e1 = 0.0 ;
		e2 = 0.0 ;
		e3 = 0.0 ;
		e4 = 0.0 ;
		e5 = 0.0 ;
		e6 = 0.0 ;
	}
};

//------------------------------------------------------------------------------------------------------------------//


double epsi[12][3][3] ;
double R = 8.314 ;
double T = 1023 ;
double dt = 8.0e-3 ;
double L = 0.239 ;  //Kinetic coefficient
int fields = 15 ;   // 12 alpha variants, 2 composition fields, 1 nucleation field. The variant fields are numbered 1-12, c_Al = no.20, c_V = no.21, and the nucleation field = no.15
const int dim = 3;
const int nx = 20;
const int ny = 20 ;
const int nz = 20;
double G_normalize = R*T ;
int steps = 70000 ;  

//--------------------------------Parameters for inducing explicit nucleation--------------------------------//
double deltag = 0.1 ;  
double kappa1 = 0.05 ; 
double kappa2 = 0.1 ; 
//-----------------------------------------------------------------------------------------------------------//


//--------------------------------Auxillary files containing helper functions--------------------------------//
#include "GradECoeff.hpp"
#include "Thermo.hpp"
//-----------------------------------------------------------------------------------------------------------//


//----------Container variables required for computing evolution of phase field and composition fields-------//
void thermo_auxillary_terms(MMSP::vector<store_type> gradient, MMSP::vector<store_type> gradientsq, double c_Al, double c_V) ;
double dGAlpha_dAl,	dGBeta_dAl, dGAlpha_dV, dGBeta_dV, del_dGAlpha_dAl, del_dGBeta_dAl, del_dGAlpha_dV, del_dGBeta_dV, delsq_dGAlpha_dAl, delsq_dGAlpha_dV, delsq_dGBeta_dAl, delsq_dGBeta_dV ;
double G_Al_alpha, G_Ti_alpha, G_V_alpha, G_Al_beta, G_Ti_beta, G_V_beta , W_Al, W_V, W_Ti ; 
const int node_total = nx*ny*nz ; 
double G_Alpha[node_total], G_Beta[node_total], W[node_total] ; 
double hphi[node_total], hphiprime[node_total], hphidoubleprime[node_total];
double gphi[node_total], gphiprime[node_total], gphidoubleprime[node_total];
const double Lx = g1(grid,0) - g0(grid,0) ;
const double Ly = g1(grid,1) - g0(grid,1) ;
const double dx = MMSP::dx(grid, 0) ;
//------------------------------------------------------------------------------------------------------------//



MMSP::grid<dim, store_type> grid(fields, 0, nx, 0, ny, 0, nz) ;


//-----------------------------------Field initialization-----------------------------------------------------//
void initialize()
{

	for(int n = 0; n < nodes(grid); n++)
	{
		MMSP::vector<int> x = position(grid, n);
		
		store_type newGrain;
		
		set(newGrain, 15) = 1.0;  //Field with the index no.15 corresponds to the nucleation field. The variable will hold a value of 1 if nucleation can happen on this node, i.e., there have not been any variant nucleated yet on this node. Since variants are introduced through a nucleation function, initially all the nodes will have a value of 1 on field 15.
		
		grid(n) = newGrain;	
		
	}
	
}

//------------------------------------------------------------------------------------------------------------//

int main()
{
	
	
string file_name ;  //For file I/O
string file_name1 ;
string file_name2 ;

	

initialize();
initialize_epsi();


//-----------------------------------Boundary conditions----------------------------------//
MMSP::b0(grid,0) = MMSP::Dirichlet;
MMSP::b1(grid,0) = MMSP::Dirichlet;

MMSP::b0(grid,1) = MMSP::Dirichlet;
MMSP::b1(grid,1) = MMSP::Dirichlet;

MMSP::b0(grid,2) = MMSP::Dirichlet;
MMSP::b1(grid,2) = MMSP::Dirichlet;
//---------------------------------------------------------------------------------------//




for(int t = 0 ; t<steps ; t++)
{
	MMSP::grid<dim, store_type > update(grid);
	
	MMSP::b0(update,0) = MMSP::Dirichlet;
	MMSP::b1(update,0) = MMSP::Dirichlet;

	MMSP::b0(update,1) = MMSP::Dirichlet;
	MMSP::b1(update,1) = MMSP::Dirichlet;

	MMSP::b0(update,2) = MMSP::Dirichlet;
	MMSP::b1(update,2) = MMSP::Dirichlet;
	
	
	int no_of_cells = (nx*ny*nz)/(5*5*5) ; //No. of cells to be checked for nucleation. 
	int cells_in_one_dir = nx/5 ; //No. of such cells in one dimension.
	
	
	
	//-----------------------------------------Nucleation Event-------------------------------------------------//
	if((t%100==0)&(t<300))     //Nucleation is only induced three times, at t = 0, 100, and 200 timesteps. 
	{
		MMSP::vector<int> x(5) ; 
		MMSP::vector<int> y(5) ; 
		MMSP::vector<int> z(5) ; 
		MMSP::vector<int> temp_pos(3); 
		for(int i = 0 ; i < cells_in_one_dir ; i++)
		{
			for(int j = 0 ; j < cells_in_one_dir ; j++)
			{
				for(int k = 0 ; k < cells_in_one_dir ; k++)
				{
					int flag = 0 ;
				
				
					x[0] = i*5 ; x[1] = i*5+1 ; x[2] = i*5+2 ; x[3] = i*5 + 3 ; x[4] = i*5+4 ;
				
					y[0] = j*5 ; y[1] = j*5+1 ; y[2] = j*5+2 ; y[3] = j*5 + 3 ; y[4] = j*5+4 ;
				
					z[0] = k*5 ; z[1] = k*5+1 ; z[2] = k*5+2 ; z[3] = k*5 + 3 ; z[4] = k*5+4 ;
					for(int p = 0; p < 5; p++)
					{
						for(int q = 0 ; q < 5 ; q++)
						{
							for(int r = 0 ; r < 5 ; r++)
							{
							
								temp_pos[0] = x[p] ; temp_pos[1] = y[q] ; temp_pos[2] = z[r] ;
								flag+= grid(temp_pos)[15] ;    //Checking if there are nodes with their field no.15 = 1
							}
						}
					}
				
					if(flag > 0) //Since a nucleation event will fill out an entire cell, either all nodes in a cell will have field no.15  = 0 or all of them will have field no.15 as 1.
					{
					
						double jstar = kappa1*exp(-kappa2*deltag);   
						double prob = 1 - exp(-jstar) ;
						double random_no = (double)(rand() % 1000) ;
						double decision_var = random_no/1000.0 ;
						if(prob > decision_var)  //Nucleation will happen if probability of nuc is > than the decision variable.
						{
							store_type newGrain;
							int random_index = (int)(rand() % 3) + 1 ;  //nucleating a random index among 1,2, and 3.
							cout<<"Nucleating index :"<<random_index<<endl;
							for(int m = 0 ; m < nodes(grid) ; m++)
							{
								set(grid(m), random_index) = 0.0 ;
							}
							for(int p = 0; p < 5; p++)
							{
								for(int q = 0 ; q < 5 ; q++)
								{
									for(int r = 0 ; r < 5 ; r++)
									{
										temp_pos[0] = x[p] ; temp_pos[1] = y[q] ; temp_pos[2] = z[r] ;
									
										if((pow((temp_pos[0]*MMSP::dx(grid,0)-x[2]*MMSP::dx(grid,0)),2) + pow((temp_pos[1]*MMSP::dx(grid,1)-y[2]*MMSP::dx(grid,0)),2) + pow((temp_pos[2]*MMSP::dx(grid,2)-z[2]*MMSP::dx(grid,0)),2)) < 3*MMSP::dx(grid,0)*MMSP::dx(grid,0))
										{
											set(grid(temp_pos), 15) = 0.0 ;   //Setting the nucleation field to 0. No more nucleation can happen on this node. 								
											set(grid(temp_pos), random_index) = 1.0; 	
										}								 
									}
								}
							}	
						}
					}
				}
			}
		}
	}
	
	//-------------------------------------------------File I/O------------------------------------------------//
	if(t%1000==0) 
	{
		std::ostringstream fileNameStream("pf_") ;
		fileNameStream<<"pf_"<<t<<".txt";
		file_name = fileNameStream.str() ;
		ofstream myfile;
		myfile.open(file_name.c_str());
		for(int n = 0 ; n < nodes(grid) ; n++)    
		{
			MMSP::vector<int> s = position(grid,n) ;
			double sum = 0 ;
			for(int h = 0 ; h < length(grid(s)) ; h++)
			{
				int hindex = MMSP::index(grid(s),h) ;
				if(hindex <= 12)
				{
					sum+=grid(s)[hindex]*grid(s)[hindex] ;
				}
			}
			
			myfile<<s[0]<<","<<s[1]<<","<<s[2]<<","<<sum<<"\n" ;   
			
			
		}
	
		myfile.close();
		
		
	}
	//---------------------------------------------------------------------------------------------------------//
	
	
	
	//----------------------------------------Updating the grid------------------------------------------------//
	for(int n = 0 ; n < nodes(grid) ; n++)
	{
		MMSP::vector<int> x = position(grid, n); 
		
		
		MMSP::vector<store_type> gradient = grad(grid, x) ;
		MMSP::vector<store_type> gradientsq = gradsq(grid, x) ;
		
		double phi_grad = 0 ;
		double phi_gradientsq = 0 ;
		for(int h = 0 ; h < length(grid(x)) ; h++)
		{
			if(h <= 12)
			{
				
				int hindex = MMSP::index(grid(x),h) ;
				for(int i = 0 ; i < dim ; i++)
				{
					phi_grad+= gradient[i][hindex] ; 
					phi_gradientsq+= gradient[i][hindex] ; 
				}
			}
		}
		
		
		W[n] = 1.1432 ;   //Depth of the double well function
				
		MMSP::sparse<phi_type> dFdp;
			
		phi_type phi_sq = 0.0; 	
		phi_type phi_sum = 0.0 ;						
		phi_type dFall = 0.0;
		phi_type dFall_no = 0 ;
					
														
		for (int j = 0; j < length(grid(n)); j++) 
		{
			int jindex = MMSP::index(grid(n), j);
			if(jindex<=12) //Only fields from 1 to 12 are counted for phi_sum
			{	
				phi_sq += grid(n)[jindex]*grid(n)[jindex];
				phi_sum += grid(n)[jindex] ;
			}
		}
				
		phi_type phi_beta = 1 - phi_sum ; //The field representing beta-matrix is a dependent one. 
			
			
		//------------Calculate anisotropic laplacian and dfdp for every variant on the node---------------//	
		for (int j = 0; j < length(grid(n)); j++)   
		{
			int jindex = MMSP::index(grid(n), j);  
			phi_type lap_aniso = 0 ;			
			MMSP::vector<int> y = position(grid, n);
			MMSP::vector<double> gradients(3) ;
			if(jindex <= 12)
			{
				for(int p = 0 ; p < dim ; p++)
				{
					for(int q = 0 ; q < dim ; q++)
					{
						lap_aniso+= gradient[p][jindex]*epsi[jindex-1][p][q]*gradient[q][jindex];  								
					}
				}
					
				double interaction_energy = 0 ;
				interaction_energy+= W[n]*grid(n)[jindex]*pow(phi_beta, 2);
				interaction_energy-= W[n]*pow(grid(n)[jindex],2)*(phi_beta);
					
				for(int k = 0 ; k < length(grid(n)); k++)
				{
					int tempindex = MMSP::index(grid(n), k); 
					if((tempindex!=jindex)&&(tempindex <=12)){ interaction_energy+= pow(W[n],2)*pow(grid(n)[tempindex],2)*grid(n)[jindex]; interaction_energy-= W[n]*pow(grid(n)[tempindex],2)*(phi_beta); }
				}
				set(dFdp, jindex) = grid(n)[jindex] * (phi_sq - grid(n)[jindex]) + interaction_energy  - (lap_aniso); 
					
				dFall += dFdp[jindex];	
				dFall_no+=1;
			}							
		}
		//---------------------------------------------------------------------------------------------------//


		//-------------------------Updating every variant associated with the node n-------------------------//
		for (int h = 0; h < length(grid(n)); h++) 
		{				
			int hindex = MMSP::index(grid(n), h);   
			if(hindex<=12)
			{    						
				store_type dpdt;
				if(dFall_no >= 1) set(dpdt, hindex) = -L * ((dFall_no+1)*dFdp[hindex] - dFall);
				else set(dpdt, hindex) = -L * ((dFall_no)*dFdp[hindex] - dFall);
				phi_type value = grid(n)[hindex] + dt * dpdt[hindex];	
				set(update(n), hindex) = value;	
			}
		}			
		//----------------------------------------------------------------------------------------------------//
									
	}
	
	
		
		swap(grid, update);


	
	
	

}

return 0 ;

}
