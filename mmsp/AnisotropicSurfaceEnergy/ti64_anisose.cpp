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
#include<deque>
#include<vector>
#include"MMSP.hpp"

double epsi[12][3][3] ;
#include "GradECoeff.hpp"

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


using namespace std;

double R = 8.314 ;
double T = 1023 ;
double dt = 3.0e-1 ;
double L = 0.239 ;
double G_normalize = R*T ; 

typedef float phi_type; 											
typedef MMSP::sparse<phi_type> store_type; 	

int variants = 10 ;
const int dim = 3;
const int nx = 10;
const int ny = 10;
const int nz = 10;
const int node_total = nx*ny*nz ; 
double G_Alpha[node_total], G_Beta[node_total], W[node_total] ; 
double hphi[node_total], hphiprime[node_total], hphidoubleprime[node_total];
double gphi[node_total], gphiprime[node_total], gphidoubleprime[node_total];

MMSP::grid<dim, store_type> grid(variants, 0, nx, 0, ny, 0, nz) ;
MMSP::grid<dim, store_type> gradgrid(variants, 0, nx, 0, ny, 0, nz) ;

const double Lx = g1(grid,0) - g0(grid,0) ;
const double Ly = g1(grid,1) - g0(grid,1) ;
const double Lz = g1(grid,2) - g0(grid,2) ;

const double dx = MMSP::dx(grid, 0) ;


void initialize()
{

	for(int n = 0; n < nodes(grid); n++)
	{
		MMSP::vector<int> x = position(grid, n);
		
		store_type newGrain;
		
		if((pow((x[0]*MMSP::dx(grid,0)-0.5*Lx),2) + pow((x[1]*MMSP::dx(grid,1)-0.5*Ly),2) + pow((x[2]*MMSP::dx(grid,2)-0.5*Lz),2)) < 3*MMSP::dx(grid,0)*MMSP::dx(grid,0))
		{
			set(newGrain, 15) = 1.0;
			
		}
			
		
			
		else
		{
			set(newGrain, 15) = 1.0;	
		}
		
		grid(n) = newGrain;		
		
	}

}

int main()
{	
	string file_name ;
	string file_name1 ;
	string file_name2 ;

	int steps = 70000 ;  	

	initialize();
	initialize_epsi();

	MMSP::b0(grid,0) = MMSP::Dirichlet;
	MMSP::b1(grid,0) = MMSP::Dirichlet;

	MMSP::b0(grid,1) = MMSP::Dirichlet;
	MMSP::b1(grid,1) = MMSP::Dirichlet;

	MMSP::b0(grid,2) = MMSP::Dirichlet;
	MMSP::b1(grid,2) = MMSP::Dirichlet;

	double lap_max ;

	for(int t = 0 ; t<steps ; t++)
	{
		MMSP::grid<dim, store_type > update(grid);
	
		MMSP::b0(update,0) = MMSP::Dirichlet;
		MMSP::b1(update,0) = MMSP::Dirichlet;

		MMSP::b0(update,1) = MMSP::Dirichlet;
		MMSP::b1(update,1) = MMSP::Dirichlet;

		MMSP::b0(update,2) = MMSP::Dirichlet;
		MMSP::b1(update,2) = MMSP::Dirichlet;
	
		if(t%100==0)
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
	
		double deltag = 0.1 ;
		double kappa1 = 0.02 ; 
		double kappa2 = 0.1 ; 
	
		int no_of_cells = (nx*ny*nz)/(3*3*3) ; 
		int cells_in_one_dir = nx/3 ; 
	
		if((t%100==0)&(t<300))   
		{
			MMSP::vector<int> x(3) ; 
			MMSP::vector<int> y(3) ; 
			MMSP::vector<int> z(3) ; 
			MMSP::vector<int> temp_pos(3); 
			for(int i = 0 ; i < cells_in_one_dir ; i++)
			{
				for(int j = 0 ; j < cells_in_one_dir ; j++)
				{
					for(int k = 0 ; k < cells_in_one_dir ; k++)
					{
						int flag = 0 ;
				
						x[0] = i*3 ; x[1] = i*3+1 ; x[2] = i*3+2 ; 
				
						y[0] = j*3 ; y[1] = j*3+1 ; y[2] = j*3+2 ; 
				
						z[0] = k*3 ; z[1] = k*3+1 ; z[2] = k*3+2 ; 
						for(int p = 0; p < 3; p++)
						{
							for(int q = 0 ; q < 3 ; q++)
							{
								for(int r = 0 ; r < 3 ; r++)
								{							
									temp_pos[0] = x[p] ; temp_pos[1] = y[q] ; temp_pos[2] = z[r] ;
									flag+= grid(temp_pos)[15] ;
								}
							}
						}
				
						if(flag > 0)
						{					
							double jstar = kappa1*exp(-kappa2*deltag);  
							double prob = 1 - exp(-jstar) ;
							double random_no = (double)(rand() % 1000) ;
							double decision_var = random_no/1000.0 ;
							if(prob > decision_var)
							{
								store_type newGrain;
								int random_index = (int)(rand() % 12) + 1 ;
								cout<<"Nucleating index :"<<random_index<<endl;
								for(int m = 0 ; m < nodes(grid) ; m++)
								{
									if(length(grid(m))==1) set(grid(m), random_index) = 0.0 ;
									else if(length(grid(m))> 1)
									{
										int temp_check = 0 ; 
										for(int h = 0 ; h < length(grid(m)) ; h++)
										{
											if(MMSP::index(grid(m), h) == random_index)
											{
												temp_check = 1 ;
												break;
											}
										}
										if(temp_check == 0) set(grid(m), random_index) = 0.0; 
									}
								}
								for(int p = 0; p < 3; p++)
								{
									for(int q = 0 ; q < 3 ; q++)
									{
										for(int r = 0 ; r < 3 ; r++)
										{
											temp_pos[0] = x[p] ; temp_pos[1] = y[q] ; temp_pos[2] = z[r] ;
									
											if((pow((temp_pos[0]*MMSP::dx(grid,0)-x[2]*MMSP::dx(grid,0)),2) + pow((temp_pos[1]*MMSP::dx(grid,1)-y[2]*MMSP::dx(grid,0)),2) + pow((temp_pos[2]*MMSP::dx(grid,2)-z[2]*MMSP::dx(grid,0)),2)) < 3*MMSP::dx(grid,0)*MMSP::dx(grid,0))
											{	
												set(grid(temp_pos), 15) = 0.0 ;   
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
	
	
		
	
	
		for(int n = 0 ; n < nodes(grid) ; n++)
		{
			MMSP::vector<int> x = position(grid, n); 		
			MMSP::vector<store_type> gradient = grad(grid, x) ;
			MMSP::vector<store_type> gradientsq = gradsq(grid, x) ;
		
			for (int j = 0; j < length(grid(n)); j++) 
			{
				int jindex = MMSP::index(grid(n), j);
				if(jindex<=12)
				{
					set(gradgrid(n),(3*(jindex-1)+1)) = gradient[0][jindex];
					set(gradgrid(n),(3*(jindex-1)+2)) = gradient[1][jindex];
					set(gradgrid(n),(3*(jindex-1)+3)) = gradient[2][jindex];
				}
			}
		
		
			double phi_grad = 0 ;
			double phi_gradientsq = 0 ;
			double lap_aniso[13];
				
			MMSP::sparse<phi_type> dFdp;
			
			phi_type phi_sq = 0.0; 	
			phi_type phi_sum = 0.0 ;						
			phi_type dFall = 0.0;
			phi_type dFall_no = 0 ;
					
														
			for (int j = 0; j < length(grid(n)); j++) 
			{
				int jindex = MMSP::index(grid(n), j);
				if(jindex<=12)
				{	
					phi_sq += grid(n)[jindex]*grid(n)[jindex];
					phi_sum += grid(n)[jindex] ;
				}
			}
				
			phi_type phi_beta = 1 - phi_sum ; 
				
			MMSP::sparse<phi_type> lap = laplacian(grid, n); 
				
			int nz = 0;
			
			for (int j = 0; j < length(grid(n)); j++)   
			{
				int jindex = MMSP::index(grid(n), j);  		
					
				MMSP::vector<int> x = position(grid, n);
				MMSP::vector<double> gradients(3) ;
					
				W[n] = 0.01432;
					
				MMSP::vector<store_type> gradgrad = grad(gradgrid, x) ;
				if(jindex <= 12)
				{
					lap_aniso[jindex] = 0.0 ;
					for(int p = 0 ; p < dim ; p++)
					{
						for(int q = 0 ; q < dim ; q++)
						{
							lap_aniso[jindex]+= 0.25*gradgrad[p][3*(jindex-1)+q+1]*epsi[jindex-1][p][q] ;
						}
					}
					
					double interaction_energy = 0 ;
					interaction_energy+= W[n]*grid(n)[jindex]*pow(phi_beta, 2);
					interaction_energy-= W[n]*pow(grid(n)[jindex],2)*(phi_beta);
					
					int check = 0;
					for(int k = 0 ; k < length(grid(n)); k++)
					{
						int tempindex = MMSP::index(grid(n), k); 
						if((tempindex!=jindex)&&(tempindex <=12))
						{
							interaction_energy+= W[n]*pow(grid(n)[tempindex],2)*grid(n)[jindex];
							interaction_energy-= W[n]*pow(grid(n)[tempindex],2)*(phi_beta); 
						}
					}
					
					double coeff ;
					if(lap[jindex]==0) coeff = 0.0 ;
					else coeff = coeff = ((lap[jindex])/(abs(lap[jindex]))) ;
					set(dFdp, jindex) = (-2.2)*(pow(grid(n)[jindex],2) - pow(grid(n)[jindex],3)) + (-1.1)*(-pow(grid(n)[jindex],2) + pow(grid(n)[jindex],3))   + interaction_energy - coeff*lap_aniso[jindex] ; 
					
					dFall += dFdp[jindex];	
					dFall_no+=1;
				}
										
			}

			for (int h = 0; h < length(grid(n)); h++) 
			{	
					
				int hindex = MMSP::index(grid(n), h);   
				if(hindex<=12)
				{    				
					
				store_type dpdt;
				if(dFall_no > 1) set(dpdt, hindex) = -L * ((dFall_no+1)*dFdp[hindex] - dFall);
				else set(dpdt, hindex) = -L * ((dFall_no+1)*dFdp[hindex] - dFall);
				phi_type value = grid(n)[hindex] + dt * dpdt[hindex];	
				set(update(n), hindex) = value;	
				}
			}
			
			set(update(n), 15) = grid(n)[15] ;
			
		}

		swap(grid, update);	

}



return 0 ;

}







