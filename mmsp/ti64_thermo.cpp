/* Coupling of thermodynamic databases to the phase field model. 
 * 
 * Dimensions: 2
 * Thermodynamic model for description of solid solution: Regular solution model
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
const int dim = 2;
const int nx = 50;
const int ny = 50 ;
double G_normalize = R*T ;
double Dalal = 0.0001 ; 
double Dalv = 0.00008 ;
double Dvv = 0.0001 ;
double Dval = 0.0008 ; 
int steps = 70000 ;
double W_prefac = 6.095 ;

//--------------------------------Auxillary files containing helper functions--------------------------------//
#include "GradECoeff.hpp"
#include "Thermo.hpp"
//-----------------------------------------------------------------------------------------------------------//


//----------Container variables required for computing evolution of phase field and composition fields-------//
void thermo_auxillary_terms(MMSP::vector<store_type> gradient, MMSP::vector<store_type> gradientsq, double c_Al, double c_V) ;
double dGAlpha_dAl,	dGBeta_dAl, dGAlpha_dV, dGBeta_dV, del_dGAlpha_dAl, del_dGBeta_dAl, del_dGAlpha_dV, del_dGBeta_dV, delsq_dGAlpha_dAl, delsq_dGAlpha_dV, delsq_dGBeta_dAl, delsq_dGBeta_dV ;
double G_Al_alpha, G_Ti_alpha, G_V_alpha, G_Al_beta, G_Ti_beta, G_V_beta , W_Al, W_V, W_Ti ; 
double c11, c22, c33, c12, c13, c23 ;
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
	
	W_Al = 1.0 ;  
	W_V =  1.0 ; 
	W_Ti = 1.5 ; 
	
	G_Al_alpha = (Al_alpha_a_3 + Al_alpha_b_3*T + Al_alpha_c_3*T*log(T) + Al_alpha_d_3*T*T) ;
	G_Al_beta = (Al_beta_a_3 + Al_beta_b_3*T + Al_beta_c_3*T*log(T) + Al_beta_d_3*T*T) ;
	
	G_Ti_alpha = (Ti_alpha_a_2 + Ti_alpha_b_2*T + Ti_alpha_c_2*T*log(T) + Ti_alpha_d_2*T*T) ;
	G_Ti_beta = (Ti_beta_a_1 + Ti_beta_b_1*T + Ti_beta_c_1*T*log(T) + Ti_beta_d_1*T*T) ;
	
	G_V_alpha = (V_alpha_a_2 + V_alpha_b_2*T + V_alpha_c_2*T*log(T) + V_alpha_d_2*T*T) ; 
	G_V_beta = (V_beta_a_2 + V_beta_b_2*T + V_beta_c_2*T*log(T) + V_beta_d_2*T*T) ;

	for(int n = 0; n < nodes(grid); n++)
	{
		MMSP::vector<int> x = position(grid, n);	
		store_type newGrain;
		
		if((pow((x[0]*MMSP::dx(grid,0)-0.5*Lx),2) + pow((x[1]*MMSP::dx(grid,1)-0.3*Ly),2)) < 5*MMSP::dx(grid,0)*MMSP::dx(grid,0))
		{
			set(newGrain, 1) = 1.0;
			set(newGrain, 20) = 0.06;
			set(newGrain, 21) = 0.04 ;
		}
				
		else
		{
			set(newGrain, 1) = 0.0;
			set(newGrain, 20) = 0.06 ;
			set(newGrain, 21) = 0.04 ;
		}
		
		grid(n) = newGrain;
		
		//-------------------------Calculation of free energy of alpha phase as a function of c_al(n) and c_v(n)--------//
		G_Alpha[n] = ((grid(n)[20]*G_Al_alpha + grid(n)[21]*G_V_alpha + (1-grid(n)[20]-grid(n)[21])*G_Ti_alpha +
					 R*T*(grid(n)[20]*log(grid(n)[20]) + grid(n)[21]*log(grid(n)[21]) + (1-grid(n)[20]-grid(n)[21])*log((1-grid(n)[20]-grid(n)[21]))) +
					 grid(n)[20]*grid(n)[21]*L0_HCP_Al_V + grid(n)[20]*(1-grid(n)[20]-grid(n)[21])*L0_HCP_Al_Ti + grid(n)[21]*(1-grid(n)[20]-grid(n)[21])*L0_HCP_Ti_V))/G_normalize ;
					 
		G_Beta[n] = (grid(n)[20]*G_Al_beta + grid(n)[21]*G_V_beta + (1-grid(n)[20]-grid(n)[21])*G_Ti_beta +
					 R*T*(grid(n)[20]*log(grid(n)[20]) + grid(n)[21]*log(grid(n)[21]) + (1-grid(n)[20]-grid(n)[21])*log((1-grid(n)[20]-grid(n)[21]))) +
					 grid(n)[20]*grid(n)[21]*L0_BCC_Al_V + grid(n)[20]*(1-grid(n)[20]-grid(n)[21])*L0_BCC_Al_Ti + grid(n)[21]*(1-grid(n)[20]-grid(n)[21])*L0_BCC_Ti_V +
					 grid(n)[20]*grid(n)[21]*(1-grid(n)[20]-grid(n)[21])*L0_BCC_Al_Ti_V)/G_normalize ;
		//--------------------------------------------------------------------------------------------------------------//
					 
		
		W[n] = W_prefac*(grid(n)[20]*W_Al + grid(n)[21]*W_V + (1-grid(n)[20]-grid(n)[21])*W_Ti) ;
		
	
	}
	
}

//------------------------------------------------------------------------------------------------------------------//

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
	
	//-------------------------------------------------File I/O------------------------------------------------//
	if(t%1000==0)
	{
		std::ostringstream fileNameStream("pf_") ;   
		fileNameStream<<"pf_"<<t<<".txt"; //Writing phi_sum 
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
			
			myfile<<s[0]<<","<<s[1]<<","<<sum<<"\n" ;    
			
			
		}
	
		myfile.close();
		
		std::ostringstream fileNameStream1 ;
		fileNameStream1<<"c_Al_"<<t<<".txt";  //Writing c_Al 
		file_name1 = fileNameStream1.str() ;
		ofstream myfile1;
		myfile1.open(file_name1.c_str());
		std::ostringstream fileNameStream2("c_V_") ;   //Writing c_V
		fileNameStream2<<"c_V_"<<t<<".txt";
		file_name2 = fileNameStream2.str() ;
		ofstream myfile2;
		myfile2.open(file_name2.c_str());
		for(int n = 0 ; n < nodes(grid) ; n++)    
		{
			MMSP::vector<int> s = position(grid,n) ;
			myfile1<<s[0]<<","<<s[1]<<","<<grid(s)[20]<<"\n" ;
			myfile2<<s[0]<<","<<s[1]<<","<<grid(s)[21]<<"\n" ;
		}
		
		std::ostringstream fileNameStream3 ;   
		fileNameStream3<<"deltachem_"<<t<<".txt";  //Writing chemical free energy driving force per volume for nucleation
		file_name3 = fileNameStream3.str() ;
		ofstream myfile3;
		myfile3.open(file_name3.c_str());
		
		
		//-----------Calculation of chemical free energy driving force per volume for nucleation------------------//
		
		for(int n = 0 ; n < nodes(grid) ; n++)
		{
			MMSP::vector<int> s = position(grid,n) ;
			double g_alpha = ((grid(n)[20]*G_Al_alpha + grid(n)[21]*G_V_alpha + (1-grid(n)[20]-grid(n)[21])*G_Ti_alpha +
					 R*T*(grid(n)[20]*log(grid(n)[20]) + grid(n)[21]*log(grid(n)[21]) + (1-grid(n)[20]-grid(n)[21])*log((1-grid(n)[20]-grid(n)[21]))) +
					 grid(n)[20]*grid(n)[21]*L0_HCP_Al_V + grid(n)[20]*(1-grid(n)[20]-grid(n)[21])*L0_HCP_Al_Ti + grid(n)[21]*(1-grid(n)[20]-grid(n)[21])*L0_HCP_Ti_V))/G_normalize ;
					 
			double g_beta = (grid(n)[20]*G_Al_beta + grid(n)[21]*G_V_beta + (1-grid(n)[20]-grid(n)[21])*G_Ti_beta +
					 R*T*(grid(n)[20]*log(grid(n)[20]) + grid(n)[21]*log(grid(n)[21]) + (1-grid(n)[20]-grid(n)[21])*log((1-grid(n)[20]-grid(n)[21]))) +
					 grid(n)[20]*grid(n)[21]*L0_BCC_Al_V + grid(n)[20]*(1-grid(n)[20]-grid(n)[21])*L0_BCC_Al_Ti + grid(n)[21]*(1-grid(n)[20]-grid(n)[21])*L0_BCC_Ti_V +
					 grid(n)[20]*grid(n)[21]*(1-grid(n)[20]-grid(n)[21])*L0_BCC_Al_Ti_V)/G_normalize ;
			myfile3<<s[0]<<","<<s[1]<<","<<(g_beta-g_alpha)<<"\n" ;
		}
		
		myfile1.close();
		myfile2.close();
		myfile3.close();
	}
	//-------------------------------------------------------------------------------------------------------//
	
	
	for(int n = 0 ; n < nodes(grid) ; n++)
	{
		MMSP::vector<int> x = position(grid, n); 
		
		//-------------------------Calculation of free energy of alpha phase as a function of c_al(n) and c_v(n)--------//
		G_Alpha[n] = ((grid(n)[20]*G_Al_alpha + grid(n)[21]*G_V_alpha + (1-grid(n)[20]-grid(n)[21])*G_Ti_alpha +
					 R*T*(grid(n)[20]*log(grid(n)[20]) + grid(n)[21]*log(grid(n)[21]) + (1-grid(n)[20]-grid(n)[21])*log((1-grid(n)[20]-grid(n)[21]))) +
					 grid(n)[20]*grid(n)[21]*L0_HCP_Al_V + grid(n)[20]*(1-grid(n)[20]-grid(n)[21])*L0_HCP_Al_Ti + grid(n)[21]*(1-grid(n)[20]-grid(n)[21])*L0_HCP_Ti_V))/G_normalize ;
					 
		G_Beta[n] = (grid(n)[20]*G_Al_beta + grid(n)[21]*G_V_beta + (1-grid(n)[20]-grid(n)[21])*G_Ti_beta +
					 R*T*(grid(n)[20]*log(grid(n)[20]) + grid(n)[21]*log(grid(n)[21]) + (1-grid(n)[20]-grid(n)[21])*log((1-grid(n)[20]-grid(n)[21]))) +
					 grid(n)[20]*grid(n)[21]*L0_BCC_Al_V + grid(n)[20]*(1-grid(n)[20]-grid(n)[21])*L0_BCC_Al_Ti + grid(n)[21]*(1-grid(n)[20]-grid(n)[21])*L0_BCC_Ti_V +
					 grid(n)[20]*grid(n)[21]*(1-grid(n)[20]-grid(n)[21])*L0_BCC_Al_Ti_V)/G_normalize ;
		//--------------------------------------------------------------------------------------------------------------//
		
		
		MMSP::vector<store_type> gradient = grad(grid, x) ;
		MMSP::vector<store_type> gradientsq = gradsq(grid, x) ;
		
		
		thermo_auxillary_terms(gradient, gradientsq, grid(n)[20], grid(n)[21]) ; // helper function for sequentially solving the fourth order cahn-hillard equation. Please refer to the file thermo.hpp for more details about this function. 

		
		double phi_grad = 0 ;
		double phi_gradientsq = 0 ;
		for(int h = 0 ; h < length(grid(x)) ; h++)
		{
			int hindex = MMSP::index(grid(x),h) ;
			if(hindex<=12)  //updating hphi, hphi', hphi'', gphi, gphi', gphi'' only if the index corresponds to a variant. 
			{
				for(int i = 0 ; i < dim ; i++)
				{
					phi_grad+= gradient[i][hindex] ; 
					phi_gradientsq+= gradient[i][hindex] ; 
				}
				hphidoubleprime[n] = 60*grid(n)[hindex]*(2*grid(n)[hindex]-1)*(grid(n)[hindex]-1) ; 
				hphi[n] = pow(grid(n)[hindex],3)*(6*pow(grid(n)[hindex],2) - 15*grid(n)[hindex] + 10) ;
				hphiprime[n] = 30*pow(grid(n)[hindex],2)*pow((grid(n)[hindex]-1),2) ;
		
				gphidoubleprime[n] = 2*(6*pow(grid(n)[hindex],2) - 6*grid(n)[hindex] + 1) ; 
				gphiprime[n] = 2*grid(n)[hindex]*(1-grid(n)[hindex])*(1-2*grid(n)[hindex]) ;
				gphi[n] = pow(grid(n)[hindex],2)*pow((1-grid(n)[hindex]),2) ;
			}
			
		}

												 
		double c_al_rhs = 2*(del_dGAlpha_dAl - del_dGBeta_dAl)*hphiprime[n]*phi_grad + delsq_dGAlpha_dAl*hphi[n] + delsq_dGBeta_dAl*(1-hphi[n]) + 
						  (dGAlpha_dAl - dGBeta_dAl)*(hphidoubleprime[n]*pow(phi_grad,2) + hphiprime[n]*phi_gradientsq) + (W_Al - W_Ti)*(gphiprime[n]*phi_gradientsq + gphidoubleprime[n]*phi_grad) ;
		double c_v_rhs = 2*(del_dGAlpha_dV - del_dGBeta_dV)*hphiprime[n]*phi_grad + delsq_dGAlpha_dV*hphi[n] + delsq_dGBeta_dV*(1-hphi[n]) + 
						  (dGAlpha_dV - dGBeta_dV)*(hphidoubleprime[n]*pow(phi_grad,2) + hphiprime[n]*phi_gradientsq) + (W_V - W_Ti)*(gphiprime[n]*phi_gradientsq + gphidoubleprime[n]*phi_grad)  ;

		//--------------------------Updating c_al and c_v-----------------------------------//
		set(update(n), 20) = grid(n)[20] + dt*(Dalal*c_al_rhs + Dalv*c_v_rhs - laplacian(grid,n)[20]) ;
		set(update(n), 21) = grid(n)[21] + dt*(Dvv*c_v_rhs + Dval*c_al_rhs -  laplacian(grid,n)[21]) ;
		//----------------------------------------------------------------------------------//
											
				
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


//Write results to file 






