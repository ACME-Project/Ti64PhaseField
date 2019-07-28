/* Introduction of strain energy into the free energy functional of the phase field model. 
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
#include "Array.h"
#include "/home/arun/Documents/mmsp-arun/fftw++-2.05/fftw++.h"
#include "Eigenstrains.hpp"
#include "strain_modules.hpp"
//-----------------------------------------------------------------------------------------------------------//


using namespace utils;
using namespace Array;
using namespace fftwpp;


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
double c_norm = 1;
double c[6][6], sigma00[12][6] ;
//------------------------------------------------------------------------------------------------------------//


MMSP::grid<3, store_type> grid(fields, 0, nx, 0, ny, 0, nz) ;



//-----------------------------------Field initialization-----------------------------------------------------//
void initialize()
{

	for(int n = 0; n < nodes(grid); n++)
	{
		MMSP::vector<int> x = position(grid, n);
		
		store_type newGrain;
		
		if((pow((x[0]*MMSP::dx(grid,0)-0.5*Lx),2) + pow((x[1]*MMSP::dx(grid,1)-0.3*Ly),2)) < 5*MMSP::dx(grid,0)*MMSP::dx(grid,0))
		{
			set(newGrain, 1) = 1.0;
		}
		
		else
		{
			set(newGrain, 1) = 0.0;
		}
		
		grid(n) = newGrain;	
		
	}
	
}

//------------------------------------------------------------------------------------------------------------//

 

int main()
{
	
initialize_alpha_eigen() ; //Refer Eigenstrains.hpp
initialize();
initialize_epsi();  //Refer GradECoeff.hpp
define_c_sigma() ;  //Refer strain_modules.hpp	

size_t align=sizeof(Complex);
	
	
//-------------------------Define fourier and inverse-fourier transforms--------------------------//
rcfft3d Forward(nx,ny,nz,f,F);    
crfft3d Backward(nx,ny,nz,HetStr,HetStr_real);
//------------------------------------------------------------------------------------------------//
	

/* array3 is a container class defined in the fftw++ library
 * I have not yet checked if the Fourier transforms works with double[][][] array 
 */
 
//---Defining two containers for each variant, for performing FFT. The 'double' stores values in real space. The 'complex' stores the values after FFT has been performed---//
array3<double> f1(nx,ny,nz,align);
array3<Complex> F1(nx,ny,nz,align);
array3<double> f2(nx,ny,nz,align);
array3<Complex> F2(nx,ny,nz,align);
array3<double> f3(nx,ny,nz,align);
array3<Complex> F3(nx,ny,nz,align);
array3<double> f4(nx,ny,nz,align);
array3<Complex> F4(nx,ny,nz,align);
array3<double> f5(nx,ny,nz,align);
array3<Complex> F5(nx,ny,nz,align);
array3<double> f6(nx,ny,nz,align);
array3<Complex> F6(nx,ny,nz,align);
array3<double> f7(nx,ny,nz,align);
array3<Complex> F7(nx,ny,nz,align);
array3<double> f8(nx,ny,nz,align);
array3<Complex> F8(nx,ny,nz,align);
array3<double> f9(nx,ny,nz,align);
array3<Complex> F9(nx,ny,nz,align);
array3<double> f10(nx,ny,nz,align);
array3<Complex> F10(nx,ny,nz,align);
array3<double> f11(nx,ny,nz,align);
array3<Complex> F11(nx,ny,nz,align);
array3<double> f12(nx,ny,nz,align);
array3<Complex> F12(nx,ny,nz,align);
//----------------------------------------------------------------------------------------------------------------------------//


//---Defining two containers for each variant, for performing inverse-FT. The 'double' stores values in real space after the inverse transform. The 'complex' stores the values in the Fourier space---//
array3<Complex> HetStr1(nx,ny,nz,align);
array3<double> HetStr_real1(nx,ny,nz,align);
array3<Complex> HetStr2(nx,ny,nz,align);
array3<double> HetStr_real2(nx,ny,nz,align);
array3<Complex> HetStr3(nx,ny,nz,align);
array3<double> HetStr_real3(nx,ny,nz,align);
array3<Complex> HetStr4(nx,ny,nz,align);
array3<double> HetStr_real4(nx,ny,nz,align);
array3<Complex> HetStr5(nx,ny,nz,align);
array3<double> HetStr_real5(nx,ny,nz,align);
array3<Complex> HetStr6(nx,ny,nz,align);
array3<double> HetStr_real6(nx,ny,nz,align);
array3<Complex> HetStr7(nx,ny,nz,align);
array3<double> HetStr_real7(nx,ny,nz,align);
array3<Complex> HetStr8(nx,ny,nz,align);
array3<double> HetStr_real8(nx,ny,nz,align);
array3<Complex> HetStr9(nx,ny,nz,align);
array3<double> HetStr_real9(nx,ny,nz,align);
array3<Complex> HetStr10(nx,ny,nz,align);
array3<double> HetStr_real10(nx,ny,nz,align);
array3<Complex> HetStr11(nx,ny,nz,align);
array3<double> HetStr_real11(nx,ny,nz,align);
array3<Complex> HetStr12(nx,ny,nz,align);
array3<double> HetStr_real12(nx,ny,nz,align);
//----------------------------------------------------------------------------------------------------------------------------//

for(int i = 0 ; i < nx ; i++)
{
	for(int j = 0; j < ny ; j++)
	{
		for(int k = 0 ; k < nz ; k++)
		{
			f1(i,j,k) = 0.0 ;
			f2(i,j,k) = 0.0 ;
			f3(i,j,k) = 0.0 ;
			f4(i,j,k) = 0.0 ;
			f5(i,j,k) = 0.0 ;
			f6(i,j,k) = 0.0 ;
			f7(i,j,k) = 0.0 ;
			f8(i,j,k) = 0.0 ;
			f9(i,j,k) = 0.0 ;
			f10(i,j,k) = 0.0 ;
			f11(i,j,k) = 0.0 ;
			f12(i,j,k) = 0.0 ;
		}
	}
}
	
string file_name ;
string file_name4 ;


//------------------------------------Initializing the grid in the Fourier space---------------------------//
double fs = nx/Lx ;
double k[nx], fk[nx] ;	
for(int i = 0 ; i < nx ; i++)
{
	k[i] = i ;
	fk[i] = k[i]*fs/nx + 0.001  ;
}
//---------------------------------------------------------------------------------------------------------//

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
	
	double strain_energy[nx*ny*nz];
	MMSP::grid<dim, store_type > update(grid);
	
	MMSP::b0(update,0) = MMSP::Dirichlet;
	MMSP::b1(update,0) = MMSP::Dirichlet;

	MMSP::b0(update,1) = MMSP::Dirichlet;
	MMSP::b1(update,1) = MMSP::Dirichlet;

	MMSP::b0(update,2) = MMSP::Dirichlet;
	MMSP::b1(update,2) = MMSP::Dirichlet;
	
	//-------------------------------------------------File I/O------------------------------------------------//
	if(t%1000==0)
	{
		std::ostringstream fileNameStream("pf_") ;  //Writing phi_sum values
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
				sum+=grid(s)[hindex]*grid(s)[hindex] ;
			}
			
			myfile<<s[0]<<","<<s[1]<<","<<s[2]<<","<<sum<<"\n" ;    
			
			
		}
	
		myfile.close();
		
		//------------Writing interaction strain energy for each of the 12 variants at node n. This energy is calculated with respect to every variant present in the grid-----//
													//-----------This module lasts from the current line + 129 lines--------------//
		for(int n = 0 ; n < nodes(grid) ; n++)
		{
			
			
			MMSP::vector<int> s = position(grid, n) ;
			for(int h = 0 ; h < length(grid(n)) ; h++)
			{
				int hindex = MMSP::index(grid(n), h) ;
				if(hindex <= 12)
				{
					if(hindex==1) f1(s[0], s[1], s[2])=grid(n)[hindex] ;
					else if(hindex==2) f2(s[0], s[1], s[2])=grid(n)[hindex] ;
					else if(hindex==3) f3(s[0], s[1], s[2])=grid(n)[hindex] ;
					else if(hindex==4) f4(s[0], s[1], s[2])=grid(n)[hindex] ;
					else if(hindex==5) f5(s[0], s[1], s[2])=grid(n)[hindex] ;
					else if(hindex==6) f6(s[0], s[1], s[2])=grid(n)[hindex] ;
					else if(hindex==7) f7(s[0], s[1], s[2])=grid(n)[hindex] ;
					else if(hindex==8) f8(s[0], s[1], s[2])=grid(n)[hindex] ;
					else if(hindex==9) f9(s[0], s[1], s[2])=grid(n)[hindex] ;
					else if(hindex==10) f10(s[0], s[1], s[2])=grid(n)[hindex] ;
					else if(hindex==11) f11(s[0], s[1], s[2])=grid(n)[hindex] ;
					else if(hindex==12) f12(s[0], s[1], s[2])=grid(n)[hindex] ;
					
				}
			}
		}
			
		Forward.fft0(f1,F1);
		Forward.fft0(f2,F2);
		Forward.fft0(f3,F3);
		Forward.fft0(f4,F4);
		Forward.fft0(f5,F5);
		Forward.fft0(f6,F6);
		Forward.fft0(f7,F7);
		Forward.fft0(f8,F8);
		Forward.fft0(f9,F9);
		Forward.fft0(f10,F10);
		Forward.fft0(f11,F11);
		Forward.fft0(f12,F12); 
		
		for(int n = 0 ; n < nodes(grid) ; n++)
		{
			MMSP::vector<int> s = position(grid, n) ;
			double k1 = fk[s[0]]*0.1 ;
			double k2 = fk[s[1]]*0.1 ;
			double k3 = fk[s[2]]*0.1 ;
			double modk = (sqrt(k1*k1 + k2*k2 + k3*k3)) ;
			if(modk==0.0) modk=0.1   ;   //Remember that variant number ranges from 1-12 everywhere in the code, whereas the index ranges from 0-11 everywhere in the code. 
			
			HetStr1(s[0], s[1], s[2]) = Bpq(k1, k2, k3, modk, 1, 1)*F1(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 1, 2)*F2(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 1, 3)*F3(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 1, 4)*F4(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 1, 5)*F5(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 1, 6)*F6(s[0], s[1], s[2]) +
			Bpq(k1, k2, k3, modk, 1, 7)*F7(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 1, 8)*F8(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 1, 9)*F9(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 1, 10)*F10(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 1, 11)*F11(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 1, 12)*F12(s[0], s[1], s[2]) ;
			
			HetStr2(s[0], s[1], s[2]) = Bpq(k1, k2, k3, modk, 2, 1)*F1(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 2, 2)*F2(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 2, 3)*F3(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 2, 4)*F4(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 2, 5)*F5(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 2, 6)*F6(s[0], s[1], s[2]) +
			Bpq(k1, k2, k3, modk, 2, 7)*F7(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 2, 8)*F8(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 2, 9)*F9(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 2, 10)*F10(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 2, 11)*F11(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 2, 12)*F12(s[0], s[1], s[2]) ;
			
			HetStr3(s[0], s[1], s[2]) = Bpq(k1, k2, k3, modk, 3, 1)*F1(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 3, 2)*F2(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 3, 3)*F3(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 3, 4)*F4(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 3, 5)*F5(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 3, 6)*F6(s[0], s[1], s[2]) +
			Bpq(k1, k2, k3, modk, 3, 7)*F7(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 3, 8)*F8(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 3, 9)*F9(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 3, 10)*F10(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 3, 11)*F11(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 3, 12)*F12(s[0], s[1], s[2]) ;
			
			HetStr4(s[0], s[1], s[2]) = Bpq(k1, k2, k3, modk, 4, 1)*F1(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 4, 2)*F2(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 4, 3)*F3(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 4, 4)*F4(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 4, 5)*F5(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 4, 6)*F6(s[0], s[1], s[2]) +
			Bpq(k1, k2, k3, modk, 4, 7)*F7(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 4, 8)*F8(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 4, 9)*F9(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 4, 10)*F10(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 4, 11)*F11(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 4, 12)*F12(s[0], s[1], s[2]) ;
			
			HetStr5(s[0], s[1], s[2]) = Bpq(k1, k2, k3, modk, 5, 1)*F1(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 5, 2)*F2(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 5, 3)*F3(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 5, 4)*F4(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 5, 5)*F5(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 5, 6)*F6(s[0], s[1], s[2]) +
			Bpq(k1, k2, k3, modk, 5, 7)*F7(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 5, 8)*F8(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 5, 9)*F9(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 5, 10)*F10(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 5, 11)*F11(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 5, 12)*F12(s[0], s[1], s[2]) ;
			
			HetStr6(s[0], s[1], s[2]) = Bpq(k1, k2, k3, modk, 6, 1)*F1(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 6, 2)*F2(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 6, 3)*F3(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 6, 4)*F4(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 6, 5)*F5(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 6, 6)*F6(s[0], s[1], s[2]) +
			Bpq(k1, k2, k3, modk, 6, 7)*F7(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 6, 8)*F8(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 6, 9)*F9(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 6, 10)*F10(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 6, 11)*F11(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 6, 12)*F12(s[0], s[1], s[2]) ;
			
			HetStr7(s[0], s[1], s[2]) = Bpq(k1, k2, k3, modk, 7, 1)*F1(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 7, 2)*F2(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 7, 3)*F3(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 7, 4)*F4(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 7, 5)*F5(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 7, 6)*F6(s[0], s[1], s[2]) +
			Bpq(k1, k2, k3, modk, 7, 7)*F7(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 7, 8)*F8(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 7, 9)*F9(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 7, 10)*F10(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 7, 11)*F11(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 7, 12)*F12(s[0], s[1], s[2]) ;
			
			HetStr8(s[0], s[1], s[2]) = Bpq(k1, k2, k3, modk, 8, 1)*F1(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 8, 2)*F2(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 8, 3)*F3(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 8, 4)*F4(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 8, 5)*F5(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 8, 6)*F6(s[0], s[1], s[2]) +
			Bpq(k1, k2, k3, modk, 8, 7)*F7(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 8, 8)*F8(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 8, 9)*F9(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 8, 10)*F10(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 8, 11)*F11(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 8, 12)*F12(s[0], s[1], s[2]) ;
			
			HetStr9(s[0], s[1], s[2]) = Bpq(k1, k2, k3, modk, 9, 1)*F1(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 9, 2)*F2(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 9, 3)*F3(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 9, 4)*F4(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 9, 5)*F5(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 9, 6)*F6(s[0], s[1], s[2]) +
			Bpq(k1, k2, k3, modk, 9, 7)*F7(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 9, 8)*F8(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 9, 9)*F9(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 9, 10)*F10(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 9, 11)*F11(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 9, 12)*F12(s[0], s[1], s[2]) ;
			
			HetStr10(s[0], s[1], s[2]) = Bpq(k1, k2, k3, modk, 10, 1)*F1(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 10, 2)*F2(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 10, 3)*F3(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 10, 4)*F4(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 10, 5)*F5(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 10, 6)*F6(s[0], s[1], s[2]) +
			Bpq(k1, k2, k3, modk, 10, 7)*F7(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 10, 8)*F8(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 10, 9)*F9(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 10, 10)*F10(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 10, 11)*F11(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 10, 12)*F12(s[0], s[1], s[2]) ;
			
			HetStr11(s[0], s[1], s[2]) = Bpq(k1, k2, k3, modk, 11, 1)*F1(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 11, 2)*F2(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 11, 3)*F3(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 11, 4)*F4(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 11, 5)*F5(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 11, 6)*F6(s[0], s[1], s[2]) +
			Bpq(k1, k2, k3, modk, 11, 7)*F7(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 11, 8)*F8(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 11, 9)*F9(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 11, 10)*F10(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 11, 11)*F11(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 11, 12)*F12(s[0], s[1], s[2]) ;
			
			HetStr12(s[0], s[1], s[2]) = Bpq(k1, k2, k3, modk, 12, 1)*F1(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 12, 2)*F2(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 12, 3)*F3(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 12, 4)*F4(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 12, 5)*F5(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 12, 6)*F6(s[0], s[1], s[2]) +
			Bpq(k1, k2, k3, modk, 12, 7)*F7(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 12, 8)*F8(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 12, 9)*F9(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 12, 10)*F10(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 12, 11)*F11(s[0], s[1], s[2]) + Bpq(k1, k2, k3, modk, 12, 12)*F12(s[0], s[1], s[2]) ;
					
				
		}
		
		
		Backward.fft0Normalized(HetStr1, HetStr_real1); 
		Backward.fft0Normalized(HetStr2, HetStr_real2); 
		Backward.fft0Normalized(HetStr3, HetStr_real3); 
		Backward.fft0Normalized(HetStr4, HetStr_real4); 
		Backward.fft0Normalized(HetStr5, HetStr_real5); 
		Backward.fft0Normalized(HetStr6, HetStr_real6); 
		Backward.fft0Normalized(HetStr7, HetStr_real7); 
		Backward.fft0Normalized(HetStr8, HetStr_real8); 
		Backward.fft0Normalized(HetStr9, HetStr_real9); 
		Backward.fft0Normalized(HetStr10, HetStr_real10); 
		Backward.fft0Normalized(HetStr11, HetStr_real11); 
		Backward.fft0Normalized(HetStr12, HetStr_real12); 
		
		std::ostringstream fileNameStream4 ;   //Check for a possible error here
		fileNameStream4<<"deltastrainmax_"<<t<<".txt";
		file_name4 = fileNameStream4.str() ;
		ofstream myfile4;
		myfile4.open(file_name4.c_str());
		
		for(int n = 0 ; n < nodes(grid) ; n++)
		 {
			MMSP::vector<int> s = position(grid, n) ; 
			double int_se_mean = (HetStr_real1(s[0], s[1], s[2]) + HetStr_real2(s[0], s[1], s[2]) + HetStr_real3(s[0], s[1], s[2]) + HetStr_real4(s[0], s[1], s[2]) + HetStr_real5(s[0], s[1], s[2]) + HetStr_real6(s[0], s[1], s[2]) + HetStr_real7(s[0], s[1], s[2]) + HetStr_real8(s[0], s[1], s[2]) + HetStr_real9(s[0], s[1], s[2]) + HetStr_real10(s[0], s[1], s[2]) + HetStr_real11(s[0], s[1], s[2]) + HetStr_real12(s[0], s[1], s[2]))/12 ;
			double int_se_std = (pow((HetStr_real1(s[0], s[1], s[2]) - int_se_mean),2) + pow((HetStr_real2(s[0], s[1], s[2]) - int_se_mean),2) + pow((HetStr_real3(s[0], s[1], s[2]) - int_se_mean),2) + pow((HetStr_real4(s[0], s[1], s[2]) - int_se_mean),2) + pow((HetStr_real5(s[0], s[1], s[2]) - int_se_mean),2) + pow((HetStr_real6(s[0], s[1], s[2]) - int_se_mean),2) + pow((HetStr_real7(s[0], s[1], s[2]) - int_se_mean),2) + pow((HetStr_real8(s[0], s[1], s[2]) - int_se_mean),2) + pow((HetStr_real9(s[0], s[1], s[2]) - int_se_mean),2) +
			pow((HetStr_real10(s[0], s[1], s[2]) - int_se_mean),2) + pow((HetStr_real11(s[0], s[1], s[2]) - int_se_mean),2) + pow((HetStr_real12(s[0], s[1], s[2]) - int_se_mean),2))/12 ;
			int_se_std = pow(int_se_std, 0.5) ;
			double int_se_min = HetStr_real1(s[0], s[1], s[2]) ;
			if(HetStr_real2(s[0], s[1], s[2]) < int_se_min) int_se_min = HetStr_real2(s[0], s[1], s[2]) ;  
			else if(HetStr_real3(s[0], s[1], s[2]) < int_se_min) int_se_min = HetStr_real3(s[0], s[1], s[2]) ; 
			else if(HetStr_real4(s[0], s[1], s[2]) < int_se_min) int_se_min = HetStr_real4(s[0], s[1], s[2]) ; 
			else if(HetStr_real5(s[0], s[1], s[2]) < int_se_min) int_se_min = HetStr_real5(s[0], s[1], s[2]) ; 
			else if(HetStr_real6(s[0], s[1], s[2]) < int_se_min) int_se_min = HetStr_real6(s[0], s[1], s[2]) ; 
			else if(HetStr_real7(s[0], s[1], s[2]) < int_se_min) int_se_min = HetStr_real7(s[0], s[1], s[2]) ; 
			else if(HetStr_real8(s[0], s[1], s[2]) < int_se_min) int_se_min = HetStr_real8(s[0], s[1], s[2]) ; 
			else if(HetStr_real9(s[0], s[1], s[2]) < int_se_min) int_se_min = HetStr_real9(s[0], s[1], s[2]) ; 
			else if(HetStr_real10(s[0], s[1], s[2]) < int_se_min) int_se_min = HetStr_real10(s[0], s[1], s[2]) ; 
			else if(HetStr_real11(s[0], s[1], s[2]) < int_se_min) int_se_min = HetStr_real11(s[0], s[1], s[2]) ; 
			else if(HetStr_real12(s[0], s[1], s[2]) < int_se_min) int_se_min = HetStr_real12(s[0], s[1], s[2]) ; 
			
			myfile4<<s[0]<<","<<s[1]<<","<<s[2]<<","<<int_se_mean<<','<<int_se_std<<','<<int_se_min<<"\n" ;
		
		}
		//----------------------------------------------------------------------------------------------------------------------------------------------------------//
			
		
		myfile4.close();
		
	}
	
	
	//-------------------------------Stress Free Transformation Strain(SFTS) energy-------------------------------------------//
	double sfts_energy =  0.0 ; 
	for(int i = 0 ; i < 6 ; i++)
	{
		for(int j = 0 ; j < 6 ; j++)
		{
			sfts_energy+= eigen_alpha[variant].e[i]*c[i][j]*eigen_alpha[variant].e[j] ;
		}
	}
	//------------------------------------------------------------------------------------------------------------------------//
	
	
	
	double phi_sum = 0 ;
	
	for(int n = 0 ; n < nodes(grid) ; n++)
	{
		MMSP::vector<int> s = position(grid, n) ;
		for(int h = 0 ; h < length(grid(n)) ; h++)
		{
			int hindex = MMSP::index(grid(n), h) ;
			phi_sum+=grid(n)[hindex] ; 
			f(s[0], s[1], s[2])=grid(n)[hindex] ;    //Copying phi values to f for FFT in the next step
		}	
		
	}
  
  
	Forward.fft0(f,F);    //FFT from f to F

  
  
	for(int n = 0 ; n < nodes(grid) ; n++)
	{
		MMSP::vector<int> x = position(grid, n) ;
		double k1 = fk[x[0]]*0.1 ;  //k1, k2, k3 are lattice points in the Fourier grid. 
		double k2 = fk[x[1]]*0.1 ;
		double k3 = fk[x[2]]*0.1 ;
		double modk = (sqrt(k1*k1 + k2*k2 + k3*k3)) ;
		if(modk==0.0) modk=0.1   ;
		HetStr(x[0],x[1],x[2]) = B(k1, k2, k3, modk, variant)*F(x[0],x[1],x[2]) ;	//Refer to strain_modules.hpp for the function B. HetStr is the heterogeneous strain component in the Fourier space.
	}

  
  
	Backward.fft0Normalized(HetStr, HetStr_real);  //Inverse Fourier transform from HetStr to HetStr_real
	
	
	for(int n = 0 ; n < nodes(grid) ; n++)
	{
		
		MMSP::vector<int> x = position(grid, n); 
		
		
		int n_index ;
		for(int h = 0 ; h < length(grid(n)) ; h++)
		{
			n_index = MMSP::index(grid(n), h) ;
			strain_energy[n] = grid(n)[n_index]*((1/(nx*ny*nz))*grid(n)[n_index]*sfts_energy - (1/(nx*ny*nz))*phi_sum*sfts_energy  - (1/(nx*ny*nz))*(HetStr_real(x[0],x[1],x[2]))) ;
			//Refer to the Strain_Energy_Doc.ipynb for description of the above equation.
	    }  
		
		MMSP::vector<store_type> gradient = grad(grid, x) ;
		MMSP::vector<store_type> gradientsq = gradsq(grid, x) ;


				
		MMSP::sparse<phi_type> dFdp;

		for (int h = 0; h < length(grid(n)); h++) 
		{	
						
			int hindex = MMSP::index(grid(n), h);   		
			phi_type phi_sq = 0.0; 							
			phi_type dFall = 0.0;
					
														
			for (int j = 0; j < length(grid(n)); j++) 
			{
				int jindex = MMSP::index(grid(n), j);
				if(jindex<=12){	phi_sq += grid(n)[jindex]*grid(n)[jindex];}
			}
					
														
			for (int j = 0; j < length(grid(n)); j++)
			{
				int jindex = MMSP::index(grid(n), j);  
				phi_type lap_aniso = 0 ;
				if(jindex <= 12)
				{
					for(int p = 0 ; p < dim ; p++)
					{
						for(int q = 0 ; q < dim ; q++)
						{
							lap_aniso+= gradient[p][jindex]*epsi[jindex-1][p][q]*gradient[q][jindex];    
						}
					}
				
				set(dFdp, jindex) = 1.0 * grid(n)[jindex] * (phi_sq - grid(n)[jindex]) + strain_energy[n] - (lap_aniso);
				dFall += dFdp[jindex];	
				}				
			}
					
			store_type dpdt;
			set(dpdt, hindex) = -L * ((variants+1)*dFdp[hindex] - dFall);
			phi_type value = grid(n)[hindex] + dt * dpdt[hindex];	
			set(update(n), hindex) = value;	
		}										
			
	}
		
		swap(grid, update);


	
	
	

}

return 0 ;

}

//Write results to file 






