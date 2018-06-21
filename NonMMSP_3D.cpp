// Symmetric Toth XMPF written with MMSP classes and no optimizations



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
#include<iostream>
#include<fstream>
#include<sstream>

double epsi1[3][3], epsi2[3][3], epsi3[3][3], epsi4[3][3], epsi5[3][3], epsi6[3][3], epsi7[3][3], epsi8[3][3], epsi9[3][3], epsi10[3][3], epsi11[3][3], epsi12[3][3];

#include "Tensor.hpp"
#include "Eigenstrains.hpp"
#include "GradECoeff.hpp"

const int dim = 3;
const int nx = 40;
const int ny = 40 ;
const int nz = 40;

const double dx = 1.0 ;
const double dy = 1.0 ;
const double dz = 1.0 ;

const double Lx = dx*nx ;
const double Ly = dy*ny ;
const double Lz = dz*nz ;

double phi[nx][ny][nz] ;
double phi_old[nx][ny][nz], dphidx[nx][ny][nz], dphidy[nx][ny][nz], dphidz[nx][ny][nz], dphidxsq[nx][ny][nz], dphidysq[nx][ny][nz], dphidzsq[nx][ny][nz], dphidxdy[nx][ny][nz], dphidzdx[nx][ny][nz], dphidydz[nx][ny][nz], dphidydx[nx][ny][nz] ;
double epsi[dim][dim] ;

using namespace std;

void initialize()
{
	for(int x = 0 ; x < nx ; x++)
	{
		for(int y = 0 ; y < ny ; y++)
		{
			for(int z = 0 ; z < nz ; z++)
			{
				
				if((pow((x*dx-0.5*Lx),2) + pow((y*dy-0.5*Ly),2) + pow((z*dz-0.5*Lz),2)) < 13*dx*dx)   //Initializes a spherical nucleus at the center of the domain
				{
					phi[x][y][z] = 1.0 ;
					phi_old[x][y][z] = 1.0 ;
				}
			
				else
				{
					phi[x][y][z] = 0.0 ;
					phi_old[x][y][z] = 0.0 ;
				}
			}
		}
	}
				
		
	
	initialize_epsi();
	
	for(int i = 0 ; i < 3 ; i++)
	{
		for(int j = 0 ; j < 3 ; j++)
		{
			epsi[i][j] = epsi1[i][j] ;
		}
	}
	
	for(int x = 1 ; x <nx-1 ; x++)
	{
		for(int y = 1 ; y < ny-1 ; y++)
		{
			for(int z = 1 ; z < nz-1 ; z++)
			{	
				dphidx[x][y][z] = (phi_old[x+1][y][z]-phi_old[x-1][y][z])/(2*dx) ;
				dphidy[x][y][z] = (phi_old[x][y+1][z]-phi_old[x][y-1][z])/(2*dy) ;	
				dphidz[x][y][z] = (phi_old[x][y][z+1]-phi_old[x][y][z-1])/(2*dz) ;
			}
		}
	}

	for(int x = 2 ; x <nx-2 ; x++)
	{
		for(int y = 2 ; y < ny-2 ; y++)
		{	
			for(int z = 1 ; z < nz-1 ; z++)
			{
				dphidxdy[x][y][z] = (dphidx[x][y+1][z] - dphidx[x][y-1][z])/(2*dy) ;
				dphidxsq[x][y][z] = (phi_old[x+1][y][z] - 2*phi_old[x][y][z] + phi_old[x-1][y][z])/(dx*dx);
				dphidysq[x][y][z] = (phi_old[x][y+1][z] - 2*phi_old[x][y][z] + phi_old[x][y-1][z])/(dy*dy);
			}
		}
	}
	
	
}


int main()
{
	
	string file_name ;
	int steps = 1000000 ;
	double dt = dx/200.0 ;
	double L = 0.1 ;
	double Cijkl[6][6] ;

	for(int i = 0 ; i < 6 ; i++)
	{
		for(int j = 0 ; j < 6 ; j++)
		{
			Cijkl[i][j] = 0 ;
		}
	}
	Cijkl[0][0] = 119.0e-05 ;
	Cijkl[0][1] = 64.0e-05 ;
	Cijkl[0][2] = 36.0e-05 ;
	Cijkl[2][2] = 0.0e-05 ;
	Cijkl[3][3] = 134.0e-05 ;
	Cijkl[1][2] = 49.0e-05 ;
	

	initialize();

	eigen_alpha[0].matrix[0][0] = -0.083  ;
	eigen_alpha[0].matrix[0][1] = 0.0095  ;
	eigen_alpha[0].matrix[0][2] = 0.0  ;
	eigen_alpha[0].matrix[1][0] = 0.0095  ;
	eigen_alpha[0].matrix[1][1] = 0.123  ;
	eigen_alpha[0].matrix[1][2] = 0.0  ;
	eigen_alpha[0].matrix[2][0] = 0.0  ;
	eigen_alpha[0].matrix[2][1] = 0.0  ;
	eigen_alpha[0].matrix[2][2] = 0.035  ;


for(int t = 0 ; t<steps ; t++)
{

	if(t%5000==0)
	{
		std::ostringstream fileNameStream("pf_") ;
		fileNameStream<<"pf_"<<t<<".txt";
		file_name = fileNameStream.str() ;
		ofstream myfile;
		myfile.open(file_name.c_str());
		for(int x = 0 ; x < nx ; x++)
		{
			for(int y = 0 ; y < ny ; y++)
			{
				for(int z = 0 ; z < nz ; z++)
				{
					myfile<<x<<","<<y<<","<<z<<","<<phi_old[x][y][z]<<"\n" ;
				}
			}
			
			//myfile<<"\n";
		}
	
		myfile.close();
	}
	
	
	for(int x = 1 ; x <nx-1 ; x++)
	{
		for(int y = 1 ; y < ny-1 ; y++)
		{
			for(int z = 1 ; z < nz-1 ; z++)
			{
				dphidx[x][y][z] = (phi_old[x+1][y][z]-phi_old[x-1][y][z])/(2*dx) ;
				dphidy[x][y][z] = (phi_old[x][y+1][z]-phi_old[x][y-1][z])/(2*dy) ;	
				dphidz[x][y][z] = (phi_old[x][y][z+1]-phi_old[x][y][z-1])/(2*dz) ;	
			}
		}
	}

	for(int x = 2 ; x <nx-2 ; x++)
	{
		for(int y = 2 ; y < ny-2 ; y++)
		{	
			for(int z = 2 ; z < nz-2 ; z++)
			{
				//dphidxdy[x][y][z] = 0.5*(phi_old[x+1][y+1][z] - phi_old[x+1][y-1][z] - phi_old[x-1][y+1][z] + phi_old[x-1][y-1][z])/(4*dy*dx)  ;
				dphidxdy[x][y][z] = 0.5*(dphidx[x][y+1][z] - dphidx[x][y-1][z])/(2*dy)  ;
				dphidydx[x][y][z] = 0.5*(dphidy[x+1][y][z] - dphidy[x-1][y][z])/(2*dx)  ;
				dphidzdx[x][y][z] = 0.5*(dphidz[x+1][y][z] - dphidz[x-1][y][z])/(2*dx)  ;
				dphidydz[x][y][z] = 0.5*(dphidy[x][y][z+1] - dphidy[x][y][z-1])/(2*dz)  ;
				dphidxsq[x][y][z] = 0.5*(phi_old[x+1][y][z] - 2*phi_old[x][y][z] + phi_old[x-1][y][z])/(dx*dx);
				dphidysq[x][y][z] = 0.5*(phi_old[x][y+1][z] - 2*phi_old[x][y][z] + phi_old[x][y-1][z])/(dy*dy);
				dphidzsq[x][y][z] = 0.5*(phi_old[x][y][z+1] - 2*phi_old[x][y][z] + phi_old[x][y][z-1])/(dz*dz);
			}
		}
	}
	
	for(int x = 2 ; x <nx-2 ; x++)
	{
		for(int y = 2 ; y < ny-2 ; y++)
		{	
			for(int z = 2 ; z < nz-2 ; z++)
			{				
				double term1 = epsi[0][0]*dphidxsq[x][y][z] + epsi[1][1]*dphidysq[x][y][z] + epsi[2][2]*dphidzsq[x][y][z] + 0.5*(epsi[0][1]*dphidxdy[x][y][z]) + 0.5*(epsi[0][1]*dphidydx[x][y][z]) + 0.5*(epsi[0][2]*dphidzdx[x][y][z]) + 0.5*(epsi[1][2]*dphidydz[x][y][z]) ;// + (epsi[1][0]*dphidydx[x][y]) ;
				double term2 = phi_old[x][y][z]*(phi_old[x][y][z]*phi_old[x][y][z] - phi_old[x][y][z]) ;	
				double rhs = term1 - term2 ; 
			
				phi[x][y][z] = phi_old[x][y][z] + L*dt*rhs ;
				if(phi[x][y][z] > 1.0)
					phi[x][y][z] = 1.0 ;
				if(phi[x][y][z] < 0.0)
					phi[x][y][z] = 0.0 ;
			
			}
		}
	}
	
	for(int x = 0 ; x <nx ; x++)
	{
		for(int y = 0 ; y < ny ; y++)
		{
			for(int z = 0 ; z < nz ; z++)
			{
				phi_old[x][y][z] = phi[x][y][z] ;
			}
		}
	}
	

}

return 0 ;

}







