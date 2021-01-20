
//------------------------Initialization of global variables/parameters-----------------------//

double M0 = 5.0e-21;
double R = 8.314 ;
double T = 1200.0 ;
double cr = 0.5;     
double V0 = 1.0e-02/60.0 ;
double epsi_sq = 1.0e-07 ; 
double G_normalize = R*T ;  //Scaling factor for terms having J/mol units
double lc = 1e-07 ;
//double dt_old = 0.44 ; 			// Dimensional dt (s)
double dt = 3e-03 ; //dt_old/tc;        //Non-dimensional dt  
double Vm = 0.00001;         //Molar volume (m^3/mol)
double w_norm = (R*T)/Vm ;       //Scaling factor for the double well height
double epsi_norm = (R*T*pow(lc,2))/Vm;    //Scaling factor for the gradient energy coeffcient
double epsi_prefac = 9.0e-07/epsi_norm ;
double W_prefac = 7.2e+03/w_norm ;   //Overall scaling for the double well depth              
int variants = 20 ;          // Number of variants : input parameter needed for defining a mmsp-grid
const int dim = 2;     //Spatial dimensions
double Dnorm = 2.4e-5*exp(-18040.0/T) ;
double Dalal = 2.4e-5*exp(-18040.0/T)/Dnorm ; //0.001 ;        // The next four lines are the components of the Onsager mobility matrix
double Dvv = 1.0e-5*exp(-17460.0/T)/Dnorm ; //0.001 ;
double L_norm = (Dnorm*Vm)/(R*T*lc*lc) ;
double T_trans = 1268.0 ;
double fp = 2.0e-03 ; 
double L_orig = fp/(T_trans-T) ; 
double L = L_orig/L_norm; 
const int nx = 50;
const int ny = 50;    
const int nz = 50;        
double W_Al = 0.5 ;       //Relative magnitudes of double well depths for the three components
double W_V =  1.0 ; 
double W_Ti = 1.0 ;
int steps = 500000 ;    //Number of simulation steps.

double kappa_c = 1.0 ;


//------------------------Definition of global containers/arrays-----------------------//


//double strain_energy[node_total][3];
double intenergy[nx][ny][nz];
const int node_total = nx*ny*nz ; 
double G_Alpha[node_total], G_Beta[node_total], W[node_total] ; 
double dGAlpha_dAl,	dGBeta_dAl, dGAlpha_dV, dGBeta_dV, del_dGAlpha_dAl, del_dGBeta_dAl, del_dGAlpha_dV, del_dGBeta_dV, delsq_dGAlpha_dAl, delsq_dGAlpha_dV, delsq_dGBeta_dAl, delsq_dGBeta_dV ;
double G_Al_alpha, G_Ti_alpha, G_V_alpha, G_Al_beta, G_Ti_beta, G_V_beta ; 
double epsi[12][3][3] ;

//---------------------Declaration of user-defined functions-------------------------//

void thermo_auxillary_terms(MMSP::vector<store_type> gradient, MMSP::vector<store_type> gradientsq, double c_Al, double c_V) ;
void customoutput(int t) ;
void calculate_strains() ;
void intstrain();
double nodesum() ;


//---------------------Declaration of grids and grid variables-----------------------//

MMSP::grid<dim, store_type> grid(variants, 0, nx, 0, ny, 0, nz) ;
MMSP::grid<dim, store_type> gradsqcal_grid(variants, 0, nx, 0, ny, 0, nz) ;
MMSP::grid<dim, store_type> gradsqcv_grid(variants, 0, nx, 0, ny, 0, nz) ;
MMSP::grid<dim, store_type> intenergies(13, 0, nx, 0, ny, 0, nz) ;
MMSP::grid<dim, store_type> selfenergies(13, 0, nx, 0, ny, 0, nz) ;
MMSP::grid<dim, store_type> nuc_grid(13, 0, nx, 0, ny, 0, nz) ;
const double Lx = g1(grid,0) - g0(grid,0) ;
const double Ly = g1(grid,1) - g0(grid,1) ;
const double Lz = g1(grid,2) - g0(grid,2) ;
const double dx = MMSP::dx(grid, 0) ;
const double dy = MMSP::dx(grid, 1) ;
const double dz = MMSP::dx(grid, 2) ;
double k[nx], fk[nx] ;

//-------------------Definition of user-defined classes---------------------------------//

class sfts    //Data structure to store the eigenstrains
{
	public:
	double e[6] ;
	
	sfts()
	{
		e[0] = 0.0 ;
		e[1] = 0.0 ;
		e[2] = 0.0 ;
		e[3] = 0.0 ;
		e[4] = 0.0 ;
		e[5] = 0.0 ;
		
	}
};

sfts eigen_alpha[12];



//----------Including the various modules-----------//

#include "Array.h"
#include "/home/arun/Documents/mmsp-arun/fftw++-2.05/fftw++.h"

//Namespaces required for fftw 
using namespace utils;
using namespace Array;
using namespace fftwpp;



#include "thermo_modules.hpp"
#include "strain_modules.hpp"



//---------Definition of user defined functions--------------//

double nodesum(MMSP::grid<dim, store_type> grid, MMSP::vector<int> s)
{
	double phisum = 0 ;
	for(int h = 0 ; h < length(grid(s)) ; h++)
	{
		int hindex = MMSP::index(grid(s),h) ;
		if((hindex <= 12))
		{
			if(grid(s)[hindex] > 0.0001)
			{
				phisum += grid(s)[hindex] ; 
			}	
		}
	}
	return phisum ;
}



void initialize_epsi()
{
	
epsi[0][0][0] = 0.0806;
epsi[0][0][1] = 0.0301 ;
epsi[0][0][2] = -0.0642 ;
epsi[0][1][0] = 0.0301 ;
epsi[0][1][1] = 0.3445;
epsi[0][1][2] = -0.0277 ;
epsi[0][2][0] = -0.0642 ;
epsi[0][2][1] = -0.0277 ;
epsi[0][2][2] = 0.3834;

epsi[1][0][0] = 0.3445;
epsi[1][0][1] = -0.1183 ;
epsi[1][0][2] = -0.0486 ;
epsi[1][1][0] = -0.1183 ;
epsi[1][1][1] = 0.0806 ; 
epsi[1][1][2] = 0.0011 ;
epsi[1][2][0] = -0.0486 ;
epsi[1][2][1] = 0.0011 ;
epsi[1][2][2] = 0.0999 ;

epsi[2][0][0] = 0.0895 ;
epsi[2][0][1] = -0.0300 ;
epsi[2][0][2] = -0.0635 ;
epsi[2][1][0] = -0.0300 ;
epsi[2][1][1] = 0.0759 ;
epsi[2][1][2] = 0.0214 ;
epsi[2][2][0] = -0.0635 ;
epsi[2][2][1] = 0.0214 ;
epsi[2][2][2] = 0.3847 ;

epsi[3][0][0] = 0.0895  ;
epsi[3][0][1] = 0.0300 ;
epsi[3][0][2] = -0.0635 ;
epsi[3][1][0] = 0.0300 ;
epsi[3][1][1] = 0.0759 ;
epsi[3][1][2] = -0.0214 ;
epsi[3][2][0] = -0.0635 ;
epsi[3][2][1] = -0.0214 ;
epsi[3][2][2] = 0.3847 ;

epsi[4][0][0] = 0.0989 ;
epsi[4][0][1] = -0.0021 ;
epsi[4][0][2] = -0.0226 ;
epsi[4][1][0] = -0.0021 ;
epsi[4][1][1] = 0.1575 ;
epsi[4][1][2] = 0.1592 ;
epsi[4][2][0] = -0.0227 ;
epsi[4][2][1] = 0.1592 ;
epsi[4][2][2] = 0.2936 ;

epsi[5][0][0] = 0.0806 ;
epsi[5][0][1] = -0.0301 ;
epsi[5][0][2] = -0.0642 ;
epsi[5][1][0] = -0.0301 ;
epsi[5][1][1] = 0.0860;
epsi[5][1][2] = 0.0277 ;
epsi[5][2][0] = -0.0642 ;
epsi[5][2][1] = 0.0277 ;
epsi[5][2][2] = 0.3834 ;

epsi[6][0][0] = 0.3240 ;
epsi[6][0][1] = -0.1376 ;
epsi[6][0][2] = -0.0411 ;
epsi[6][1][0] = -0.1376 ;
epsi[6][1][1] = 0.1322 ;
epsi[6][1][2] = -0.0016 ;
epsi[6][2][0] = -0.0411 ;
epsi[6][2][1] = -0.0016 ;
epsi[6][2][2] = 0.0938 ;

epsi[7][0][0] = 0.0938 ;
epsi[7][0][1] = 0.0016 ;
epsi[7][0][2] = -0.0411 ;
epsi[7][1][0] = 0.0016 ;
epsi[7][1][1] = 0.1322 ;
epsi[7][1][2] = 0.1376 ;
epsi[7][2][0] = -0.0411 ;
epsi[7][2][1] = 0.1376 ;
epsi[7][2][2] = 0.3240 ;

epsi[8][0][0] = 0.0758 ;
epsi[8][0][1] = -0.0259 ;
epsi[8][0][2] = -0.0278 ;
epsi[8][1][0] = -0.0259 ;
epsi[8][1][1] = 0.0771 ;
epsi[8][1][2] = 0.0101 ;
epsi[8][2][0] = -0.0278 ;
epsi[8][2][1] = 0.0101 ;
epsi[8][2][2] = 0.3971 ;

epsi[9][0][0] = 0.3456 ;
epsi[9][0][1] = -0.1251 ;
epsi[9][0][2] = -0.0102 ;
epsi[9][1][0] = -0.1251 ;
epsi[9][1][1] = 0.1123 ;
epsi[9][1][2] = -0.0155 ;
epsi[9][2][0] = -0.0102 ;
epsi[9][2][1] = -0.0155 ;
epsi[9][2][2] = 0.0121 ;

epsi[10][0][0] = 0.0758 ;
epsi[10][0][1] = 0.0259 ;
epsi[10][0][2] = -0.0278 ;
epsi[10][1][0] = 0.0259 ;
epsi[10][1][1] = 0.0771 ;
epsi[10][1][2] = -0.0101 ;
epsi[10][2][0] = -0.0278 ;
epsi[10][2][1] = -0.0101 ;
epsi[10][2][2] = 0.3971 ;

epsi[11][0][0] = 0.0863 ;
epsi[11][0][1] = 0.0177 ;
epsi[11][0][2] = -0.0254 ;
epsi[11][1][0] = 0.0177 ;
epsi[11][1][1] = 0.0819 ;
epsi[11][1][2] = 0.0729 ;
epsi[11][2][0] = -0.0254 ;
epsi[11][2][1] = 0.0729 ;
epsi[11][2][2] = 0.3818 ;

}


void initialize_alpha_eigen()
{

eigen_alpha[0].e[0] = -0.083  ;
eigen_alpha[0].e[5] = 0.0095  ;
eigen_alpha[0].e[4] = 0.0  ;
eigen_alpha[0].e[1] = 0.123  ;
eigen_alpha[0].e[3] = 0.0  ;
eigen_alpha[0].e[2] = 0.035  ;

eigen_alpha[1].e[0] = -0.083  ;
eigen_alpha[1].e[5] = 0.0  ;
eigen_alpha[1].e[4] = 0.0095  ;
eigen_alpha[1].e[1] = 0.035  ;
eigen_alpha[1].e[3] = 0.0  ;
eigen_alpha[1].e[2] = 0.123  ;

eigen_alpha[2].e[0] = 0.079  ;
eigen_alpha[2].e[5] = -0.0359  ;
eigen_alpha[2].e[4] = -0.0264  ;
eigen_alpha[2].e[1] = 0.0047  ;
eigen_alpha[2].e[3] = 0.0810  ;
eigen_alpha[2].e[2] = -0.0087  ;

eigen_alpha[3].e[0] = -0.079  ;
eigen_alpha[3].e[5] = 0.0359  ;
eigen_alpha[3].e[4] = 0.0264  ;
eigen_alpha[3].e[1] = 0.0047  ;
eigen_alpha[3].e[3] = 0.0810  ;
eigen_alpha[3].e[2] = -0.0087  ;

eigen_alpha[4].e[0] = 0.079  ;
eigen_alpha[4].e[5] = -0.0359  ;
eigen_alpha[4].e[4] = 0.0264  ;
eigen_alpha[4].e[1] = 0.0047  ;
eigen_alpha[4].e[3] = -0.0810  ;
eigen_alpha[4].e[2] = -0.0087  ;

eigen_alpha[5].e[0] = 0.079  ;
eigen_alpha[5].e[5] = 0.0359  ;
eigen_alpha[5].e[4] = -0.0264  ;
eigen_alpha[5].e[1] = 0.0047  ;
eigen_alpha[5].e[3] = -0.0810  ;
eigen_alpha[5].e[2] = -0.0087  ;

eigen_alpha[6].e[0] = -0.083  ;
eigen_alpha[6].e[5] = -0.0095  ;
eigen_alpha[6].e[4] = 0.0  ;
eigen_alpha[6].e[1] = 0.123  ;
eigen_alpha[6].e[3] = 0.0  ;
eigen_alpha[6].e[2] = 0.035  ;

eigen_alpha[7].e[0] = -0.083  ;
eigen_alpha[7].e[5] = 0.0  ;
eigen_alpha[7].e[4] = -0.0095  ;
eigen_alpha[7].e[1] = 0.035  ;
eigen_alpha[7].e[3] = 0.0  ;
eigen_alpha[7].e[2] = 0.123  ;

eigen_alpha[8].e[0] = 0.079  ;
eigen_alpha[8].e[5] = -0.0264  ;
eigen_alpha[8].e[4] = -0.0359  ;
eigen_alpha[8].e[1] = -0.0087  ;
eigen_alpha[8].e[3] = 0.0810  ;
eigen_alpha[8].e[2] = 0.0047  ;

eigen_alpha[9].e[0] = 0.079  ;
eigen_alpha[9].e[5] = 0.0264  ;
eigen_alpha[9].e[4] = 0.0359  ;
eigen_alpha[9].e[1] = -0.0087  ;
eigen_alpha[9].e[3] = 0.0810  ;
eigen_alpha[9].e[2] = 0.0047  ;

eigen_alpha[10].e[0] = 0.079  ;
eigen_alpha[10].e[5] = -0.0264  ;
eigen_alpha[10].e[4] = 0.0359  ;
eigen_alpha[10].e[1] = -0.0087  ;
eigen_alpha[10].e[3] = -0.0810  ;
eigen_alpha[10].e[2] = 0.0047  ;

eigen_alpha[11].e[0] = 0.079  ; 
eigen_alpha[11].e[5] = 0.0264  ; 
eigen_alpha[11].e[4] = -0.0359  ; 
eigen_alpha[11].e[1] = -0.0087  ;
eigen_alpha[11].e[3] = -0.0810  ; 
eigen_alpha[11].e[2] = 0.0047  ; 
}




