
//------------------------Initialization of global variables/parameters-----------------------//

double R = 8.314 ;
double T = 1023 ;
double T_orig = 1023;
double cr = 10.0;        //Cooling rate (K/s)
double L_orig = 0.139 ;  //Phase field mobility
double G_normalize = R*T ;  //Scaling factor for terms having J/mol units
double lold = 9.677*pow(10, -8) ;  //Grid size (m)
double lc = 9.677*pow(10,-8) ;    //Characteristic length scale, or capillary length
double dx_nd = lold/lc ;          // Non-dimensional grid size
double tc = 4.4 ;               //Characteristic time scale (s)
double dt_old = 0.44 ; 			// Dimensional dt (s)
double dt = dt_old/tc;        //Non-dimensional dt
double Dalal = 0.0001 ;        // The next four lines are the components of the Onsager mobility matrix
double Dalv = -0.00008 ;
double Dvv = 0.0001 ;
double Dval = -0.00008 ;
double Vm = 0.00001;         //Molar volume (m^3/mol)
double w_norm = (R*T)/Vm ;       //Scaling factor for the double well height
double epsi_norm = (R*T*pow(lc,2))/Vm;    //Scaling factor for the gradient energy coeffcient
int variants = 20 ;          // Number of variants : input parameter needed for defining a mmsp-grid
const int dim = 3;     //Spatial dimensions
const int nx = 20;
const int ny = 20;
const int nz = 20;
double W_prefac = 6.0 ;   //Overall scaling for the double well depth              
double W_Al = 1.0 ;       //Relative magnitudes of double well depths for the three components
double W_V =  1.0 ; 
double W_Ti = 1.5 ;
double scaling = 1/(2*3.14*3.14) ;        //Scaling factor, needed for Fourier Transformation integrals. 
int steps = 70000 ;    //Number of simulation steps.

//------------------------Definition of global containers/arrays-----------------------//

const int node_total = nx*ny*nz ; 
double G_Alpha[node_total], G_Beta[node_total], W[node_total] ; 
double hphi[node_total], hphiprime[node_total], hphidoubleprime[node_total];
double gphi[node_total], gphiprime[node_total], gphidoubleprime[node_total];
double dGAlpha_dAl,	dGBeta_dAl, dGAlpha_dV, dGBeta_dV, del_dGAlpha_dAl, del_dGBeta_dAl, del_dGAlpha_dV, del_dGBeta_dV, delsq_dGAlpha_dAl, delsq_dGAlpha_dV, delsq_dGBeta_dAl, delsq_dGBeta_dV ;
double G_Al_alpha, G_Ti_alpha, G_V_alpha, G_Al_beta, G_Ti_beta, G_V_beta ; 
double fs = nx/Lx;
double k[nx], fk[nx] ;
double epsi[12][3][3] ;


//-----------------------Controlling the number and type of output files----------------//

// If a particular file is not needed, change the value to 0
bool output_pfsq = 1 ;
bool output_cal = 1 ;
bool output_cv = 1 ;
bool output_chemnuc = 1 ;
bool output_strainintnuc = 1 ;

//---------------------Definition of functions-------------------------//
void thermo_auxillary_terms(MMSP::vector<store_type> gradient, MMSP::vector<store_type> gradientsq, double c_Al, double c_V) ;
void customoutput(MMSP::grid<dim, store_type> grid, int t) ;
void calculate_strains(MMSP::grid<dim, store_type> grid) ;
double nodesum(MMSP::grid<dim, store_type> grid, MMSP::vector<int> s) ;


//---------------------Definition of grids and grid variables-----------------------//
MMSP::grid<dim, store_type> grid(variants, 0, nx, 0, ny, 0, nz) ;
MMSP::grid<dim, store_type> gradsqcal_grid(variants, 0, nx, 0, ny, 0, nz) ;
MMSP::grid<dim, store_type> gradsqcv_grid(variants, 0, nx, 0, ny, 0, nz) ;
MMSP::grid<dim, store_type> gradxgrid(variants, 0, nx, 0, ny, 0, nz) ;
MMSP::grid<dim, store_type> gradygrid(variants, 0, nx, 0, ny, 0, nz) ;
MMSP::grid<dim, store_type> gradzgrid(variants, 0, nx, 0, ny, 0, nz) ;
const double Lx = g1(grid,0) - g0(grid,0) ;
const double Ly = g1(grid,1) - g0(grid,1) ;
const double Lz = g1(grid,2) - g0(grid,2) ;
const double dx = MMSP::dx(grid, 0) ;
const double dy = MMSP::dx(grid, 1) ;
const double dz = MMSP::dx(grid, 2) ;

//-------------------Definition of user-defined classes---------------------------------//
class sfts
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


//--------------Definition of variables from the array class as defined by the fftw module, needed for the FFT operations----------//

#include "Array.h"
#include "/home/arun/Documents/mmsp-arun/fftw++-2.05/fftw++.h"

using namespace utils;
using namespace Array;
using namespace fftwpp;


size_t align=sizeof(Complex);
  
					//Container arrays for the 12 phase fields to be used in a regular simulation timestep. 'fi' for the real space, and 'Fi' for the Fourier space. 
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
	
					//Container arrays for the 12 phase fields to be used in a nucleation timestep. 'finuc' for the real space, and 'Finuc' for the Fourier space. 
array3<double> f1nuc(nx,ny,nz,align);
array3<Complex> F1nuc(nx,ny,nz,align);
array3<double> f2nuc(nx,ny,nz,align);
array3<Complex> F2nuc(nx,ny,nz,align);
array3<double> f3nuc(nx,ny,nz,align);
array3<Complex> F3nuc(nx,ny,nz,align);
array3<double> f4nuc(nx,ny,nz,align);
array3<Complex> F4nuc(nx,ny,nz,align);
array3<double> f5nuc(nx,ny,nz,align);
array3<Complex> F5nuc(nx,ny,nz,align);

					//Container arrays for the strain energy densities for each of the 12 variants. 'dfdstri' for the Fourier space, and 'dfdstr_reali' for the real space. 
array3<Complex> dfdstr1(nx,ny,nz,align);
array3<double> dfdstr_real1(nx,ny,nz,align);
array3<Complex> dfdstr2(nx,ny,nz,align);
array3<double> dfdstr_real2(nx,ny,nz,align);
array3<Complex> dfdstr3(nx,ny,nz,align);
array3<double> dfdstr_real3(nx,ny,nz,align);
array3<Complex> dfdstr4(nx,ny,nz,align);
array3<double> dfdstr_real4(nx,ny,nz,align);
array3<Complex> dfdstr5(nx,ny,nz,align);
array3<double> dfdstr_real5(nx,ny,nz,align);
array3<Complex> dfdstr6(nx,ny,nz,align);
array3<double> dfdstr_real6(nx,ny,nz,align);
array3<Complex> dfdstr7(nx,ny,nz,align);
array3<double> dfdstr_real7(nx,ny,nz,align);
array3<Complex> dfdstr8(nx,ny,nz,align);
array3<double> dfdstr_real8(nx,ny,nz,align);
array3<Complex> dfdstr9(nx,ny,nz,align);
array3<double> dfdstr_real9(nx,ny,nz,align);
array3<Complex> dfdstr10(nx,ny,nz,align);
array3<double> dfdstr_real10(nx,ny,nz,align);
array3<Complex> dfdstr11(nx,ny,nz,align);
array3<double> dfdstr_real11(nx,ny,nz,align);
array3<Complex> dfdstr12(nx,ny,nz,align);
array3<double> dfdstr_real12(nx,ny,nz,align);
	
					//Container arrays for the elastic interaction energies for each of the 12 variants, to be used in a nucleation timestep. 'elinti' for the Fourier space, and 'elint_reali' for the real space. 
array3<Complex> elint1(nx,ny,nz,align);
array3<double> elint_real1(nx,ny,nz,align);
array3<Complex> elint2(nx,ny,nz,align);
array3<double> elint_real2(nx,ny,nz,align);
array3<Complex> elint3(nx,ny,nz,align);
array3<double> elint_real3(nx,ny,nz,align);
array3<Complex> elint4(nx,ny,nz,align);
array3<double> elint_real4(nx,ny,nz,align);
array3<Complex> elint5(nx,ny,nz,align);
array3<double> elint_real5(nx,ny,nz,align);
	
					//Arrays to calculate and store total strain energy in the system as a function of the preexisting microstructure. 
array3<Complex> TotStr(nx,ny,nz,align);
array3<double> TotStr_real(nx,ny,nz,align);

					
rcfft3d Forward1(nx,ny,nz,f1,F1);
rcfft3d Forward2(nx,ny,nz,f2,F2);
rcfft3d Forward3(nx,ny,nz,f3,F3);
rcfft3d Forward4(nx,ny,nz,f4,F4);
rcfft3d Forward5(nx,ny,nz,f5,F5);
rcfft3d Forward6(nx,ny,nz,f6,F6);
rcfft3d Forward7(nx,ny,nz,f7,F7);
rcfft3d Forward8(nx,ny,nz,f8,F8);
rcfft3d Forward9(nx,ny,nz,f9,F9);
rcfft3d Forward10(nx,ny,nz,f10,F10);
rcfft3d Forward11(nx,ny,nz,f11,F11);
rcfft3d Forward12(nx,ny,nz,f12,F12);
	
crfft3d Backward1(nx,ny,nz,dfdstr1,dfdstr_real1);
crfft3d Backward2(nx,ny,nz,dfdstr2,dfdstr_real2);
crfft3d Backward3(nx,ny,nz,dfdstr3,dfdstr_real3);
crfft3d Backward4(nx,ny,nz,dfdstr4,dfdstr_real4);
crfft3d Backward5(nx,ny,nz,dfdstr5,dfdstr_real5);
crfft3d Backward6(nx,ny,nz,dfdstr6,dfdstr_real6);
crfft3d Backward7(nx,ny,nz,dfdstr7,dfdstr_real7);
crfft3d Backward8(nx,ny,nz,dfdstr8,dfdstr_real8);
crfft3d Backward9(nx,ny,nz,dfdstr9,dfdstr_real9);
crfft3d Backward10(nx,ny,nz,dfdstr10,dfdstr_real10);
crfft3d Backward11(nx,ny,nz,dfdstr11,dfdstr_real11);
crfft3d Backward12(nx,ny,nz,dfdstr12,dfdstr_real12);
crfft3d Backward13(nx,ny,nz,TotStr,TotStr_real);
	
rcfft3d Forward1n(nx,ny,nz,f1nuc,F1nuc);
rcfft3d Forward2n(nx,ny,nz,f2nuc,F2nuc);
rcfft3d Forward3n(nx,ny,nz,f3nuc,F3nuc);
rcfft3d Forward4n(nx,ny,nz,f4nuc,F4nuc);
rcfft3d Forward5n(nx,ny,nz,f5nuc,F5nuc);
	
crfft3d Backward1n(nx,ny,nz,elint1,elint_real1);
crfft3d Backward2n(nx,ny,nz,elint2,elint_real2);
crfft3d Backward3n(nx,ny,nz,elint3,elint_real3);
crfft3d Backward4n(nx,ny,nz,elint4,elint_real4);
crfft3d Backward5n(nx,ny,nz,elint5,elint_real5);

//----------Including the various modules-----------//


#include "thermo_modules.hpp"
#include "strain_modules.hpp"
#include "outputgeneration.hpp"



//---------User defined functions--------------//

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
	
epsi[0][0][0] = 0.0806 ;
epsi[0][0][1] = 0.0301 ;
epsi[0][0][2] = -0.0642 ;
epsi[0][1][0] = 0.0301 ;
epsi[0][1][1] = 0.0806 ;
epsi[0][1][2] = -0.0277 ;
epsi[0][2][0] = -0.0642 ;
epsi[0][2][1] = -0.0277 ;
epsi[0][2][2] = 0.3834 ;

epsi[1][0][0] = 0.3445 ;
epsi[1][0][1] = -0.1183 ;
epsi[1][0][2] = -0.0486 ;
epsi[1][1][0] = -0.1183 ;
epsi[1][1][1] = 0.1056 ;
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

epsi[3][0][0] = 0.0895 ;
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
epsi[5][1][1] = 0.0860 ;
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
eigen_alpha[11].e[2] = 0.0047  ; }




