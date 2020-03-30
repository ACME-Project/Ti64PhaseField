
//------------------------Initialization of global variables/parameters-----------------------//

double M0 = 5.0e-21;
double R = 8.314 ;
double T = 1200.0 ;
double cr = 0.1;        //Cooling rate (K/s)
double G_normalize = R*T ;  //Scaling factor for terms having J/mol units
double tc = 4.4 ;               //Characteristic time scale (s)
double lc = 1e-07 ;
double dt = 0.01 ;        //Non-dimensional dt
double Vm = 0.00001;         //Molar volume (m^3/mol)
double w_norm = (R*T)/Vm ;       //Scaling factor for the double well height
double epsi_norm = (R*T*pow(lc,2))/Vm;    //Scaling factor for the gradient energy coeffcient
double W_prefac = 7.2e+06/w_norm ;   //Overall scaling for the double well depth 
double epsi_prefac = 9.0e-07/epsi_norm ;
double D = 2.4e-05*exp(-18040/T) ;             
double L_norm = (D*Vm)/(R*T*lc*lc) ;
double v_alpha =  3.34e-06 ;
double heat = 4.17e+8 ;
double T_trans = 1268.0 ;
double L_orig = v_alpha*T_trans/((T_trans - T)*heat*5e-07) ;
double L = L_orig/L_norm;
int variants = 20 ;          // Number of variants : input parameter needed for defining a mmsp-grid
std::vector<int> variant_ids ;
const int dim = 3;     //Spatial dimensions
double Dalal = 2.02e-03 ;        // The next four lines are the components of the Onsager mobility matrix
double Dvv = 2.02e-03 ;
const int nx = 24;
const int ny = 24;
const int nz = 24;
unsigned int nzp=nz/2+1;
double W_Al = 1.0 ;       //Relative magnitudes of double well depths for the three components
double W_V =  1.0 ; 
double W_Ti = 1.0 ;
long steps = 5e+8 ;    //Number of simulation steps.

double kappa_c = 0.1 ; //1.76 ;

//------------------------Definition of global containers/arrays-----------------------//

const int node_total = nx*ny*nz; 
double G_Alpha[node_total], G_Beta[node_total], W[node_total] ; 
double hphi[node_total][12], hphiprime[node_total][12], hphidoubleprime[node_total][12];
double gphi[node_total][12], gphiprime[node_total][12], gphidoubleprime[node_total][12];
double strain_energy[node_total][3];
double intenergy[nx][ny];
double dGAlpha_dAl,	dGBeta_dAl, dGAlpha_dV, dGBeta_dV, del_dGAlpha_dAl, del_dGBeta_dAl, del_dGAlpha_dV, del_dGBeta_dV, delsq_dGAlpha_dAl, delsq_dGAlpha_dV, delsq_dGBeta_dAl, delsq_dGBeta_dV ;
double G_Al_alpha, G_Ti_alpha, G_V_alpha, G_Al_beta, G_Ti_beta, G_V_beta ; 
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
void customoutput(int t) ;
void calculate_strains() ;
double nodesum() ;


//---------------------Definition of grids and grid variables-----------------------//
MMSP::grid<dim, store_type> grid(variants, 0, nx, 0, ny, 0, nz) ;
MMSP::grid<dim, store_type> update(variants, 0, nx, 0, ny, 0, nz) ;
MMSP::grid<dim, store_type> gradsqcal_grid(variants, 0, nx, 0, ny, 0, nz) ;
MMSP::grid<dim, store_type> gradsqcv_grid(variants, 0, nx, 0, ny, 0, nz) ;

const double Lx = g1(grid,0) - g0(grid,0) ;
const double Ly = g1(grid,1) - g0(grid,1) ;
const double Lz = g1(grid,2) - g0(grid,2) ;
const double dx = MMSP::dx(grid, 0) ;
const double dy = MMSP::dx(grid, 1) ;
const double dz = MMSP::dx(grid, 2) ;
double fs = nx/Lx;
double k[nx], fk[nx] ;


size_t align=sizeof(Complex);

array3<double> f1(nx,ny,nz,align);
array3<Complex> F1(nx,ny,nzp,align);

array3<double> f2(nx,ny,nz,align);
array3<Complex> F2(nx,ny,nzp,align);

array3<double> f3(nx,ny,nz,align);
array3<Complex> F3(nx,ny,nzp,align);

array3<double> f4(nx,ny,nz,align);
array3<Complex> F4(nx,ny,nzp,align);

array3<double> f5(nx,ny,nz,align);
array3<Complex> F5(nx,ny,nzp,align);

array3<Complex> dfdstr1(nx,ny,nzp,align);
array3<double> dfdstr_real1(nx,ny,nz,align);

array3<Complex> dfdstr2(nx,ny,nzp,align);
array3<double> dfdstr_real2(nx,ny,nz,align);

array3<Complex> dfdstr3(nx,ny,nzp,align);
array3<double> dfdstr_real3(nx,ny,nz,align);

array3<Complex> dfdstr4(nx,ny,nzp,align);
array3<double> dfdstr_real4(nx,ny,nz,align);

array3<Complex> dfdstr5(nx,ny,nzp,align);
array3<double> dfdstr_real5(nx,ny,nz,align);

rcfft3d Forward1(nx,ny,nz,f1,F1);
crfft3d Backward1(nx,ny,nz,dfdstr1,dfdstr_real1);

rcfft3d Forward2(nx,ny,nz,f2,F2);
crfft3d Backward2(nx,ny,nz,dfdstr2,dfdstr_real2);

rcfft3d Forward3(nx,ny,nz,f3,F3);
crfft3d Backward3(nx,ny,nz,dfdstr3,dfdstr_real3);

rcfft3d Forward4(nx,ny,nz,f4,F4);
crfft3d Backward4(nx,ny,nz,dfdstr4,dfdstr_real4);

rcfft3d Forward5(nx,ny,nz,f5,F5);
crfft3d Backward5(nx,ny,nz,dfdstr5,dfdstr_real5);


sfts eigen_alpha[12];
double c[6][6], sigma00[12][6] ;
double c_norm = 1;

#include "grid_dep_funcs.h"






