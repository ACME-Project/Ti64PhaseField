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
#include <complex>
#include <valarray>
typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

#include "global.h"
//#include "grid_dep_vars.hpp"
#include "grid_ind_funcs.h"

using namespace std;


double gdiff(double cal, double cv, double T)
{
    double cal_alpha = cal ;
    double cal_beta =  cal ;
    double cv_alpha = cv ;
    double cv_beta = cv ;
    
    G_Al_alpha = (Al_alpha_a_3 + Al_alpha_b_3*T+ Al_alpha_c_3*T*log(T) + Al_alpha_d_3*T*T) ;
    G_Al_beta = (Al_beta_a_3 + Al_beta_b_3*T+ Al_beta_c_3*T*log(T) + Al_beta_d_3*T*T) ;
	
	G_Ti_alpha = (Ti_alpha_a_2 + Ti_alpha_b_2*T+ Ti_alpha_c_2*T*log(T) + Ti_alpha_d_2*T*T) ;
	G_Ti_beta = (Ti_beta_a_1 + Ti_beta_b_1*T+ Ti_beta_c_1*T*log(T) + Ti_beta_d_1*T*T) ;
	
	G_V_alpha = (V_alpha_a_2 + V_alpha_b_2*T+ V_alpha_c_2*T*log(T) + V_alpha_d_2*T*T) ; //(21.96 - T/50)
	G_V_beta = (V_beta_a_2 + V_beta_b_2*T+ V_beta_c_2*T*log(T) + V_beta_d_2*T*T) ;
    
    double alpha_alv = cal*cv*(L0_HCP_Al_V) ;
    double alpha_alti = cal*(1-cal-cv)*(L0_HCP_Al_Ti + L1_HCP_Al_Ti*(2*cal + cv - 1)) ;
    double alpha_vti = cv*(1-cal-cv)*(L0_HCP_Ti_V) ;
                                
    double beta_alv = cal*cv*(L0_BCC_Al_V) ;
    double beta_alti = cal*(1-cal-cv)*(L0_BCC_Al_Ti) ;
    double beta_vti = cv*(1-cal-cv)*(L0_BCC_Ti_V) ;
    double beta_alvti = cal*cv*(1-cal-cv)*L0_BCC_Al_Ti_V ;                        

    double galpha = (cal_alpha*G_Al_alpha + cv_alpha*G_V_alpha + (1-cal_alpha-cv_alpha)*G_Ti_alpha +
             R*T*(cal_alpha*log(cal_alpha) + cv_alpha*log(cv_alpha) + (1-cal_alpha-cv_alpha)*log((1-cal_alpha-cv_alpha))) + alpha_alv + alpha_alti + alpha_vti)/G_normalize ;
            

    double gbeta = (cal_beta*G_Al_beta + cv_beta*G_V_beta + (1-cal_beta-cv_beta)*G_Ti_beta +
             R*T*(cal_beta*log(cal_beta) + cv_beta*log(cv_beta) +
                  (1-cal_beta-cv_beta)*log((1-cal_beta-cv_beta))) 
             + beta_alv + beta_alti + beta_vti + beta_alvti)/G_normalize ;
             
    return (gbeta-galpha) ; 
    
}

double q_alpha_al(double cal, double cv, double T)
{
    double ans = cal*q_alpha_al_al + cv*q_alpha_al_v + (1-cv-cal)*q_alpha_al_ti + cal*cv*a_alpha_al_al_v + cal*(1-cal-cv)*a_alpha_al_al_ti + cv*(1-cal-cv)*(a_alpha_al_ti_v) ;
    return ans ;
}

double q_alpha_v(double cal, double cv, double T)
{
    double ans = cal*q_alpha_v_al + cv*q_alpha_v_v + (1-cv-cal)*q_alpha_v_ti + cal*cv*a_alpha_v_al_v + cal*(1-cal-cv)*a_alpha_v_al_ti + cv*(1-cal-cv)*(a_alpha_v_ti_v) ;
    return ans ;
}

double q_beta_al(double cal, double cv, double T)
{
    double ans = cal*q_beta_al_al + cv*q_beta_al_v + (1-cv-cal)*q_beta_al_ti + cal*cv*a_beta_al_al_v + cal*(1-cal-cv)*(a_beta_al_al_ti + (2*cal + cv - 1)*a1_beta_al_al_ti)  ;
    return ans ;
}

double q_beta_v(double cal, double cv, double T)
{
    double ans = cal*q_beta_v_al + cv*q_beta_v_v + (1-cv-cal)*q_beta_v_ti + cal*cv*a_beta_v_al_v + cv*(1-cal-cv)*(a_beta_v_ti_v + (2*cv + cal - 1)*a1_beta_v_ti_v) ;
    return ans ;
}


double Malal_alpha(double cal, double cv, double T)
{
    double ans = cal*(M0/(R*T))*exp(-q_alpha_al(cal,cv,T)/(R*T)) ;
    return ans ;
}
    
double Malal_beta(double cal, double cv, double T)
{
    double ans = cal*(M0/(R*T))*exp(-q_beta_al(cal,cv,T)/(R*T)) ;
    return ans ;
}

double Mvv_alpha(double cal, double cv, double T)
{
    double ans = cv*(M0/(R*T))*exp(-q_alpha_v(cal,cv,T)/(R*T)) ;
    return ans ;
}

double Mvv_beta(double cal, double cv, double T)
{
    double ans = cv*(M0/(R*T))*exp(-q_beta_v(cal,cv,T)/(R*T)) ;
    return ans ;
}

double d_Malal_alpha_dcal(double cal, double cv, double T) 
{
    double ans = (M0/(R*T))*exp(-q_alpha_al(cal,cv,T)/(R*T)) - (cal*M0/(pow((R*T),2)))*exp(-q_alpha_al(cal,cv,T)/(R*T))*(q_alpha_al_al/G_normalize - q_alpha_al_ti/G_normalize + cv*a_alpha_al_al_v/G_normalize + (1-2*cal-cv)*a_alpha_al_al_ti/G_normalize - cv*a_alpha_al_ti_v/G_normalize) ; 
    return ans ;
}

double d_Malal_alpha_dcv(double cal, double cv, double T) 
{
    double ans = -(cal*M0/(pow((R*T),2)))*exp(-q_alpha_al(cal,cv,T)/(R*T))*(q_alpha_al_v/G_normalize - q_alpha_al_ti/G_normalize + cal*a_alpha_al_al_v/G_normalize + (1-2*cv-cal)*a_alpha_al_ti_v/G_normalize - cal*a_alpha_al_al_ti/G_normalize) ; 
    return ans ;
}

double d_Malal_beta_dcal(double cal, double cv, double T) 
{
    double ans = (M0/(R*T))*exp(-q_beta_al(cal,cv,T)/(R*T)) - (cal*M0/(pow((R*T),2)))*exp(-q_beta_al(cal,cv,T)/(R*T))*(q_beta_al_al/G_normalize - q_beta_al_ti/G_normalize + cv*a_beta_al_al_v/G_normalize + (1-2*cal-cv)*a_beta_al_al_ti/G_normalize - cv*a_beta_al_ti_v/G_normalize) ; 
    return ans ;
}

double d_Malal_beta_dcv(double cal, double cv, double T) 
{
    double ans = -(cal*M0/(pow((R*T),2)))*exp(-q_beta_al(cal,cv,T)/(R*T))*(q_beta_al_v/G_normalize - q_beta_al_ti/G_normalize + cal*a_beta_al_al_v/G_normalize + (1-2*cv-cal)*a_beta_al_ti_v/G_normalize - cal*a_beta_al_al_ti/G_normalize) ; 
    return ans ;
}


double d_Mvv_alpha_dcal(double cal, double cv, double T) 
{
    double ans = -(cv*M0/(pow((R*T),2)))*exp(-q_alpha_v(cal,cv,T)/(R*T))*(q_alpha_v_al/G_normalize - q_alpha_v_ti/G_normalize + cv*a_alpha_v_al_v/G_normalize + (1-2*cal-cv)*a_alpha_v_al_ti/G_normalize - cv*a_alpha_v_ti_v/G_normalize) ; 
    return ans ;
}

double d_Mvv_alpha_dcv(double cal, double cv, double T) 
{
    double ans = (M0/(R*T))*exp(-q_alpha_v(cal,cv,T)/(R*T)) - (cv*M0/(pow((R*T),2)))*exp(-q_alpha_v(cal,cv,T)/(R*T))*(q_alpha_v_v/G_normalize - q_alpha_v_ti/G_normalize + cal*a_alpha_v_al_v/G_normalize + (1-2*cv-cal)*a_alpha_v_ti_v/G_normalize - cal*a_alpha_v_al_ti/G_normalize) ; 
    return ans ;
}

double d_Mvv_beta_dcal(double cal, double cv, double T) 
{
    double ans = -(cv*M0/(pow((R*T),2)))*exp(-q_beta_v(cal,cv,T)/(R*T))*(q_beta_v_al/G_normalize - q_beta_v_ti/G_normalize + cv*a_beta_v_al_v/G_normalize + (1-2*cal-cv)*a_beta_v_al_ti/G_normalize - cv*a_beta_v_ti_v/G_normalize) ; 
    return ans ;
}

double d_Mvv_beta_dcv(double cal, double cv, double T) 
{
    double ans = (M0/(R*T))*exp(-q_beta_v(cal,cv,T)/(R*T)) - (cv*M0/(pow((R*T),2)))*exp(-q_beta_v(cal,cv,T)/(R*T))*(q_beta_v_v/G_normalize - q_beta_v_ti/G_normalize + cal*a_beta_v_al_v/G_normalize + (1-2*cv-cal)*a_beta_v_ti_v/G_normalize - cal*a_beta_v_al_ti/G_normalize) ; 
    return ans ;
}

void initialize_epsi()
{
	
epsi[0][0][0] = 0.0806/0.3445 ;  //0.0806 ;
epsi[0][0][1] = 0.0301 ;
epsi[0][0][2] = -0.0642 ;
epsi[0][1][0] = 0.0301 ;
epsi[0][1][1] = 0.3445/0.3445 ; //0.1056 ;
epsi[0][1][2] = -0.0277 ;
epsi[0][2][0] = -0.0642 ;
epsi[0][2][1] = -0.0277 ;
epsi[0][2][2] = 0.3834 ;

epsi[1][0][0] = 0.3445/0.3445 ; //0.3445 ;
epsi[1][0][1] = -0.1183 ;
epsi[1][0][2] = -0.0486 ;
epsi[1][1][0] = -0.1183 ;
epsi[1][1][1] = 0.0806/0.3445 ; 
epsi[1][1][2] = 0.0011 ;
epsi[1][2][0] = -0.0486 ;
epsi[1][2][1] = 0.0011 ;
epsi[1][2][2] = 0.0999 ;

epsi[2][0][0] = 0.0895/0.0895 ;
epsi[2][0][1] = -0.0300 ;
epsi[2][0][2] = -0.0635 ;
epsi[2][1][0] = -0.0300 ;
epsi[2][1][1] = 0.0759/0.0895 ;
epsi[2][1][2] = 0.0214 ;
epsi[2][2][0] = -0.0635 ;
epsi[2][2][1] = 0.0214 ;
epsi[2][2][2] = 0.3847 ;

epsi[3][0][0] = 0.0895/0.0895  ;
epsi[3][0][1] = 0.0300 ;
epsi[3][0][2] = -0.0635 ;
epsi[3][1][0] = 0.0300 ;
epsi[3][1][1] = 0.0759/0.0895  ;
epsi[3][1][2] = -0.0214 ;
epsi[3][2][0] = -0.0635 ;
epsi[3][2][1] = -0.0214 ;
epsi[3][2][2] = 0.3847 ;

epsi[4][0][0] = 0.0989/0.1575 ;
epsi[4][0][1] = -0.0021 ;
epsi[4][0][2] = -0.0226 ;
epsi[4][1][0] = -0.0021 ;
epsi[4][1][1] = 0.1575/0.1575 ;
epsi[4][1][2] = 0.1592 ;
epsi[4][2][0] = -0.0227 ;
epsi[4][2][1] = 0.1592 ;
epsi[4][2][2] = 0.2936 ;

epsi[5][0][0] = 0.0806/0.0860 ;
epsi[5][0][1] = -0.0301 ;
epsi[5][0][2] = -0.0642 ;
epsi[5][1][0] = -0.0301 ;
epsi[5][1][1] = 0.0860/0.0860 ;
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


double omega_inv_ij(double* cnew, double* k ,  double modk)
{
	double sum = 0 ;
	for(int i = 0 ; i < 9 ; i++)
	{
		sum += cnew[i]*k[i] ;
		//cout<<cnew[i]<<" "<<k[i]<<endl;
	}

	
	return sum ; 
}

void define_c_sigma()
{

	for(int i = 0 ; i < 6 ; i++)
	{
		for(int j = 0 ; j < 6 ; j++)
		{
			c[i][j] = 0.0 ;   
		}
	}
    double scale = 100.0 ;
    c[0][0] = 175.0/scale;
    c[1][1] = 175.0/scale;
    c[2][2] = 220.0/scale;
    c[0][1] = 88.7/scale ;
    c[0][2] = 62.3/scale ;
    c[1][2] = 62.3/scale ;
    c[3][3] = 62.2/scale ;
    c[4][4] = 62.2/scale ;
    c[5][5] = (c[0][0]-c[0][1])/2.0 ;
    c[5][5] /= scale ;



	for(int variant = 0 ; variant < 12 ; variant++)
	{
		for(int i = 0 ; i < 6 ; i++)
		{
			sigma00[variant][i] = 0.0 ;
			for(int j = 0 ; j < 6 ; j++)
			{
				sigma00[variant][i] += c[i][j]*eigen_alpha[variant].e[j] ;
			}
		}			
	} 
}

double Bpq(double k1, double k2, double k3, double modk, int p, int q) {
	
	double c11[9] = {c[0][0], c[0][5], c[0][4], c[5][0], c[5][5], c[5][4], c[4][0], c[4][5], c[4][4]} ; 
	double c12[9] = {c[0][5], c[0][1], c[0][3], c[5][5], c[5][1], c[5][3], c[4][5], c[4][1], c[4][3]} ; 
	double c13[9] = {c[0][4], c[0][3], c[0][2], c[5][4], c[5][3], c[5][2], c[4][4], c[4][3], c[4][2]} ; 
	double c21[9] = {c[5][0], c[5][5], c[5][4], c[1][0], c[1][5], c[1][4], c[3][0], c[3][5], c[3][4]} ; 
	double c22[9] = {c[5][5], c[5][1], c[5][3], c[1][5], c[1][1], c[1][3], c[3][5], c[3][1], c[3][3]} ; 
	double c23[9] = {c[5][4], c[5][3], c[5][2], c[1][4], c[1][3], c[1][2], c[3][4], c[3][3], c[3][2]} ; 
	double c31[9] = {c[4][0], c[4][5], c[4][4], c[3][0], c[3][5], c[3][4], c[2][0], c[2][5], c[2][4]} ; 
	double c32[9] = {c[4][5], c[4][1], c[4][3], c[3][5], c[3][1], c[3][3], c[2][5], c[2][1], c[2][3]} ; 
	double c33[9] = {c[4][4], c[4][3], c[4][2], c[3][4], c[3][3], c[3][2], c[2][4], c[2][3], c[2][2]} ; 
	double k[9] = {k1*k1, k1*k2, k1*k3, k2*k1, k2*k2, k2*k3, k3*k1, k3*k2, k3*k3} ;

	double omega_inv[3][3] = {{omega_inv_ij(c11, k, modk), omega_inv_ij(c12, k, modk), omega_inv_ij(c13, k, modk)}, {omega_inv_ij(c21, k, modk), omega_inv_ij(c22, k, modk), omega_inv_ij(c23, k, modk)},
						{omega_inv_ij(c31, k, modk), omega_inv_ij(c32, k, modk), omega_inv_ij(c33, k, modk)}} ;


	double det_omega_inv = omega_inv[0][0]*(omega_inv[1][1]*omega_inv[2][2] - omega_inv[1][2]*omega_inv[2][1]) - omega_inv[0][1]*(omega_inv[1][0]*omega_inv[2][2] - omega_inv[1][2]*omega_inv[2][0]) +
					   omega_inv[0][2]*(omega_inv[1][0]*omega_inv[2][1] - omega_inv[1][1]*omega_inv[2][0]) + 1e-05  ;
					   
	double omega[3][3] ; 

	for(int i = 0 ; i < 2 ; i++)
	{
		for(int j = 0 ; j < 2 ; j++)
		{
			omega[i][j] =  pow(-1, (i+j))*(1/det_omega_inv)*(omega_inv[(j+1)%3][(i+1)%3]*omega_inv[(j+2)%3][(i+2)%3] - omega_inv[(j+1)%3][(i+2)%3]*omega_inv[(j+2)%3][(i+1)%3]) ;  //Inverse of omega_inv
		}
	
	}
	double k1l = sigma00[q-1][0]*k1/modk + sigma00[q-1][5]*k2/modk + sigma00[q-1][5]*k3/modk ; 
	double k2l = sigma00[q-1][5]*k1/modk + sigma00[q-1][1]*k2/modk + sigma00[q-1][3]*k3/modk ; 
	double k3l = sigma00[q-1][4]*k1/modk + sigma00[q-1][3]*k2/modk + sigma00[q-1][2]*k3/modk ; 
	
	double j1kl = omega[0][0]*k1l + omega[0][1]*k2l + omega[0][2]*k3l ;
	double j2kl = omega[1][0]*k1l + omega[1][1]*k2l + omega[1][2]*k3l ;
	double j3kl = omega[2][0]*k1l + omega[2][1]*k2l + omega[2][2]*k3l ;
	
	double i1jkl = (k1/modk)*(sigma00[p-1][0]*j1kl + sigma00[p-1][5]*j2kl + sigma00[p-1][4]*j3kl) ;
	double i2jkl = (k2/modk)*(sigma00[p-1][5]*j1kl + sigma00[p-1][1]*j2kl + sigma00[p-1][3]*j3kl) ;
	double i3jkl = (k3/modk)*(sigma00[p-1][4]*j1kl + sigma00[p-1][3]*j2kl + sigma00[p-1][2]*j3kl) ;
	
	double sfts_sum = 0 ;
	for(int i = 0 ; i < 6 ; i++)
	{
		for(int j = 0 ; j < 6 ; j++)
		{
			sfts_sum+= eigen_alpha[p-1].e[i]*c[i][j]*eigen_alpha[q-1].e[j] ;
		}
	}
			
    
    return (sfts_sum - (i1jkl + i2jkl + i3jkl)) ; 
}

double B(double k1, double k2, double k3, double modk, int variant) {
	double c11[9] = {c[0][0], c[0][5], c[0][4], c[5][0], c[5][5], c[5][4], c[4][0], c[4][5], c[4][4]} ; 
	double c12[9] = {c[0][5], c[0][1], c[0][3], c[5][5], c[5][1], c[5][3], c[4][5], c[4][1], c[4][3]} ; 
	double c13[9] = {c[0][4], c[0][3], c[0][2], c[5][4], c[5][3], c[5][2], c[4][4], c[4][3], c[4][2]} ; 
	double c21[9] = {c[5][0], c[5][5], c[5][4], c[1][0], c[1][5], c[1][4], c[3][0], c[3][5], c[3][4]} ; 
	double c22[9] = {c[5][5], c[5][1], c[5][3], c[1][5], c[1][1], c[1][3], c[3][5], c[3][1], c[3][3]} ; 
	double c23[9] = {c[5][4], c[5][3], c[5][2], c[1][4], c[1][3], c[1][2], c[3][4], c[3][3], c[3][2]} ; 
	double c31[9] = {c[4][0], c[4][5], c[4][4], c[3][0], c[3][5], c[3][4], c[2][0], c[2][5], c[2][4]} ; 
	double c32[9] = {c[4][5], c[4][1], c[4][3], c[3][5], c[3][1], c[3][3], c[2][5], c[2][1], c[2][3]} ; 
	double c33[9] = {c[4][4], c[4][3], c[4][2], c[3][4], c[3][3], c[3][2], c[2][4], c[2][3], c[2][2]} ; 
	double k[9] = {k1*k1, k1*k2, k1*k3, k2*k1, k2*k2, k2*k3, k3*k1, k3*k2, k3*k3} ;


	double omega_inv[3][3] = {{omega_inv_ij(c11, k, modk), omega_inv_ij(c12, k, modk), omega_inv_ij(c13, k, modk)}, {omega_inv_ij(c21, k, modk), omega_inv_ij(c22, k, modk), omega_inv_ij(c23, k, modk)},
						{omega_inv_ij(c31, k, modk), omega_inv_ij(c32, k, modk), omega_inv_ij(c33, k, modk)}} ;


	double det_omega_inv = omega_inv[0][0]*(omega_inv[1][1]*omega_inv[2][2] - omega_inv[1][2]*omega_inv[2][1]) - omega_inv[0][1]*(omega_inv[1][0]*omega_inv[2][2] - omega_inv[1][2]*omega_inv[2][0]) +
					   omega_inv[0][2]*(omega_inv[1][0]*omega_inv[2][1] - omega_inv[1][1]*omega_inv[2][0]) + 1e-05  ;
	if(det_omega_inv==0) cout<<"Found zero det"<<endl ;
	
					   
	double omega[3][3] ; 

	for(int i = 0 ; i < 2 ; i++)
	{
		for(int j = 0 ; j < 2 ; j++)
		{
			omega[i][j] =  pow(-1, (i+j))*(1/det_omega_inv)*(omega_inv[(j+1)%3][(i+1)%3]*omega_inv[(j+2)%3][(i+2)%3] - omega_inv[(j+1)%3][(i+2)%3]*omega_inv[(j+2)%3][(i+1)%3]) ;  
		}
	
	}
	double k1l = sigma00[variant][0]*k1/modk + sigma00[variant][5]*k2/modk + sigma00[variant][5]*k3/modk ; 
	double k2l = sigma00[variant][5]*k1/modk + sigma00[variant][1]*k2/modk + sigma00[variant][3]*k3/modk ; 
	double k3l = sigma00[variant][4]*k1/modk + sigma00[variant][3]*k2/modk + sigma00[variant][2]*k3/modk ; 
	
	double j1kl = omega[0][0]*k1l + omega[0][1]*k2l + omega[0][2]*k3l ;
	double j2kl = omega[1][0]*k1l + omega[1][1]*k2l + omega[1][2]*k3l ;
	double j3kl = omega[2][0]*k1l + omega[2][1]*k2l + omega[2][2]*k3l ;
	
	double i1jkl = (k1/modk)*(sigma00[variant][0]*j1kl + sigma00[variant][5]*j2kl + sigma00[variant][4]*j3kl) ;
	double i2jkl = (k2/modk)*(sigma00[variant][5]*j1kl + sigma00[variant][1]*j2kl + sigma00[variant][3]*j3kl) ;
	double i3jkl = (k3/modk)*(sigma00[variant][4]*j1kl + sigma00[variant][3]*j2kl + sigma00[variant][2]*j3kl) ;
	
	double sfts_sum = 0 ;
	for(int i = 0 ; i < 6 ; i++)
	{
		for(int j = 0 ; j < 6 ; j++)
		{
			sfts_sum+= eigen_alpha[variant-1].e[i]*c[i][j]*eigen_alpha[variant-1].e[j] ;
		}
	}
    
    return (sfts_sum - (i1jkl + i2jkl + i3jkl)) ; 
} 

void update_T_dependent_vars(double T)
{
    T= T- dt*cr ; 
    
    G_Al_alpha = (Al_alpha_a_3 + Al_alpha_b_3*T+ Al_alpha_c_3*T*log(T) + Al_alpha_d_3*T*T) ;
    G_Al_beta = (Al_beta_a_3 + Al_beta_b_3*T+ Al_beta_c_3*T*log(T) + Al_beta_d_3*T*T) ;
    G_Ti_alpha = (Ti_alpha_a_2 + Ti_alpha_b_2*T+ Ti_alpha_c_2*T*log(T) + Ti_alpha_d_2*T*T) ;
    G_Ti_beta = (Ti_beta_a_1 + Ti_beta_b_1*T+ Ti_beta_c_1*T*log(T) + Ti_beta_d_1*T*T) ;
    G_V_alpha = (V_alpha_a_2 + V_alpha_b_2*T+ V_alpha_c_2*T*log(T) + V_alpha_d_2*T*T) ; //(21.96 - T/50)
    G_V_beta = (V_beta_a_2 + V_beta_b_2*T+ V_beta_c_2*T*log(T) + V_beta_d_2*T*T) ;	

    w_norm = (R*T)/Vm ;       //Scaling factor for the well height
    epsi_norm = (R*T*pow(lc,2))/Vm;    //Scaling factor for the gradient energy coeffcient
    epsi_prefac = 9.0e-07/epsi_norm ;
    W_prefac = 7.2e+06/w_norm ;   //Overall scaling for the well depth    
    cout<<W_prefac<<endl;
    D = 2.4e-05*exp(-18040/T) ;          
    L_norm = (D*Vm)/(R*T*lc*lc) ;
    L_orig = v_alpha*T_trans/(6*(T_trans - T)*heat*5e-07) ;
    L = L_orig/L_norm;
}

