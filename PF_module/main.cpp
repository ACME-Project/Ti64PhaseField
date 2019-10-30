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


using namespace std;

typedef float phi_type; 											
typedef MMSP::sparse<phi_type> store_type; 	

#include "definitions.hpp"



void initialize()
{
	
	G_Al_alpha = (Al_alpha_a_3 + Al_alpha_b_3*T+ Al_alpha_c_3*T*log(T) + Al_alpha_d_3*T*T) ;
	G_Al_beta = (Al_beta_a_3 + Al_beta_b_3*T+ Al_beta_c_3*T*log(T) + Al_beta_d_3*T*T) ;
	
	G_Ti_alpha = (Ti_alpha_a_2 + Ti_alpha_b_2*T+ Ti_alpha_c_2*T*log(T) + Ti_alpha_d_2*T*T) ;
	G_Ti_beta = (Ti_beta_a_1 + Ti_beta_b_1*T+ Ti_beta_c_1*T*log(T) + Ti_beta_d_1*T*T) ;
	
	G_V_alpha = (V_alpha_a_2 + V_alpha_b_2*T+ V_alpha_c_2*T*log(T) + V_alpha_d_2*T*T) ; //(21.96 - T/50)
	G_V_beta = (V_beta_a_2 + V_beta_b_2*T+ V_beta_c_2*T*log(T) + V_beta_d_2*T*T) ;	

	for(int n = 0; n < nodes(grid); n++)
	{
		MMSP::vector<int> x = position(grid, n);
		
		set(gradsqcal_grid(n), 1) = 0.0 ;
		set(gradsqcv_grid(n), 1) = 0.0 ;
		set(gradxgrid(n), 1) = 0.0 ;
		set(gradygrid(n), 1) = 0.0 ;
		set(gradzgrid(n), 1) = 0.0 ;
		set(gradxgrid(n), 2) = 0.0 ;
		set(gradygrid(n), 2) = 0.0 ;
		set(gradzgrid(n), 2) = 0.0 ;
		set(gradxgrid(n), 3) = 0.0 ;
		set(gradygrid(n), 3) = 0.0 ;
		set(gradzgrid(n), 3) = 0.0 ;
		set(gradxgrid(n), 4) = 0.0 ;
		set(gradygrid(n), 4) = 0.0 ;
		set(gradzgrid(n), 4) = 0.0 ;
		set(gradxgrid(n), 5) = 0.0 ;
		set(gradygrid(n), 5) = 0.0 ;
		set(gradzgrid(n), 5) = 0.0 ;
		
		store_type newGrain;
		
		set(newGrain, 1) = 0.0;
		set(newGrain, 2) = 0.0;
		set(newGrain, 3) = 0.0;
		set(newGrain, 4) = 0.0;
		set(newGrain, 5) = 0.0;
		set(newGrain, 20) = 0.06;
		set(newGrain, 21) = 0.04 ;
		
		grid(n) = newGrain;		
		
	}

}



int main()
{	
	initialize();
	initialize_epsi();
	initialize_alpha_eigen() ;
	define_c_sigma() ;

	MMSP::b0(grid,0) = MMSP::Dirichlet;
	MMSP::b1(grid,0) = MMSP::Dirichlet;

	MMSP::b0(grid,1) = MMSP::Dirichlet;
	MMSP::b1(grid,1) = MMSP::Dirichlet;

	MMSP::b0(grid,2) = MMSP::Dirichlet;
	MMSP::b1(grid,2) = MMSP::Dirichlet;

 
	for(int i = 0 ; i < nx ; i++)
	{
		k[i] = i ;
		fk[i] = 2*3.14*k[i]*fs/nx + 0.001  ;
	}
	
	fftw::maxthreads=get_max_threads();
	

	for(int t = 0 ; t<steps ; t++)
	{

		if((T-dt*tc*cr) > 100.0) 
		{ T= T- dt*tc*cr ; }
		
		
		if(cr!=0.0)
		{
			G_Al_alpha = (Al_alpha_a_3 + Al_alpha_b_3*T+ Al_alpha_c_3*T*log(T) + Al_alpha_d_3*T*T) ;
			G_Al_beta = (Al_beta_a_3 + Al_beta_b_3*T+ Al_beta_c_3*T*log(T) + Al_beta_d_3*T*T) ;
	
			G_Ti_alpha = (Ti_alpha_a_2 + Ti_alpha_b_2*T+ Ti_alpha_c_2*T*log(T) + Ti_alpha_d_2*T*T) ;
			G_Ti_beta = (Ti_beta_a_1 + Ti_beta_b_1*T+ Ti_beta_c_1*T*log(T) + Ti_beta_d_1*T*T) ;
	
			G_V_alpha = (V_alpha_a_2 + V_alpha_b_2*T+ V_alpha_c_2*T*log(T) + V_alpha_d_2*T*T) ; //(21.96 - T/50)
			G_V_beta = (V_beta_a_2 + V_beta_b_2*T+ V_beta_c_2*T*log(T) + V_beta_d_2*T*T) ;
		}

		MMSP::grid<dim, store_type > update(grid);
	
		MMSP::b0(update,0) = MMSP::Dirichlet;
		MMSP::b1(update,0) = MMSP::Dirichlet;

		MMSP::b0(update,1) = MMSP::Dirichlet;
		MMSP::b1(update,1) = MMSP::Dirichlet;

		MMSP::b0(update,2) = MMSP::Dirichlet;
		MMSP::b1(update,2) = MMSP::Dirichlet;
	
		if(t%100==0)
		{
			customoutput(grid, t) ; 		
		}
	
		double deltag = 0.1 ;
		double kappa1 = 0.001 ; 
		double kappa2 = 0.01 ; 
	
		int no_of_cells = (nx*ny*nz)/(3*3*3) ; 
		int cells_in_one_dir = nx/3 ; 
	
		if((t%100==0)&(t<300))   
		// This is not physical.
		{
			for(int n = 0 ; n < nodes(grid) ; n++)
			{
					double jstar = kappa1*exp(-kappa2*deltag);  
					double prob = 1 - exp(-jstar) ;
					double random_no = (double)(rand() % 1000) ;
					double decision_var = random_no/1000.0 ;
					if(prob > decision_var)
					{
						store_type newGrain;
						int random_index = (int)(rand() % 5) + 1 ;
						cout<<"Nucleating index :"<<random_index<<endl;
						
					    MMSP::vector<int> s = position(grid, n) ;
					    MMSP::vector<int> s1 = position(grid, n) ;
					    
					    for(int i = -1 ; i <= 1 ; i++)
					    {
							s[0] += i ;
							for(int j = -1 ; j <= 1 ; j++)
							{
								s[1] += j ;
								for(int k = -1 ; k <= 1 ; k++)
								{
									s[2] += k ;
									if((pow((s[0]*MMSP::dx(grid,0)-s1[0]*MMSP::dx(grid,0)),2) + pow((s[1]*MMSP::dx(grid,1)-s1[1]*MMSP::dx(grid,0)),2) + pow((s[2]*MMSP::dx(grid,2)-s1[2]*MMSP::dx(grid,0)),2)) < 3*MMSP::dx(grid,0)*MMSP::dx(grid,0))
									{
										if(nodesum(grid, s)==0.0)
										{ 
											set(grid(s), random_index) = 1.0; 
										}
									}
									s[2] -= k ;
								}
								s[1] -= j ;
							}
							s[0] -= i ;
						}
					}

						
			}
		}
		
		calculate_strains(grid) ;
	
		
		for(int n = 0 ; n<nodes(grid) ; n++)
		{
			MMSP::vector<int> x = position(grid, n); 
			MMSP::vector<store_type> gradientsqtemp = gradsq(grid, x) ;
			double temp1 = 0.0 ;
			double temp2 = 0.0 ;
			for(int i = 0 ; i < dim ; i++)
			{
					temp1 += gradientsqtemp[i][20] ;
					temp1 += gradientsqtemp[i][21] ;
			}
			set(gradsqcal_grid(n),1) = temp1 ;
			set(gradsqcv_grid(n), 1) = temp2 ;
		
		}
	
	
	
		for(int n = 0 ; n < nodes(grid) ; n++)
		{
			MMSP::vector<int> x = position(grid, n); 
			
			double strain_energy[13];
			strain_energy[0] = 0.0 ;
			strain_energy[1] = (1.0/(double)(nx*ny*nz))*(dfdstr_real1(x[0],x[1],x[2])) ; 
			strain_energy[2] = (1.0/(double)(nx*ny*nz))*(dfdstr_real2(x[0],x[1],x[2])) ; 
			strain_energy[3] = (1.0/(double)(nx*ny*nz))*(dfdstr_real3(x[0],x[1],x[2])) ; 
			strain_energy[4] = (1.0/(double)(nx*ny*nz))*(dfdstr_real4(x[0],x[1],x[2])) ; 
			strain_energy[5] = (1.0/(double)(nx*ny*nz))*(dfdstr_real5(x[0],x[1],x[2])) ; 
			
			
			G_Alpha[n] = ((grid(n)[20]*G_Al_alpha + grid(n)[21]*G_V_alpha + (1-grid(n)[20]-grid(n)[21])*G_Ti_alpha +
					 R*T*(grid(n)[20]*log(grid(n)[20]) + grid(n)[21]*log(grid(n)[21]) + (1-grid(n)[20]-grid(n)[21])*log((1-grid(n)[20]-grid(n)[21]))) +
					 grid(n)[20]*grid(n)[21]*L0_HCP_Al_V + grid(n)[20]*(1-grid(n)[20]-grid(n)[21])*L0_HCP_Al_Ti + grid(n)[21]*(1-grid(n)[20]-grid(n)[21])*L0_HCP_Ti_V))/G_normalize ;
					 
			G_Beta[n] = (grid(n)[20]*G_Al_beta + grid(n)[21]*G_V_beta + (1-grid(n)[20]-grid(n)[21])*G_Ti_beta +
					 R*T*(grid(n)[20]*log(grid(n)[20]) + grid(n)[21]*log(grid(n)[21]) + (1-grid(n)[20]-grid(n)[21])*log((1-grid(n)[20]-grid(n)[21]))) +
					 grid(n)[20]*grid(n)[21]*L0_BCC_Al_V + grid(n)[20]*(1-grid(n)[20]-grid(n)[21])*L0_BCC_Al_Ti + grid(n)[21]*(1-grid(n)[20]-grid(n)[21])*L0_BCC_Ti_V +
					 grid(n)[20]*grid(n)[21]*(1-grid(n)[20]-grid(n)[21])*L0_BCC_Al_Ti_V)/G_normalize ;
					 
			W[n] = W_prefac*(grid(n)[20]*W_Al + grid(n)[21]*W_V + (1-grid(n)[20]-grid(n)[21])*W_Ti) ;
		
		
			MMSP::vector<store_type> gradient = grad(grid, x) ;
			MMSP::vector<store_type> gradientsq = gradsq(grid, x) ;
			MMSP::sparse<phi_type> dFdp;				
			MMSP::sparse<phi_type> lap = laplacian(grid, n); 
		
			thermo_auxillary_terms(gradient, gradientsq, grid(n)[20], grid(n)[21]) ;

			MMSP::vector<store_type> gradcalp4 = gradsq(gradsqcal_grid, x) ;
			MMSP::vector<store_type> gradcvp4 = gradsq(gradsqcv_grid, x) ;
		
			double gradp4_cal = 0 ;
			double gradp4_cv  = 0 ;
		
			for(int i = 0 ; i < dim ; i++)
			{
				gradp4_cal += gradcalp4[i][1]  ;
				gradp4_cv += gradcvp4[i][1] ;
			}
				
			phi_type phi_grad = 0 ;
			phi_type phi_gradientsq = 0 ;
			phi_type phi_sq = 0.0; 	
			phi_type phi_sum = 0.0 ;						
			phi_type dFall = 0.0;
			phi_type dFall_no = 0 ;
			for(int h = 0 ; h < length(grid(x)) ; h++)
			{
				int hindex = MMSP::index(grid(x),h) ;
				if((hindex <= 12)&&(grid(n)[hindex]>1e-02))
				{
				
					for(int i = 0 ; i < dim ; i++)
					{
						phi_grad+= gradient[i][hindex] ; 
						phi_gradientsq+= gradientsq[i][hindex] ; 
						phi_sq += grid(n)[jindex]*grid(n)[jindex];
						phi_sum += grid(n)[jindex] ;
					}
				}
			}
			
			phi_beta = 1 - phi_sum ; 
		
			hphidoubleprime[n] = 0.0 ; 
			hphi[n] = 0.0 ;
			hphiprime[n] = 0.0 ;
		
			gphidoubleprime[n] = 0.0 ; 
			gphiprime[n] = 0.0 ;
			gphi[n] = 0.0 ;
		
			for(int h = 0 ; h < length(grid(x)) ; h++)
			{
				int hindex = MMSP::index(grid(x),h) ;
				if((hindex <= 12))
				{
					hphidoubleprime[n] +=  2*pow(grid(n)[hindex],1) - 3*pow(grid(n)[hindex],2) ; 
					hphi[n] += pow(grid(n)[hindex],3)/3 - pow(grid(n)[hindex],4)/4 ;
					hphiprime[n] += pow(grid(n)[hindex],2) - pow(grid(n)[hindex],3) ;
				
					gphidoubleprime[n] += (2*(6*pow(grid(n)[hindex],2) - 6*grid(n)[hindex] + 1)) ; 
					gphiprime[n] += (2*grid(n)[hindex]*(1-grid(n)[hindex])*(1-2*grid(n)[hindex])) ;
					gphi[n] += pow(grid(n)[hindex],2)*pow((1-grid(n)[hindex]),2) ;
				}
			}
													 
			double c_al_rhs = 2*(del_dGAlpha_dAl - del_dGBeta_dAl)*hphiprime[n]*phi_grad + delsq_dGAlpha_dAl*hphi[n] + delsq_dGBeta_dAl*(1-hphi[n]) + 
						  (dGAlpha_dAl - dGBeta_dAl)*(hphidoubleprime[n]*pow(phi_grad,2) + hphiprime[n]*phi_gradientsq) + (W_Al - W_Ti)*(gphiprime[n]*phi_gradientsq + gphidoubleprime[n]*phi_grad) ;
			double c_v_rhs = 2*(del_dGAlpha_dV - del_dGBeta_dV)*hphiprime[n]*phi_grad + delsq_dGAlpha_dV*hphi[n] + delsq_dGBeta_dV*(1-hphi[n]) + 
						  (dGAlpha_dV - dGBeta_dV)*(hphidoubleprime[n]*pow(phi_grad,2) + hphiprime[n]*phi_gradientsq) + (W_V - W_Ti)*(gphiprime[n]*phi_gradientsq + gphidoubleprime[n]*phi_grad)  ;
						  

			
			set(update(n), 20) = grid(n)[20] + dt*(Dalal*c_al_rhs + Dalv*c_v_rhs - 0.00076*gradp4_cal) ;
			set(update(n), 21) = grid(n)[21] + dt*(Dvv*c_v_rhs + Dval*c_al_rhs -  0.00076*gradp4_cv) ; 
			
			
			double lap_aniso[13];				
			
			for (int j = 0; j < length(grid(n)); j++)   
			{
				int jindex = MMSP::index(grid(n), j);  		
					
				MMSP::vector<int> x = position(grid, n);
					
				W[n] = 0.01432 ; 
				if(jindex <= 12)
				{
					lap_aniso[jindex] = 0.0 ;
					for(int p = 0 ; p < dim ; p++)
					{
						for(int q = 0 ; q < dim ; q++)
						{
							lap_aniso[jindex]+= (6.8*pow(10,-6)/epsi_norm)*(gradient[p][jindex])*epsi[jindex-1][p][q]*(gradient[q][jindex]);   //
						}
					}
					
					double coeff ;
					if(lap[jindex]==0) coeff = 0.0 ;
					set(dFdp, jindex) = (G_Alpha[n] - G_Beta[n])*(30*pow(grid(n)[jindex],2)*pow((grid(n)[jindex]-1),2))  + 0.01*(2*grid(n)[jindex]*(1-grid(n)[jindex])*(1-2*grid(n)[jindex]))  + strain_energy[jindex]  - coeff*(lap_aniso[jindex]); 
					
					dFall += dFdp[jindex];	
					dFall_no+=1;
				}
										
			}
			
			double L = L_orig*(T/T_orig) ;

			for (int h = 0; h < length(grid(n)); h++) 
			{	
					
				int hindex = MMSP::index(grid(n), h);   
				if(hindex<=12)
				{
					store_type dpdt;
					set(dpdt, hindex) = -L * ((dFall_no+1)*dFdp[hindex] - dFall);
					set(update(n), hindex) = grid(n)[hindex] + dt * dpdt[hindex];	
				}
			}
			
		}

		swap(grid, update);	

}



return 0 ;

}







