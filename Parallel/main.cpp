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
#include <time.h>
#include<deque>
#include<mpi.h>     //Relevant include for parallel implementation.
#include <stdio.h>
#include <stdlib.h>
#include<vector>
#include <complex>
#include <valarray>
 
#include"MMSP.hpp"


using namespace std;

typedef float phi_type; 											
typedef MMSP::sparse<phi_type> store_type; 	

#include "definitions_v2.hpp"


int main(int argc, char* argv[])
{	
    
    //Note: A separate initialization function, like in the serial mode, has not been written because grid needs to be 'defined' after MPI_Init and MPI_Init needs to be inside main. 
    
    initialize_epsi();
	initialize_alpha_eigen() ;
	//define_c_sigma() ;
	
	MPI_Init(&argc, &argv);
	
	MMSP::grid<dim, store_type> grid(variants, 0, nx, 0, ny) ;
	MMSP::grid<dim, store_type> gradsqcal_grid(variants, 0, nx, 0, ny) ;
	MMSP::grid<dim, store_type> gradsqcv_grid(variants, 0, nx, 0, ny) ;
    MMSP::grid<dim, store_type> intenergies(13, 0, nx, 0, ny) ;
    std::vector<int> variant_ids ;
	
	int nx_local = x1(grid) - x0(grid);
	int ny_local = y1(grid) - y0(grid);

	const double Lx = g1(grid,0) - g0(grid,0) ;
	const double Ly = g1(grid,1) - g0(grid,1) ;
	
	const double dx = MMSP::dx(grid, 0) ;
	const double dy = MMSP::dx(grid, 1) ;
	
    //Confirming that domain has been split between processes.
	int world_size, world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    cout<<"Hello world from rank:"<< world_rank<<", out of processors: "<<world_size<<" with local limits as: "<<x0(grid)<<","<<x1(grid)<<" "<<y0(grid)<<","<<y1(grid)<<endl;
	const int node_total = nx_local*ny_local ; 
	double G_Alpha[node_total], G_Beta[node_total], W[node_total] ; 

	//-------Setting the boundary conditions for various grids-------//
    
	MMSP::b0(grid,0) = MMSP::parallel;
	MMSP::b1(grid,0) = MMSP::parallel;

	MMSP::b0(grid,1) = MMSP::parallel;
	MMSP::b1(grid,1) = MMSP::parallel;
	
	MMSP::b0(gradsqcal_grid,0) = MMSP::parallel;
	MMSP::b1(gradsqcal_grid,0) = MMSP::parallel;

	MMSP::b0(gradsqcal_grid,1) = MMSP::parallel;
	MMSP::b1(gradsqcal_grid,1) = MMSP::parallel;

	MMSP::b0(gradsqcv_grid,0) = MMSP::parallel;
	MMSP::b1(gradsqcv_grid,0) = MMSP::parallel;

	MMSP::b0(gradsqcv_grid,1) = MMSP::parallel;
	MMSP::b1(gradsqcv_grid,1) = MMSP::parallel;
    
    
	//Initialization of thermodynamic terms.
    G_Al_alpha = (Al_alpha_a_3 + Al_alpha_b_3*T+ Al_alpha_c_3*T*log(T) + Al_alpha_d_3*T*T) ;
	G_Al_beta = (Al_beta_a_3 + Al_beta_b_3*T+ Al_beta_c_3*T*log(T) + Al_beta_d_3*T*T) ;
	
	G_Ti_alpha = (Ti_alpha_a_2 + Ti_alpha_b_2*T+ Ti_alpha_c_2*T*log(T) + Ti_alpha_d_2*T*T) ;
	G_Ti_beta = (Ti_beta_a_1 + Ti_beta_b_1*T+ Ti_beta_c_1*T*log(T) + Ti_beta_d_1*T*T) ;
	
	G_V_alpha = (V_alpha_a_2 + V_alpha_b_2*T+ V_alpha_c_2*T*log(T) + V_alpha_d_2*T*T) ; 
	G_V_beta = (V_beta_a_2 + V_beta_b_2*T+ V_beta_c_2*T*log(T) + V_beta_d_2*T*T) ;	

    //Grid initialization. 
	for(int n = 0; n < nodes(grid); n++)
	{
		MMSP::vector<int> x = position(grid, n);
		store_type newGrain;
		MMSP::set(gradsqcal_grid(n), 1) = 0.0 ;
		MMSP::set(gradsqcv_grid(n), 1) = 0.0 ;
        MMSP::set(gradsqcal_grid(n), 2) = 0.0 ;
		MMSP::set(gradsqcv_grid(n), 2) = 0.0 ;
        
        //store_type newGrain;
        MMSP::set(newGrain, 1) = 0.0;
        MMSP::set(newGrain, 3) = 0.0;
        MMSP::set(newGrain, 4) = 0.0;
        MMSP::set(newGrain, 5) = 0.0;
        MMSP::set(newGrain, 15) = 0.0 ;
		MMSP::set(newGrain, 20) = 0.10 ;
		MMSP::set(newGrain, 21) = 0.036 ;

        
        if(x[0]==1 & x[1] == 1) MMSP::set(newGrain, 2) = 1.0 ;  
        else if(x[0]==2 & x[1] == (nx/50)*8) MMSP::set(newGrain, 2) = 1.0 ; 
        else if(x[0]==2 & x[1] == 12*(nx/50)) MMSP::set(newGrain, 2) = 1.0 ; 
        else if(x[0]==2 & x[1] == 16*(nx/50)) MMSP::set(newGrain, 2) = 1.0 ; 
        else if(x[0]==2 & x[1] == 20*(nx/50)) MMSP::set(newGrain, 2) = 1.0 ; 
        else if(x[0]==2 & x[1] == 24*(nx/50)) MMSP::set(newGrain, 2) = 1.0 ; 
        else if(x[0]==2 & x[1] == 32*(nx/50)) MMSP::set(newGrain, 2) = 1.0 ; 
        else if(x[0]==2 & x[1] == 36*(nx/50)) MMSP::set(newGrain, 2) = 1.0 ; 
        else MMSP::set(newGrain, 2) = 0.0 ; 

		
		MMSP::set(gradsqcal_grid(n),1) = 0.0 ;
        MMSP::set(gradsqcv_grid(n), 1) = 0.0 ;
        MMSP::set(gradsqcal_grid(n),2) = 0.0 ;
        MMSP::set(gradsqcv_grid(n), 2) = 0.0 ;
        
		grid(n) = newGrain;
		
	}

	
    srand(time(NULL));
    
    for(int t = 0 ; t<steps ; t++)
	{
	
        if(T <= 700.0)
        {
            return 0 ;
        }

        // Updating temperature dependent variables.
		if(cr > 0.0 and (T>700))    //Implementing the variation of Temperature based on the given cooling rate
		{ 
            T= T- dt*cr ; 
            G_Al_alpha = (Al_alpha_a_3 + Al_alpha_b_3*T+ Al_alpha_c_3*T*log(T) + Al_alpha_d_3*T*T) ;
            G_Al_beta = (Al_beta_a_3 + Al_beta_b_3*T+ Al_beta_c_3*T*log(T) + Al_beta_d_3*T*T) ;
	
            G_Ti_alpha = (Ti_alpha_a_2 + Ti_alpha_b_2*T+ Ti_alpha_c_2*T*log(T) + Ti_alpha_d_2*T*T) ;
            G_Ti_beta = (Ti_beta_a_1 + Ti_beta_b_1*T+ Ti_beta_c_1*T*log(T) + Ti_beta_d_1*T*T) ;
	
            G_V_alpha = (V_alpha_a_2 + V_alpha_b_2*T+ V_alpha_c_2*T*log(T) + V_alpha_d_2*T*T) ; 
            G_V_beta = (V_beta_a_2 + V_beta_b_2*T+ V_beta_c_2*T*log(T) + V_beta_d_2*T*T) ;	
            w_norm = (R*T)/Vm ;       //Scaling factor for the well height
            epsi_norm = (R*T*pow(lc,2))/Vm;    //Scaling factor for the gradient energy coeffcient
            epsi_prefac = 9.0e-07/epsi_norm ;
            W_prefac = 7.2e+05/w_norm ;   //Overall scaling for the well depth  
            Dnorm = 2.4e-5*exp(-18040.0/T) ; //0.001 ;        // The next four lines are the components of the Onsager mobility matrix
            Dalal = 2.4e-5*exp(-18040.0/T)/(100*Dnorm) ; 
            Dvv = (1.0e-5*exp(-17460.0/T))/(100*Dnorm) ; 
            L_norm = (Dnorm*Vm)/(R*T*lc*lc) ;
            L_orig = fp/(T_trans-T) ; //v_alpha*T_trans/((T_trans - T)*heat*5e-07) ;
            L = L_orig/L_norm; //exp(a-b/T) ; //  
        }        
       

        for(int n = 0 ; n<nodes(grid) ; n++)
        {
			MMSP::vector<int> x = position(grid, n); 
			MMSP::vector<store_type> gradientsqtemp = gradsq(grid, x) ;
			double temp1 = 0.0 ;
			double temp2 = 0.0 ;
            MMSP::set(gradsqcal_grid(n),1) = gradientsqtemp[0][20] ;
            MMSP::set(gradsqcal_grid(n),2) = gradientsqtemp[1][20] ;
            
			MMSP::set(gradsqcv_grid(n),1) = gradientsqtemp[0][21] ;
            MMSP::set(gradsqcv_grid(n),2) = gradientsqtemp[1][21] ;	
		}
            
        // Generating output
        if(t%500==0)
		{
			std::string file_name = "output_" + to_string(t) + ".dat" ;
			char* temp_array = new char[100 + 1];
			strcpy(temp_array, file_name.c_str()) ; ;
			MMSP::output(grid, temp_array);		
			delete [] temp_array;
		}
        
        if(t%6000==0)    //Example implementation of boundary nucleation
        {
           double deltag = 0.0001 ;
           double kappa1 = 0.001 ; 
           double kappa2 = 0.001 ; 
            
            for(int n = 0 ; n < nodes(grid) ; n++)
            {
                MMSP::vector<int> x = position(grid, n);
                
                if((x[0]!=1 & x[0]!=nx-2) & (x[1]!=1 & x[1]!=ny-2)) continue ;
                
                double phi_sum = 0.0 ;
                for(int h = 0 ; h < length(grid(n)) ; h++)
                {
                    int hindex = MMSP::index(grid(n),h) ;
                    if(hindex <= 12)
                    {
                        phi_sum+=grid(n)[hindex] ;
                    }
                }
                if(phi_sum > 0.1) continue ;
                double prob1 ; 
                double df = gdiff(grid(n)[20], grid(n)[21], T) ;
                double jstar = kappa1*exp(-kappa2*deltag) ;
                double prob = 1 - exp(-jstar) ;
                double random_no = (double)(rand() % 1000) ;
                double decision_var = random_no/1000.0 ;
                if(prob > decision_var)
                {
                    int nuc_index ;
                    if(x[0]==1 || x[0]==nx-2) nuc_index = 2 ; //(int)(rand()%2) + 1 ;
                    else if(x[1]==1 || x[1]==ny-2) nuc_index = 1 ; 
                    cout<<"Nucleating index: "<<nuc_index<<" at position "<<x[0]<<" "<<x[1]<<endl;
                    variant_ids.push_back(nuc_index);
                    MMSP::set(grid(n), nuc_index) = 1.0; 
                    MMSP::set(grid(n), 15) = 1.0 ;
                    if(variant_ids.size() > 75) break ;
                }
            }
        }
        
       if(t%100==0 & t >= 200 & t <= 600)  //In-grain nucleation
        {
           double kappa1 = 0.01 ; 
           double kappa2 = (16*3.14*(pow(0.02,3)))/(3*1.38e-23*T) ;
            
            for(int n = 0 ; n < nodes(grid) ; n++)
            {
                MMSP::vector<int> x = position(grid, n);
                
                if((x[0]==1 || x[0]==nx-2) || (x[1]==1 || x[1]==ny-2)) continue ;
                
                double phi_sum = 0.0 ;
                for(int h = 0 ; h < length(grid(n)) ; h++)
                {
                    int hindex = MMSP::index(grid(n),h) ;
                    if(hindex <= 12)
                    {
                        phi_sum+=grid(n)[hindex] ;
                    }
                }
                if(phi_sum > 0.1) continue ;
                double df = (gdiff(grid(n)[20], grid(n)[21], T)*8.314*T)/1.0e-5 ;
                double jstar = kappa1*exp(-kappa2/(pow(df,2))) ;
                double prob = 1 - exp(-jstar) ;
                double random_no = (double)(rand() % 1000) ;
                double decision_var = random_no/1000.0 ;
                if(prob > decision_var)
                {
                    int nuc_index ;
                    nuc_index = (int)(rand()%5) + 1 ;
                    cout<<"Nucleating index: "<<nuc_index<<" at position "<<x[0]<<" "<<x[1]<<endl;
                    variant_ids.push_back(nuc_index);
                    MMSP::set(grid(n), nuc_index) = 1.0; 
                }
            }
        }

        // Update loop. In the absence of documentation, please refer to the main file in the serial folder. The following implementation is identical to what is found there. 
		for(int n = 0 ; n < nodes(grid) ; n++)
		{
			MMSP::vector<int> x = position(grid, n); 
			
			G_Alpha[n] = ((grid(n)[20]*G_Al_alpha + grid(n)[21]*G_V_alpha + (1-grid(n)[20]-grid(n)[21])*G_Ti_alpha +
					 R*T*(grid(n)[20]*log(grid(n)[20]) + grid(n)[21]*log(grid(n)[21]) + (1-grid(n)[20]-grid(n)[21])*log((1-grid(n)[20]-grid(n)[21]))) +
					 grid(n)[20]*grid(n)[21]*L0_HCP_Al_V + grid(n)[20]*(1-grid(n)[20]-grid(n)[21])*(L0_HCP_Al_Ti + L1_HCP_Al_Ti*(2*grid(n)[20] + grid(n)[21] - 1)) + grid(n)[21]*(1-grid(n)[20]-grid(n)[21])*L0_HCP_Ti_V))/G_normalize ;
					 
			G_Beta[n] = (grid(n)[20]*G_Al_beta + grid(n)[21]*G_V_beta + (1-grid(n)[20]-grid(n)[21])*G_Ti_beta +
					 R*T*(grid(n)[20]*log(grid(n)[20]) + grid(n)[21]*log(grid(n)[21]) + (1-grid(n)[20]-grid(n)[21])*log((1-grid(n)[20]-grid(n)[21]))) +
					 grid(n)[20]*grid(n)[21]*L0_BCC_Al_V + grid(n)[20]*(1-grid(n)[20]-grid(n)[21])*L0_BCC_Al_Ti + grid(n)[21]*(1-grid(n)[20]-grid(n)[21])*L0_BCC_Ti_V +
					 grid(n)[20]*grid(n)[21]*(1-grid(n)[20]-grid(n)[21])*L0_BCC_Al_Ti_V)/G_normalize ;
					 
            W[n] = W_prefac*((0.10-grid(n)[20])*50*W_Al + (grid(n)[21]-0.036)*50*W_V) ;
			MMSP::vector<store_type> gradient = grad(grid, x) ;
			MMSP::vector<store_type> gradientsq = gradsq(grid, x) ;
            		
			thermo_auxillary_terms(gradient, gradientsq, grid(n)[20], grid(n)[21]) ;

			MMSP::vector<store_type> gradcalp4 = gradsq(gradsqcal_grid, x) ;	
			MMSP::vector<store_type> gradcvp4 = gradsq(gradsqcv_grid, x) ;
		
            double gradp4_cal = gradcalp4[0][1] + gradcalp4[1][2] ;
            double gradp4_cv  = gradcvp4[0][1] + gradcvp4[1][2] ;
		
            double hphi[12], hphiprime[12], hphidoubleprime[12];
            double gphi[12], gphiprime[12], gphidoubleprime[12];
		
			for(int h = 0 ; h < length(grid(x)) ; h++)
			{
				int hindex = MMSP::index(grid(x),h) ;
				if(hindex <= 12)
				{
					hphidoubleprime[hindex-1] = 60*grid(n)[hindex]*(2*grid(n)[hindex]-1)*(grid(n)[hindex]-1) ; 
					hphi[hindex-1] = pow(grid(n)[hindex],3)*(6*pow(grid(n)[hindex],2) - 15*grid(n)[hindex] + 10) ;
					hphiprime[hindex-1] = 30*pow(grid(n)[hindex],2)*pow((grid(n)[hindex]-1),2) ;
		
					gphidoubleprime[hindex-1] = 2*(6*pow(grid(n)[hindex],2) - 6*grid(n)[hindex] + 1) ; 
					gphiprime[hindex-1] = 2*grid(n)[hindex]*(1-grid(n)[hindex])*(1-2*grid(n)[hindex]) ;
					gphi[hindex-1] = pow(grid(n)[hindex],2)*pow((1-grid(n)[hindex]),2) ;
				}
			}
            
            double cal = grid(n)[20] ;
            double cv = grid(n)[21] ;
            double grad_cal = gradient[0][20] + gradient[1][20] ;
            double grad_cv = gradient[0][21] + gradient[1][21] ;
            
            double hphidoubleprimesum = 0;
            double hphiprimesum1 = 0;
            double hphiprimesum2 = 0;
            double hphisum = 0;
            double gphidoubleprimesum = 0;
            double gphiprimesum = 0;
            double gphisum = 0;
	
            for(int h = 0 ; h < length(grid(x)) ; h++)
            {
                int hindex = MMSP::index(grid(x),h) ;
                if(hindex <= 12)
                {
                    double delphi = gradient[0][hindex] + gradient[1][hindex] ;
                    double delphisq = gradientsq[0][hindex] + gradientsq[1][hindex]  ;
                    hphidoubleprimesum += hphidoubleprime[hindex-1]*pow(delphi,2) ;
                    hphiprimesum1 += hphiprime[hindex-1]*delphisq;
                    hphiprimesum2 += hphiprime[hindex-1]*delphi ;
                    hphisum += hphi[hindex-1] ;
				
                    gphiprimesum += gphiprime[hindex-1]*delphisq ;
                    gphidoubleprimesum += gphidoubleprime[hindex-1]*delphi ;
                }
            }
               
            
            double c_al_rhs = 2*(del_dGAlpha_dAl - del_dGBeta_dAl)*hphiprimesum2 + delsq_dGAlpha_dAl*hphisum + delsq_dGBeta_dAl*(1-hphisum) + 
					  (dGAlpha_dAl - dGBeta_dAl)*(hphidoubleprimesum + hphiprimesum1) + (W_Al - W_Ti)*(gphiprimesum + gphidoubleprimesum) ;
            double c_v_rhs = 2*(del_dGAlpha_dV - del_dGBeta_dV)*hphiprimesum2 + delsq_dGAlpha_dV*hphisum + delsq_dGBeta_dV*(1-hphisum) + 
					  (dGAlpha_dV - dGBeta_dV)*(hphidoubleprimesum + hphiprimesum1) + (W_V - W_Ti)*(gphiprimesum + gphidoubleprimesum)  ;


            if((grid(n)[20] + dt*(Dalal*(c_al_rhs)- kappa_c*gradp4_cal)) < 0.005)
            {
                MMSP::set(grid(n), 20) = 0.005 ;
            }
            else
            {
                MMSP::set(grid(n), 20) = grid(n)[20] + dt*(Dalal*(c_al_rhs)- kappa_c*gradp4_cal) ; 
            }
            if((grid(n)[21] + dt*(Dvv*(c_v_rhs) - kappa_c*gradp4_cv)) < 0.005)
            {
                MMSP::set(grid(n), 21) = 0.005 ;
            }
            else
            {
                MMSP::set(grid(n), 21) = grid(n)[21] + dt*(Dvv*(c_v_rhs) - kappa_c*gradp4_cv) ; 
            }
		
    

			
			MMSP::sparse<phi_type> dFdp;
			phi_type dFall = 0.0;		
			phi_type phi_sum = 0.0 ;
			
			for (int j = 0; j < length(grid(n)); j++) 
			{
				int jindex = MMSP::index(grid(n), j);
				if(jindex<=12)
				{		
					phi_sum += grid(n)[jindex] ;
				}
			}
				
			phi_type phi_beta = 1 - phi_sum ; 
			
			for (int j = 0; j < length(grid(n)); j++)   
			{
				int jindex = MMSP::index(grid(n), j);  		
				MMSP::vector<int> x = position(grid, n);
                W[n] = W_prefac*((0.10-grid(n)[20])*W_Al + (grid(n)[21]-0.036)*W_V) ;
				if(jindex <= 12)
				{
					double lap_aniso = 0.0  ;
                    for(int i = 0 ; i < dim ; i++)
                    {
                        lap_aniso += epsi_prefac*gradientsq[i][jindex]*epsi[jindex-1][i][i] ;
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
							interaction_energy+= 10.0*pow(grid(n)[tempindex],2)*grid(n)[jindex];
							interaction_energy-= W[n]*pow(grid(n)[tempindex],2)*(phi_beta); 
						}
					}
					double df = gdiff(grid(n)[20], grid(n)[21], T);
					MMSP::set(dFdp, jindex) = -(df)*hphiprime[jindex-1] + interaction_energy - lap_aniso ;    //+ strain_energy[jindex] 
					dFall += dFdp[jindex];	
				}
			}

			for (int h = 0; h < length(grid(n)); h++) 
			{	
					
				int hindex = MMSP::index(grid(n), h);   
				if(hindex<=12)
				{
					store_type dpdt;
					MMSP::set(dpdt, hindex) = -L * (dFdp[hindex]);
					MMSP::set(grid(n), hindex) = grid(n)[hindex] + dt * dpdt[hindex];	
				}
			}				
		}

	//ghostswap(update);
	//swap(grid, update);	
    ghostswap(grid);

}

#ifdef MPI_VERSION
	MPI_Finalize();
#endif
return 0 ;

}







