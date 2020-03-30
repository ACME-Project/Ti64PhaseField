void initialize()
{
	
	G_Al_alpha = (Al_alpha_a_3 + Al_alpha_b_3*T+ Al_alpha_c_3*T*log(T) + Al_alpha_d_3*T*T) ;
	G_Al_beta = (Al_beta_a_3 + Al_beta_b_3*T+ Al_beta_c_3*T*log(T) + Al_beta_d_3*T*T) ;
	
	G_Ti_alpha = (Ti_alpha_a_2 + Ti_alpha_b_2*T+ Ti_alpha_c_2*T*log(T) + Ti_alpha_d_2*T*T) ;
	G_Ti_beta = (Ti_beta_a_1 + Ti_beta_b_1*T+ Ti_beta_c_1*T*log(T) + Ti_beta_d_1*T*T) ;
	
	G_V_alpha = (V_alpha_a_2 + V_alpha_b_2*T+ V_alpha_c_2*T*log(T) + V_alpha_d_2*T*T) ; //(21.96 - T/50)
	G_V_beta = (V_beta_a_2 + V_beta_b_2*T+ V_beta_c_2*T*log(T) + V_beta_d_2*T*T) ;	

    variant_ids.push_back(2) ;
    variant_ids.push_back(2) ;
    variant_ids.push_back(2) ;
    variant_ids.push_back(2) ;
    variant_ids.push_back(2) ;
    variant_ids.push_back(2) ;
    variant_ids.push_back(2) ;
    variant_ids.push_back(2) ;
	for(int n = 0; n < nodes(grid); n++)
	{
		MMSP::vector<int> x = position(grid, n);
        //cout<<x[0]<<" "<<x[1]<<endl;
        
        store_type newGrain;
        
        
        if(pow((x[0]*MMSP::dx(grid,0) - 0.98*Lx),2) + pow((x[1]*MMSP::dx(grid,1) - 0.8*Ly),2) + pow((x[2]*MMSP::dx(grid,2) - 0.5*Lz),2) <= 8*MMSP::dx(grid,0)*MMSP::dx(grid,0))
		{
            set(newGrain, 1) = 0.0;
            //set(newGrain, 3) = 0.0;
			set(newGrain, 2) = 1.0;
            set(newGrain, 3) = 0.0;
            set(newGrain, 4) = 0.0;
            set(newGrain, 5) = 0.0;
			set(newGrain, 20) = 0.10;
			set(newGrain, 21) = 0.036 ;
            
		}
        
        /*else if(pow((x[0]*MMSP::dx(grid,0) - 0.98*Lx),2) + pow((x[1]*MMSP::dx(grid,1) - 0.7*Ly),2) + pow((x[2]*MMSP::dx(grid,2) - 0.5*Lz),2) <= 8*MMSP::dx(grid,0)*MMSP::dx(grid,0))
		{
            set(newGrain, 1) = 0.0;
            //set(newGrain, 3) = 0.0;
			set(newGrain, 2) = 1.0;
            set(newGrain, 3) = 0.0;
            set(newGrain, 4) = 0.0;
            set(newGrain, 5) = 0.0;
			set(newGrain, 20) = 0.10;
			set(newGrain, 21) = 0.036 ;
            //variant_ids.push_back(2) ;
		}*/
        
        else if(pow((x[0]*MMSP::dx(grid,0) - 0.98*Lx),2) + pow((x[1]*MMSP::dx(grid,1) - 0.6*Ly),2) + pow((x[2]*MMSP::dx(grid,2) - 0.5*Lz),2) <= 8*MMSP::dx(grid,0)*MMSP::dx(grid,0))
		{
            set(newGrain, 1) = 0.0;
            //set(newGrain, 3) = 0.0;
			set(newGrain, 2) = 1.0;
            set(newGrain, 3) = 0.0;
            set(newGrain, 4) = 0.0;
            set(newGrain, 5) = 0.0;
			set(newGrain, 20) = 0.10;
			set(newGrain, 21) = 0.036 ;
            //variant_ids.push_back(2) ;
		}
        
        /*else if(pow((x[0]*MMSP::dx(grid,0) - 0.98*Lx),2) + pow((x[1]*MMSP::dx(grid,1) - 0.5*Ly),2) + pow((x[2]*MMSP::dx(grid,2) - 0.5*Lz),2) <= 8*MMSP::dx(grid,0)*MMSP::dx(grid,0))
		{
            set(newGrain, 1) = 0.0;
            //set(newGrain, 3) = 0.0;
			set(newGrain, 2) = 1.0;
            set(newGrain, 3) = 0.0;
            set(newGrain, 4) = 0.0;
            set(newGrain, 5) = 0.0;
			set(newGrain, 20) = 0.10;
			set(newGrain, 21) = 0.036 ;
            //variant_ids.push_back(2) ;
		}*/
        
        else if(pow((x[0]*MMSP::dx(grid,0) - 0.98*Lx),2) + pow((x[1]*MMSP::dx(grid,1) - 0.4*Ly),2) + pow((x[2]*MMSP::dx(grid,2) - 0.5*Lz),2) <= 8*MMSP::dx(grid,0)*MMSP::dx(grid,0))
		{
            set(newGrain, 1) = 0.0;
            //set(newGrain, 3) = 0.0;
			set(newGrain, 2) = 1.0;
            set(newGrain, 3) = 0.0;
            set(newGrain, 4) = 0.0;
            set(newGrain, 5) = 0.0;
			set(newGrain, 20) = 0.10;
			set(newGrain, 21) = 0.036 ;
            //variant_ids.push_back(2) ;
		}
        
        /*else if(pow((x[0]*MMSP::dx(grid,0) - 0.98*Lx),2) + pow((x[1]*MMSP::dx(grid,1) - 0.3*Ly),2) + pow((x[2]*MMSP::dx(grid,2) - 0.5*Lz),2) <= 8*MMSP::dx(grid,0)*MMSP::dx(grid,0))
		{
            set(newGrain, 1) = 0.0;
            //set(newGrain, 3) = 0.0;
			set(newGrain, 2) = 1.0;
            set(newGrain, 3) = 0.0;
            set(newGrain, 4) = 0.0;
            set(newGrain, 5) = 0.0;
			set(newGrain, 20) = 0.10;
			set(newGrain, 21) = 0.036 ;
            //variant_ids.push_back(2) ;
		}*/
        
        else if(pow((x[0]*MMSP::dx(grid,0) - 0.98*Lx),2) + pow((x[1]*MMSP::dx(grid,1) - 0.2*Ly),2) + pow((x[2]*MMSP::dx(grid,2) - 0.5*Lz),2) <= 8*MMSP::dx(grid,0)*MMSP::dx(grid,0))
		{
            set(newGrain, 1) = 0.0;
            //set(newGrain, 3) = 0.0;
			set(newGrain, 2) = 1.0;
            set(newGrain, 3) = 0.0;
            set(newGrain, 4) = 0.0;
            set(newGrain, 5) = 0.0;
			set(newGrain, 20) = 0.10;
			set(newGrain, 21) = 0.036 ;
            //variant_ids.push_back(2) ;
		}
        
        /*else if(pow((x[0]*MMSP::dx(grid,0) - 0.98*Lx),2) + pow((x[1]*MMSP::dx(grid,1) - 0.1*Ly),2) + pow((x[2]*MMSP::dx(grid,2) - 0.5*Lz),2) <= 8*MMSP::dx(grid,0)*MMSP::dx(grid,0))
		{
            set(newGrain, 1) = 0.0;
            //set(newGrain, 3) = 0.0;
			set(newGrain, 2) = 1.0;
            set(newGrain, 3) = 0.0;
            set(newGrain, 4) = 0.0;
            set(newGrain, 5) = 0.0;
			set(newGrain, 20) = 0.10;
			set(newGrain, 21) = 0.036 ;
            //variant_ids.push_back(2) ;
		}*/
        
        else if(pow((x[0]*MMSP::dx(grid,0) - 0.03*Lx),2) + pow((x[1]*MMSP::dx(grid,1) - 0.8*Ly),2) + pow((x[2]*MMSP::dx(grid,2) - 0.5*Lz),2) <= 8*MMSP::dx(grid,0)*MMSP::dx(grid,0))
		{
            set(newGrain, 1) = 0.0;
            //set(newGrain, 3) = 0.0;
			set(newGrain, 2) = 1.0;
            set(newGrain, 3) = 0.0;
            set(newGrain, 4) = 0.0;
            set(newGrain, 5) = 0.0;
			set(newGrain, 20) = 0.10;
			set(newGrain, 21) = 0.036 ;
            //variant_ids.push_back(2) ;
		}
        
        /*else if(pow((x[0]*MMSP::dx(grid,0) - 0.03*Lx),2) + pow((x[1]*MMSP::dx(grid,1) - 0.7*Ly),2) + pow((x[2]*MMSP::dx(grid,2) - 0.5*Lz),2) <= 8*MMSP::dx(grid,0)*MMSP::dx(grid,0))
		{
            set(newGrain, 1) = 0.0;
            //set(newGrain, 3) = 0.0;
			set(newGrain, 2) = 1.0;
            set(newGrain, 3) = 0.0;
            set(newGrain, 4) = 0.0;
            set(newGrain, 5) = 0.0;
			set(newGrain, 20) = 0.10;
			set(newGrain, 21) = 0.036 ;
            //variant_ids.push_back(2) ;
		}*/
        
        else if(pow((x[0]*MMSP::dx(grid,0) - 0.03*Lx),2) + pow((x[1]*MMSP::dx(grid,1) - 0.6*Ly),2) + pow((x[2]*MMSP::dx(grid,2) - 0.5*Lz),2) <= 8*MMSP::dx(grid,0)*MMSP::dx(grid,0))
		{
            set(newGrain, 1) = 0.0;
            //set(newGrain, 3) = 0.0;
			set(newGrain, 2) = 1.0;
            set(newGrain, 3) = 0.0;
            set(newGrain, 4) = 0.0;
            set(newGrain, 5) = 0.0;
			set(newGrain, 20) = 0.10;
			set(newGrain, 21) = 0.036 ;
            //variant_ids.push_back(2) ;
		}
        
                
        /*else if(pow((x[0]*MMSP::dx(grid,0) - 0.03*Lx),2) + pow((x[1]*MMSP::dx(grid,1) - 0.5*Ly),2) + pow((x[2]*MMSP::dx(grid,2) - 0.5*Lz),2) <= 8*MMSP::dx(grid,0)*MMSP::dx(grid,0))
		{
            set(newGrain, 1) = 0.0;
            //set(newGrain, 3) = 0.0;
			set(newGrain, 2) = 1.0;
            set(newGrain, 3) = 0.0;
            set(newGrain, 4) = 0.0;
            set(newGrain, 5) = 0.0;
			set(newGrain, 20) = 0.10;
			set(newGrain, 21) = 0.036 ;
            //variant_ids.push_back(2) ;
		}*/
        
        else if(pow((x[0]*MMSP::dx(grid,0) - 0.03*Lx),2) + pow((x[1]*MMSP::dx(grid,1) - 0.4*Ly),2) + pow((x[2]*MMSP::dx(grid,2) - 0.5*Lz),2) <= 8*MMSP::dx(grid,0)*MMSP::dx(grid,0))
		{
            set(newGrain, 1) = 0.0;
            //set(newGrain, 3) = 0.0;
			set(newGrain, 2) = 1.0;
            set(newGrain, 3) = 0.0;
            set(newGrain, 4) = 0.0;
            set(newGrain, 5) = 0.0;
			set(newGrain, 20) = 0.10;
			set(newGrain, 21) = 0.036 ;
            //variant_ids.push_back(2) ;
		}
        
                
        /*else if(pow((x[0]*MMSP::dx(grid,0) - 0.03*Lx),2) + pow((x[1]*MMSP::dx(grid,1) - 0.3*Ly),2) + pow((x[2]*MMSP::dx(grid,2) - 0.5*Lz),2) <= 8*MMSP::dx(grid,0)*MMSP::dx(grid,0))
		{
            set(newGrain, 1) = 0.0;
            //set(newGrain, 3) = 0.0;
			set(newGrain, 2) = 1.0;
            set(newGrain, 3) = 0.0;
            set(newGrain, 4) = 0.0;
            set(newGrain, 5) = 0.0;
			set(newGrain, 20) = 0.10;
			set(newGrain, 21) = 0.036 ;
            //variant_ids.push_back(2) ;
		}*/
        
        
        else if(pow((x[0]*MMSP::dx(grid,0) - 0.03*Lx),2) + pow((x[1]*MMSP::dx(grid,1) - 0.2*Ly),2) + pow((x[2]*MMSP::dx(grid,2) - 0.5*Lz),2) <= 8*MMSP::dx(grid,0)*MMSP::dx(grid,0))
		{
            set(newGrain, 1) = 0.0;
            //set(newGrain, 3) = 0.0;
			set(newGrain, 2) = 1.0;
            set(newGrain, 3) = 0.0;
            set(newGrain, 4) = 0.0;
            set(newGrain, 5) = 0.0;
			set(newGrain, 20) = 0.10;
			set(newGrain, 21) = 0.036 ;
            //variant_ids.push_back(2) ;
		}   
        
                
        /*else if(pow((x[0]*MMSP::dx(grid,0) - 0.03*Lx),2) + pow((x[1]*MMSP::dx(grid,1) - 0.1*Ly),2) + pow((x[2]*MMSP::dx(grid,2) - 0.5*Lz),2) <= 8*MMSP::dx(grid,0)*MMSP::dx(grid,0))
		{
            set(newGrain, 1) = 0.0;
            //set(newGrain, 3) = 0.0;
			set(newGrain, 2) = 1.0;
            set(newGrain, 3) = 0.0;
            set(newGrain, 4) = 0.0;
            set(newGrain, 5) = 0.0;
			set(newGrain, 20) = 0.10;
			set(newGrain, 21) = 0.036 ;
            //variant_ids.push_back(2) ;
		}*/
        
        
			
		else
		{
            set(newGrain, 1) = 0.0;
            //set(newGrain, 3) = 0.0;
			set(newGrain, 2) = 0.0;
            set(newGrain, 3) = 0.0;
            set(newGrain, 4) = 0.0;
            set(newGrain, 5) = 0.0;
			set(newGrain, 20) = 0.10 ;
			set(newGrain, 21) = 0.036 ;
            //variant_ids.push_back(2) ;
		}
			
		
		set(gradsqcal_grid(n),1) = 0.0 ;
        set(gradsqcv_grid(n), 1) = 0.0 ;
        set(gradsqcal_grid(n),2) = 0.0 ;
        set(gradsqcv_grid(n), 2) = 0.0 ;
        set(gradsqcal_grid(n),3) = 0.0 ;
        set(gradsqcv_grid(n), 3) = 0.0 ;
		grid(n) = newGrain;
		
	}
    
    for(int nuc_t = 0 ; nuc_t < 500 ; nuc_t++)
    {
        int nuc_index = 2 ;
        for(int n = 0; n < nodes(grid) ; n++)
        {
            MMSP::vector<int> temp_pos = position(grid, n);
            double old_phi = grid(temp_pos)[nuc_index] ;
            double df = gdiff(grid(n)[20], grid(n)[21], T) ;
            double W = W_prefac*((0.10-grid(temp_pos)[20])*20*W_Al + (grid(temp_pos)[21]-0.036)*20*W_V) ;
            double hphiprime = 30*pow(old_phi,2)*pow((old_phi-1),2) ;
            MMSP::vector<store_type> gradient = grad(grid, temp_pos) ;
            MMSP::vector<store_type> gradientsq = gradsq(grid, temp_pos) ;
            double lap_aniso = 0.0  ;
            for(int i = 0 ; i < dim ; i++)
            {
                lap_aniso += epsi_prefac*gradientsq[i][nuc_index]*epsi[nuc_index-1][i][i] ;
            }
            double interaction_energy = 0 ;
            interaction_energy+= W*grid(temp_pos)[nuc_index]*pow((1-grid(temp_pos)[nuc_index]), 2);
            interaction_energy-= W*pow(grid(temp_pos)[nuc_index],2)*((1-grid(temp_pos)[nuc_index]));
            double dfdp = -(df)*hphiprime + interaction_energy - lap_aniso ;    
            double dpdt = -L*10 * (dfdp);
            set(grid(temp_pos), nuc_index) = old_phi + dt * dpdt;	
            
        }
    }

}

void init_fourier()
{
    for(int i = 0 ; i < nx ; i++)
	{
		fk[i] = (2.0*3.14*i)/(double)nx  ;
	}
}


void thermo_auxillary_terms(MMSP::vector<store_type> gradient, MMSP::vector<store_type> gradientsq, double cal, double cv)
{
	double grad_cal = gradient[0][20] + gradient[1][20] + gradient[2][20] ;
	double grad_cv = gradient[0][21]  + gradient[1][21]  + gradient[2][21] ;
	double gradsq_cal = gradientsq[0][20]  + gradientsq[1][20] + gradientsq[2][20] ;
	double gradsq_cv = gradientsq[0][21] + gradientsq[1][21] + gradientsq[2][21] ;                                
                                
    dGAlpha_dAl = (G_Al_alpha/G_normalize - G_Ti_alpha/G_normalize) + log(cal) - log((1-cal-cv)) + cv*L0_HCP_Al_V/G_normalize +  (1-2*cal-cv)*L0_HCP_Al_Ti/G_normalize - cv*L0_HCP_Ti_V/G_normalize ;
    dGBeta_dAl = (G_Al_beta/G_normalize - G_Ti_beta/G_normalize) + log(cal) - log((1-cal-cv)) + cv*L0_BCC_Al_V/G_normalize +
						    (1-2*cal-cv)*L0_BCC_Al_Ti/G_normalize - cv*L0_BCC_Ti_V/G_normalize + (cv*(1-cal-cv) - cal*cv)*L0_BCC_Al_Ti_V/G_normalize ;
    dGAlpha_dV = (G_V_alpha/G_normalize - G_Ti_alpha/G_normalize) + log(cv) - log((1-cal-cv)) + cal*L0_HCP_Al_V/G_normalize +
							(1-cal-2*cv)*L0_HCP_Ti_V/G_normalize - cal*L0_HCP_Al_Ti/G_normalize ;	
    dGBeta_dV = (G_V_beta/G_normalize - G_Ti_beta/G_normalize) + log(cv) - log((1-cal-cv)) + cal*L0_BCC_Al_V/G_normalize +
							(1-cal-2*cv)*L0_BCC_Ti_V/G_normalize - cal*L0_BCC_Al_Ti/G_normalize + (cal*(1-cal-cv) - cal*cv)*L0_BCC_Al_Ti_V/G_normalize ;
    del_dGAlpha_dAl = grad_cal*(1.0/cal + 1/(1-cal-cv) - 2*L0_HCP_Al_Ti/G_normalize) + grad_cv*(1/(1-cal-cv) +
							L0_HCP_Al_V/G_normalize - L0_HCP_Al_Ti/G_normalize -L0_HCP_Ti_V/G_normalize) ;
    del_dGBeta_dAl = grad_cal*(1.0/cal + 1/(1-cal-cv) - 2*L0_BCC_Al_Ti/G_normalize - 2*cv*L0_BCC_Al_Ti_V/G_normalize) + 
								grad_cv*(1/(1-cal-cv) + L0_HCP_Al_V/G_normalize - L0_HCP_Al_Ti/G_normalize -L0_HCP_Ti_V/G_normalize + (1-2*cal-2*cv)*L0_BCC_Al_Ti_V/G_normalize) ;	
    del_dGAlpha_dV = grad_cv*(1.0/cv + 1/(1-cal-cv) - 2*L0_HCP_Ti_V/G_normalize) + grad_cal*(1/(1-cal-cv) + 
							    L0_HCP_Al_V/G_normalize - L0_HCP_Al_Ti/G_normalize - L0_HCP_Ti_V/G_normalize) ;	
    del_dGBeta_dV = 	grad_cv*(1.0/cv + 1/(1-cal-cv) -2*L0_BCC_Ti_V/G_normalize - 2*cal*L0_BCC_Al_Ti_V/G_normalize) + 
							    grad_cal*(1/(1-cal-cv) + L0_BCC_Al_V/G_normalize - L0_BCC_Al_Ti/G_normalize - L0_BCC_Ti_V/G_normalize + (1-2*cal-2*cv)*L0_BCC_Al_Ti_V/G_normalize) ;
							    
    delsq_dGAlpha_dAl = gradsq_cal*(1.0/cal + 1/(1-cal-cv) - 2*L0_HCP_Al_Ti/G_normalize) + gradsq_cv*(1/(1-cal-cv) +
								L0_HCP_Al_V/G_normalize - L0_HCP_Al_Ti/G_normalize -L0_HCP_Ti_V/G_normalize) + grad_cal*(-grad_cal/pow(cal,2) +
								(grad_cal + grad_cv)/pow((1-cal-cv),2)) + grad_cv*((grad_cal + grad_cv)/pow((1-cal-cv),2)) ;
								
    delsq_dGAlpha_dV = gradsq_cv*(1.0/cv + 1/(1-cal-cv) - 2*L0_HCP_Ti_V/G_normalize) + gradsq_cal*(1/(1-cal-cv) +
									L0_HCP_Al_V/G_normalize - L0_HCP_Al_Ti/G_normalize - L0_HCP_Ti_V/G_normalize) + grad_cv*(-grad_cv/pow(cv,2) +
									(grad_cal + grad_cv)/pow((1-cal-cv),2)) + grad_cal*((grad_cal + grad_cv)/pow((1-cal-cv),2)) ;								
								
    delsq_dGBeta_dAl = gradsq_cal*(1.0/cal + 1/(1-cal-cv) - 2*L0_BCC_Al_Ti/G_normalize - 2*cv*L0_BCC_Al_Ti_V/G_normalize) +
								  gradsq_cv*(1/(1-cal-cv) + L0_HCP_Al_V/G_normalize - L0_HCP_Al_Ti/G_normalize -L0_HCP_Ti_V/G_normalize +
								   (1-2*cal-2*cv)*L0_BCC_Al_Ti_V/G_normalize) + grad_cal*(-grad_cal/pow(cal,2) + (grad_cal + grad_cv)/pow((1-cal-cv),2)
								   - 2*grad_cv*L0_BCC_Al_Ti_V/G_normalize) + grad_cv*((grad_cal + grad_cv)/pow((1-cal-cv),2) + (-2*grad_cal*grad_cv)*L0_BCC_Al_Ti_V/G_normalize) ;
									
    delsq_dGBeta_dV = gradsq_cv*(1.0/cv + 1/(1-cal-cv) -2*L0_BCC_Ti_V/G_normalize - 2*cal*L0_BCC_Al_Ti_V/G_normalize) +
								 gradsq_cal*(1/(1-cal-cv) + L0_BCC_Al_V/G_normalize - L0_BCC_Al_Ti/G_normalize - L0_BCC_Ti_V/G_normalize + 
								(1-2*cal-2*cv)*L0_BCC_Al_Ti_V/G_normalize) + grad_cv*(-grad_cv/pow(cv,2) + (grad_cal + grad_cv)/pow((1-cal-cv),2) -
								 2*grad_cal*L0_BCC_Al_Ti_V/G_normalize) + grad_cal*((grad_cal + grad_cv)/pow((1-cal-cv),2) + (-2*grad_cal*grad_cv)*L0_BCC_Al_Ti_V/G_normalize) ;
}

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

void calculate_strains()
{
    for(int n = 0 ; n < nodes(grid) ; n++)
    {
        MMSP::vector<int> x = position(grid, n) ;
        f1(x[0], x[1], x[2])=grid(n)[1] ;
        f2(x[0], x[1], x[2])=grid(n)[2] ; 
        f3(x[0], x[1], x[2])=grid(n)[3] ; 
        f4(x[0], x[1], x[2])=grid(n)[4] ; 
        f5(x[0], x[1], x[2])=grid(n)[5] ; 
    }
        
    Forward1.fft0(f1,F1);  
    Forward2.fft0(f2,F2);
    Forward3.fft0(f3,F3);  
    Forward4.fft0(f4,F4);
    Forward5.fft0(f5,F5);
        
    for(int i = 0 ; i < nx ; i++)
    {
        for(int j = 0 ; j < ny ; j++)
        {
            for(int k = 0 ; k < nzp ; k++)
            {
                double k1 = fk[i] ;
                double k2 = fk[j] ;
                double k3 = fk[k] ;
                double modk = (sqrt(k1*k1 + k2*k2 + k3*k3)) ;
                if(modk==0.0) 
                {
                    dfdstr1(i,j,k) = 0.0  ;
                    dfdstr2(i,j,k) = 0.0  ;
                    dfdstr3(i,j,k) = 0.0  ;
                    dfdstr4(i,j,k) = 0.0  ;
                    dfdstr5(i,j,k) = 0.0  ;
                }
                else 
                {
                    dfdstr1(i,j,k) = Bpq(k1, k2, k3, modk, 1, 1)*F1(i,j,k) + Bpq(k1, k2, k3, modk, 1, 2)*F2(i,j,k) + Bpq(k1, k2, k3, modk, 1, 3)*F3(i,j,k) + Bpq(k1, k2, k3, modk, 1, 4)*F4(i,j,k) + Bpq(k1, k2, k3, modk, 1, 5)*F5(i,j,k) ;
                    dfdstr2(i,j,k) = Bpq(k1, k2, k3, modk, 2, 1)*F1(i,j,k) + Bpq(k1, k2, k3, modk, 2, 2)*F2(i,j,k) + Bpq(k1, k2, k3, modk, 2, 3)*F3(i,j,k) + Bpq(k1, k2, k3, modk, 2, 4)*F4(i,j,k) + Bpq(k1, k2, k3, modk, 2, 5)*F5(i,j,k) ;
                    dfdstr3(i,j,k) = Bpq(k1, k2, k3, modk, 3, 1)*F1(i,j,k) + Bpq(k1, k2, k3, modk, 3, 2)*F2(i,j,k) + Bpq(k1, k2, k3, modk, 3, 3)*F3(i,j,k) + Bpq(k1, k2, k3, modk, 3, 4)*F4(i,j,k) + Bpq(k1, k2, k3, modk, 3, 5)*F5(i,j,k) ;
                    dfdstr4(i,j,k) = Bpq(k1, k2, k3, modk, 4, 1)*F1(i,j,k) + Bpq(k1, k2, k3, modk, 4, 2)*F2(i,j,k) + Bpq(k1, k2, k3, modk, 4, 3)*F3(i,j,k) + Bpq(k1, k2, k3, modk, 4, 4)*F4(i,j,k) + Bpq(k1, k2, k3, modk, 4, 5)*F5(i,j,k) ;
                    dfdstr5(i,j,k) = Bpq(k1, k2, k3, modk, 5, 1)*F1(i,j,k) + Bpq(k1, k2, k3, modk, 5, 2)*F2(i,j,k) + Bpq(k1, k2, k3, modk, 5, 3)*F3(i,j,k) + Bpq(k1, k2, k3, modk, 5, 4)*F4(i,j,k) + Bpq(k1, k2, k3, modk, 5, 5)*F5(i,j,k) ;
                }
             
            }
        }
    }
        
    Backward1.fft0Normalized(dfdstr1, dfdstr_real1); 
    Backward2.fft0Normalized(dfdstr2, dfdstr_real2); 
    Backward3.fft0Normalized(dfdstr3, dfdstr_real3); 
    Backward4.fft0Normalized(dfdstr4, dfdstr_real4); 
    Backward5.fft0Normalized(dfdstr5, dfdstr_real5); 
}

void calc_gradsq_c()
{
    for(int n = 0 ; n<nodes(grid) ; n++)
    {
        MMSP::vector<int> x = position(grid, n); 
        MMSP::vector<store_type> gradientsqtemp = gradsq(grid, x) ;
        double temp1 = 0.0 ;
        double temp2 = 0.0 ;
        set(gradsqcal_grid(n),1) = gradientsqtemp[0][20] ;
        set(gradsqcal_grid(n),2) = gradientsqtemp[1][20] ;
        set(gradsqcal_grid(n),3) = gradientsqtemp[2][20] ;
            
        set(gradsqcv_grid(n),1) = gradientsqtemp[0][21] ;
        set(gradsqcv_grid(n),2) = gradientsqtemp[1][21] ;		
        set(gradsqcv_grid(n),3) = gradientsqtemp[2][20] ;
    }
}

void write_output(int t)
{
    std::string file_name = "output_" + to_string(t) + ".txt" ;
            ofstream myfile;
            myfile.open(file_name.c_str());
            for(int n = 0 ; n < nodes(grid) ; n++)    //for n < nodes(grid) ; 
            {
                MMSP::vector<int> s = position(grid,n) ;
                MMSP::vector<store_type> gradient = grad(grid, s) ;
                MMSP::vector<store_type> gradientsq = gradsq(grid, s) ;
                double lap_aniso = 1.5*(gradientsq[0][2]*epsi[1][0][0] + gradientsq[1][2]*epsi[1][1][1] + gradientsq[2][2]*epsi[1][2][2]) ; 
                W[n] = W_prefac*((0.10-grid(n)[20])*40*W_Al + (grid(n)[21]-0.036)*40*W_V) ;
                MMSP::vector<store_type> gradcalp4 = gradsq(gradsqcal_grid, s) ;	
                MMSP::vector<store_type> gradcvp4 = gradsq(gradsqcv_grid, s) ;
                double ss = 2*grid(n)[20]/0.10 + 2*0.036/grid(n)[21] ;
                L = (ss*L_orig)/L_norm;
                double gradp4_cal = gradcalp4[0][1] + gradcalp4[1][2] + gradcalp4[2][3] ; ;
                double gradp4_cv  = gradcvp4[0][1] + gradcvp4[1][2] + gradcvp4[2][3] ; ;
                double sum = 0 ;
                for(int h = 0 ; h < length(grid(s)) ; h++)
                {
                    int hindex = MMSP::index(grid(s),h) ;
                    if(hindex <= 12)
                    {
                        sum+=grid(s)[hindex]*grid(s)[hindex] ;
                    }
                }
                myfile<<s[0]<<","<<s[1]<<","<<s[2]<<","<<sum<<","<<grid(s)[20]<<","<<grid(s)[21]<<","<<dfdstr_real1[s[0]][s[1]][s[2]]<<","<<dfdstr_real2[s[0]][s[1]][s[2]]<<","<<dfdstr_real3[s[0]][s[1]][s[2]]<<","<<dfdstr_real4[s[0]][s[1]][s[2]]<<","<<dfdstr_real5[s[0]][s[1]][s[2]]<<"\n" ; 
            }
            myfile.close();
}

void nucleation()
{
            double deltag = 0.1 ;
            double kappa1 = 0.05 ; 
            double kappa2 = 0.01 ; 
            const int cell_length = 4 ;
            int no_of_cells = (nx*ny)/(cell_length*cell_length*cell_length) ; 
            int cells_in_one_dir = nx/cell_length ; 
            MMSP::vector<int> x(cell_length) ; 
            MMSP::vector<int> y(cell_length) ; 
            MMSP::vector<int> z(cell_length) ; 
        
            for(int i = 0 ; i < cells_in_one_dir ; i++)
            {
                for(int j = 0 ; j < cells_in_one_dir ; j++)
                {
                    for(int k = 0 ; k < cells_in_one_dir ; k++)
                    {
                        double flag = 0.0 ; double mid_x, mid_y, mid_z ;
                    
                        for(int l = 0 ; l < cell_length ; l++)
                        {
                            x[l] = i*cell_length + l ;
                            y[l] = j*cell_length + l ;
                            z[l] = k*cell_length + l ;
                        }
                        mid_x = x[cell_length/2] ;
                        mid_y = y[cell_length/2] ;
                        mid_z = z[cell_length/2] ;
                
                        MMSP::vector<int> temp_pos(3); 
                        double cal_avg = 0.0 ;
                        double cv_avg = 0.0 ;
                        for(int p = 0; p < cell_length; p++)
                        {
                            for(int q = 0 ; q < cell_length ; q++)
                            {
                                for(int r = 0 ; r < cell_length ; r++)
                                {
                                    temp_pos[0] = x[p] ; temp_pos[1] = y[q] ; temp_pos[2] = z[r] ;
                                                         
                                    double sum = 0.0 ;
                                    for(int h = 0 ; h < length(grid(temp_pos)) ; h++)
                                    {
                                        int hindex = MMSP::index(grid(temp_pos),h) ;
                                        if(hindex <= 12)
                                        {
                                            sum+=grid(temp_pos)[hindex] ;
                                        }
                                    }
                                    flag+= (1.0 - sum) ;                           
                        
                                    cal_avg += grid(temp_pos)[20]/1000.0 ; 
                                    cv_avg += grid(temp_pos)[21]/1000.0 ; 
                                }
                            }
                        }
                     //cout<<flag<<" "<<endl;
                        if(flag >= pow(cell_length,3)-20)
                        {    				
                    
                            double prob1 ; 
                            double df = gdiff(cal_avg, cv_avg, T) ;
                            double jstar = kappa1*exp(-kappa2*deltag) ;
                            double prob = 1 - exp(-jstar) ;
                            double random_no = (double)(rand() % 1000) ;
                            double decision_var = random_no/1000.0 ;
                            if(prob > decision_var)
                            {
                                store_type newGrain;
                                int nuc_index = (int)(rand()%2) + 1 ;
                                cout<<"Nucleating index: "<<nuc_index<<endl;
                                variant_ids.push_back(nuc_index);
                                for(int p = 0; p < cell_length; p++)
                                {
                                    for(int q = 0 ; q < cell_length ; q++)
                                    {
                                        for(int r = 0 ; r < cell_length ; r++)
                                        {
                                            temp_pos[0] = x[p] ; temp_pos[1] = y[q] ; temp_pos[2] = z[r] ;
                                            temp_pos[0] = x[p] ; temp_pos[1] = y[q] ;
                                            if(pow((x[p]*MMSP::dx(grid,0) - mid_x*MMSP::dx(grid,0)),2) + pow((y[q]*MMSP::dx(grid,1) - mid_y*MMSP::dx(grid,1)),2) + pow((z[r]*MMSP::dx(grid,2) - mid_z*MMSP::dx(grid,2)),2) <= 8*MMSP::dx(grid,0)*MMSP::dx(grid,0))
                                            {
                                                set(grid(temp_pos), nuc_index) = 1.0; 
                                            }
                                            else
                                            {
                                                set(grid(temp_pos), nuc_index) = 0.0;
                                            }
                                        }
                                    }
                                }
                                for(int nuc_t = 0 ; nuc_t < 500 ; nuc_t++)
                                {
                                    for(int p = 0; p < cell_length; p++)
                                    {
                                        for(int q = 0 ; q < cell_length ; q++)
                                        {
                                            for(int r = 0 ; r < cell_length ; r++)
                                            {
                                                temp_pos[0] = x[p] ; temp_pos[1] = y[q] ; temp_pos[2] = z[r] ;
                                                double old_phi = grid(temp_pos)[nuc_index] ;
                                                double W = W_prefac*((0.10-grid(temp_pos)[20])*20*W_Al + (grid(temp_pos)[21]-0.036)*20*W_V) ;
                                                double hphiprime = 30*pow(old_phi,2)*pow((old_phi-1),2) ;
                                                MMSP::vector<store_type> gradient = grad(grid, temp_pos) ;
                                                MMSP::vector<store_type> gradientsq = gradsq(grid, temp_pos) ;
                                                double lap_aniso = epsi_prefac*(gradientsq[0][nuc_index]*epsi[nuc_index-1][0][0] + gradientsq[1][nuc_index]*epsi[nuc_index-1][1][1] + gradientsq[2][nuc_index]*epsi[nuc_index-1][2][2])  ;
                                                double interaction_energy = 0 ;
                                                interaction_energy+= W*grid(temp_pos)[nuc_index]*pow((1-grid(temp_pos)[nuc_index]), 2);
                                                interaction_energy-= W*pow(grid(temp_pos)[nuc_index],2)*((1-grid(temp_pos)[nuc_index]));
                                                double dfdp = -(df)*hphiprime + interaction_energy - lap_aniso ;    
                                                double dpdt = -L*10 * (dfdp);
                                                set(grid(temp_pos), nuc_index) = old_phi + dt * dpdt;	
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

void update_node_n()
{
        for(int n = 0 ; n < nodes(grid) ; n++)
		{
			MMSP::vector<int> x = position(grid, n); 
            
            double strain_energy[13];
			strain_energy[0] = 0.0 ;
			strain_energy[1] = (double)(1.0/(nx*ny*nz))*(dfdstr_real1[x[0]][x[1]][x[2]]) ; 
            strain_energy[2] = (double)(1.0/(nx*ny*nz))*(dfdstr_real2[x[0]][x[1]][x[2]]) ; 
            strain_energy[3] = (double)(1.0/(nx*ny*nz))*(dfdstr_real3[x[0]][x[1]][x[2]]) ; 
            strain_energy[4] = (double)(1.0/(nx*ny*nz))*(dfdstr_real4[x[0]][x[1]][x[2]]) ; 
            strain_energy[5] = (double)(1.0/(nx*ny*nz))*(dfdstr_real5[x[0]][x[1]][x[2]]) ; 
			
			G_Alpha[n] = ((grid(n)[20]*G_Al_alpha + grid(n)[21]*G_V_alpha + (1-grid(n)[20]-grid(n)[21])*G_Ti_alpha +
					 R*T*(grid(n)[20]*log(grid(n)[20]) + grid(n)[21]*log(grid(n)[21]) + (1-grid(n)[20]-grid(n)[21])*log((1-grid(n)[20]-grid(n)[21]))) +
					 grid(n)[20]*grid(n)[21]*L0_HCP_Al_V + grid(n)[20]*(1-grid(n)[20]-grid(n)[21])*(L0_HCP_Al_Ti + L1_HCP_Al_Ti*(2*grid(n)[20] + grid(n)[21] - 1)) + grid(n)[21]*(1-grid(n)[20]-grid(n)[21])*L0_HCP_Ti_V))/G_normalize ;
					 
			G_Beta[n] = (grid(n)[20]*G_Al_beta + grid(n)[21]*G_V_beta + (1-grid(n)[20]-grid(n)[21])*G_Ti_beta +
					 R*T*(grid(n)[20]*log(grid(n)[20]) + grid(n)[21]*log(grid(n)[21]) + (1-grid(n)[20]-grid(n)[21])*log((1-grid(n)[20]-grid(n)[21]))) +
					 grid(n)[20]*grid(n)[21]*L0_BCC_Al_V + grid(n)[20]*(1-grid(n)[20]-grid(n)[21])*L0_BCC_Al_Ti + grid(n)[21]*(1-grid(n)[20]-grid(n)[21])*L0_BCC_Ti_V +
					 grid(n)[20]*grid(n)[21]*(1-grid(n)[20]-grid(n)[21])*L0_BCC_Al_Ti_V)/G_normalize ;
					 
			//W[n] = W_prefac*(grid(n)[20]*W_Al + grid(n)[21]*W_V + (1-grid(n)[20]-grid(n)[21])*W_Ti) ;
            W[n] = W_prefac*((0.10-grid(n)[20])*50*W_Al + (grid(n)[21]-0.036)*50*W_V) ;
			MMSP::vector<store_type> gradient = grad(grid, x) ;
			MMSP::vector<store_type> gradientsq = gradsq(grid, x) ;
            		
			thermo_auxillary_terms(gradient, gradientsq, grid(n)[20], grid(n)[21]) ;

			MMSP::vector<store_type> gradcalp4 = gradsq(gradsqcal_grid, x) ;	
			MMSP::vector<store_type> gradcvp4 = gradsq(gradsqcv_grid, x) ;
		
            double gradp4_cal = gradcalp4[0][1] + gradcalp4[1][2] + gradcalp4[2][3]  ;
            double gradp4_cv  = gradcvp4[0][1] + gradcvp4[1][2] + gradcvp4[2][3] ;
		
		
			for(int h = 0 ; h < length(grid(x)) ; h++)
			{
				int hindex = MMSP::index(grid(x),h) ;
				if(hindex <= 12)
				{
					hphidoubleprime[n][hindex-1] = 60*grid(n)[hindex]*(2*grid(n)[hindex]-1)*(grid(n)[hindex]-1) ; 
					hphi[n][hindex-1] = pow(grid(n)[hindex],3)*(6*pow(grid(n)[hindex],2) - 15*grid(n)[hindex] + 10) ;
					hphiprime[n][hindex-1] = 30*pow(grid(n)[hindex],2)*pow((grid(n)[hindex]-1),2) ;
		
					gphidoubleprime[n][hindex-1] = 2*(6*pow(grid(n)[hindex],2) - 6*grid(n)[hindex] + 1) ; 
					gphiprime[n][hindex-1] = 2*grid(n)[hindex]*(1-grid(n)[hindex])*(1-2*grid(n)[hindex]) ;
					gphi[n][hindex-1] = pow(grid(n)[hindex],2)*pow((1-grid(n)[hindex]),2) ;
				
				
				}
			}
            
            double cal = grid(n)[20] ;
            double cv = grid(n)[21] ;
            double grad_cal = gradient[0][20] + gradient[1][20] + gradient[2][20] ;
            double grad_cv = gradient[0][21] + gradient[1][21] + gradient[2][21] ;
            
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
                    double delphi = gradient[0][hindex] + gradient[1][hindex] + gradient[2][hindex]  ;
                    double delphisq = gradientsq[0][hindex] + gradientsq[1][hindex] + gradientsq[2][hindex]  ;
                    hphidoubleprimesum += hphidoubleprime[n][hindex-1]*pow(delphi,2) ;
                    hphiprimesum1 += hphiprime[n][hindex-1]*delphisq;
                    hphiprimesum2 += hphiprime[n][hindex-1]*delphi ;
                    hphisum += hphi[n][hindex-1] ;
				
                    gphiprimesum += gphiprime[n][hindex-1]*delphisq ;
                    gphidoubleprimesum += gphidoubleprime[n][hindex-1]*delphi ;
                }
            }
               
            
            double c_al_rhs = 2*(del_dGAlpha_dAl - del_dGBeta_dAl)*hphiprimesum2 + delsq_dGAlpha_dAl*hphisum + delsq_dGBeta_dAl*(1-hphisum) + 
					  (dGAlpha_dAl - dGBeta_dAl)*(hphidoubleprimesum + hphiprimesum1) + (W_Al - W_Ti)*(gphiprimesum + gphidoubleprimesum) ;
            double c_v_rhs = 2*(del_dGAlpha_dV - del_dGBeta_dV)*hphiprimesum2 + delsq_dGAlpha_dV*hphisum + delsq_dGBeta_dV*(1-hphisum) + 
					  (dGAlpha_dV - dGBeta_dV)*(hphidoubleprimesum + hphiprimesum1) + (W_V - W_Ti)*(gphiprimesum + gphidoubleprimesum)  ;


            if((grid(n)[20] + dt*(Dalal*(c_al_rhs)- kappa_c*gradp4_cal)) < 0.005)
            {
                set(update(n), 20) = 0.005 ;
            }
            else
            {
                set(update(n), 20) = grid(n)[20] + dt*(Dalal*(c_al_rhs)- kappa_c*gradp4_cal) ; 
            }
            if((grid(n)[21] + dt*(Dvv*(c_v_rhs) - kappa_c*gradp4_cv)) < 0.005)
            {
                set(update(n), 21) = 0.005 ;
            }
            else
            {
                set(update(n), 21) = grid(n)[21] + dt*(Dvv*(c_v_rhs) - kappa_c*gradp4_cv) ; 
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
				W[n] = W_prefac*((0.10-grid(n)[20])*40*W_Al + (grid(n)[21]-0.036)*40*W_V) ;
                //cout<<W[n]<<endl; ; 
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
					set(dFdp, jindex) = -(df)*hphiprime[n][jindex-1] + interaction_energy + strain_energy[jindex]  - lap_aniso ;    //
					dFall += dFdp[jindex];	
				}
			}
            double ss = 2*grid(n)[20]/0.10 + 2*0.036/grid(n)[21] ;
            L = (ss*L_orig)/L_norm;
			for (int h = 0; h < length(grid(n)); h++) 
			{	
					
				int hindex = MMSP::index(grid(n), h);   
				if(hindex<=12)
				{
					store_type dpdt;
					set(dpdt, hindex) = -L * (dFdp[hindex]);
					set(update(n), hindex) = grid(n)[hindex] + dt * dpdt[hindex];	
				}
			}			
		}
}
