
void customoutput(MMSP::grid<dim, store_type> grid, int t)
{
	
	if(output_pfsq == 1)
	{
		string file_name = "pf_" + to_string(t) + ".txt" ;
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
			
			myfile<<s[0]<<","<<s[1]<<","<<s[2]<<","<<sum<<"\n" ;    
			
		}
	
		myfile.close();
	}
	
	if(output_cal == 1)
	{
		string file_name1 = "c_Al_" + to_string(t) + ".txt" ;
		ofstream myfile1;
		myfile1.open(file_name1.c_str());
		for(int n = 0 ; n < nodes(grid) ; n++)    //for n < nodes(grid) ; 
		{
			MMSP::vector<int> s = position(grid,n) ;
			myfile1<<s[0]<<","<<s[1]<<","<<s[2]<<","<<grid(s)[20]<<"\n" ;
		}
		myfile1.close();
	}
	
	if(output_cv==1)
	{	
		string file_name2 = "c_V_" + to_string(t) + ".txt" ;
		ofstream myfile2;
		myfile2.open(file_name2.c_str());
		for(int n = 0 ; n < nodes(grid) ; n++)    //for n < nodes(grid) ; 
		{
			MMSP::vector<int> s = position(grid,n) ;
			myfile2<<s[0]<<","<<s[1]<<","<<s[2]<<","<<grid(s)[21]<<"\n" ;
		}
		myfile2.close();
	}	
	
	if(output_chemnuc == 1)
	{
		string file_name3 = "deltachem_" + to_string(t) + ".txt" ;
		ofstream myfile3;
		myfile3.open(file_name3.c_str());
	
		
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
					 
			myfile3<<s[0]<<","<<s[1]<<","<<s[2]<<","<<(g_beta-g_alpha)<<"\n" ;
		}
		myfile3.close();
	}
	
	if(output_strainintnuc == 1)
	{
		string file_name4 = "strainint_" + to_string(t) + ".txt" ;
		ofstream myfile4;
		myfile4.open(file_name4.c_str());
		
		for(int n = 0 ; n < nodes(grid) ; n++)
		{
			MMSP::vector<int> s = position(grid, n) ;
			for(int h = 0 ; h < length(grid(n)) ; h++)
			{
				int hindex = MMSP::index(grid(n), h) ;
				f1nuc(s[0], s[1], s[2])=grid(n)[1] ;   
				f2nuc(s[0], s[1], s[2])=grid(n)[2] ; 
				f3nuc(s[0], s[1], s[2])=grid(n)[3] ; 
				f4nuc(s[0], s[1], s[2])=grid(n)[4] ; 
				f5nuc(s[0], s[1], s[2])=grid(n)[5] ; 
			}	
		
		}
	

  
  
		Forward1n.fft0(f1nuc,F1nuc);    
		Forward2n.fft0(f2nuc,F2nuc); 
		Forward3n.fft0(f3nuc,F3nuc); 
		Forward4n.fft0(f4nuc,F4nuc); 
		Forward5n.fft0(f5nuc,F5nuc);
		
		
		
		for(int n = 0 ; n < nodes(grid) ; n++)
		{
			MMSP::vector<int> x = position(grid, n) ;
			double k1 = fk[x[0]]*0.1 ;
			double k2 = fk[x[1]]*0.1 ;
			double k3 = fk[x[2]]*0.1 ;
			double modk = (sqrt(k1*k1 + k2*k2 + k3*k3)) ;
			if(modk==0.0) modk=0.1   ;
			elint1(x[0],x[1],x[2]) = scaling*(Bpq(k1, k2, k3, modk, 1, 1)*F1nuc(x[0],x[1],x[2])*F1nuc(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 1, 2)*F1nuc(x[0],x[1],x[2])*F2nuc(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 1, 3)*F1nuc(x[0],x[1],x[2])*F3nuc(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 1, 4)*F1nuc(x[0],x[1],x[2])*F4nuc(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 1, 5)*F1nuc(x[0],x[1],x[2])*F5nuc(x[0],x[1],x[2]))  ;
			elint2(x[0],x[1],x[2]) = scaling*(Bpq(k1, k2, k3, modk, 2, 1)*F2nuc(x[0],x[1],x[2])*F1nuc(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 2, 2)*F2nuc(x[0],x[1],x[2])*F2nuc(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 2, 3)*F2nuc(x[0],x[1],x[2])*F3nuc(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 2, 4)*F2nuc(x[0],x[1],x[2])*F4nuc(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 2, 5)*F2nuc(x[0],x[1],x[2])*F5nuc(x[0],x[1],x[2])) ;
			elint3(x[0],x[1],x[2]) = scaling*(Bpq(k1, k2, k3, modk, 3, 1)*F3nuc(x[0],x[1],x[2])*F1nuc(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 3, 2)*F3nuc(x[0],x[1],x[2])*F2nuc(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 3, 3)*F3nuc(x[0],x[1],x[2])*F3nuc(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 3, 4)*F3nuc(x[0],x[1],x[2])*F4nuc(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 3, 5)*F3nuc(x[0],x[1],x[2])*F5nuc(x[0],x[1],x[2]))  ;
			elint4(x[0],x[1],x[2]) = scaling*(Bpq(k1, k2, k3, modk, 4, 1)*F4nuc(x[0],x[1],x[2])*F1nuc(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 4, 2)*F4nuc(x[0],x[1],x[2])*F2nuc(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 4, 3)*F4nuc(x[0],x[1],x[2])*F3nuc(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 4, 4)*F4nuc(x[0],x[1],x[2])*F4nuc(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 4, 5)*F4nuc(x[0],x[1],x[2])*F5nuc(x[0],x[1],x[2]))  ;
			elint5(x[0],x[1],x[2]) = scaling*(Bpq(k1, k2, k3, modk, 5, 1)*F5nuc(x[0],x[1],x[2])*F1nuc(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 5, 2)*F5nuc(x[0],x[1],x[2])*F2nuc(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 5, 3)*F5nuc(x[0],x[1],x[2])*F3nuc(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 5, 4)*F5nuc(x[0],x[1],x[2])*F4nuc(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 5, 5)*F5nuc(x[0],x[1],x[2])*F5nuc(x[0],x[1],x[2]))  ;
		}
		
		Backward1n.fft0Normalized(elint1, elint_real1);  
		Backward2n.fft0Normalized(elint2, elint_real2); 
		Backward3n.fft0Normalized(elint3, elint_real3); 
		Backward4n.fft0Normalized(elint4, elint_real4); 
		Backward5n.fft0Normalized(elint5, elint_real5); 
		
		for(int n = 0 ; n < nodes(grid) ; n++)
		{
			MMSP::vector<int> s = position(grid, n) ; 
			
			myfile4<<s[0]<<","<<s[1]<<","<<s[2]<<","<<elint_real1(s[0], s[1], s[2])<<','<<elint_real2(s[0], s[1], s[2])<<','<<elint_real3(s[0], s[1], s[2])<<','<<elint_real4(s[0], s[1], s[2])<<','<<elint_real5(s[0], s[1], s[2])<<"\n" ;
		
		}
		myfile4.close();
	}
}
	
	
		
			
			
			

