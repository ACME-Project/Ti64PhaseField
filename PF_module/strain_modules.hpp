using namespace std ;


/*double * c_subset_for_omega_inv(int i1, int i2, int i3, int j1, int j2, int j3)
{
	double c_temp[9] = {c[i1][j1], c[i1][j2], c[i1][j3], c[i2][j1], c[i2][j2], c[i2][j3], c[i3][j1], c[i3][j2], c[i3][j3]} ;
	//cout<<c[i1][j1]<<" "<<c[i1][j2]<<" "<<c[i2][j2]<<endl;
	return c_temp ;
}*/

double c[6][6], sigma00[12][6] ;
double c_norm = 1;

double omega_inv_ij(double* cnew, double* k ,  double modk)
{
	double sum = 0 ;
	for(int i = 0 ; i < 9 ; i++)
	{
		sum += cnew[i]*k[i] ;
		//cout<<cnew[i]<<" "<<k[i]<<endl;
	}
	
	//cout<<sum<<endl ;
	
	return sum ; 
}

void define_c_sigma()
{

	for(int i = 0 ; i < 6 ; i++)
	{
		for(int j = 0 ; j < 6 ; j++)
		{
			c[i][j] = 1*(0.1 + (0.1*i + 0.1*j)/2) ;   //Edit later
		}
	}
	
	
/*c[0][0] = 119.0e-01 /c_norm ;
c[1][1] = 64.0e-01/c_norm ;
c[2][2] = 36.0e-01/c_norm ;
c[0][1] = c[1][0] = 0.0e-03/c_norm ;
c[1][2] = c[2][1] = 134.0e-01/c_norm ;
c[0][2] = c[2][0] = 49.0e-01/c_norm ;*/


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


double B(double k1, double k2, double k3, double modk, int variant) {
	
	//cout<<"In B variant"<<endl;
	
	double c11[9] = {c[0][0], c[0][5], c[0][4], c[5][0], c[5][5], c[5][4], c[4][0], c[4][5], c[4][4]} ; 
	//cout<<c11[0]<<" "<<c11[1]<<" "<<c11[2]<<endl;
	double c12[9] = {c[0][5], c[0][1], c[0][3], c[5][5], c[5][1], c[5][3], c[4][5], c[4][1], c[4][3]} ; 
	double c13[9] = {c[0][4], c[0][3], c[0][2], c[5][4], c[5][3], c[5][2], c[4][4], c[4][3], c[4][2]} ; 
	double c21[9] = {c[5][0], c[5][5], c[5][4], c[1][0], c[1][5], c[1][4], c[3][0], c[3][5], c[3][4]} ; 
	double c22[9] = {c[5][5], c[5][1], c[5][3], c[1][5], c[1][1], c[1][3], c[3][5], c[3][1], c[3][3]} ; 
	double c23[9] = {c[5][4], c[5][3], c[5][2], c[1][4], c[1][3], c[1][2], c[3][4], c[3][3], c[3][2]} ; 
	double c31[9] = {c[4][0], c[4][5], c[4][4], c[3][0], c[3][5], c[3][4], c[2][0], c[2][5], c[2][4]} ; 
	double c32[9] = {c[4][5], c[4][1], c[4][3], c[3][5], c[3][1], c[3][3], c[2][5], c[2][1], c[2][3]} ; 
	double c33[9] = {c[4][4], c[4][3], c[4][2], c[3][4], c[3][3], c[3][2], c[2][4], c[2][3], c[2][2]} ; 
	double k[9] = {k1*k1, k1*k2, k1*k3, k2*k1, k2*k2, k2*k3, k3*k1, k3*k2, k3*k3} ;
	
	//cout<<k[0]<<" "<<k[7]<<endl;

	double omega_inv[3][3] = {{omega_inv_ij(c11, k, modk), omega_inv_ij(c12, k, modk), omega_inv_ij(c13, k, modk)}, {omega_inv_ij(c21, k, modk), omega_inv_ij(c22, k, modk), omega_inv_ij(c23, k, modk)},
						{omega_inv_ij(c31, k, modk), omega_inv_ij(c32, k, modk), omega_inv_ij(c33, k, modk)}} ;

	/*for(int i = 0 ; i < 3 ; i++) 
	{
		for(int j = 0 ; j < 3 ; j++)
		{
			if(omega_inv[i][j]==0) cout<<"Found zero omega_inv"<<endl ;
		}
	}*/
    //cout<<omega_inv[0][0]<<" "<<omega_inv[0][1]<<" "<<omega_inv[0][2]<<" "<<omega_inv[1][0]<<" "<<omega_inv[1][1]<<" "<<omega_inv[1][2]<<" "<<omega_inv[2][0]<<" "<<omega_inv[2][1]<<" "<<omega_inv[2][2]<<" "<<endl ;
	double det_omega_inv = omega_inv[0][0]*(omega_inv[1][1]*omega_inv[2][2] - omega_inv[1][2]*omega_inv[2][1]) - omega_inv[0][1]*(omega_inv[1][0]*omega_inv[2][2] - omega_inv[1][2]*omega_inv[2][0]) +
					   omega_inv[0][2]*(omega_inv[1][0]*omega_inv[2][1] - omega_inv[1][1]*omega_inv[2][0]) + 1e-05  ;
	if(det_omega_inv==0) cout<<"Found zero det"<<endl ;
	
					   
	double omega[3][3] ; 

	for(int i = 0 ; i < 2 ; i++)
	{
		for(int j = 0 ; j < 2 ; j++)
		{
			omega[i][j] =  pow(-1, (i+j))*(1/det_omega_inv)*(omega_inv[(j+1)%3][(i+1)%3]*omega_inv[(j+2)%3][(i+2)%3] - omega_inv[(j+1)%3][(i+2)%3]*omega_inv[(j+2)%3][(i+1)%3]) ;  //Inverse of omega_inv
			//cout<<omega_inv[0][0]<<" "<<det_omega_inv<<" "<<omega[i][j]<<endl;
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
    
    return (sfts_sum - (i1jkl + i2jkl + i3jkl)) ; }
    
double Bpq(double k1, double k2, double k3, double modk, int p, int q) {
	
	double c11[9] = {c[0][0], c[0][5], c[0][4], c[5][0], c[5][5], c[5][4], c[4][0], c[4][5], c[4][4]} ; 
	//cout<<c11[0]<<" "<<c11[1]<<" "<<c11[2]<<endl;
	double c12[9] = {c[0][5], c[0][1], c[0][3], c[5][5], c[5][1], c[5][3], c[4][5], c[4][1], c[4][3]} ; 
	double c13[9] = {c[0][4], c[0][3], c[0][2], c[5][4], c[5][3], c[5][2], c[4][4], c[4][3], c[4][2]} ; 
	double c21[9] = {c[5][0], c[5][5], c[5][4], c[1][0], c[1][5], c[1][4], c[3][0], c[3][5], c[3][4]} ; 
	double c22[9] = {c[5][5], c[5][1], c[5][3], c[1][5], c[1][1], c[1][3], c[3][5], c[3][1], c[3][3]} ; 
	double c23[9] = {c[5][4], c[5][3], c[5][2], c[1][4], c[1][3], c[1][2], c[3][4], c[3][3], c[3][2]} ; 
	double c31[9] = {c[4][0], c[4][5], c[4][4], c[3][0], c[3][5], c[3][4], c[2][0], c[2][5], c[2][4]} ; 
	double c32[9] = {c[4][5], c[4][1], c[4][3], c[3][5], c[3][1], c[3][3], c[2][5], c[2][1], c[2][3]} ; 
	double c33[9] = {c[4][4], c[4][3], c[4][2], c[3][4], c[3][3], c[3][2], c[2][4], c[2][3], c[2][2]} ; 
	//double k[9] = {k1*k1, k1*k2, k1*k3, k2*k1, k2*k2, k2*k3, k3*k1, k3*k2, k3*k3} ;
	double k[9] = {k1*k1, k1*k2, k1*k3, k2*k1, k2*k2, k2*k3, k3*k1, k3*k2, k3*k3} ;

	double omega_inv[3][3] = {{omega_inv_ij(c11, k, modk), omega_inv_ij(c12, k, modk), omega_inv_ij(c13, k, modk)}, {omega_inv_ij(c21, k, modk), omega_inv_ij(c22, k, modk), omega_inv_ij(c23, k, modk)},
						{omega_inv_ij(c31, k, modk), omega_inv_ij(c32, k, modk), omega_inv_ij(c33, k, modk)}} ;

	/*for(int i = 0 ; i < 3 ; i++) 
	{
		for(int j = 0 ; j < 3 ; j++)
		{
			if(omega_inv[i][j]==0) cout<<"Found zero omega_inv"<<endl ;
			//cout<<omega_inv[i][j]<<endl;
		}
	}*/

	double det_omega_inv = omega_inv[0][0]*(omega_inv[1][1]*omega_inv[2][2] - omega_inv[1][2]*omega_inv[2][1]) - omega_inv[0][1]*(omega_inv[1][0]*omega_inv[2][2] - omega_inv[1][2]*omega_inv[2][0]) +
					   omega_inv[0][2]*(omega_inv[1][0]*omega_inv[2][1] - omega_inv[1][1]*omega_inv[2][0]) + 1e-05  ;
	if(det_omega_inv==0) cout<<"Found zero det"<<endl ;
	//cout<<det_omega_inv<<endl;
					   
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
			
    
    return (sfts_sum - (i1jkl + i2jkl + i3jkl)) ; } 
    
void calculate_strains(MMSP::grid<dim, store_type> grid) 
{
	for(int n = 0 ; n < nodes(grid) ; n++)
		{
			MMSP::vector<int> s = position(grid, n) ;
			for(int h = 0 ; h < length(grid(n)) ; h++)
			{
				int hindex = MMSP::index(grid(n), h) ;
				f1(s[0], s[1], s[2])=grid(n)[1] ;   
				f2(s[0], s[1], s[2])=grid(n)[2] ; 
				f3(s[0], s[1], s[2])=grid(n)[3] ; 
				f4(s[0], s[1], s[2])=grid(n)[4] ; 
				f5(s[0], s[1], s[2])=grid(n)[5] ; 
				/*f6(s[0], s[1], s[2])=grid(n)[6] ; 
				f7(s[0], s[1], s[2])=grid(n)[7] ; 
				f8(s[0], s[1], s[2])=grid(n)[8] ; 
				f9(s[0], s[1], s[2])=grid(n)[9] ; 
				f10(s[0], s[1], s[2])=grid(n)[10] ; 
				f11(s[0], s[1], s[2])=grid(n)[11] ; 
				f12(s[0], s[1], s[2])=grid(n)[12] ;*/ 
			}	
		
		}
	
	


  
  
		Forward1.fft0(f1,F1);    
		Forward2.fft0(f2,F2); 
		Forward3.fft0(f3,F3); 
		Forward4.fft0(f4,F4); 
		Forward5.fft0(f5,F5); 
		/*Forward6.fft0(f6,F6); 
		Forward7.fft0(f7,F7); 
		Forward8.fft0(f8,F8); 
		Forward9.fft0(f9,F9); 
		Forward10.fft0(f10,F10); 
		Forward11.fft0(f11,F11); 
		Forward12.fft0(f12,F12); */
    

  
  
		for(int n = 0 ; n < nodes(grid) ; n++)
		{
			MMSP::vector<int> x = position(grid, n) ;
			double k1 = fk[x[0]]*0.1 ;
			double k2 = fk[x[1]]*0.1 ;
			double k3 = fk[x[2]]*0.1 ;
			double modk = (sqrt(k1*k1 + k2*k2 + k3*k3)) ;
			if(modk==0.0) modk=0.1   ;
			dfdstr1(x[0],x[1],x[2]) = Bpq(k1, k2, k3, modk, 1, 1)*F1(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 1, 2)*F2(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 1, 3)*F3(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 1, 4)*F4(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 1, 5)*F5(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 1, 6)*F6(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 1, 7)*F7(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 1, 8)*F8(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 1, 9)*F9(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 1, 10)*F10(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 1, 11)*F11(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 1, 12)*F12(x[0],x[1],x[2]) ;
			dfdstr2(x[0],x[1],x[2]) = Bpq(k1, k2, k3, modk, 2, 1)*F1(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 2, 2)*F2(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 2, 3)*F3(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 2, 4)*F4(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 2, 5)*F5(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 2, 6)*F6(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 2, 7)*F7(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 2, 8)*F8(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 2, 9)*F9(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 2, 10)*F10(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 2, 11)*F11(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 2, 12)*F12(x[0],x[1],x[2]) ;
			dfdstr3(x[0],x[1],x[2]) = Bpq(k1, k2, k3, modk, 3, 1)*F1(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 3, 2)*F2(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 3, 3)*F3(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 3, 4)*F4(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 3, 5)*F5(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 3, 6)*F6(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 3, 7)*F7(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 3, 8)*F8(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 3, 9)*F9(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 3, 10)*F10(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 3, 11)*F11(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 3, 12)*F12(x[0],x[1],x[2]) ;
			dfdstr4(x[0],x[1],x[2]) = Bpq(k1, k2, k3, modk, 4, 1)*F1(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 4, 2)*F2(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 4, 3)*F3(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 4, 4)*F4(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 4, 5)*F5(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 4, 6)*F6(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 4, 7)*F7(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 4, 8)*F8(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 4, 9)*F9(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 4, 10)*F10(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 4, 11)*F11(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 4, 12)*F12(x[0],x[1],x[2]) ;
			dfdstr5(x[0],x[1],x[2]) = Bpq(k1, k2, k3, modk, 5, 1)*F1(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 5, 2)*F2(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 5, 3)*F3(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 5, 4)*F4(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 5, 5)*F5(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 5, 6)*F6(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 5, 7)*F7(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 5, 8)*F8(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 5, 9)*F9(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 5, 10)*F10(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 5, 11)*F11(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 5, 12)*F12(x[0],x[1],x[2]) ;
			
			/*dfdstr6(x[0],x[1],x[2]) = Bpq(k1, k2, k3, modk, 6, 1)*F1(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 6, 2)*F2(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 6, 3)*F3(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 6, 4)*F4(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 6, 5)*F5(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 6, 6)*F6(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 6, 7)*F7(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 6, 8)*F8(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 6, 9)*F9(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 6, 10)*F10(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 6, 11)*F11(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 6, 12)*F12(x[0],x[1],x[2]) ;
			dfdstr7(x[0],x[1],x[2]) = Bpq(k1, k2, k3, modk, 7, 1)*F1(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 7, 2)*F2(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 7, 3)*F3(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 7, 4)*F4(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 7, 5)*F5(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 7, 6)*F6(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 7, 7)*F7(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 7, 8)*F8(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 7, 9)*F9(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 7, 10)*F10(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 7, 11)*F11(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 7, 12)*F12(x[0],x[1],x[2]) ;
			dfdstr8(x[0],x[1],x[2]) = Bpq(k1, k2, k3, modk, 8, 1)*F1(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 8, 2)*F2(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 8, 3)*F3(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 8, 4)*F4(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 8, 5)*F5(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 8, 6)*F6(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 8, 7)*F7(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 8, 8)*F8(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 8, 9)*F9(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 8, 10)*F10(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 8, 11)*F11(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 8, 12)*F12(x[0],x[1],x[2]) ;
			dfdstr9(x[0],x[1],x[2]) = Bpq(k1, k2, k3, modk, 9, 1)*F1(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 9, 2)*F2(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 9, 3)*F3(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 9, 4)*F4(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 9, 5)*F5(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 9, 6)*F6(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 9, 7)*F7(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 9, 8)*F8(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 9, 9)*F9(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 9, 10)*F10(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 9, 11)*F11(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 9, 12)*F12(x[0],x[1],x[2]) ;
			dfdstr10(x[0],x[1],x[2]) = Bpq(k1, k2, k3, modk, 10, 1)*F1(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 10, 2)*F2(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 10, 3)*F3(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 10, 4)*F4(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 10, 5)*F5(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 10, 6)*F6(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 10, 7)*F7(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 10, 8)*F8(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 10, 9)*F9(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 10, 10)*F10(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 10, 11)*F11(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 10, 12)*F12(x[0],x[1],x[2]) ;	
			dfdstr11(x[0],x[1],x[2]) = Bpq(k1, k2, k3, modk, 11, 1)*F1(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 11, 2)*F2(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 11, 3)*F3(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 11, 4)*F4(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 11, 5)*F5(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 11, 6)*F6(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 11, 7)*F7(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 11, 8)*F8(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 11, 9)*F9(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 11, 10)*F10(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 11, 11)*F11(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 11, 12)*F12(x[0],x[1],x[2]) ;
			dfdstr12(x[0],x[1],x[2]) = Bpq(k1, k2, k3, modk, 12, 1)*F1(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 12, 2)*F2(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 12, 3)*F3(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 12, 4)*F4(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 12, 5)*F5(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 12, 6)*F6(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 12, 7)*F7(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 12, 8)*F8(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 12, 9)*F9(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 12, 10)*F10(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 12, 11)*F11(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 12, 12)*F12(x[0],x[1],x[2]) ;*/
		
			TotStr(x[0],x[1],x[2]) = Bpq(k1, k2, k3, modk, 1, 1)*F1(x[0],x[1],x[2])*F1(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 1, 2)*F1(x[0],x[1],x[2])*F2(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 1, 3)*F1(x[0],x[1],x[2])*F3(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 2, 2)*F2(x[0],x[1],x[2])*F2(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 2, 3)*F2(x[0],x[1],x[2])*F3(x[0],x[1],x[2]) + Bpq(k1, k2, k3, modk, 3, 3)*F3(x[0],x[1],x[2])*F3(x[0],x[1],x[2]);
			
		}

  
  
		Backward1.fft0Normalized(dfdstr1, dfdstr_real1);  
		Backward2.fft0Normalized(dfdstr2, dfdstr_real2); 
		Backward3.fft0Normalized(dfdstr3, dfdstr_real3); 
		Backward4.fft0Normalized(dfdstr4, dfdstr_real4); 
		Backward5.fft0Normalized(dfdstr5, dfdstr_real5); 
		/*Backward6.fft0Normalized(dfdstr6, dfdstr_real6); 
		Backward7.fft0Normalized(dfdstr7, dfdstr_real7); 
		Backward8.fft0Normalized(dfdstr8, dfdstr_real8); 
		Backward9.fft0Normalized(dfdstr9, dfdstr_real9); 
		Backward10.fft0Normalized(dfdstr10, dfdstr_real10); 
		Backward11.fft0Normalized(dfdstr11, dfdstr_real11); 
		Backward12.fft0Normalized(dfdstr12, dfdstr_real12); */
	
		Backward13.fft0Normalized(TotStr, TotStr_real); 
	}
