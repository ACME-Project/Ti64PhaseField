using namespace std ;
double c[6][6], sigma00[12][6] ;
double c_norm = 1;

double omega_inv_ij(double* cnew, double* k ,  double modk)
{
	double sum = 0 ;
	for(int i = 0 ; i < 9 ; i++)
	{
		sum += cnew[i]*k[i] ;
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



double Bpq(double k1, double k2, double k3, double modk, int p, int q) 
{  	
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
    

