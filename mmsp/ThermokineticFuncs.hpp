#include "ThermokineticDefs.hpp"

double Q_V_Alpha(double cV, double cAl) 
{
	return (Q_HCP_V_Al*cAl + Q_HCP_V_V*cV + Q_HCP_V_Ti*(1-cAl-cV)) ; // + cAl*cV*(A_HCP_V_Al_V_0 + A_HCP_V_Al_V_1*(cAl-cV)) + cAl*(1-cAl-cV)*(A_HCP_V_Al_Ti_0 + A_HCP_V_Al_Ti_1*(2*cAl+cV-1)) +cV*(1-cAl-cV)*(A_HCP_V_V_Ti_0 + A_HCP_V_V_Ti_1*(2*cAl+cV-1))) ;
}

double Q_V_Beta(double cV, double cAl){
	return (Q_BCC_V_Al*cAl + Q_BCC_V_V*cV + Q_BCC_V_Ti*(1-cAl-cV)) ; //+ cAl*cV*(A_BCC_V_Al_V_0 + A_BCC_V_Al_V_1*(cAl-cV)) + cAl*(1-cAl-cV)*(A_BCC_V_Al_Ti_0 + A_BCC_V_Al_Ti_1*(2*cAl+cV-1)) +cV*(1-cAl-cV)*(A_BCC_V_V_Ti_0 + A_BCC_V_V_Ti_1*(2*cAl+cV-1))) ; 
}
			
double Q_Al_Alpha(double cV, double cAl) {
	return (Q_HCP_Al_Al*cAl + Q_HCP_Al_V*cV + Q_HCP_Al_Ti*(1-cAl-cV)) ; //  + cAl*cV*(A_HCP_Al_Al_V_0 + A_HCP_Al_Al_V_1*(cAl-cV)) + cAl*(1-cAl-cV)*(A_HCP_Al_Al_Ti_0 + A_HCP_Al_Al_Ti_1*(2*cAl+cV-1)) +cV*(1-cAl-cV)*(A_HCP_Al_V_Ti_0 + A_HCP_Al_V_Ti_1*(2*cAl+cV-1))) ; 
}
			
double Q_Al_Beta(double cV, double cAl) {
	return (Q_BCC_Al_Al*cAl + Q_BCC_Al_V*cV + Q_BCC_Al_Ti*(1-cAl-cV)) ; //+ cAl*cV*(A_BCC_Al_Al_V_0 + A_BCC_Al_Al_V_1*(cAl-cV)) + cAl*(1-cAl-cV)*(A_BCC_Al_Al_Ti_0 + A_BCC_Al_Al_Ti_1*(2*cAl+cV-1)) +cV*(1-cAl-cV)*(A_BCC_Al_V_Ti_0 + A_BCC_Al_V_Ti_1*(2*cAl+cV-1))) ; 
}

double Q_Ti_Alpha(double cV, double cAl) {
	return (Q_HCP_Al_Al*cAl + Q_HCP_Al_V*cV + Q_HCP_Al_Ti*(1-cAl-cV)) ; //  + cAl*cV*(A_HCP_Al_Al_V_0 + A_HCP_Al_Al_V_1*(cAl-cV)) + cAl*(1-cAl-cV)*(A_HCP_Al_Al_Ti_0 + A_HCP_Al_Al_Ti_1*(2*cAl+cV-1)) +cV*(1-cAl-cV)*(A_HCP_Al_V_Ti_0 + A_HCP_Al_V_Ti_1*(2*cAl+cV-1))) ; 
}
			
double Q_Ti_Beta(double cV, double cAl) {
	return (Q_BCC_Al_Al*cAl + Q_BCC_Al_V*cV + Q_BCC_Al_Ti*(1-cAl-cV)) ; //+ cAl*cV*(A_BCC_Al_Al_V_0 + A_BCC_Al_Al_V_1*(cAl-cV)) + cAl*(1-cAl-cV)*(A_BCC_Al_Al_Ti_0 + A_BCC_Al_Al_Ti_1*(2*cAl+cV-1)) +cV*(1-cAl-cV)*(A_BCC_Al_V_Ti_0 + A_BCC_Al_V_Ti_1*(2*cAl+cV-1))) ; 
}

double M_Al_Alpha(double cV, double cAl){
	return (M0/(R*T))*exp(-Q_Al_Alpha(cV, cAl)/(R*T)) ;
}

double M_Al_Beta(double cV, double cAl){
	return (M0/(R*T))*exp(-Q_Al_Beta(cV, cAl)/(R*T)) ;
}

double M_V_Alpha(double cV, double cAl){
	return (M0/(R*T))*exp(-Q_V_Alpha(cV, cAl)/(R*T)) ;
}

double M_V_Beta(double cV, double cAl){
	return (M0/(R*T))*exp(-Q_V_Beta(cV, cAl)/(R*T)) ;
}

double M_Al(double cV, double cAl){
	return 	((M0/(R*T))*(exp(-Q_Al_Alpha(cV,cAl)/(R*T)))) + ((M0/(R*T))*(exp(-Q_Al_Beta(cV,cAl)/(R*T)))) - pow(((M0/(R*T))*(exp(-Q_Al_Alpha(cV,cAl)/(R*T)))), eta) - pow(((M0/(R*T))*(exp(-Q_Al_Beta(cV,cAl)/(R*T)))), (1-eta)) ;
}

double M_V(double cV, double cAl){
	return 	((M0/(R*T))*(exp(-Q_V_Alpha(cV,cAl)/(R*T)))) + ((M0/(R*T))*(exp(-Q_V_Beta(cV,cAl)/(R*T)))) - pow(((M0/(R*T))*(exp(-Q_V_Alpha(cV,cAl)/(R*T)))), eta) - pow(((M0/(R*T))*(exp(-Q_V_Beta(cV,cAl)/(R*T)))), (1-eta)) ;
}

double M_Ti(double cV, double cAl){
	return 	((M0/(R*T))*(exp(-Q_Ti_Alpha(cV,cAl)/(R*T)))) + ((M0/(R*T))*(exp(-Q_Ti_Beta(cV,cAl)/(R*T)))) - pow(((M0/(R*T))*(exp(-Q_Ti_Alpha(cV,cAl)/(R*T)))), eta) - pow(((M0/(R*T))*(exp(-Q_Ti_Beta(cV,cAl)/(R*T)))), (1-eta)) ;
}
			
//double Q_Ti_Alpha(double cV, double cAl)

//double Q_Ti_Beta(double cV, double cAl)

double d_M_V_Alpha_d_C_V(double cV, double cAl){
	return (M0/(pow((R*T),2)))*(exp(-Q_V_Alpha(cV,cAl)/(R*T)))*(-Q_HCP_V_V + Q_HCP_V_Ti) ; } // - cAl*A_HCP_V_Al_V_0 - (1-cAl-cV)*A_HCP_V_Al_Ti_0 + cV*A_HCP_V_Al_Ti_0 + cAl*A_HCP_V_Al_Ti_0) ; }
	
double d_M_V_Alpha_d_C_Al(double cV, double cAl){
	return (M0/(pow((R*T),2)))*(exp(-Q_V_Alpha(cV,cAl)/(R*T)))*(-Q_HCP_V_Al + Q_HCP_V_Ti) ; }// cV*A_HCP_V_Ti_V_0 - (1-cAl-cV)*A_HCP_V_Al_Ti_0 + cAl*A_HCP_V_Al_Ti_0 - cV*A_HCP_V_Al_V_0) ; }

double d_M_V_Beta_d_C_V(double cV, double cAl){
	return (M0/(pow((R*T),2)))*(exp(-Q_V_Beta(cV,cAl)/(R*T)))*(-Q_BCC_V_V + Q_BCC_V_Ti); } // - cAl*A_BCC_V_Al_V_0 - (1-cAl-cV)*A_BCC_V_Al_Ti_0 + cV*A_BCC_V_Al_Ti_0 + cAl*A_BCC_V_Al_Ti_0) ; }

double d_M_V_Beta_d_C_Al(double cV, double cAl){
	return (M0/(pow((R*T),2)))*(exp(-Q_V_Beta(cV,cAl)/(R*T)))*(-Q_BCC_V_Al + Q_BCC_V_Ti); } // + cV*A_BCC_V_Ti_V_0 - (1-cAl-cV)*A_BCC_V_Al_Ti_0 + cAl*A_BCC_V_Al_Ti_0 - cV*A_BCC_V_Al_V_0) ; }
	
	
double d_M_Al_Alpha_d_C_V(double cV, double cAl) {
	return (M0/(pow((R*T),2)))*(exp(-Q_Al_Alpha(cV,cAl)/(R*T)))*(-Q_HCP_Al_V + Q_HCP_Al_Ti);}// - (1-cV-cAl)*(A_HCP_Al_V_Ti_0) + cV*A_HCP_Al_V_Ti_0 + cAl*A_HCP_Al_Al_Ti_0 - cAl*A_HCP_Al_Al_V_0) ; }
	
double d_M_Al_Alpha_d_C_Al(double cV, double cAl) {
	return (M0/(pow((R*T),2)))*(exp(-Q_Al_Alpha(cV,cAl)/(R*T)))*(-Q_HCP_Al_Al + Q_HCP_Al_Ti);}//- (1-cV-cAl)*(A_HCP_Al_Al_Ti_0) + cV*A_HCP_Al_V_Ti_0 + cAl*A_HCP_Al_Al_Ti_0 - cV*A_HCP_Al_Al_V_0) ; }
	
double d_M_Al_Beta_d_C_V(double cV, double cAl)  {
	return (M0/(pow((R*T),2)))*(exp(-Q_Al_Alpha(cV,cAl)/(R*T)))*(-Q_BCC_Al_V + Q_BCC_Al_Ti);}// - (1-cV-cAl)*(A_BCC_Al_V_Ti_0) + cV*A_BCC_Al_V_Ti_0 + cAl*A_BCC_Al_Al_Ti_0 - cAl*A_BCC_Al_Al_V_0) ; }
	
double d_M_Al_Beta_d_C_Al(double cV, double cAl) {
	return (M0/(pow((R*T),2)))*(exp(-Q_Al_Alpha(cV,cAl)/(R*T)))*(-Q_BCC_Al_Al + Q_BCC_Al_Ti);}// - (1-cV-cAl)*(A_BCC_Al_Al_Ti_0) + cV*A_BCC_Al_V_Ti_0 + cAl*A_BCC_Al_Al_Ti_0 - cV*A_BCC_Al_Al_V_0) ; }
	
//double d_M_Ti_Alpha_d_C_V(double cV, double cAl) 
	

//double d_M_Ti_Alpha_d_C_Al(double cV, double cAl)

//double d_M_Ti_Beta_d_C_V(double cV, double cAl)

//double d_M_Ti_Beta_d_C_Al(double cV, double cAl)
	
double d_M_V_d_C_V(double cV, double cAl){
	double coeff1 = eta*pow(M_V_Alpha(cV,cAl),(eta-1))*d_M_V_Alpha_d_C_V(cV,cAl)*pow(M_V_Beta(cV, cAl), (1-eta));
	double coeff2 = pow(M_V_Alpha(cV,cAl),(eta))*d_M_V_Beta_d_C_V(cV,cAl)*(1-eta)*pow(M_V_Beta(cV, cAl), (-eta));
	//double coeff1 = eta*pow( ((M0/(pow((R*T),2)))*(exp(-Q_V_Alpha(cV,cAl)/(R*T)))) , (eta-1))*d_M_V_Alpha_d_C_V*pow( ((M0/(pow((R*T),2)))*(exp(-Q_V_Beta(cV,cAl)/(R*T)))) , (1-eta)) ;
	//double coeff2 = (1-eta)*pow( ((M0/(pow((R*T),2)))*(exp(-Q_V_Beta(cV,cAl)/(R*T)))) , (-eta))*d_M_V_Beta_d_C_V*pow( ((M0/(pow((R*T),2)))*(exp(-Q_V_Alpha(cV,cAl)/(R*T)))) , eta) ;
	return d_M_V_Alpha_d_C_V(cV,cAl) + d_M_V_Beta_d_C_V(cV, cAl) - coeff1 - coeff2 ; }
	
double d_M_V_d_C_Al(double cV, double cAl){
	double coeff1 = eta*pow(M_V_Alpha(cV,cAl),(eta-1))*d_M_V_Alpha_d_C_Al(cV,cAl)*pow(M_V_Beta(cV, cAl), (1-eta));
	double coeff2 = pow(M_V_Alpha(cV,cAl),(eta))*d_M_V_Beta_d_C_Al(cV,cAl)*(1-eta)*pow(M_V_Beta(cV, cAl), (-eta));
	//double coeff1 = eta*pow( ((M0/(pow((R*T),2)))*(exp(-Q_V_Alpha(cV,cAl)/(R*T)))) , (eta-1))*d_M_V_Alpha_d_C_Al*pow( ((M0/(pow((R*T),2)))*(exp(-Q_V_Beta(cV,cAl)/(R*T)))) , (1-eta)) ;
	//double coeff2 = (1-eta)*pow( ((M0/(pow((R*T),2)))*(exp(-Q_V_Beta(cV,cAl)/(R*T)))) , (-eta))*d_M_V_Beta_d_C_Al*pow( ((M0/(pow((R*T),2)))*(exp(-Q_V_Alpha(cV,cAl)/(R*T)))) , eta) ;
	return d_M_V_Alpha_d_C_Al(cV,cAl) + d_M_V_Beta_d_C_Al(cV, cAl) - coeff1 - coeff2 ; }
	
double d_M_Al_d_C_V(double cV, double cAl){
	double coeff1 = eta*pow(M_Al_Alpha(cV,cAl),(eta-1))*d_M_Al_Alpha_d_C_V(cV,cAl)*pow(M_Al_Beta(cV, cAl), (1-eta));
	double coeff2 = pow(M_Al_Alpha(cV,cAl),(eta))*d_M_Al_Beta_d_C_V(cV,cAl)*(1-eta)*pow(M_Al_Beta(cV, cAl), (-eta));
	//double coeff1 = eta*pow( ((M0/(pow((R*T),2)))*(exp(-Q_Al_Alpha(cV,cAl)/(R*T)))) , (eta-1))*d_M_Al_Alpha_d_C_V*pow( ((M0/(pow((R*T),2)))*(exp(-Q_Al_Beta(cV,cAl)/(R*T)))) , (1-eta)) ;
	//double coeff2 = (1-eta)*pow( ((M0/(pow((R*T),2)))*(exp(-Q_Al_Beta(cV,cAl)/(R*T)))) , (-eta))*d_M_Al_Beta_d_C_V*pow( ((M0/(pow((R*T),2)))*(exp(-Q_Al_Alpha(cV,cAl)/(R*T)))) , eta) ;
	return d_M_Al_Alpha_d_C_V(cV,cAl) + d_M_Al_Beta_d_C_V(cV, cAl) - coeff1 - coeff2 ; }
	
double d_M_Al_d_C_Al(double cV, double cAl) {
	double coeff1 = eta*pow(M_Al_Alpha(cV,cAl),(eta-1))*d_M_Al_Alpha_d_C_Al(cV,cAl)*pow(M_Al_Beta(cV, cAl), (1-eta));
	double coeff2 = pow(M_Al_Alpha(cV,cAl),(eta))*d_M_Al_Beta_d_C_Al(cV,cAl)*(1-eta)*pow(M_Al_Beta(cV, cAl), (-eta));
	//double coeff1 = eta*pow( ((M0/(pow((R*T),2)))*(exp(-Q_Al_Alpha(cV,cAl)/(R*T)))) , (eta-1))*d_M_Al_Alpha_d_C_Al*pow( ((M0/(pow((R*T),2)))*(exp(-Q_Al_Beta(cV,cAl)/(R*T)))) , (1-eta)) ;
	//double coeff2 = (1-eta)*pow( ((M0/(pow((R*T),2)))*(exp(-Q_Al_Beta(cV,cAl)/(R*T)))) , (-eta))*d_M_Al_Beta_d_C_Al*pow( ((M0/(pow((R*T),2)))*(exp(-Q_V_Alpha(cV,cAl)/(R*T)))) , eta) ;
	return d_M_Al_Alpha_d_C_Al(cV,cAl) + d_M_Al_Beta_d_C_Al(cV, cAl) - coeff1 - coeff2 ; }
	
double d_M_Ti_d_C_Al(double cV, double cAl){
		double coeff1 = eta*pow(M_Al_Alpha(cV,cAl),(eta-1))*d_M_Al_Alpha_d_C_V(cV,cAl)*pow(M_Al_Beta(cV, cAl), (1-eta));
	double coeff2 = pow(M_Al_Alpha(cV,cAl),(eta))*d_M_Al_Beta_d_C_V(cV,cAl)*(1-eta)*pow(M_Al_Beta(cV, cAl), (-eta));
	//double coeff1 = eta*pow( ((M0/(pow((R*T),2)))*(exp(-Q_Al_Alpha(cV,cAl)/(R*T)))) , (eta-1))*d_M_Al_Alpha_d_C_V*pow( ((M0/(pow((R*T),2)))*(exp(-Q_Al_Beta(cV,cAl)/(R*T)))) , (1-eta)) ;
	//double coeff2 = (1-eta)*pow( ((M0/(pow((R*T),2)))*(exp(-Q_Al_Beta(cV,cAl)/(R*T)))) , (-eta))*d_M_Al_Beta_d_C_V*pow( ((M0/(pow((R*T),2)))*(exp(-Q_Al_Alpha(cV,cAl)/(R*T)))) , eta) ;
	return d_M_Al_Alpha_d_C_V(cV,cAl) + d_M_Al_Beta_d_C_V(cV, cAl) - coeff1 - coeff2 ; }
	

double d_M_Ti_d_C_V(double cV, double cAl){
		double coeff1 = eta*pow(M_Al_Alpha(cV,cAl),(eta-1))*d_M_Al_Alpha_d_C_V(cV,cAl)*pow(M_Al_Beta(cV, cAl), (1-eta));
	double coeff2 = pow(M_Al_Alpha(cV,cAl),(eta))*d_M_Al_Beta_d_C_V(cV,cAl)*(1-eta)*pow(M_Al_Beta(cV, cAl), (-eta));
	//double coeff1 = eta*pow( ((M0/(pow((R*T),2)))*(exp(-Q_Al_Alpha(cV,cAl)/(R*T)))) , (eta-1))*d_M_Al_Alpha_d_C_V*pow( ((M0/(pow((R*T),2)))*(exp(-Q_Al_Beta(cV,cAl)/(R*T)))) , (1-eta)) ;
	//double coeff2 = (1-eta)*pow( ((M0/(pow((R*T),2)))*(exp(-Q_Al_Beta(cV,cAl)/(R*T)))) , (-eta))*d_M_Al_Beta_d_C_V*pow( ((M0/(pow((R*T),2)))*(exp(-Q_Al_Alpha(cV,cAl)/(R*T)))) , eta) ;
	return d_M_Al_Alpha_d_C_V(cV,cAl) + d_M_Al_Beta_d_C_V(cV, cAl) - coeff1 - coeff2 ; }
	

	
double dMAlAl_dcV(double cV, double cAl) {
	return  cAl*cAl*M_V(cV, cAl) - cAl*cAl*M_Ti(cV,cAl) + pow((1-cAl),2)*cAl*d_M_Al_d_C_V(cV,cAl) + cAl*cAl*cV*d_M_V_d_C_V(cV,cAl) + cAl*cAl*(1-cAl-cV)*d_M_Ti_d_C_V(cV,cAl) ; 
}
    
double dMAlAl_dcAl(double cV, double cAl) {
	return (-2*(1-cAl)*cAl*M_Al(cV,cAl) + pow((1-cAl),2)*M_Al(cV, cAl) + 2*cAl*cV*M_V(cV, cAl) + 2*cAl*(1-cAl-cV)*M_Ti(cV,cAl) - cAl*cAl*M_Ti(cV,cAl) +
			pow((1-cAl),2)*cAl*d_M_Al_d_C_Al(cV,cAl) + cAl*cAl*cV*d_M_V_d_C_Al(cV,cAl) + cAl*cAl*(1-cAl-cV)*d_M_Ti_d_C_Al(cV,cAl)) ; 
}

double dMVV_dcV(double cV, double cAl) {
	return (2*cV*cAl*M_Al(cV,cAl) + 2*(1-cV)*cV*M_V(cV,cAl) + pow((1-cV),2)*M_V(cV,cAl) + 2*cV*(1-cV-cAl)*M_Ti(cV,cAl) - cV*cV*M_Ti(cV,cAl) +
			cV*cV*cAl*d_M_Al_d_C_V(cV,cAl) + cV*pow((1-cV),2)*d_M_V_d_C_V(cV,cAl) + cV*cV*(1-cAl-cV)*d_M_Ti_d_C_V(cV,cAl)) ; 
}
	
	
double dMVV_dcAl(double cV, double cAl) {
	return  cV*cV*M_Al(cV, cAl)  - cV*cV*M_Ti(cV,cAl) + cV*cV*cAl*d_M_Al_d_C_Al(cV,cAl) + pow((1-cV),2)*cV*d_M_V_d_C_Al(cV,cAl) + cV*cV*(1-cV-cAl)*d_M_Ti_d_C_Al(cV,cAl)  ;
}
 
double dMAlV_dcV(double cV, double cAl) {
	return (-(1-cAl)*cAl*M_Al(cV,cAl) - (cAl)*(1-cV)*M_V(cV,cAl) + (cAl)*(1)*cV*M_V(cV,cAl) + cAl*(1-cAl-cV)*M_Ti(cV,cAl) - cAl*cV*M_Ti(cV,cAl) +
		   (1-cAl)*cAl*(-cV)*d_M_Al_d_C_V(cV,cAl) + (-cAl)*(1-cV)*cV*d_M_V_d_C_V(cV,cAl) + cV*cAl*(1-cV-cAl)*d_M_Ti_d_C_V(cV,cAl)) ;
}
double dMAlV_dcAl(double cV, double cAl) {
	return (cAl*cV*M_Al(cV,cAl) - (1-cAl)*cV*M_Al(cV,cAl) - (1-cV)*cV*M_V(cV,cAl) + cV*(1-cAl-cV)*M_Ti(cV,cAl) - cV*cAl*M_Ti(cV,cAl) +
		   (1-cAl)*cAl*(-cV)*d_M_Al_d_C_Al(cV,cAl) + (-cAl)*(1-cV)*cV*d_M_V_d_C_Al(cV,cAl) + cV*cAl*(1-cV-cAl)*d_M_Ti_d_C_V(cV,cAl)) ;
}

		
