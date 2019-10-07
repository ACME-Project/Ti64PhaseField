
//#include "MersenneTwister.h"
//#include "Tensor.hpp"
//#include "Tensor.hpp"

//double eigen_alpha[12][3][3];

using namespace std ;



sfts eigen_alpha[12];


void initialize_alpha_eigen()
{

eigen_alpha[0].e[0] = -4*0.083  ;
eigen_alpha[0].e[5] = 0.0 ; //0.0095  ;
eigen_alpha[0].e[4] = 0.0  ;
//eigen_alpha[0].matrix[1][0] = 0.0095  ;
eigen_alpha[0].e[1] = 4*0.123  ;
eigen_alpha[0].e[3] = 0.0  ;
//eigen_alpha[0].matrix[2][0] = 0.0  ;
//eigen_alpha[0].matrix[2][1] = 0.0  ;
eigen_alpha[0].e[2] = 4*0.035  ;

eigen_alpha[1].e[0] = -2*0.083  ;
eigen_alpha[1].e[5] = 0.0  ;
eigen_alpha[1].e[4] = 2*0.0095  ;
//eigen_alpha[1].matrix[1][0] = 0.0  ;
eigen_alpha[1].e[1] = 2*0.035  ;
eigen_alpha[1].e[3] = 0.0  ;
//eigen_alpha[1].matrix[2][0] = 0.0095  ;
//eigen_alpha[1].matrix[2][1] = 0.0  ;
eigen_alpha[1].e[2] = 2*0.123  ;

eigen_alpha[2].e[0] = 0.079  ;
eigen_alpha[2].e[5] = -0.0359  ;
eigen_alpha[2].e[4] = -0.0264  ;
//eigen_alpha[2].matrix[1][0] = -0.0359  ;
eigen_alpha[2].e[1] = 0.0047  ;
eigen_alpha[2].e[3] = 0.0810  ;
//eigen_alpha[2].matrix[2][0] = -0.0264  ;
//eigen_alpha[2].matrix[2][1] = 0.0810  ;
eigen_alpha[2].e[2] = -0.0087  ;

eigen_alpha[3].e[0] = -0.079  ;
eigen_alpha[3].e[5] = 0.0359  ;
eigen_alpha[3].e[4] = 0.0264  ;
//eigen_alpha[3].matrix[1][0] = 0.0359  ;
eigen_alpha[3].e[1] = 0.0047  ;
eigen_alpha[3].e[3] = 0.0810  ;
//eigen_alpha[3].matrix[2][0] = 0.0264  ;
//eigen_alpha[3].matrix[2][1] = 0.0810  ;
eigen_alpha[3].e[2] = -0.0087  ;

eigen_alpha[4].e[0] = 0.079  ;
eigen_alpha[4].e[5] = -0.0359  ;
eigen_alpha[4].e[4] = 0.0264  ;
//eigen_alpha[4].matrix[1][0] = -0.0359  ;
eigen_alpha[4].e[1] = 0.0047  ;
eigen_alpha[4].e[3] = -0.0810  ;
//eigen_alpha[4].matrix[2][0] = 0.0264  ;
//eigen_alpha[4].matrix[2][1] = -0.0810  ;
eigen_alpha[4].e[2] = -0.0087  ;

eigen_alpha[5].e[0] = 0.079  ;
eigen_alpha[5].e[5] = 0.0359  ;
eigen_alpha[5].e[4] = -0.0264  ;
//eigen_alpha[5].matrix[1][0] = 0.0359  ;
eigen_alpha[5].e[1] = 0.0047  ;
eigen_alpha[5].e[3] = -0.0810  ;
//eigen_alpha[5].matrix[2][0] = -0.0264  ;
//eigen_alpha[5].matrix[2][1] = -0.0810  ;
eigen_alpha[5].e[2] = -0.0087  ;

eigen_alpha[6].e[0] = -0.083  ;
eigen_alpha[6].e[5] = -0.0095  ;
eigen_alpha[6].e[4] = 0.0  ;
//eigen_alpha[6].matrix[1][0] = -0.0095  ;
eigen_alpha[6].e[1] = 0.123  ;
eigen_alpha[6].e[3] = 0.0  ;
//eigen_alpha[6].matrix[2][0] = 0.0  ;
//eigen_alpha[6].matrix[2][1] = 0.0  ;
eigen_alpha[6].e[2] = 0.035  ;

eigen_alpha[7].e[0] = -2*0.083  ;
eigen_alpha[7].e[5] = 0.0  ;
eigen_alpha[7].e[4] = -2*0.0095  ;
//eigen_alpha[7].matrix[1][0] = 0.0  ;
eigen_alpha[7].e[1] = 2*0.035  ;
eigen_alpha[7].e[3] = 0.0  ;
//eigen_alpha[7].matrix[2][0] = -2*0.0095  ;
//eigen_alpha[7].matrix[2][1] = 0.0  ;
eigen_alpha[7].e[2] = 2*0.123  ;

eigen_alpha[8].e[0] = 0.079  ;
eigen_alpha[8].e[5] = -0.0264  ;
eigen_alpha[8].e[4] = -0.0359  ;
//eigen_alpha[8].matrix[1][0] = -0.0264  ;
eigen_alpha[8].e[1] = -0.0087  ;
eigen_alpha[8].e[3] = 0.0810  ;
//eigen_alpha[8].matrix[2][0] = -0.0359 ;
//eigen_alpha[8].matrix[2][1] = 0.0810  ;
eigen_alpha[8].e[2] = 0.0047  ;

eigen_alpha[9].e[0] = 2*0.079  ;
eigen_alpha[9].e[5] = 2*0.0264  ;
eigen_alpha[9].e[4] = 2*0.0359  ;
//eigen_alpha[9].matrix[1][0] = 0.0264  ;
eigen_alpha[9].e[1] = 2*-0.0087  ;
eigen_alpha[9].e[3] = 2*0.0810  ;
//eigen_alpha[9].matrix[2][0] = 0.0359 ;
//eigen_alpha[9].matrix[2][1] = 0.0810  ;
eigen_alpha[9].e[2] = 2*0.0047  ;

eigen_alpha[10].e[0] = 2*0.079  ;
eigen_alpha[10].e[5] = -2*0.0264  ;
eigen_alpha[10].e[4] = 2*0.0359  ;
//eigen_alpha[10].matrix[1][0] = -2*0.0264  ;
eigen_alpha[10].e[1] = -2*0.0087  ;
eigen_alpha[10].e[3] = -2*0.0810  ;
//eigen_alpha[10].matrix[2][0] = 2*0.0359 ;
//eigen_alpha[10].matrix[2][1] = -2*0.0810  ;
eigen_alpha[10].e[2] = 2*0.0047  ;

eigen_alpha[11].e[0] = 2*0.079  ; 
eigen_alpha[11].e[5] = 2*0.0264  ; 
eigen_alpha[11].e[4] = 2*-0.0359  ; 
 eigen_alpha[11].e[1] = 2*-0.0087  ;
 eigen_alpha[11].e[3] = 2*-0.0810  ; 
//eigen_alpha[11].matrix[2][0] = -0.0359 ; 
//eigen_alpha[11].matrix[2][1] = -0.0810  ; 
eigen_alpha[11].e[2] = 2*0.0047  ; }
