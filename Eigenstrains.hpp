
//#include "MersenneTwister.h"
//#include "Tensor.hpp"
//#include "Tensor.hpp"

//double eigen_alpha[12][3][3];

using namespace std ;



tensor eigen_alpha(3)[12];


void initialize_alpha_eigen()
{

eigen_alpha[0].matrix[0][0] = -0.083  ;
eigen_alpha[0].matrix[0][1] = 0.0095  ;
eigen_alpha[0].matrix[0][2] = 0.0  ;
eigen_alpha[0].matrix[1][0] = 0.0095  ;
eigen_alpha[0].matrix[1][1] = 0.123  ;
eigen_alpha[0].matrix[1][2] = 0.0  ;
eigen_alpha[0].matrix[2][0] = 0.0  ;
eigen_alpha[0].matrix[2][1] = 0.0  ;
eigen_alpha[0].matrix[2][2] = 0.035  ;

eigen_alpha[1].matrix[0][0] = -0.083  ;
eigen_alpha[1].matrix[0][1] = 0.0  ;
eigen_alpha[1].matrix[0][2] = 0.0095  ;
eigen_alpha[1].matrix[1][0] = 0.0  ;
eigen_alpha[1].matrix[1][1] = 0.035  ;
eigen_alpha[1].matrix[1][2] = 0.0  ;
eigen_alpha[1].matrix[2][0] = 0.0095  ;
eigen_alpha[1].matrix[2][1] = 0.0  ;
eigen_alpha[1].matrix[2][2] = 0.123  ;

eigen_alpha[2].matrix[0][0] = 0.079  ;
eigen_alpha[2].matrix[0][1] = -0.0359  ;
eigen_alpha[2].matrix[0][2] = -0.0264  ;
eigen_alpha[2].matrix[1][0] = -0.0359  ;
eigen_alpha[2].matrix[1][1] = 0.0047  ;
eigen_alpha[2].matrix[1][2] = 0.0810  ;
eigen_alpha[2].matrix[2][0] = -0.0264  ;
eigen_alpha[2].matrix[2][1] = 0.0810  ;
eigen_alpha[2].matrix[2][2] = -0.0087  ;

eigen_alpha[3].matrix[0][0] = -0.079  ;
eigen_alpha[3].matrix[0][1] = 0.0359  ;
eigen_alpha[3].matrix[0][2] = 0.0264  ;
eigen_alpha[3].matrix[1][0] = 0.0359  ;
eigen_alpha[3].matrix[1][1] = 0.0047  ;
eigen_alpha[3].matrix[1][2] = 0.0810  ;
eigen_alpha[3].matrix[2][0] = 0.0264  ;
eigen_alpha[3].matrix[2][1] = 0.0810  ;
eigen_alpha[3].matrix[2][2] = -0.0087  ;

eigen_alpha[4].matrix[0][0] = 0.079  ;
eigen_alpha[4].matrix[0][1] = -0.0359  ;
eigen_alpha[4].matrix[0][2] = 0.0264  ;
eigen_alpha[4].matrix[1][0] = -0.0359  ;
eigen_alpha[4].matrix[1][1] = 0.0047  ;
eigen_alpha[4].matrix[1][2] = -0.0810  ;
eigen_alpha[4].matrix[2][0] = 0.0264  ;
eigen_alpha[4].matrix[2][1] = -0.0810  ;
eigen_alpha[4].matrix[2][2] = -0.0087  ;

eigen_alpha[5].matrix[0][0] = 0.079  ;
eigen_alpha[5].matrix[0][1] = 0.0359  ;
eigen_alpha[5].matrix[0][2] = -0.0264  ;
eigen_alpha[5].matrix[1][0] = 0.0359  ;
eigen_alpha[5].matrix[1][1] = 0.0047  ;
eigen_alpha[5].matrix[1][2] = -0.0810  ;
eigen_alpha[5].matrix[2][0] = -0.0264  ;
eigen_alpha[5].matrix[2][1] = -0.0810  ;
eigen_alpha[5].matrix[2][2] = -0.0087  ;

eigen_alpha[6].matrix[0][0] = -0.083  ;
eigen_alpha[6].matrix[0][1] = -0.0095  ;
eigen_alpha[6].matrix[0][2] = 0.0  ;
eigen_alpha[6].matrix[1][0] = -0.0095  ;
eigen_alpha[6].matrix[1][1] = 0.123  ;
eigen_alpha[6].matrix[1][2] = 0.0  ;
eigen_alpha[6].matrix[2][0] = 0.0  ;
eigen_alpha[6].matrix[2][1] = 0.0  ;
eigen_alpha[6].matrix[2][2] = 0.035  ;

eigen_alpha[7].matrix[0][0] = -0.083  ;
eigen_alpha[7].matrix[0][1] = 0.0  ;
eigen_alpha[7].matrix[0][2] = -0.0095  ;
eigen_alpha[7].matrix[1][0] = 0.0  ;
eigen_alpha[7].matrix[1][1] = 0.035  ;
eigen_alpha[7].matrix[1][2] = 0.0  ;
eigen_alpha[7].matrix[2][0] = -0.0095  ;
eigen_alpha[7].matrix[2][1] = 0.0  ;
eigen_alpha[7].matrix[2][2] = 0.123  ;

eigen_alpha[8].matrix[0][0] = 0.079  ;
eigen_alpha[8].matrix[0][1] = -0.0264  ;
eigen_alpha[8].matrix[0][2] = -0.0359  ;
eigen_alpha[8].matrix[1][0] = -0.0264  ;
eigen_alpha[8].matrix[1][1] = -0.0087  ;
eigen_alpha[8].matrix[1][2] = 0.0810  ;
eigen_alpha[8].matrix[2][0] = -0.0359 ;
eigen_alpha[8].matrix[2][1] = 0.0810  ;
eigen_alpha[8].matrix[2][2] = 0.0047  ;

eigen_alpha[9].matrix[0][0] = 0.079  ;
eigen_alpha[9].matrix[0][1] = 0.0264  ;
eigen_alpha[9].matrix[0][2] = 0.0359  ;
eigen_alpha[9].matrix[1][0] = 0.0264  ;
eigen_alpha[9].matrix[1][1] = -0.0087  ;
eigen_alpha[9].matrix[1][2] = 0.0810  ;
eigen_alpha[9].matrix[2][0] = 0.0359 ;
eigen_alpha[9].matrix[2][1] = 0.0810  ;
eigen_alpha[9].matrix[2][2] = 0.0047  ;

eigen_alpha[10].matrix[0][0] = 0.079  ;
eigen_alpha[10].matrix[0][1] = -0.0264  ;
eigen_alpha[10].matrix[0][2] = 0.0359  ;
eigen_alpha[10].matrix[1][0] = -0.0264  ;
eigen_alpha[10].matrix[1][1] = -0.0087  ;
eigen_alpha[10].matrix[1][2] = -0.0810  ;
eigen_alpha[10].matrix[2][0] = 0.0359 ;
eigen_alpha[10].matrix[2][1] = -0.0810  ;
eigen_alpha[10].matrix[2][2] = 0.0047  ;

eigen_alpha[11].matrix[0][0] = 0.079  ;
eigen_alpha[11].matrix[0][1] = 0.0264  ;
eigen_alpha[11].matrix[0][2] = -0.0359  ;
eigen_alpha[11].matrix[1][0] = 0.0264  ;
eigen_alpha[11].matrix[1][1] = -0.0087  ;
eigen_alpha[11].matrix[1][2] = -0.0810  ;
eigen_alpha[11].matrix[2][0] = -0.0359 ;
eigen_alpha[11].matrix[2][1] = -0.0810  ;
eigen_alpha[11].matrix[2][2] = 0.0047  ;
}
