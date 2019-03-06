//
//  cal_BTDB.hpp
//  hybrid_fem_bie
//
//  Created by Max on 2/6/18.
//
//

#ifndef cal_BTDB_hpp
#define cal_BTDB_hpp

#include <stdio.h>
#include <Eigen/Eigen>
using namespace Eigen;
void cal_BTDB (MatrixXd coord,double E, double nu,MatrixXd &B_T_D_B);


#endif /* cal_BTDB_hpp */
