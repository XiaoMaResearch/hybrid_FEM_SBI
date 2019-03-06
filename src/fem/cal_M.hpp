//
//  cal_M.hpp
//  hybrid_fem_bie
//
//  Created by Max on 2/10/18.
//
//

#ifndef cal_M_hpp
#define cal_M_hpp

#include <stdio.h>
#include <Eigen/Eigen>
#include <Eigen/KroneckerProduct>

using namespace Eigen;
void cal_M(MatrixXd coord, double density, MatrixXd &M_el);


#endif /* cal_M_hpp */
