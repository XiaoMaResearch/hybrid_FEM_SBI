//
//  cal_ke.hpp
//  hybrid_fem_bie
//
//  Created by Max on 2/7/18.
//
//

#ifndef cal_ke_hpp
#define cal_ke_hpp

#include <stdio.h>
#include <Eigen/Eigen>

using namespace Eigen;

void cal_ke(MatrixXd coord ,double E, double nu, MatrixXd &ke);

#endif /* cal_ke_hpp */
