//
//  cal_T_abc.hpp
//  hybrid_fem_bie
//
//  Created by Max on 2/7/18.
//
//

#ifndef cal_T_abc_hpp
#define cal_T_abc_hpp

#include <stdio.h>
#include <Eigen/Eigen>

using namespace Eigen;

void cal_T_abc(MatrixXd coord ,double density, MatrixXd &T_abc, double v_p, double v_s);

#endif /* cal_T_abc_hpp */
