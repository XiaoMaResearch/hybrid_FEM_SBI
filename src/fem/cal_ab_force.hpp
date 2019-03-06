//
//  cal_ab_force.hpp
//  hybrid_fem_bie
//
//  Created by Max on 2/7/18.
//
//

#ifndef cal_ab_force_hpp
#define cal_ab_force_hpp

#include <stdio.h>
#include <iostream>
#include <Eigen/Eigen>
#include "maplocal.hpp"
#include "mapglobal.hpp"

using namespace Eigen;
extern double time_fem;
void cal_ab_force(MatrixXd &T_abc,VectorXd &F_total, VectorXi &ab_surf_index, VectorXd &v_n, int Ndofn);
#endif /* cal_ab_force_hpp */
