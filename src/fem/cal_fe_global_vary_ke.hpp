//
//  cal_fe_global_vary_ke.hpp
//  hybrid_fem_bie
//
//  Created by Max on 2/7/18.
//
//

#ifndef cal_fe_global_vary_ke_hpp
#define cal_fe_global_vary_ke_hpp

#include <stdio.h>
#include <iostream>
#include <Eigen/Eigen>
#include "maplocal.hpp"
#include "mapglobal.hpp"

using namespace Eigen;
extern double time_fem;
void cal_fe_global_vary_ke(int n_nodes, int n_el, MatrixXi &index_store, double q, VectorXd &u_n, VectorXd &v_n, int Ndofn,
                           std::vector<MatrixXd> &ke, VectorXd &fe_global);

#endif /* cal_fe_global_vary_ke_hpp */
