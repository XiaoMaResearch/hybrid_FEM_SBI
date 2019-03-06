//
//  cal_fe_global_const_ke_bimat.hpp
//  hybrid_fem_bie
//
//  Created by Max on 2/7/18.
//
//

#ifndef cal_fe_global_const_ke_bimat_hpp
#define cal_fe_global_const_ke_bimat_hpp

#include <stdio.h>
#include <iostream>
#include <Eigen/Eigen>
#include "maplocal.hpp"
#include "mapglobal.hpp"

using namespace Eigen;

typedef Eigen::Matrix<int, -1, -1,RowMajor> MatrixXi_rm;

void cal_fe_global_const_ke_bimat(MatrixXd &Node, MatrixXi_rm &Element,int n_nodes, int n_el, MatrixXi &index_store, double q, VectorXd &u_n, VectorXd &v_n, int Ndofn, MatrixXd ke_1, MatrixXd ke_2, VectorXd &fe_global);

#endif /* cal_fe_global_const_ke_bimat_hpp */
