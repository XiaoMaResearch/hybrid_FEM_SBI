//
//  cal_fe_global_plastic.hpp
//  hybrid_fem_bie
//
//  Created by Max on 3/21/18.
//
//

#ifndef cal_fe_global_plastic_hpp
#define cal_fe_global_plastic_hpp

#include <stdio.h>
#include <Eigen/Eigen>
#include "mapglobal.hpp"
#include "maplocal.hpp"
#include "cal_fe_int_el.hpp"
using namespace Eigen;

void cal_fe_global_plastic(int n_el, MatrixXi &index_store,MatrixXd &ke, std::vector<MatrixXd> &B_mat, double detJ, double E, double nu, double q, VectorXd &u_n, VectorXd &v_n, int Ndofn, VectorXd &fe_global, std::vector<MatrixXd> &strain_n_store,std::vector<MatrixXd> &stress_n_store,double sxx_initial,double syy_initial, double sxy_initial, double cohes, double blkfric, std::vector<double> &eq_ep_out,double &E_p);

#endif /* cal_fe_global_plastic_hpp */
