//
//  cal_fe_int_el.hpp
//  hybrid_fem_bie
//
//  Created by Max on 3/21/18.
//
//

#ifndef cal_fe_int_el_hpp
#define cal_fe_int_el_hpp

#include <stdio.h>
#include <Eigen/Eigen>
#include "Getplastic.hpp"
using namespace Eigen;

void cal_fe_int_el(double E, double nu, std::vector<Eigen::MatrixXd> &B_mat, double &detJ,VectorXd &u_n_local, VectorXd &fe_int_el,
                   MatrixXd &strain_n, MatrixXd &stress_n,double sxx_initial,
                   double syy_initial, double sxy_initial, double cohes, double blkfric, double &eq_ep, double &E_p);



#endif /* cal_fe_int_el_hpp */
