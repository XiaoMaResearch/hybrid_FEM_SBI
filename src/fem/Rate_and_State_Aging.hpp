//
//  Rate_and_State_Aging.hpp
//  hybrid_fem_bie
//
//  Created by Max on 2/10/18.
//
//

#ifndef Rate_and_State_Aging_hpp
#define Rate_and_State_Aging_hpp

#include <stdio.h>
#include <iostream>
#include <Eigen/Eigen>
#include <math.h>
#include "maplocal.hpp"
#include "NR_solve_eq10.hpp"
using namespace Eigen;

void Rate_and_State_Aging(VectorXd &M_global_vec, ArrayXi &top_surf_index, ArrayXi &bot_surf_index, VectorXd &fe_global, double dt, double dx,double dy, int nx, VectorXd &delt_u_n, VectorXd &delt_v_n,VectorXd &T_0, VectorXd &tau_s, VectorXd a, VectorXd b, double L, double V_0, double f_0, int Ndofn, double M, VectorXd &F_fault, VectorXd &T_c, VectorXd &theta_n, VectorXd &theta_dot_n);

//void Rate_and_State_Aging(VectorXd &M_global_vec, ArrayXi &top_surf_index, ArrayXi &bot_surf_index, VectorXd &fe_global, double dt, double dx,double dy, int nx, ArrayXd &delt_u_n, ArrayXd &delt_v_n, ArrayXd &T_0, ArrayXd &tau_s, VectorXd a, VectorXd b, double L, double V_0, double f_0, int Ndofn, double M, VectorXd &F_fault, ArrayXd &T_c, VectorXd &theta_n, VectorXd &theta_dot_n);

#endif /* Rate_and_State_Aging_hpp */
