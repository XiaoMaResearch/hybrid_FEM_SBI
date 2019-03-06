//
//  Slip_Weakening_lumpM_with_buffer.hpp
//  hybrid_fem_bie
//
//  Created by Max on 2/8/18.
//
//

#ifndef Slip_Weakening_lumpM_with_buffer_hpp
#define Slip_Weakening_lumpM_with_buffer_hpp

#include <stdio.h>
#include <iostream>
#include <Eigen/Eigen>
#include "maplocal.hpp"

void Slip_Weakening_lumpM_with_buffer(VectorXd &M_global_vec, ArrayXi &top_surf_index, ArrayXi &bot_surf_index, VectorXd &fe_global, double dt, double dx,double dy, int nx, ArrayXd &delt_v_n, ArrayXd &delt_u_n, ArrayXd &T_0, ArrayXd &tau_s, ArrayXd &mu_s, double mu_d, double Dc, int Ndofn, VectorXd &M_pos, VectorXd &M_neg, VectorXd &F_fault, ArrayXd &T_c, double fault_angle);
#endif /* Slip_Weakening_lumpM_with_buffer_hpp */
