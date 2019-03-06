//
//  BIE_correct.hpp
//  hybrid_fem_bie
//
//  Created by Max on 2/6/18.
//
//

#ifndef BIE_correct_hpp
#define BIE_correct_hpp

#include <stdio.h>
#include <Eigen/Eigen>
#include "get_disp_from_BIE.hpp"
using namespace Eigen;

void BIE_correct(const Array<int,Dynamic,1> &BIE_top_surf_index, const Array<int,Dynamic,1> &BIE_bot_surf_index, const VectorXd &fe_global, int Ndofn, int nx, const double dx,InfiniteBoundary &BIE_inf_top, InfiniteBoundary &BIE_inf_bot, VectorXd &u_n, VectorXd &v_n);

#endif /* BIE_correct_hpp */
