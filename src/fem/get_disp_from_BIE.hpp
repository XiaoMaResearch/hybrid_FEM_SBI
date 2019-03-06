//
//  get_disp_from_BIE.hpp
//  hybrid_fem_bie
//
//  Created by Max on 2/8/18.
//
//

#ifndef get_disp_from_BIE_hpp
#define get_disp_from_BIE_hpp

#include <stdio.h>
#include <iostream>
#include <Eigen/Eigen>
#include "maplocal.hpp"
#include "infinite_boundary.hh"


using namespace Eigen;
void get_disp_from_BIE(ArrayXi &BIE_top_surf_index, ArrayXi &BIE_bot_surf_index, VectorXd &fe_global, int nx, double dx, int Ndofn, InfiniteBoundary & BIE_inf_top, InfiniteBoundary &BIE_inf_bot, NodalField *BIE_top_disp_x_ptr, NodalField *BIE_top_disp_y_ptr, NodalField *BIE_top_vel_x_ptr, NodalField *BIE_top_vel_y_ptr, NodalField *BIE_bot_disp_x_ptr, NodalField *BIE_bot_disp_y_ptr, NodalField *BIE_bot_vel_x_ptr, NodalField *BIE_bot_vel_y_ptr);

#endif /* get_disp_from_BIE_hpp */
