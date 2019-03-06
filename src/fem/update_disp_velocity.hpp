//
//  update_disp_velocity.hpp
//  hybrid_fem_bie
//
//  Created by Max on 2/11/18.
//
//

#ifndef update_disp_velocity_hpp
#define update_disp_velocity_hpp

#include <stdio.h>
#include <Eigen/Eigen>
#include "mapglobal.hpp"

using namespace Eigen;

void update_disp_velocity(VectorXd &u_n, VectorXd &v_n, ArrayXi &top_surf_index, ArrayXi &bot_surf_index, int Ndofn, int nx, VectorXd &delt_u_n, VectorXd &delt_v_n);

#endif /* update_disp_velocity_hpp */
