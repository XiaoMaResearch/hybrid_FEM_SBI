//
//  cal_slip_sliprate.hpp
//  hybrid_fem_bie
//
//  Created by Max on 2/8/18.
//
//

#ifndef cal_slip_sliprate_hpp
#define cal_slip_sliprate_hpp

#include <stdio.h>
#include <Eigen/Eigen>
#include "maplocal.hpp"

using namespace Eigen;
void cal_slip_slip_rate(VectorXd &u_n, VectorXd &v_n, ArrayXi &top_surf_index, ArrayXi &bot_surf_index, int Ndofn, int nx, ArrayXd &delt_u_n, ArrayXd &delt_v_n);

#endif /* cal_slip_sliprate_hpp */
