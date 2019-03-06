//
//  time_advance.hpp
//  hybrid_fem_bie
//
//  Created by Max on 2/8/18.
//
//

#ifndef time_advance_hpp
#define time_advance_hpp

#include <stdio.h>
#include <Eigen/Eigen>
using namespace Eigen;

void time_advance(VectorXd &u_n, VectorXd &v_n, VectorXd &F_total, VectorXd &M_global_vec, double dt);
#endif /* time_advance_hpp */
