//
//  time_advance.cpp
//  hybrid_fem_bie
//
//  Created by Max on 2/8/18.
//
//

#include "time_advance.hpp"

void time_advance(VectorXd &u_n, VectorXd &v_n, VectorXd &F_total, VectorXd &M_global_vec, double dt)
{
    // Central differecne time integration
    VectorXd  a_n = F_total.cwiseQuotient(M_global_vec);
    VectorXd v_new = v_n+dt*a_n;
    VectorXd u_new = u_n+dt*v_new;
    v_n = v_new;
    u_n = u_new;
}
