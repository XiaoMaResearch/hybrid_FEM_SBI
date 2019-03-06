//
//  update_disp_velocity.cpp
//  hybrid_fem_bie
//
//  Created by Max on 2/11/18.
//
//

#include "update_disp_velocity.hpp"

//void update_disp_velocity(VectorXd &u_n, VectorXd &v_n, ArrayXi &top_surf_index, ArrayXi &bot_surf_index, int Ndofn, int nx, ArrayXd &delt_u_n, ArrayXd &delt_v_n)

void update_disp_velocity(VectorXd &u_n, VectorXd &v_n, ArrayXi &top_surf_index, ArrayXi &bot_surf_index, int Ndofn, int nx, VectorXd &delt_u_n, VectorXd &delt_v_n)
{
    for (int i=0;i<nx+1;i++)
    {
        u_n(top_surf_index(2*i)) = delt_u_n(2*i)/2.0;
        u_n(bot_surf_index(2*i)) = -delt_u_n(2*i)/2.0;
        v_n(top_surf_index(2*i)) = delt_v_n(2*i)/2.0;
        v_n(bot_surf_index(2*i)) = -delt_v_n(2*i)/2.0;
    }
}
