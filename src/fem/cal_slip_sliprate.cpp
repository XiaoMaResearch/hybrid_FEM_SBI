//
//  cal_slip_sliprate.cpp
//  hybrid_fem_bie
//
//  Created by Max on 2/8/18.
//
//

#include "cal_slip_sliprate.hpp"

//void cal_slip_slip_rate(VectorXd &u_n, VectorXd &v_n, ArrayXi &top_surf_index, ArrayXi &bot_surf_index, int Ndofn, int nx, VectorXd &delt_u_n, VectorXd &delt_v_n)
void cal_slip_slip_rate(VectorXd &u_n, VectorXd &v_n, ArrayXi &top_surf_index, ArrayXi &bot_surf_index, int Ndofn, int nx, ArrayXd &delt_u_n, ArrayXd &delt_v_n)

{
    VectorXd u_n_fault_pos = VectorXd::Zero(Ndofn*(nx+1),1);
    VectorXd v_n_fault_pos = VectorXd::Zero(Ndofn*(nx+1),1);
    VectorXd u_n_fault_neg = VectorXd::Zero(Ndofn*(nx+1),1);
    VectorXd v_n_fault_neg = VectorXd::Zero(Ndofn*(nx+1),1);
    // Get the local displacement on the fault
    maplocal(top_surf_index,u_n,u_n_fault_pos);
    maplocal(top_surf_index,v_n,v_n_fault_pos);
    maplocal(bot_surf_index,u_n,u_n_fault_neg);
    maplocal(bot_surf_index,v_n,v_n_fault_neg);
    // Update Slip and slip rate
    delt_u_n = u_n_fault_pos-u_n_fault_neg;
    delt_v_n = v_n_fault_pos-v_n_fault_neg;
}
