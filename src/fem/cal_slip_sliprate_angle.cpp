//
//  cal_slip_sliprate_angle.cpp
//  hybrid_fem_bie
//
//  Created by Max on 2/8/18.
//
//

#include "cal_slip_sliprate_angle.hpp"

void cal_slip_slip_rate_angle(VectorXd &u_n, VectorXd &v_n, ArrayXi &top_surf_index, ArrayXi &bot_surf_index, int Ndofn, int nx, ArrayXd &delt_u_n, ArrayXd &delt_v_n,double fault_angle)

{
    double c = std::cos(fault_angle);
    double s = std::sin(fault_angle);
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
    // Rotate delt_u_n delt_v_n cal in global coordianteds to Local fault coordiantes using R
    for (int i=0; i<delt_u_n.size()/Ndofn; i++)
    {
    
        delt_u_n.data()[2*i] = c*delt_u_n.data()[2*i] + s*delt_u_n.data()[2*i+1];
        delt_u_n.data()[2*i+1] = -s*delt_u_n.data()[2*i] + c*delt_u_n.data()[2*i+1];
        
        delt_v_n.data()[2*i] = c*delt_v_n.data()[2*i] + s*delt_v_n.data()[2*i+1];
        delt_v_n.data()[2*i+1] = -s*delt_v_n.data()[2*i] + c*delt_v_n.data()[2*i+1];
        
     
        
    }
    
}
