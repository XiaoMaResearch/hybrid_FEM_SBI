//
//  get_disp_from_BIE.cpp
//  hybrid_fem_bie
//
//  Created by Max on 2/8/18.
//
//

#include "get_disp_from_BIE.hpp"

void get_disp_from_BIE(ArrayXi &BIE_top_surf_index, ArrayXi &BIE_bot_surf_index, VectorXd &fe_global, int nx, double dx, int Ndofn, InfiniteBoundary & BIE_inf_top, InfiniteBoundary &BIE_inf_bot, NodalField *BIE_top_disp_x_ptr, NodalField *BIE_top_disp_y_ptr, NodalField *BIE_top_vel_x_ptr, NodalField *BIE_top_vel_y_ptr, NodalField *BIE_bot_disp_x_ptr, NodalField *BIE_bot_disp_y_ptr, NodalField *BIE_bot_vel_x_ptr, NodalField *BIE_bot_vel_y_ptr)
{
    VectorXd BIE_top_surf_force = VectorXd::Zero(Ndofn*(nx+1),1);
    VectorXd BIE_bot_surf_force = VectorXd::Zero(Ndofn*(nx+1),1);
    maplocal(BIE_top_surf_index,fe_global,BIE_top_surf_force);
    maplocal(BIE_bot_surf_index,fe_global,BIE_bot_surf_force);
    // Get the BIE surf force
    Map<VectorXd,0,InnerStride<2> > BIE_top_surf_force_x(BIE_top_surf_force.data(), BIE_top_surf_force.size()/2);
    Map<VectorXd,0,InnerStride<2> > BIE_top_surf_force_y(BIE_top_surf_force.data()+1, BIE_top_surf_force.size()/2);
    Map<VectorXd,0,InnerStride<2> > BIE_bot_surf_force_x(BIE_bot_surf_force.data(), BIE_bot_surf_force.size()/2);
    Map<VectorXd,0,InnerStride<2> > BIE_bot_surf_force_y(BIE_bot_surf_force.data()+1, BIE_bot_surf_force.size()/2);
    BIE_top_surf_force_x = -1.0/dx*BIE_top_surf_force_x;
    BIE_top_surf_force_y = -1.0/dx*BIE_top_surf_force_y;
    BIE_bot_surf_force_x = -1.0/dx*BIE_bot_surf_force_x;
    BIE_bot_surf_force_y = -1.0/dx*BIE_bot_surf_force_y;
    // Creat the NodalField pointer point to the BIE surf force. [TOP]
    NodalField *BIE_top_force_x_ptr = new NodalField(nx+1);
    NodalField *BIE_top_force_y_ptr = new NodalField(nx+1);
    BIE_top_force_x_ptr->setValuesTo(BIE_top_surf_force_x.data());
    BIE_top_force_y_ptr->setValuesTo(BIE_top_surf_force_y.data());
    // Creat the NodalField pointer point to the BIE surf force. [BOT]
    NodalField *BIE_bot_force_x_ptr = new NodalField(nx+1);
    NodalField *BIE_bot_force_y_ptr = new NodalField(nx+1);
    BIE_bot_force_x_ptr->setValuesTo(BIE_bot_surf_force_x.data());
    BIE_bot_force_y_ptr->setValuesTo(BIE_bot_surf_force_y.data());
    // Get the surf displacement velocity from BIE
    BIE_inf_top.getDisplacement(BIE_top_force_x_ptr,BIE_top_force_y_ptr,BIE_top_disp_x_ptr,BIE_top_disp_y_ptr,BIE_top_vel_x_ptr,BIE_top_vel_y_ptr);
    BIE_inf_bot.getDisplacement(BIE_bot_force_x_ptr,BIE_bot_force_y_ptr,BIE_bot_disp_x_ptr,BIE_bot_disp_y_ptr,BIE_bot_vel_x_ptr,BIE_bot_vel_y_ptr);

}
