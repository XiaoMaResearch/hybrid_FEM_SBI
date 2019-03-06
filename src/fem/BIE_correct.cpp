//
//  BIE_correct.cpp
//  hybrid_fem_bie
//
//  Created by Max on 2/6/18.
//
//

#include "BIE_correct.hpp"
using namespace std;

void BIE_correct(const Array<int,Dynamic,1> &BIE_top_surf_index, const Array<int,Dynamic,1> &BIE_bot_surf_index, const VectorXd &fe_global, int Ndofn, int nx, const double dx,InfiniteBoundary &BIE_inf_top, InfiniteBoundary &BIE_inf_bot, VectorXd &u_n, VectorXd &v_n)
{
    VectorXd BIE_top_surf_force = VectorXd::Zero(Ndofn*(nx+1),1);
    VectorXd BIE_bot_surf_force = VectorXd::Zero(Ndofn*(nx+1),1);
    maplocal(BIE_top_surf_index,fe_global,BIE_top_surf_force);
    maplocal(BIE_bot_surf_index,fe_global,BIE_bot_surf_force);
    // Get the traction from force
    BIE_top_surf_force /= -dx;
    BIE_bot_surf_force /= dx;
    // Set the force pointer to be passed into the BIE call
    NodalField *BIE_top_force_x_ptr = new NodalField(nx+1);
    NodalField *BIE_top_force_y_ptr = new NodalField(nx+1);
    NodalField *BIE_bot_force_x_ptr = new NodalField(nx+1);
    NodalField *BIE_bot_force_y_ptr = new NodalField(nx+1);
    BIE_top_force_x_ptr->setValuesTo(BIE_top_surf_force.data(),Ndofn);
    BIE_top_force_y_ptr->setValuesTo(BIE_top_surf_force.data()+1,Ndofn);
    BIE_bot_force_x_ptr->setValuesTo(BIE_bot_surf_force.data(),Ndofn);
    BIE_bot_force_y_ptr->setValuesTo(BIE_bot_surf_force.data()+1,Ndofn);
    // Pointers to for the disp and velocity from BIE
    NodalField *BIE_top_disp_x_ptr = new NodalField(nx+1);
    NodalField *BIE_top_disp_y_ptr = new NodalField(nx+1);
    NodalField *BIE_top_vel_x_ptr  = new NodalField(nx+1);
    NodalField *BIE_top_vel_y_ptr  = new NodalField(nx+1);
    NodalField *BIE_bot_disp_x_ptr = new NodalField(nx+1);
    NodalField *BIE_bot_disp_y_ptr = new NodalField(nx+1);
    NodalField *BIE_bot_vel_x_ptr  = new NodalField(nx+1);
    NodalField *BIE_bot_vel_y_ptr  = new NodalField(nx+1);
    // Call BIE getdisplacement function
    BIE_inf_top.getDisplacement(BIE_top_force_x_ptr,BIE_top_force_y_ptr,BIE_top_disp_x_ptr,BIE_top_disp_y_ptr,BIE_top_vel_x_ptr,BIE_top_vel_y_ptr);
    BIE_inf_bot.getDisplacement(BIE_bot_force_x_ptr,BIE_bot_force_y_ptr,BIE_bot_disp_x_ptr,BIE_bot_disp_y_ptr,BIE_bot_vel_x_ptr,BIE_bot_vel_y_ptr);
    
    // Correct the displacement and velocity on the BIE top and bot surf.
    for (int i=0;i<nx+1;i++)
    {
        // [BIE-TOP]
        u_n(BIE_top_surf_index(Ndofn*i))  = BIE_top_disp_x_ptr->storage()[i];
        u_n(BIE_top_surf_index(Ndofn*i+1))= BIE_top_disp_y_ptr->storage()[i];
        v_n(BIE_top_surf_index(Ndofn*i))  = BIE_top_vel_x_ptr->storage()[i];
        v_n(BIE_top_surf_index(Ndofn*i+1))= BIE_top_vel_y_ptr->storage()[i];
        // [BIE-BOT]
        u_n(BIE_bot_surf_index(Ndofn*i))  = BIE_bot_disp_x_ptr->storage()[i];
        u_n(BIE_bot_surf_index(Ndofn*i+1))= BIE_bot_disp_y_ptr->storage()[i];
        v_n(BIE_bot_surf_index(Ndofn*i))  = BIE_bot_vel_x_ptr->storage()[i];
        v_n(BIE_bot_surf_index(Ndofn*i+1))= BIE_bot_vel_y_ptr->storage()[i];
    }
    
    delete BIE_top_disp_x_ptr;
    delete BIE_top_disp_y_ptr;
    delete BIE_top_vel_x_ptr;
    delete BIE_top_vel_y_ptr;
    delete BIE_bot_disp_x_ptr;
    delete BIE_bot_disp_y_ptr;
    delete BIE_bot_vel_x_ptr;
    delete BIE_bot_vel_y_ptr;
    delete BIE_top_force_x_ptr;
    delete BIE_top_force_y_ptr;
    delete BIE_bot_force_x_ptr;
    delete BIE_bot_force_y_ptr;
    
    
}

