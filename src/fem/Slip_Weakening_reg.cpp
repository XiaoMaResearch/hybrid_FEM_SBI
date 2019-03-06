//
//  Slip_Weakening_reg.cpp
//  hybrid_fem_bie
//
//  Created by Max on 2/8/18.
//
//

#include "Slip_Weakening_reg.hpp"

void Slip_Weakening_reg(VectorXd &M_global_vec, ArrayXi &top_surf_index, ArrayXi &bot_surf_index, VectorXd &fe_global, double dt, double dx,double dy, int nx, ArrayXd &delt_v_n, ArrayXd &delt_u_n, ArrayXd &T_0, ArrayXd &tau_s, double mu_s, double mu_d, double Dc, int Ndofn, double M, VectorXd &F_fault, ArrayXd &T_c)
{
    // TSN formulation Day and Lapusta 2005
    ArrayXd mu_sw = ArrayXd::Zero(nx+1, 1);
    // Calculate the friction strenght based on the normal stress at t-dt
    for (int k =0;k<nx+1;k++)
    {
        if (std::abs(delt_u_n(2*k))<Dc)
        {
            mu_sw(k) = mu_s-(mu_s-mu_d)*std::abs(delt_u_n(2*k))/Dc;
        }
        else
        {
            mu_sw(k) = mu_d;
        }
    }
//    for (int k=0; k<nx+1;k++)
//    {
//        tau_s(k) = mu_sw(k)*(-T_c(2*k+1));
//    }
    VectorXd T = VectorXd::Zero((nx+1)*Ndofn,1);
    // Mesh vectors on the two sideds of the fault (Not general)
    VectorXd M_pos = VectorXd::Zero(nx+1,1);
    M_pos(0) = M/4.0;
    M_pos.segment(1,nx)=M/2*VectorXd::Ones(nx,1);
    M_pos(nx)= M/4.0;
    VectorXd M_neg=M_pos;
    // Area of the element edge on the fault
    double area = dx;
    // Elastic restoring force
    VectorXd R_pos = VectorXd::Zero(Ndofn*(nx+1),1);
    VectorXd R_neg = VectorXd::Zero(Ndofn*(nx+1),1);
    maplocal(top_surf_index, -fe_global, R_pos);
    maplocal(bot_surf_index, -fe_global, R_neg);
    // Sticking force
    for (int i=0; i<T.size()/2; i++)
    {
        T(2*i) =(1.0/dt*M_neg(i)*M_pos(i)*delt_v_n(2*i)+(M_neg(i)*R_pos(2*i)-M_pos(i)*R_neg(2*i)))/(area*(M_neg(i)+M_pos(i)))+T_0(2*i);
        T(2*i+1) = (1.0/dt*M_neg(i)*M_pos(i)*(delt_v_n(2*i+1)+1.0/dt*delt_u_n(2*i+1))+(M_neg(i)*R_pos(2*i+1)-M_pos(i)*R_neg(2*i+1)))/(area*(M_neg(i)+M_pos(i)))+T_0(2*i+1);
    }
    // Check if the normal traction is tension (positive)
    for (int k=0; k<nx+1; k++)
    {
        if(T(2*k+1)<0.0)
        {
            T_c(2*k+1) = T(2*k+1);
        }
        else
        {
            T_c(2*k+1) = 0.0;
        }
    }
    // Calculate the tau_s_dot with tau_s(using normals stress at (t-dt) and mu_sw*normal stress at (t))
    double V_star = 1.0;
    double L = 0.1;
    ArrayXd tau_s_dot = ArrayXd::Zero(nx+1,1);
    for (int k =0;k<nx+1;k++)
    {
        tau_s_dot(k) = -(std::abs(delt_u_n(2*k))+V_star)/L*(tau_s(k)-mu_sw(k)*(-T_c(2*k+1)));
       // if (std::abs(tau_s_dot(k))<1.0e-8)
        if (std::abs(tau_s(k)-mu_sw(k)*(-T_c(2*k+1)))>tau_s(k)*0.01)
        {
            tau_s(k) = mu_sw(k)*(-T_c(2*k+1));
        }
        else
        {
            tau_s(k) += tau_s_dot(k)*dt;
        }
    }
    
    // Adding the sign of the sticking force to the friction force
    for (int k =0;k<nx+1;k++)
    {
        if (std::abs(T(2*k))<=tau_s(k))
        {
            T_c(2*k) = T(2*k);
        }
        else
        {
            T_c(2*k) = tau_s(k)*T(2*k)/std::abs(T(2*k));
        }
    }
    F_fault = area*(T_c-T_0);

    
   }
