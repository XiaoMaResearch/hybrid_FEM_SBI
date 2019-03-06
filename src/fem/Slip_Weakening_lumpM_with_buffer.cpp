//
//  Slip_Weakening_lumpM_with_buffer.cpp
//  hybrid_fem_bie
//
//  Created by Max on 2/8/18.
//
//

#include "Slip_Weakening_lumpM_with_buffer.hpp"

void Slip_Weakening_lumpM_with_buffer(VectorXd &M_global_vec, ArrayXi &top_surf_index, ArrayXi &bot_surf_index, VectorXd &fe_global, double dt, double dx,double dy, int nx, ArrayXd &delt_v_n, ArrayXd &delt_u_n, ArrayXd &T_0, ArrayXd &tau_s, ArrayXd &mu_s, double mu_d, double Dc, int Ndofn, VectorXd &M_pos, VectorXd &M_neg, VectorXd &F_fault, ArrayXd &T_c, double fault_angle)
{
    
    double c = std::cos(fault_angle);
    double s = std::sin(fault_angle);
//    MatrixXd Rot= MatrixXd::Zero(2, 2);
//    Rot << c , s ,
//          -s , c;
    // TSN formulation Day and Lapusta 2005
    ArrayXd mu_sw = ArrayXd::Zero(nx+1, 1);
    ArrayXd T_s_dot = ArrayXd::Zero(nx+1,1);
    VectorXd T = VectorXd::Zero((nx+1)*Ndofn,1);
    // Area of the element edge on the fault
    double area = dx;
    // Elastic restoring force
    VectorXd R_pos = VectorXd::Zero(Ndofn*(nx+1),1);
    VectorXd R_neg = VectorXd::Zero(Ndofn*(nx+1),1);
    
    maplocal(top_surf_index, -fe_global, R_pos);
    maplocal(bot_surf_index, -fe_global, R_neg);
    
    // Rotation R_pos and R_neg to local fault coordinates using R
    for (int i =0 ; i<R_pos.size()/Ndofn;i++)
    {
        R_pos.data()[2*i] = c*R_pos.data()[2*i] + s*R_pos.data()[2*i+1];
        R_pos.data()[2*i+1] = -s*R_pos.data()[2*i] + c*R_pos.data()[2*i+1];
        
        R_neg.data()[2*i] = c*R_neg.data()[2*i] + s*R_neg.data()[2*i+1];
        R_neg.data()[2*i+1] = -s*R_neg.data()[2*i] + c*R_neg.data()[2*i+1];
    }
    // Sticking force
    for (int i=0; i<T.size()/2; i++)
    {
//        T(2*i) =(1.0/dt*M_neg(i)*M_pos(i)*delt_v_n(2*i)+(M_neg(i)*R_pos(2*i)-M_pos(i)*R_neg(2*i)))/(area*(M_neg(i)+M_pos(i)))+T_0(2*i);
//        T(2*i+1) = (1.0/dt*M_neg(i)*M_pos(i)*(delt_v_n(2*i+1)+1.0/dt*delt_u_n(2*i+1))+(M_neg(i)*R_pos(2*i+1)-M_pos(i)*R_neg(2*i+1)))/(area*(M_neg(i)+M_pos(i)))+T_0(2*i+1);
        
        T(2*i) =(1.0/dt*M_neg(2*i)*M_pos(2*i)*delt_v_n(2*i)+(M_neg(2*i)*R_pos(2*i)-M_pos(2*i)*R_neg(2*i)))/(area*(M_neg(2*i)+M_pos(2*i)))+T_0(2*i);
        T(2*i+1) = (1.0/dt*M_neg(2*i)*M_pos(2*i)*(delt_v_n(2*i+1)+1.0/dt*delt_u_n(2*i+1))+(M_neg(2*i)*R_pos(2*i+1)-M_pos(2*i)*R_neg(2*i+1)))/(area*(M_neg(2*i)+M_pos(2*i)))+T_0(2*i+1);
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
    // Calculate the friction force
    for (int k =0;k<nx+1;k++)
    {
        if (std::abs(delt_u_n(2*k))<Dc)
        {
            tau_s(k) = (mu_s(k)-(mu_s(k)-mu_d)*std::abs(delt_u_n(2*k))/Dc)*(-T_c(2*k+1));
        }
        else
        {
            tau_s(k) = mu_d*(-T_c(2*k+1));
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
    // Rotate F_fault to global cooridantes x y   using R'
    for (int i = 0; i<F_fault.size()/Ndofn; i++)
    {
        F_fault.data()[2*i] = c*F_fault.data()[2*i] - s*F_fault.data()[2*i+1];
        F_fault.data()[2*i+1] = s*F_fault.data()[2*i] + c*F_fault.data()[2*i+1];
    }

    
   }
