//
//  Rate_and_State_Aging.cpp
//  hybrid_fem_bie
//
//  Created by Max on 2/10/18.
//
//

#include "Rate_and_State_Aging.hpp"

//function [F_fault,T_c,delt_u_n,delt_v_n,theta_n,theta_dot_n] = Rate_and_State_Aging(M_global_vec,top_surf_index,bot_surf_index,fe_global,Rot,angle,dt,dx,nx,delt_v_n,delt_u_n,T_0,Ndofn,theta_n,theta_dot_n,a,b,L,V_0,f_0)


void Rate_and_State_Aging(VectorXd &M_global_vec, ArrayXi &top_surf_index, ArrayXi &bot_surf_index, VectorXd &fe_global, double dt, double dx,double dy, int nx, VectorXd &delt_u_n, VectorXd &delt_v_n,VectorXd &T_0, VectorXd &tau_s, VectorXd a, VectorXd b, double L, double V_0, double f_0, int Ndofn, double M, VectorXd &F_fault, VectorXd &T_c, VectorXd &theta_n, VectorXd &theta_dot_n)

//void Rate_and_State_Aging(VectorXd &M_global_vec, ArrayXi &top_surf_index, ArrayXi &bot_surf_index, VectorXd &fe_global, double dt, double dx,double dy, int nx, ArrayXd &delt_u_n, ArrayXd &delt_v_n, ArrayXd &T_0, ArrayXd &tau_s, VectorXd a, VectorXd b, double L, double V_0, double f_0, int Ndofn, double M, VectorXd &F_fault, ArrayXd &T_c, VectorXd &theta_n, VectorXd &theta_dot_n)

{
    // TSN formulation Day and Lapusta 2005
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
        if(T(2*k+1)<=0.0)
        {
            T_c(2*k+1) = T(2*k+1);
        }
        else
        {
            T_c(2*k+1) = 0.0;
        }
    }
//    theta_pre = theta_n + dt/2 * theta_dot_n;
//    for k = 1: nx+1
//        x0 = delt_v_n(2*k-1);
//    [x_solve] =  NR_solve_eq10(x0,dt,theta_pre(k),delt_v_n(2*k-1),M_pos(k),area,R_pos(2*k-1),R_neg(2*k-1),-T_c(2*k),T_0(2*k-1),a(k),b,L,V_0,f_0);
//    V_pre(k) = x_solve ;
//    
    VectorXd theta_pre = VectorXd::Zero(nx+1,1);
    VectorXd V_pre = VectorXd::Zero(nx+1,1);
    VectorXd V_cor = VectorXd::Zero(nx+1,1);
    // Rate and state
    // Predictor
    theta_pre = theta_n + dt/2 * theta_dot_n;
    for (int k=0; k<nx+1;k++)
    {
        double x;
        x = delt_v_n(2*k);
        NR_solve_eq10(dt, theta_pre(k), delt_v_n(2*k),M_pos(k), area, R_pos(2*k), R_neg(2*k), -T_c(2*k+1),T_0(2*k), a(k), b(k), L, V_0, f_0, x);
        V_pre(k) = x;
    }
    int iter = 1.0;
    double err = 1.0 ;
    int max_iter = 500;
    double TOL= 1e-8;
    // Corrector
    VectorXd delt_u_n_x = VectorXd::Zero(nx+1, 1);
    VectorXd delt_v_n_x = VectorXd::Zero(nx+1, 1);
    VectorXd V_t_cor = VectorXd::Zero(nx+1,1);
    VectorXd theta_cor = VectorXd::Zero(nx+1,1);
    for(int i=0;i<nx+1;i++)
    {
        delt_u_n_x(i) = delt_u_n(2*i);
        delt_v_n_x(i) = delt_v_n(2*i);
    }
    while (err>TOL&& iter<max_iter)
   // while (iter<max_iter)
    {
        VectorXd V_t_pre = 0.5*(V_pre+delt_v_n_x);
        VectorXd theta_t_dot_pre = VectorXd::Ones(nx+1,1)-theta_pre.cwiseProduct(V_t_pre)/L;
        theta_cor = theta_n + 0.5*(theta_dot_n+theta_t_dot_pre)*dt/2;
        for (int k = 0; k<nx+1;k++)
        {
            double x;
            x = delt_v_n(2*k);
            NR_solve_eq10(dt, theta_cor(k), delt_v_n(2*k),M_pos(k), area, R_pos(2*k), R_neg(2*k), -T_c(2*k+1),T_0(2*k), a(k), b(k), L, V_0, f_0, x);
            V_cor(k) = x;
        }
        V_t_cor = 0.5*(V_cor+delt_v_n_x);
        VectorXd temp_max_V_t = V_t_cor.cwiseMax(V_t_pre);
        VectorXd temp_max_theta_t = theta_cor.cwiseMax(theta_pre);
//        err_V_t=max(abs((V_t_cor-V_t_pre)/temp_max_V_t));
        double err_V_T = ((V_t_cor-V_t_pre).cwiseQuotient(temp_max_V_t)).cwiseAbs().maxCoeff();
        double err_theta_t = ((theta_cor-theta_pre).cwiseQuotient(temp_max_theta_t)).cwiseAbs().maxCoeff();
        err = std::max(err_V_T,err_theta_t);
        
        std::cout<<"err="<<err<<std::endl;
        iter = iter+1;
        V_pre = V_cor;
        theta_pre = theta_cor;
    }
    std::cout<<"*************** iter = \n"<<iter<<std::endl;

 //   theta_t_dot = 1.0 - V_t_cor.*theta_cor/L;
 //   theta_new = theta_n+ theta_t_dot*dt ;
 //   theta_dot_new = 1.0- V_cor(:).*theta_new/L;
  //  delt_u_new_x  = delt_u_n(1:2:end) + + V_cor(:) *dt;
//    delt_u_n(1:2:end) = delt_u_new_x;
//    delt_v_n(1:2:end) = V_cor(:);
//    theta_n = theta_new ;
//    theta_dot_n = theta_dot_new;
    
    VectorXd theta_t_dot = VectorXd::Ones(nx+1,1)-V_t_cor.cwiseProduct(theta_cor)/L;
    VectorXd theta_new = theta_n + theta_t_dot*dt;
    VectorXd theta_dot_new = VectorXd::Ones(nx+1,1)-V_cor.cwiseProduct(theta_new)/L;
    VectorXd delt_u_new_x = delt_u_n_x + V_cor*dt;
    for(int i=0;i<nx+1;i++)
    {
        delt_u_n(2*i) = delt_u_new_x(i);
        delt_v_n(2*i) = V_cor(i);
    }
    theta_n = theta_new ;
    theta_dot_n = theta_dot_new;
    for(int i=0;i<nx+1;i++)
    {
        T_c(2*i) = a(i)*(-T_c(2*i+1))*asinh(delt_v_n(2*i)/(2*V_0)*exp((f_0+b(i)*log(V_0*theta_n(i)/L))/a(i)));
    }
    F_fault = area*(T_c-T_0);
    
}
