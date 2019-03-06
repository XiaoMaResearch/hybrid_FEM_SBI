//
//  NR_solve_eq10.cpp
//  hybrid_fem_bie
//
//  Created by Max on 2/10/18.
//
//

#include "NR_solve_eq10.hpp"

//[x_solve] =  NR_solve_eq10(x0,dt,theta_pre(k),delt_v_n(2*k-1),M_pos(k),area,R_pos(2*k-1),R_neg(2*k-1),-T_c(2*k),T_0(2*k-1),a(k),b,L,V_0,f_0);

//function [x]= NR_solve_eq10(x,dt,theta,V_n,M,area,Rx_plus,Rx_minus,sigma,Tx_0,a,b,L,V_0,f_0)

void NR_solve_eq10(double dt, double theta, double V_n, double M, double area, double Rx_plus, double Rx_minus, double sigma, double Tx_0, double a, double b, double L, double V_0, double f_0,double &x)
{
    double t_fric = a* sigma *asinh(0.5*(x+V_n)/(2*V_0)*exp((f_0+b*log(V_0*theta/L))/a));
    double Tx = t_fric;
    
    double F = x - V_n - dt/M*(Rx_plus - Rx_minus - 2*area*Tx+2*area*Tx_0);
    double J = (a*area*dt*sigma*exp((f_0 + b*log((V_0*theta)/L))/a))/(2*M*V_0*pow((exp((2*f_0 + 2*b*log((V_0*theta)/L))/a)*pow(V_n/2 + x/2,2))/(4*pow(V_0,2)) + 1,0.5)) + 1;
    double TOL = 1e-6;
    int iter_max = 100;
    double err = 1.0;
    int iter = 0;
    while (err>TOL&&iter<iter_max)
    {
        double x0 = x ;
        t_fric = a* sigma *asinh(0.5*(x+V_n)/(2*V_0)*exp((f_0+b*log(V_0*theta/L))/a));
        Tx = t_fric;
        F = x - V_n - dt/M*(Rx_plus - Rx_minus - 2*area*Tx+2*area*Tx_0);
        J = (a*area*dt*sigma*exp((f_0 + b*log((V_0*theta)/L))/a))/(2*M*V_0*pow((exp((2*f_0 + 2*b*log((V_0*theta)/L))/a)*pow(V_n/2 + x/2,2))/(4*pow(V_0,2)) + 1,0.5)) + 1;
        x  = x - F/J; 
        err =fabs(x -x0)/fabs(x);
        iter = iter +1 ;
    }
}
