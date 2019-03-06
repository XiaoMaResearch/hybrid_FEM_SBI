//
//  cal_T_abc.cpp
//  hybrid_fem_bie
//
//  Created by Max on 2/7/18.
//
//

#include "cal_T_abc.hpp"

void cal_T_abc(MatrixXd coord ,double density, MatrixXd &T_abc, double v_s, double v_p)
{
    T_abc = MatrixXd::Zero(4, 4);
    VectorXd xi_gp=VectorXd::Zero(2,1);
    VectorXd eta_gp=VectorXd::Zero(2,1);
    VectorXd w_gp=VectorXd::Zero(2,1);
    xi_gp<<-1.0/sqrt(3),1.0/sqrt(3);
    eta_gp<<-1.0/sqrt(3),1.0/sqrt(3);
    w_gp<<1.0,1.0;
    for (int i=0; i<2;i++)
    {
        double xi = xi_gp(i);
        double wgp = w_gp(i);
        // Jacobian Matrix (2x2)
        MatrixXd dNdxi=MatrixXd::Zero(1,2);
        dNdxi << 0.5, -0.5;
        MatrixXd Jac=MatrixXd::Zero(2,2);
        Jac = dNdxi*coord;
        double detJ = Jac.norm();
        MatrixXd dNdx = MatrixXd::Zero(2,4);
        MatrixXd N_shape = MatrixXd::Zero(2,4);
        double N1 = 0.5*(1+xi);
        double N2 = 0.5*(1-xi);
        N_shape<< N1, 0.0 , N2 , 0,
                   0 , N1 , 0 , N2 ;
        T_abc += density*N_shape.transpose()*N_shape*detJ*wgp;
    }
    T_abc.row(0) *= v_p;
    T_abc.row(1) *= v_s;
    T_abc.row(2) *= v_p;
    T_abc.row(3) *= v_s;
}


