//
//  cal_M.cpp
//  hybrid_fem_bie
//
//  Created by Max on 2/10/18.
//
//

#include "cal_M.hpp"


void cal_M(MatrixXd coord, double density, MatrixXd &M_el)
{
    MatrixXd M = MatrixXd::Zero(4,4);
    VectorXd xi_gp=VectorXd::Zero(2,1);
    VectorXd eta_gp=VectorXd::Zero(2,1);
    VectorXd w_gp=VectorXd::Zero(2,1);
    xi_gp<<-1.0/sqrt(3),1.0/sqrt(3);
    eta_gp<<-1.0/sqrt(3),1.0/sqrt(3);
    w_gp<<1.0,1.0;
    for (int i=0; i<2;i++)
    {
        for (int j=0;j<2;j++)
        {
            double xi = xi_gp(i);
            double eta = eta_gp(j);
            double wgp = w_gp(j);
            // Jacobian Matrix (2x2)
            MatrixXd dNdxi=MatrixXd::Zero(2,4);
            dNdxi << eta-1.0, 1.0-eta, 1.0+eta, -eta-1.0,
                     xi-1.0 , -xi-1.0, xi+1.0, 1.0-xi;
            dNdxi=dNdxi/4.0;
            MatrixXd Jac=MatrixXd::Zero(2,2);
            Jac = dNdxi*coord;
            double detJ = Jac.determinant();
            MatrixXd dNdx = MatrixXd::Zero(2,4);
            dNdx = Jac.inverse()*dNdxi;
            VectorXd N(4);
            N << (1-xi)*(1-eta),(1+xi)*(1-eta),(1+xi)*(1+eta),(1-xi)*(1+eta);
            N /= 4.0;
            M +=(density*(N*N.transpose())*detJ*wgp);

        }
    }
    MatrixXd eye = MatrixXd::Identity(2, 2);
    M_el = kroneckerProduct(M,eye);
    
}

