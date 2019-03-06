//
//  cal_Bmat.cpp
//  hybrid_fem_bie
//
//  Created by Max on 3/21/18.
//
//

#include "cal_Bmat.hpp"

void cal_Bmat(MatrixXd coord ,double E, double nu, std::vector<Eigen::MatrixXd> &B_mat, double &detJ)
{
    MatrixXd D=MatrixXd::Zero(3,3);
    D << 1-nu, nu , 0,
    nu,   1-nu, 0,
    0,   0 ,  (1-2*nu)/2;
    D = E/((1+nu)*(1-2*nu))*D;
    VectorXd xi_gp=VectorXd::Zero(2,1);
    VectorXd eta_gp=VectorXd::Zero(2,1);
    VectorXd w_gp=VectorXd::Zero(2,1);
    xi_gp<<-1.0/sqrt(3),1.0/sqrt(3);
    eta_gp<<-1.0/sqrt(3),1.0/sqrt(3);
    w_gp<<1.0,1.0;
    int k = 0;
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
            detJ = Jac.determinant();
            MatrixXd dNdx = MatrixXd::Zero(2,4);
            dNdx = Jac.inverse()*dNdxi;
            MatrixXd B= MatrixXd::Zero(3,8);
            B<< dNdx(0,0),         0, dNdx(0,1),         0, dNdx(0,2),         0,  dNdx(0,3),         0,
                0, dNdx(1,0),         0, dNdx(1,1),        0 , dNdx(1,2),          0,  dNdx(1,3),
                dNdx(1,0), dNdx(0,0), dNdx(1,1), dNdx(0,1), dNdx(1,2), dNdx(0,2),  dNdx(1,3),  dNdx(0,3);
            
            B_mat[k] = B;
            k= k+1;
            
        }
    }
}


