//
//  cal_BTDB.cpp
//  hybrid_fem_bie
//
//  Created by Max on 2/6/18.
//
//

#include "cal_BTDB.hpp"
void cal_BTDB (MatrixXd coord ,double E, double nu, MatrixXd &B_T_D_B)
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
            //MatrixXd Jac=MatrixXd::Zero(2,2);
            MatrixXd Jac=MatrixXd::Zero(2,2);
            Jac = dNdxi*coord;
            //cout << "dNdxi:\n" << dNdxi << "\ncoord:\n" << coord << "\n\n";
            //cout<<Jac<<"\n\n";
            // Here is the problem
            //  Jac << 0.5,0,
            //         0 , 0.5;
            
            
            double detJ = Jac.determinant();
            MatrixXd dNdx = MatrixXd::Zero(2,4);
            dNdx = Jac.inverse()*dNdxi;
            MatrixXd B= MatrixXd::Zero(3,8);
            B<< dNdx(0,0),         0, dNdx(0,1),         0, dNdx(0,2),         0,  dNdx(0,3),         0,
            0, dNdx(1,0),         0, dNdx(1,1),        0 , dNdx(1,2),          0,  dNdx(1,3),
            dNdx(1,0), dNdx(0,0), dNdx(1,1), dNdx(0,1), dNdx(1,2), dNdx(0,2),  dNdx(1,3),  dNdx(0,3);
            B_T_D_B.block(k,0,8,8)= B.transpose()*D*B*detJ*wgp;
            
            k=k+8;
            
            
        }
    }
}

