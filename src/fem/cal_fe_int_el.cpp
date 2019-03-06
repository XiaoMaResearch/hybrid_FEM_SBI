//
//  cal_fe_int_el.cpp
//  hybrid_fem_bie
//
//  Created by Max on 3/21/18.
//
//

#include "cal_fe_int_el.hpp"
#include <iostream>

void cal_fe_int_el(double E, double nu, std::vector<Eigen::MatrixXd> &B_mat, double &detJ,VectorXd &u_n_local, VectorXd &fe_int_el,
                   MatrixXd &strain_n, MatrixXd &stress_n,double sxx_initial,
                   double syy_initial, double sxy_initial, double cohes, double blkfric, double &eq_ep, double &E_p)

{
    
    std::vector<MatrixXd> de_p(4);
   // VectorXd de_p = VectorXd::Zero(3, 1);

    
    for (int i=0;i<4;i++)
    {
        VectorXd strain = B_mat[i]*u_n_local;
       // VectorXd stress = D*strain;
        VectorXd stress = VectorXd::Zero(3,1);
        VectorXd strain_temp = strain_n.col(i);
        VectorXd stress_temp = stress_n.col(i);
        Getplastic(E, nu, strain, stress, stress_temp, strain_temp,sxx_initial,syy_initial,sxy_initial,cohes,blkfric,de_p[i]);
        fe_int_el = fe_int_el+B_mat[i].transpose()*stress*detJ;
        //fe_int_el = fe_int_el+B_mat[i].transpose()*D*B_mat[i]*u_n_local*detJ;
        strain_n.col(i) = strain;
        stress_n.col(i) = stress;
    }
    // Calculating the equivalent plastic strain
    for (int i=0;i<4;i++)
    {
       // std::cout<<de_p[i].rows()<<de_p[i].cols()<<std::endl;
        double em = (de_p[i](0,0)+de_p[i](1,0))/3.0;
        double exx = de_p[i](0,0)-em;
        double eyy = de_p[i](1,0)-em;
        double exy = de_p[i](2,0);
        
        double Y2 = 0.5*(exx*exx+eyy*eyy+2*exy*exy);
        eq_ep += 1.0/4.0*(std::pow(4.0/3.0*Y2,0.5));
        
        E_p += ((stress_n(0,i)+sxx_initial)*de_p[i](0,0)+(stress_n(1,i)+syy_initial)*de_p[i](1,0)+
                        (stress_n(2,i)+sxy_initial)*de_p[i](2,0))*detJ;
        //std::cout<<eq_ep<<std::endl;
    }
    
    
}
