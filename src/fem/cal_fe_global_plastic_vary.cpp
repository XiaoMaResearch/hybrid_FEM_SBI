//
//  cal_fe_global_plastic_vary.cpp
//  hybrid_fem_bie
//
//  Created by Max on 3/21/18.
//
//

#include "cal_fe_global_plastic_vary.hpp"


void cal_fe_global_plastic_vary(int n_el, MatrixXi &index_store,std::vector<MatrixXd> &ke,std::vector<std::vector<MatrixXd>> &B_mat,
                                std::vector<double> &detJ, double E, double nu, double q, VectorXd &u_n, VectorXd &v_n, int Ndofn, VectorXd &fe_global, std::vector<MatrixXd> &strain_n_store,std::vector<MatrixXd> &stress_n_store,double sxx_initial,double syy_initial, double sxy_initial, double cohes, double blkfric,std::vector<double> &eq_ep_out,double &E_p)
{
    for (int i=0;i<n_el;i++)
    {
        VectorXd u_n_local=VectorXd::Zero(8,1) ;
        VectorXd v_n_local=VectorXd::Zero(8,1) ;
        ArrayXi index = index_store.col(i);
        //VectorXi index = index_store.col(i);
        maplocal(index,u_n,u_n_local);
        maplocal(index,v_n,v_n_local);
        VectorXd fe_int = VectorXd::Zero(8, 1);
        double E_p_temp=0.0;
        cal_fe_int_el(E, nu, B_mat[i], detJ[i], u_n_local, fe_int,strain_n_store[i],stress_n_store[i],sxx_initial,syy_initial,sxy_initial,cohes,blkfric,eq_ep_out[i],E_p_temp);
        //eq_ep_out.push_back(eq_ep);
        fe_int  = fe_int + ke[i]*q*v_n_local;
        //VectorXd fe_int = ke*(u_n_local+q*v_n_local);
        mapglobal(index,fe_global,fe_int);
        E_p +=E_p_temp;
    }

}
