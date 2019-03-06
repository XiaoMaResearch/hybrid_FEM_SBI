//
//  cal_fe_global_const_ke_bimat.cpp
//  hybrid_fem_bie
//
//  Created by Max on 2/7/18.
//
//

#include "cal_fe_global_const_ke_bimat.hpp"


void cal_fe_global_const_ke_bimat(MatrixXd &Node, MatrixXi_rm &Element,int n_nodes, int n_el, MatrixXi &index_store, double q, VectorXd &u_n, VectorXd &v_n, int Ndofn, MatrixXd ke_1, MatrixXd ke_2, VectorXd &fe_global)
{
    
    
    for (int i=0;i<n_el;i++)
    {
        MatrixXd coord = MatrixXd::Zero(4,2);
        VectorXi Element_0= Element.row(i);
        double y1 =Node.row(Element_0(0))(1);
        double y2 = Node.row(Element_0(2))(1);
        double y_avg = 0.5*(y1+y2);
        MatrixXd ke = MatrixXd::Zero(8,8);
        if (y_avg>0.0)
        {
            //   ke = ke_background;
            ke = ke_2;
            
        }
        else
        {
            //   ke = ke_inclusion;
            ke = ke_1;
            
        }
        VectorXd u_n_local=VectorXd::Zero(8,1) ;
        VectorXd v_n_local=VectorXd::Zero(8,1) ;
        ArrayXi index = index_store.col(i);
        //VectorXi index = index_store.col(i);
        maplocal(index,u_n,u_n_local);
        maplocal(index,v_n,v_n_local);
        VectorXd fe_int = ke*(u_n_local+q*v_n_local);
        mapglobal(index,fe_global,fe_int);
    }
}
