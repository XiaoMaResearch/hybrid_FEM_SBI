//
//  cal_M_global_vec.cpp
//  hybrid_fem_bie
//
//  Created by Max on 2/10/18.
//
//
#include <iostream>
#include "cal_M_global_vec.hpp"
void cal_M_global_vec(MatrixXd &Node, MatrixXi &Element, double density, MatrixXi &index_store, int Ndofn, VectorXd &M_global_vec)
{
    int n_el = Element.rows();
    int Nnel = Element.cols();
    MatrixXd M_el = MatrixXd::Zero(Ndofn*Nnel,Ndofn*Nnel);
    for (int i=0; i<n_el; i++)
    {
        MatrixXd coord = MatrixXd::Zero(Nnel,Ndofn);
        VectorXi Element_0= Element.row(i);
        coord.row(0) = Node.row(Element_0(0));
        coord.row(1) = Node.row(Element_0(1));
        coord.row(2) = Node.row(Element_0(2));
        coord.row(3) = Node.row(Element_0(3));
        cal_M(coord, density, M_el);
        VectorXd M_el_vec =M_el.rowwise().sum();
        mapglobal(index_store.col(i),M_global_vec,M_el_vec);

    }
}
