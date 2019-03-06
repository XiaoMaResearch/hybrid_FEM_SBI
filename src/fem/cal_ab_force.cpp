//
//  cal_ab_force.cpp
//  hybrid_fem_bie
//
//  Created by Max on 2/7/18.
//
//
#include "cal_ab_force.hpp"
//#include <omp.h>
//extern double time_fem; // declaration (not definition)

void cal_ab_force(MatrixXd &T_abc,VectorXd &F_total, VectorXi &ab_surf_index, VectorXd &v_n, int Ndofn)
{
 //   #pragma omp parallel for
    for (int i=0;i<ab_surf_index.size()/2-1;i++)
    {
        VectorXd v_n_local=VectorXd::Zero(4,1) ;
        ArrayXi index = ab_surf_index.segment(2*i, 4);
        maplocal(index,v_n,v_n_local);
        VectorXd fe_abc = -T_abc*v_n_local;
        mapglobal(index,F_total,fe_abc);
    }
}
