//
//  cal_fe_global_const_ke.cpp
//  hybrid_fem_bie
//
//  Created by Max on 2/7/18.
//
//
#include "cal_fe_global_const_ke.hpp"
//#include <omp.h>
//extern double time_fem; // declaration (not definition)

void cal_fe_global_const_ke(int n_nodes, int n_el, MatrixXi &index_store, double q, VectorXd &u_n, VectorXd &v_n, int Ndofn, MatrixXd ke, VectorXd &fe_global)
{  
  //int tid;
    //int nthreads;
  //  double start = omp_get_wtime();
    #pragma omp parallel for
    for (int i=0;i<n_el;i++)
    {
        VectorXd u_n_local=VectorXd::Zero(8,1) ;
        VectorXd v_n_local=VectorXd::Zero(8,1) ;
        ArrayXi index = index_store.col(i);
        //VectorXi index = index_store.col(i);
        maplocal(index,u_n,u_n_local);
        maplocal(index,v_n,v_n_local);
        VectorXd fe_int = ke*(u_n_local+q*v_n_local);
        mapglobal(index,fe_global,fe_int);
       // tid = omp_get_thread_num(); 
       /* printf("wrting from thread = %d\n", tid);
        if (tid==0)
	{
	   nthreads = omp_get_num_threads();
        printf("Number of threads = %d\n", nthreads);
	}
       */
    }   
  // double end = omp_get_wtime();


 //   time_fem +=(end-start);
}
