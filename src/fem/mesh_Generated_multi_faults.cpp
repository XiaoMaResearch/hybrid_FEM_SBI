//
//  mesh_Generated_multi_faults.cpp
//  hybrid_fem_bie
//
//  Created by Max on 2/13/18.
//
//

#include "mesh_Generated_multi_faults.hpp"
//#include <omp.h>
#include <numeric>
#include <vector>
void mesh_Generated_multi_faults(int x_min, int x_max ,int y_min, int y_max, double dx , double dy, int nx, int ny, MatrixXd &Nodes ,MatrixXi_rm &Element, ArrayXi &BIE_top_surf_nodes, ArrayXi &BIE_bot_surf_nodes, MatrixXd &fault_pos, std::vector<std::vector<int>> &fault_nodes)
{
    Meshrectangular(x_min, x_max ,y_min, y_max, dx, dx, nx, ny, Nodes, Element);
    std::cout<< "Done"<<std::endl;    
// Getting fault nodes
    //std::vector<int> fault_surf_nodes_new;
    double TOL = 1e-13;
    //#pragma omp parallel for
    for (int j = 0 ; j<fault_pos.rows(); j++)
    {
     double x_L = fault_pos(j,1);
     double x_R = fault_pos(j,2);
     int    n_fault = (x_R-x_L)/dx+1;
     double   y = fault_pos(j,0);
     std::vector<int> fault_temp(n_fault);
     std::iota(fault_temp.begin(),fault_temp.end(),(y-y_min)/dy*(nx+1)+(x_L-x_min)/dx);
   //  std::vector<int> fault_nodes_temp
    // fault_nodes[2*j].insert
     auto it = fault_nodes[2*j].begin();
     fault_nodes[2*j].insert(it,fault_temp.begin(),fault_temp.end());
     std::cout<<"Done_finding_fault_one"<<std::endl;
     update_mesh_topology(Nodes,Element,fault_nodes[2*j],fault_nodes[2*j+1]);
     std::cout<<"Done_adding_fault_one"<<std::endl;
    }
    std::cout<<"Done_finding_nodes"<<std::endl;
    // Getting the BIE boundary nodes
    int m = 0 , n = 0;
    for(int i=0; i<Nodes.rows(); i++)
    {
        if (std::abs(Nodes(i,1)-y_max)<TOL)
        {
            BIE_top_surf_nodes(m) = i ;
            m+=1;
        };
        if (std::abs(Nodes(i,1)-y_min)<TOL)
        {
            BIE_bot_surf_nodes(n) = i ;
            n+=1;
        };
    }
}
