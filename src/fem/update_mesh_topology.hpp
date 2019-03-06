//
//  update_mesh_topology.hpp
//  hybrid_fem_bie
//
//  Created by Max on 2/6/18.
//
//

#ifndef update_mesh_topology_hpp
#define update_mesh_topology_hpp

#include <stdio.h>
#include <Eigen/Eigen>
#include "Meshrectangular.hpp"
#include "share_header.hpp"

void update_mesh_topology(MatrixXd &Nodes, MatrixXi_rm &Element,std::vector<int> &fault_surf_nodes, std::vector<int> &fault_surf_nodes_new);
#endif /* update_mesh_topology_hpp */
