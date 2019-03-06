//
//  mesh_Generated_multi_faults.hpp
//  hybrid_fem_bie
//
//  Created by Max on 2/13/18.
//
//

#ifndef mesh_Generated_multi_faults_hpp
#define mesh_Generated_multi_faults_hpp

#include <stdio.h>
#include <iostream>
#include <Eigen/Eigen>
#include "Meshrectangular.hpp"
#include "update_mesh_topology.hpp"
#include "share_header.hpp"
using namespace Eigen;

void mesh_Generated_multi_faults(int x_min, int x_max ,int y_min, int y_max, double dx , double dy, int nx, int ny, MatrixXd &Nodes ,MatrixXi_rm &Element, ArrayXi &BIE_top_surf_nodes, ArrayXi &BIE_bot_surf_nodes, MatrixXd &fault_pos, std::vector<std::vector<int>> &fault_nodes);



#endif /* mesh_Generated_multi_faults_hpp */
