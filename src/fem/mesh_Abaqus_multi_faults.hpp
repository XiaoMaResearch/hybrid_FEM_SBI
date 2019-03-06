//
//  mesh_Abaqus_multi_faults.hpp
//  hybrid_fem_bie
//
//  Created by Max on 2/13/18.
//
//

#ifndef mesh_Abaqus_multi_faults_hpp
#define mesh_Abaqus_multi_faults_hpp

#include <stdio.h>
#include <iostream>
#include <Eigen/Eigen>
#include "Meshrectangular.hpp"
#include "update_mesh_topology.hpp"
#include "share_header.hpp"
#include "Abaqus_read.hpp"
using namespace Eigen;

void mesh_Abaqus_multi_faults(std::string &fname,MatrixXd &Node ,MatrixXi_rm &Element,
                              std::vector<int> &BIE_top, std::vector<int> &BIE_bot,
                              std::vector<std::vector<int>> &fault_nodes, double num_faults);
#endif /* mesh_Abaqus_multi_faults_hpp */
