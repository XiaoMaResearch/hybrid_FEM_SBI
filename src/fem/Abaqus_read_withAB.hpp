//
//  Abaqus_read_withAB.hpp
//  hybrid_fem_bie
//
//  Created by Max on 3/23/18.
//
//

#ifndef Abaqus_read_withAB_hpp
#define Abaqus_read_withAB_hpp

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <Eigen/Eigen>
using namespace Eigen;
typedef Eigen::Matrix<int, -1, -1,Eigen::RowMajor> MatrixXi_rm;

void Abaqus_read_withAB(std::string &fname, MatrixXd &Nodes, MatrixXi_rm &Element, std::vector<int> &fault_main, std::vector<std::vector<int>> &fault_nodes_org,std::vector<int> &BIE_top, std::vector<int> &BIE_bot, std::vector<int> &AB_left, std::vector<int> &AB_right);


#endif /* Abaqus_read_withAB_hpp */
