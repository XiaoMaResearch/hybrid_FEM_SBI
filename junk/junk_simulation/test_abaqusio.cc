// Read Abaqus input file
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <Eigen/Eigen>
#include "Abaqus_read.hpp"
#include "mesh_Abaqus_multi_faults.hpp"
#include "share_header.hpp"
using namespace Eigen;
//void Abaqus_read(std::string & fname);

int main (int argc, char *argv[]) 
{
    std::string filename = "abaqus_cpp_5degree_100m_10km_4km.inp";
    std::cout<<filename<<std::endl;
    double x_min = -30e3;
    double x_max = 30e3;
    double y_min = -0.5e3;
    double y_max = 0.5e3;
    int dim = 2.0;
    double dx = 3.125;
    double dy = 3.125;
    int nx_el = (x_max-x_min)/dx;
    int ny = (y_max-y_min)/dy;
    //ArrayXi BIE_top_surf_nodes = ArrayXi::Zero((nx_el+1),1);
    //ArrayXi BIE_bot_surf_nodes = ArrayXi::Zero((nx_el+1),1);
    std::vector<int> BIE_top, BIE_bot;
    MatrixXd Nodes;
    MatrixXi_rm Element;
    std::vector<std::vector<int>> fault_nodes(2);
    mesh_Abaqus_multi_faults(filename, Nodes, Element, BIE_top, BIE_bot,fault_nodes);
    std::cout<<Nodes<<std::endl;
    std::cout<<"Node size=" << Nodes.rows()<<std::endl;
    std::cout<<Element<<std::endl;
    std::cout<<"Element size="<<Element.rows()<<std::endl;
    return 1;
}
