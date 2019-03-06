//
//  mesh_Abaqus_multi_faults_withAB.cpp
// With Absorbing the bc 
//  hybrid_fem_bie
//
//  Created by Max on 2/13/18.
//
//

#include "mesh_Abaqus_multi_faults_withAB.hpp"
#include "Eiglab.h"
//#include <omp.h>
#include <numeric>
#include <vector>
void mesh_Abaqus_multi_faults_withAB(std::string &fname,MatrixXd &Node ,MatrixXi_rm &Element,
                                     std::vector<int> &BIE_top, std::vector<int> &BIE_bot,
                                     std::vector<int> &AB_left, std::vector<int> &AB_right,
                              std::vector<std::vector<int>> &fault_nodes, double num_faults)
{
    std::vector<std::vector<int>> fault_nodes_org(num_faults-1);
    std::vector<int> fault_main;
    //Abaqus_read(fname, Node, Element,fault_main, fault_nodes_org,BIE_top,BIE_bot);
    Abaqus_read_withAB(fname, Node, Element, fault_main, fault_nodes_org, BIE_top, BIE_bot, AB_left,AB_right);
    std::cout<< "Done"<<std::endl;
    //  std::vector<int> fault_nodes_temp
    // fault_nodes[2*j].insert
    auto it = fault_nodes[0].begin();
    fault_nodes[0].insert(it,fault_main.begin(),fault_main.end());
    update_mesh_topology(Node,Element,fault_nodes[0],fault_nodes[1]);
    int k=0;
    for (int i=1;i<fault_nodes.size()/2;i++)
    {
        it = fault_nodes[2*i].begin();
        fault_nodes[2*i].insert(it,fault_nodes_org[k].begin(),fault_nodes_org[k].end());
       // update_mesh_topology(Node,Element,fault_nodes[2],fault_nodes[3]);
        update_mesh_topology(Node,Element,fault_nodes[2*i],fault_nodes[2*i+1]);
        k=k+1;

    }
    std::cout<<"Done_finding_nodes"<<std::endl;
    // Sorting the BIE nodes to make sure it is align as left to right in the array
    //TOP
    VectorXd x_dummy = VectorXd::Zero(BIE_top.size(), 1) ;
    for (int i=0; i<BIE_top.size(); i++)
    {
        x_dummy(i) = Node(BIE_top[i],0);
    }
   // std::cout<<x_dummy<<std::endl;
    ArrayXi index = EigLab::sort(x_dummy, 1);
    std::vector<int> BIE_top_temp = BIE_top;
    for (int i =0; i<BIE_top.size();i++)
    {
        BIE_top[i]= BIE_top_temp[index(i)];
    }
    // BOT
    x_dummy = VectorXd::Zero(BIE_bot.size(), 1) ;
    for (int i=0; i<BIE_bot.size(); i++)
    {
        x_dummy(i) = Node(BIE_bot[i],0);
    }
    // std::cout<<x_dummy<<std::endl;
    index = EigLab::sort(x_dummy, 1);
    std::vector<int> BIE_bot_temp = BIE_bot;
    for (int i =0; i<BIE_bot.size();i++)
    {
        BIE_bot[i]= BIE_bot_temp[index(i)];
    }
    
    // AB_left
    x_dummy = VectorXd::Zero(AB_left.size(), 1) ;
    for (int i=0; i<AB_left.size(); i++)
    {
        x_dummy(i) = Node(AB_left[i],1);
    }
    index = EigLab::sort(x_dummy, 1);
    std::vector<int> AB_left_temp = AB_left;
    for (int i =0; i<AB_left_temp.size();i++)
    {
        AB_left[i]= AB_left_temp[index(i)];
    }
    
    // AB_right
    x_dummy = VectorXd::Zero(AB_right.size(), 1) ;
    for (int i=0; i<AB_right.size(); i++)
    {
        x_dummy(i) = Node(AB_right[i],1);
    }
    index = EigLab::sort(x_dummy, 1);
    std::vector<int> AB_right_temp = AB_right;
    for (int i =0; i<AB_right_temp.size();i++)
    {
        AB_right[i]= AB_right_temp[index(i)];
    }
    
    
    
    
 // sort the fault nodes from left to right in x
    for(int i = 0; i<2*num_faults; i++)
    {
        VectorXd x_dummy_temp = VectorXd::Zero(fault_nodes[i].size(),1);

        for (int k=0; k<fault_nodes[i].size(); k++)
        {
            x_dummy_temp(k) = Node(fault_nodes[i][k],0);
        }
        index = EigLab::sort(x_dummy_temp, 1);
        std::vector<int> temp = fault_nodes[i];
        for (int j = 0;j<fault_nodes[i].size();j++)
        {
            fault_nodes[i][j]= temp[index(j)];
        }
    }
    
   }
