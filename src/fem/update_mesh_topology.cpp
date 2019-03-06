//
//  update_mesh_topology.cpp
//  hybrid_fem_bie
//
//  Created by Max on 2/6/18.
//
//

#include "update_mesh_topology.hpp"
#include <iostream>
//#include <omp.h>
void update_mesh_topology(MatrixXd &Nodes, MatrixXi_rm &Element,std::vector<int> &fault_surf_nodes, std::vector<int> &fault_surf_nodes_new)
{
    // number of nodes
    int n_nodes = Nodes.rows();
    // number of elements
    int n_el = Element.rows();
    //
    std::vector<Eigen::Array<int, -1, 2> > NodesElem(n_nodes);
    // 
   // double start = omp_get_wtime();
    for(int i=0; i<n_el; i++)
    {
        for (int j=0; j<4; j++)
        {
            NodesElem[Element(i,j)].conservativeResize(NodesElem[Element(i,j)].rows()+1,2);
            NodesElem[Element(i,j)].row(NodesElem[Element(i,j)].rows()-1) << i, j;
        }
    }
   // double end = omp_get_wtime();
   // std::cout<<"Loop_1="<<end-start<<std::endl;
    // Set the new fault_surf_nodes_new numbers (start from n_nodes)
    //start = omp_get_wtime();
    for (int i=0; i<fault_surf_nodes.size(); i++)
    {
       // fault_surf_nodes_new[i] = n_nodes+i;
        fault_surf_nodes_new.push_back(n_nodes+i);
    }
   // end = omp_get_wtime();
   // std::cout<<"loop_2="<<end-start<<std::endl;
    // Update the Nodes by append the duplcated nodes coordiations at the end
    //start = omp_get_wtime();
    Nodes.conservativeResize(n_nodes+fault_surf_nodes_new.size(),2);
    for (int i=0; i<fault_surf_nodes_new.size();i++)
    {
        Nodes.row(n_nodes+i) = Nodes.row(fault_surf_nodes[i]);
    }
    //end = omp_get_wtime();
    //std::cout<<"loop_3="<<end-start<<std::endl;
    // Element_logical , nodes that are on the fault is 1, otherwise it is -1
    MatrixXi Elements_logcial = -MatrixXi::Ones(Element.rows(),4);
    // VectorXd Fault_Elements = VectorXd::Zero(4*fault_surf_nodes.rows(),1);
    std::vector<int> Fault_Elements;
    Fault_Elements.reserve(4*fault_surf_nodes.size());
    //start = omp_get_wtime();
    for (int i = 0; i<fault_surf_nodes.size();i++)
    {
        for (int j = 0; j<NodesElem[fault_surf_nodes[i]].rows(); j++)
        {
            Elements_logcial(NodesElem[fault_surf_nodes[i]](j,0),NodesElem[fault_surf_nodes[i]](j,1)) = 1.0;
            Fault_Elements.push_back(NodesElem[fault_surf_nodes[i]](j,0));

            
        }
    }
    //end = omp_get_wtime();
    //std::cout<<"loop_4="<<end-start<<std::endl;
    // Find the unique Elements
    std::sort(Fault_Elements.begin(),Fault_Elements.end());
    auto last = std::unique(Fault_Elements.begin(), Fault_Elements.end());
    Fault_Elements.erase(last,Fault_Elements.end());
    // Calculating the normal direction of the elemnts on the faults
    std::vector<int> Element_fault_new(Fault_Elements.size()/2);
    int k = 0;
    //start = omp_get_wtime();
    for (auto el:Fault_Elements)
    {
        // Loop over fault boundary elements and get vector tangent to fault
        VectorXi temp = Elements_logcial.row(el);
        VectorXd vec(2);
        if (temp.sum() == -2)
        {
            continue;
        }
        else if (temp(0)==1 && temp(temp.size()-1)==1)
            //Special case where fault nodes are first and last of element
        {
            vec = Nodes.row(Element(el,0))-Nodes.row(Element(el,Element.cols()-1));
        }
        else
            // Normal case where fault nodes are consecutive nodes in the element
        {
            for (int j=0; j<temp.size(); j++)
            {
                if(temp(j)==1)
                {
                    vec = Nodes.row(Element(el,j+1))-Nodes.row(Element(el,j));
                    break;
                }
            }
        }
        // Check direction, check if the tangential direction projection in [1.0 , 0.0] is postive or negative
        // If negative , then that means the elements are in the negative side of the fault we define positive as the downward direction.
        VectorXd n_hat(2);
        n_hat << 1.0, 0.0;
        //std::cout<< temp_dir<<std::endl;
        if (vec.dot(n_hat)<0)
        {
            Element_fault_new[k] = el;
            k=k+1;
        }
        // Normal direction points downward --nothing to do
    }
    //end = omp_get_wtime();
    //std::cout<<"loop_5="<<end-start<<std::endl;
    VectorXi Element_logical_array = VectorXi::Zero(n_el, 1);
    //start = omp_get_wtime();
    for (int i=0; i<Element_fault_new.size(); i++)
    {
        Element_logical_array(Element_fault_new[i]) = 1;
    }
    //end = omp_get_wtime();
    //std::cout<<"loop_6="<<end-start<<std::endl;
    //start = omp_get_wtime();
    for (int i =0; i<fault_surf_nodes_new.size(); i++)
    {
        for (int j=0; j<NodesElem[fault_surf_nodes[i]].rows(); j++)
        {
            int el = NodesElem[fault_surf_nodes[i]](j,0);
            int nd = NodesElem[fault_surf_nodes[i]](j,1);
            if (Element_logical_array(el)==1)
            {
                Element(el,nd) = fault_surf_nodes_new[i];
            }
        }
    } 
    //end = omp_get_wtime();
    //std::cout<<"loop_7="<<end-start<<std::endl;
}
