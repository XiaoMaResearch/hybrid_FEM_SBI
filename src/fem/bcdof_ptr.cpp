//
//  bcdof_ptr.cpp
//  hybrid_fem_bie
//
//  Created by Max on 2/14/18.
//
//

#include "bcdof_ptr.hpp"

void bcdof_ptr(std::vector<int> &node, int dim, int *ptr)
{
    for (int i =0; i<node.size(); i++)
    {
        ptr[2*i] = dim*node[i];
        ptr[2*i+1] = dim*node[i]+1;
    }
    
}

