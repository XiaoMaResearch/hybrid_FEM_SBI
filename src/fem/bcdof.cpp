//
//  bcdof.cpp
//  hybrid_fem_bie
//
//  Created by Max on 2/5/18.
//
//

#include "bcdof.hpp"
void bcdof(VectorXi node, int dim, VectorXi &index)
{
    for (int i =0; i<node.size(); i++)
    {
        index(2*i) = dim*node(i);
        index(2*i+1) = dim*node(i)+1;
    }
    
}

