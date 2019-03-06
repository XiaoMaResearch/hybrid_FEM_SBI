//
//  maplocal.cpp
//  hybrid_fem_bie
//
//  Created by Max on 2/6/18.
//
//

#include "maplocal.hpp"
void maplocal(const ArrayXi &index, const VectorXd &u_global , VectorXd & u_local)
{
    
    for (int i=0;i<u_local.size();i++)
    {
        u_local(i)=u_global(index(i));
    }
}
