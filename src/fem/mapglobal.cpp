//
//  mapglobal.cpp
//  hybrid_fem_bie
//
//  Created by Max on 2/6/18.
//
//

#include "mapglobal.hpp"
void mapglobal(const ArrayXi &index, VectorXd &u_global , const VectorXd &u_local)
{
    for (int i=0;i<u_local.size();i++)
    {
        //  u_local(i)=u_global(index(i));
        u_global(index(i))=u_global(index(i))+u_local(i);
    }
}
