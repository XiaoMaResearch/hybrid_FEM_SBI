//
//  maplocal.hpp
//  hybrid_fem_bie
//
//  Created by Max on 2/6/18.
//
//

#ifndef maplocal_hpp
#define maplocal_hpp

#include <stdio.h>
#include <Eigen/Eigen>
using namespace Eigen;
void maplocal(const ArrayXi &index, const VectorXd &u_global , VectorXd &u_local);


#endif /* maplocal_hpp */
