//
//  bcdof.hpp
//  hybrid_fem_bie
//
//  Created by Max on 2/5/18.
//
//

#ifndef bcdof_hpp
#define bcdof_hpp

#include <stdio.h>
#include <Eigen/Eigen>
using namespace Eigen;
void bcdof(VectorXi node, int dim, VectorXi &index);

#endif /* bcdof_hpp */
