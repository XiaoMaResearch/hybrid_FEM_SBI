//
//  bcdof_ptr.hpp
//  hybrid_fem_bie
//
//  Created by Max on 2/14/18.
//
//

#ifndef bcdof_ptr_hpp
#define bcdof_ptr_hpp

#include <stdio.h>
#include <Eigen/Eigen>
using namespace Eigen;

void bcdof_ptr(std::vector<int> &node, int dim, int *ptr);

#endif /* bcdof_ptr_hpp */
