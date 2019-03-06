//
//  mapglobal.hpp
//  hybrid_fem_bie
//
//  Created by Max on 2/6/18.
//
//

#ifndef mapglobal_hpp
#define mapglobal_hpp

#include <stdio.h>
#include <Eigen/Eigen>
using namespace Eigen;
void mapglobal(const ArrayXi &index, VectorXd &u_global ,const VectorXd &u_local);

#endif /* mapglobal_hpp */
