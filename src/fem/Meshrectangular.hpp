//
//  Meshrectangular.hpp
//  hybrid_fem_bie
//
//  Created by Max on 2/5/18.
//
//

#ifndef Meshrectangular_hpp
#define Meshrectangular_hpp

#include <stdio.h>
#include <Eigen/Eigen>
#include "share_header.hpp"
using namespace Eigen;

void Meshrectangular(double x_min, double x_max, double y_min, double y_max, double dx ,double dy ,int nx, int ny , MatrixXd &Node, MatrixXi_rm &Element);


#endif /* Meshrectangular_hpp */
