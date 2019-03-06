//
//  Getplastic.hpp
//  hybrid_fem_bie
//
//  Created by Max on 3/20/18.
//
//

#ifndef Getplastic_hpp
#define Getplastic_hpp

#include <stdio.h>
#include <Eigen/Eigen>

using namespace Eigen;

void Getplastic(double E, double nu , VectorXd &strain, VectorXd &stress , VectorXd &stress_n, VectorXd &strain_n,double sxx_initial,
                double syy_initial, double sxy_initial, double cohes, double blkfric, MatrixXd &de_p);
#endif /* Getplastic_hpp */
