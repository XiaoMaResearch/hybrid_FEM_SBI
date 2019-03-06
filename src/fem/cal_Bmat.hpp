//
//  cal_Bmat.hpp
//  hybrid_fem_bie
//
//  Created by Max on 3/21/18.
//
//

#ifndef cal_Bmat_hpp
#define cal_Bmat_hpp

#include <stdio.h>
#include <Eigen/Eigen>

using namespace Eigen;

void cal_Bmat(MatrixXd coord ,double E, double nu, std::vector<Eigen::MatrixXd> &B_mat, double &detJ);


#endif /* cal_Bmat_hpp */
