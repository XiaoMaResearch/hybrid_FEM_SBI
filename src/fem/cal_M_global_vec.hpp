//
//  cal_M_global_vec.hpp
//  hybrid_fem_bie
//
//  Created by Max on 2/10/18.
//
//

#ifndef cal_M_global_vec_hpp
#define cal_M_global_vec_hpp

#include <stdio.h>
#include <Eigen/Eigen>
#include "cal_M.hpp"
#include "mapglobal.hpp"

using namespace Eigen;
void cal_M_global_vec(MatrixXd &Node, MatrixXd &Element, double density, MatrixXd &index_store, int Ndofn, VectorXd &M_global_vec);


#endif /* cal_M_global_vec_hpp */
