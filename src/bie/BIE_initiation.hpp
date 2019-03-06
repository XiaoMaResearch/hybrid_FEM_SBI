//
//  BIE_initiation.hpp
//  hybrid_fem_bie
//
//  Created by Max on 2/8/18.
//
//

#ifndef BIE_initiation_hpp
#define BIE_initiation_hpp

#include <stdio.h>
#include "material.hh"
#include "precomputed_kernel.hh"
#include "bimat_interface.hh"
#include "infinite_boundary.hh"

void BIE_initiation(double E, double nu, double density, double x_max, double x_min, int nx, double dt, InfiniteBoundary &BIE_inf_top, InfiniteBoundary &BIE_inf_bot);


#endif /* BIE_initiation_hpp */
