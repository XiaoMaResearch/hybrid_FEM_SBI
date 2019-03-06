//
//  BIE_initiation.cpp
//  hybrid_fem_bie
//
//  Created by Max on 2/8/18.
//
//

#include "BIE_initiation.hpp"

void BIE_initiation(double E, double nu, double density, double x_max, double x_min, int nx, double dt, InfiniteBoundary &BIE_inf_top, InfiniteBoundary &BIE_inf_bot)
{
//    // Setting up the material property for the BIE code
//    Material BIE_top_mat = Material(E,nu,density);
//    Material BIE_bot_mat = Material(E,nu,density);
//    double length = x_max-x_min;
//    // infinte bc BIE call infinite_boundary.cc
//    PrecomputedKernel h11("kernels/nu_.25_h11.dat");
//    PrecomputedKernel h12("kernels/nu_.25_k12.dat");
//    PrecomputedKernel h22("kernels/nu_.25_h22.dat");
//    InfiniteBoundary BIE_inf_top(length,nx+1,1.0,&BIE_top_mat,&h11,&h12,&h22);
//    InfiniteBoundary BIE_inf_bot(length,nx+1,-1.0,&BIE_bot_mat,&h11,&h12,&h22);
//    // BIE setting time step
//    BIE_inf_top.setTimeStep(dt);
//    BIE_inf_bot.setTimeStep(dt);
//    // BIE initialization
//    BIE_inf_top.init();
//    BIE_inf_bot.init();
}


