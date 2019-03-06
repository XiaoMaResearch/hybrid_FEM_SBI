//
//  Getplastic.cpp
//  hybrid_fem_bie
//
//  Created by Max on 3/20/18.
//
//

#include "Getplastic.hpp"
#include <iostream>

void Getplastic(double E, double nu , VectorXd &strain, VectorXd &stress , VectorXd &stress_n, VectorXd &strain_n,double sxx_initial,
                double syy_initial, double sxy_initial, double cohes, double blkfric, MatrixXd &de_p)
{
//    double cohes = 0.0;
//    double blkfric = 0.6;
    double angfric = std::atan(blkfric);
    
    VectorXd strain_inc = strain-strain_n;
    double exx = strain_inc(0);
    double eyy = strain_inc(1);
    double exy = strain_inc(2);
    
    double sxx = stress_n(0);
    double syy = stress_n(1);
    double sxy = stress_n(2);
    
//    double density = 2670;
//    double vp = 6000.0 ;
//    double vs = 3464.0;
//    double mu = density*(vs*vs);
//    double lambda = density*(vp*vp-2.0*(vs*vs));
    double lambda = E*nu/((1+nu)*(1-2*nu));
    double mu = E/(2.0*(1+nu));
    
    double etrace = exx+eyy;
    sxx = sxx + lambda*etrace+2.0*mu*exx;
    syy = syy + lambda*etrace+2.0*mu*eyy;
    sxy = sxy + 2.0*mu*exy/2.0;

    sxx = sxx + sxx_initial;
    syy = syy + syy_initial;
    sxy = sxy + sxy_initial;
//    
   // double szz = 0.5*(sxx+syy);
    double szz = nu*(sxx+syy);

    double sm = (sxx+syy+szz)/3.0;
    // Find stress divatoric components
    double sdxx = sxx - sm;
    double sdyy = syy - sm;
    double sdzz = szz - sm;
    double sdxy = sxy;
    // Secondar invariant of stress
    double secinv = 0.5*(sdxx*sdxx+sdyy*sdyy+sdzz*sdzz)+sdxy*sdxy;
    // scalar measure of shear stress
    double tau = std::sqrt(secinv);
    double taulim = cohes*std::cos(angfric)-(sm)*std::sin(angfric);
   // std::cout<<"tau="<<tau<<"taulim="<<taulim<<std::endl;
   // std::cout<<"sm="<<sm<<"angfric="<<angfric<<"cohes="<<cohes<<std::endl;
    taulim = std::max(0.0, taulim);
    
//   std::cout<<"tau="<<tau<<"taulim="<<taulim<<std::endl;
  //  VectorXd de_p = VectorXd::Zero(3, 1);
  //  std::cout<<"tau="<<tau<<std::endl;
  //  std::cout<<"tau_lim="<<taulim<<std::endl;
    de_p = MatrixXd::Zero(3, 1);
//
    if (tau>taulim)
    {
        double yldfac = taulim/tau;
        sxx = sdxx*yldfac + sm;
        syy = sdyy*yldfac + sm;
        szz = sdzz*yldfac + sm;
        sxy = sdxy*yldfac;
        
        // Adding the calculation of the plastic strain rate
        double dsxx = sxx-stress_n(0);
        double dsyy = syy-stress_n(1);
        double dsxy = sxy-stress_n(2);
        double dszz = szz-nu*(stress_n(0)+stress_n(1));
        
        double de_xx  = (dsxx - nu*(dsyy+dszz))/E;
        double de_yy  = (dsyy - nu*(dsxx+dszz))/E;
        double de_xy  = (dszz - nu*(dsxx+dsyy))/E;
        
        de_p(0) = exx-de_xx;
        de_p(1) = eyy-de_yy;
        de_p(2) = exy-de_xy;
        
        //
        
        
        sxx = sxx-sxx_initial;
        syy = syy-syy_initial;
        sxy = sxy-sxy_initial;
        
        
    }
    else
    {
        sxx = sxx-sxx_initial;
        syy = syy-syy_initial;
        sxy = sxy-sxy_initial;
    }
    
    stress(0) = sxx;
    stress(1) = syy;
    stress(2) = sxy;
}
