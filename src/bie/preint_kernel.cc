/**
 * @file preint_kernel.cc
 *
 * @author David S. Kammer <kammer@cornell.edu>
 *
 * @date creation: Mon Oct 09 2017
 * @date last modification: Mon Oct 09 2017
 *
 * @brief 
 *
 * @section LICENSE
 *
 * MIT License
 * Copyright (c) 2017 David S. Kammer 
 *
 */
#include "preint_kernel.hh"
//#include <omp.h>
/* -------------------------------------------------------------------------- */
PreintKernel::PreintKernel(const Kernel * kernel) :
  kernel(kernel) {
  
}

/* -------------------------------------------------------------------------- */
PreintKernel::~PreintKernel() {

}

/* -------------------------------------------------------------------------- */
void PreintKernel::preintegrate(double time_factor, double time_step) {
  
  // time for truncation of this mode
  // Tcut = tcut * q * cs  <-> tcut = Tcut / q / cs
  double trunc = this->kernel->getTruncation() / time_factor;

  int nb_integration_int = (int)(trunc / time_step);
  this->preintegrated_kernel.resize(nb_integration_int);

  // compute trapezoidal integral over time step
  double k_i = this->kernel->at(0.);
  for (unsigned int i=0; i<nb_integration_int; ++i) {
    double k_ii = this->kernel->at( (i+1)*time_step * time_factor); // k after i-th step
    this->preintegrated_kernel[i] = 0.5 * (k_i + k_ii) * time_step * time_factor;
    k_i = k_ii;
  }

}

/* -------------------------------------------------------------------------- */
void PreintKernel::multiplyBy(double factor) {

  unsigned int nb_integration_int = this->preintegrated_kernel.size();
  for (unsigned int i=0; i<nb_integration_int; ++i) {
    this->preintegrated_kernel[i] *= factor;
  }
  
}

/* -------------------------------------------------------------------------- */
std::complex<double> PreintKernel::convolve(const LimitedHistory * U_r,
					    const LimitedHistory * U_i) {

  unsigned int nb_U = U_r->getNbHistoryPoints();

  double real = 0.;
  double imag = 0.;
//  #pragma omp parallel for
  for (unsigned int i=0; i<nb_U; ++i) {

    double K_cum = this->preintegrated_kernel[i]; 

    real += K_cum * U_r->at(i);
    imag += K_cum * U_i->at(i);
  }

  return {real, imag};
}
