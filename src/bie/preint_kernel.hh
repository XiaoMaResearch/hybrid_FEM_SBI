/**
 * @file preint_kernel.hh
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
/* -------------------------------------------------------------------------- */
#ifndef __PREINT_KERNEL_H__
#define __PREINT_KERNEL_H__
/* -------------------------------------------------------------------------- */

#include <vector>
#include <complex>

#include "kernel.hh"
#include "limited_history.hh"

class PreintKernel {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  PreintKernel(const Kernel * kernel);
  virtual ~PreintKernel();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  // compute pre-integration of kernel (only kernel: no other factors)
  // e.g., integral of H11(q*c_s*t) with time_factor q*c_s
  void preintegrate(double time_factor, 
		    double time_step);

  // muliply preintegrated kernel by factor
  void multiplyBy(double factor);

  // compute convolution of kernel with a history
  std::complex<double> convolve(const LimitedHistory * U_r,
				const LimitedHistory * U_i);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  unsigned int getSize() const { return this->preintegrated_kernel.size(); };
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  const Kernel * kernel;
  
  std::vector<double> preintegrated_kernel;

};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */


#endif /* __PREINT_KERNEL_H__ */
