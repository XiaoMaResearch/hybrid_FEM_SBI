/**
 * @file precomputed_kernel.hh
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
#ifndef __PRECOMPUTED_KERNEL_H__
#define __PRECOMPUTED_KERNEL_H__
/* -------------------------------------------------------------------------- */
#include <fstream>
#include <iostream>
#include <sstream>

#include <vector>

#include "kernel.hh"

class PrecomputedKernel : public Kernel {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  PrecomputedKernel(const std::string & fname);
  virtual ~PrecomputedKernel();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  // read kernel from file generated before with fortran script
  void readKernelFromFile(const std::string & fname);

  // Templated binary reader
  template<class T>
  void readBinary(std::ifstream & file, T & val);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  // get value of Kernel at given T (by interpolation)
  virtual double at(double T) const;

  // get T for truncation
  virtual double getTruncation() const { return this->trunc; };

  // get number of nodes
  unsigned int getSize() const { return this->values.size(); };

  // get direct access to values (only used to debug)
  std::vector<double> & getDirectAccess() { return this->values; };
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  // nodal field
  std::vector<double> values;

  // truncation of kernel
  double trunc;

  // delta time of kernel entries
  double delta_t;

};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
template<>
inline void PrecomputedKernel::readBinary<std::vector<double> >(std::ifstream & file, 
								std::vector<double> & val) {
  
  int nval = val.size();
  double * tmp_values = new double[2*nval];

  file.read((char*) tmp_values, 2*nval*sizeof(double));
  for (int i = 0; i < nval; ++i) {
    val[i] = tmp_values[2*i+1];
  }
  
  delete[] tmp_values;
}


#endif /* __PRECOMPUTED_KERNEL_H__ */
