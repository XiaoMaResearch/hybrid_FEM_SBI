/**
 * @file precomputed_kernel.cc
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
#include "precomputed_kernel.hh"

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

/* -------------------------------------------------------------------------- */
PrecomputedKernel::PrecomputedKernel(const std::string & fname) {
  this->readKernelFromFile(fname);
}

/* -------------------------------------------------------------------------- */
PrecomputedKernel::~PrecomputedKernel() {
}

/* -------------------------------------------------------------------------- */
double PrecomputedKernel::at(double T) const {
  
  int index = (int)(T/this->delta_t);

  if (index >= this->values.size()) {
    std::cout << "Try to get value of PrecomputedKernel for T=" << T 
	      << ", which is beyond the PrecomputedKernel's range" << std::endl;
    throw T;
  }

  double dk = this->values[index+1] - this->values[index];
  double dT = T - index*this->delta_t;

  return this->values[index] +  dk * dT/this->delta_t;
}


/* -------------------------------------------------------------------------- */
/* 
   This function reads Kernels with the following format
   It contains:
   - not needed information (1 int)
   - trucation (1 double)
   - number of entries in kernel (1 int)
   - time step of kernel entries (1 double)
   - not needed information (cd/cs and nu: 2 double)
   - kernel entries
 */
/* -------------------------------------------------------------------------- */
void PrecomputedKernel::readKernelFromFile(const std::string & fname) {
  
  std::ifstream file(fname, std::ios::binary);

  if (!file.is_open()){
    std::cout << "Cannot find file named: " << fname << std::endl;
    return;
  }

  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << "Read Kernel from fortran-generated binary file: " << fname << std::endl;
  
  double dtmp;
  int itmp;

  this->readBinary(file, itmp);
  std::cout << "  tmp   = " << itmp << std::endl;
  this->readBinary(file, this->trunc);
  std::cout << "  trunc = " << this->trunc << std::endl;
  int size;
  this->readBinary(file, size);
  std::cout << "  size  = " << size << std::endl;
  this->readBinary(file, this->delta_t);
  std::cout << "  delta = " << this->delta_t << std::endl;

  this->readBinary(file, dtmp);
  std::cout << "  tmp   = " << dtmp << std::endl;
  this->readBinary(file, dtmp);
  std::cout << "  tmp   = " << dtmp << std::endl;

  this->values.resize(size);
  this->readBinary(file, this->values);
  std::cout << "-----------------------------------------------" << std::endl << std::endl;

  /*
  std::ofstream kout;
  kout.open("kout.txt");

  int ksize = values.size();
  for (int i=0; i<ksize; ++i) {
    kout << values[i] << "\n";
  }
  kout.close();
  */
  
}

/* -------------------------------------------------------------------------- */
template<class T>
void PrecomputedKernel::readBinary(std::ifstream & file, T & val) {
  char * val_char = reinterpret_cast<char *>(&val);
  file.read((char*)val_char, sizeof(T));
}
