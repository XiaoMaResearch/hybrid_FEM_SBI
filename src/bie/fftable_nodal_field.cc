/**
 * @file fftable_nodal_field.cc
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
#include "fftable_nodal_field.hh"
#include <cstring>

/* -------------------------------------------------------------------------- */
FFTableNodalField::FFTableNodalField(unsigned int nb_nodes) :
  NodalField(nb_nodes) {
  
  this->nb_fft = this->nb_nodes / 2 + 1;

  this->freq_dom_field = new fftw_complex[this->nb_fft];
  memset(this->freq_dom_field, 0., this->nb_fft*sizeof(fftw_complex));

  this->forward_plan = fftw_plan_dft_r2c_1d(this->nb_nodes,
					    this->field,
					    this->freq_dom_field,
					    FFTW_MEASURE);

  this->backward_plan = fftw_plan_dft_c2r_1d(this->nb_nodes,
					     this->freq_dom_field,
					     this->field,
					     FFTW_MEASURE);
  
}

/* -------------------------------------------------------------------------- */
FFTableNodalField::~FFTableNodalField() {
  
  fftw_destroy_plan(this->forward_plan);
  fftw_destroy_plan(this->backward_plan);
  //fftw_cleanup(); // not needed because we don't know if all plans are destroyed

  delete[] this->freq_dom_field;

}

/* -------------------------------------------------------------------------- */
void FFTableNodalField::forwardFFT() {

  fftw_execute(this->forward_plan);

}

/* -------------------------------------------------------------------------- */
void FFTableNodalField::backwardFFT() {

  fftw_execute(this->backward_plan);

  for (unsigned int i=0; i<this->nb_nodes; ++i) {
    this->field[i] /= double(this->nb_nodes);
  }

}
