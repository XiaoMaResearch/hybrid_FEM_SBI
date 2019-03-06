/**
 * @file fftable_nodal_field.hh
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
#ifndef __FFTABLE_NODAL_FIELD_H__
#define __FFTABLE_NODAL_FIELD_H__
/* -------------------------------------------------------------------------- */
#include <fftw3.h>

#include "nodal_field.hh"
/* -------------------------------------------------------------------------- */

class FFTableNodalField : public NodalField {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  FFTableNodalField(unsigned int nb_nodes);
  virtual ~FFTableNodalField();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void forwardFFT();
  void backwardFFT();

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  // get number of points in frequency domain
  unsigned int getNbFFT() { return this->nb_fft; };

  // get one value of frequency domain
  inline fftw_complex & fd(unsigned int f);

  // get access directly to frequency domain
  // WARNING: convert it to double (assuming that fftw_complex is double[2]
  inline fftw_complex * fd_storage() { return this->freq_dom_field; };
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  // number of points in frequency domain
  unsigned int nb_fft;

  // values in frequency domain in complex form
  fftw_complex * freq_dom_field;

  // the forward and backward plan
  fftw_plan forward_plan;
  fftw_plan backward_plan;

};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
inline fftw_complex & FFTableNodalField::fd(unsigned int f) {
  return this->freq_dom_field[f];
};

#endif /* __FFTABLE_NODAL_FIELD_H__ */
