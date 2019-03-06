/**
 * @file half_space.cc
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
#include "half_space.hh"

#include <iostream>
#include <cmath>
//#include <omp.h>
extern double time_bie;
/* -------------------------------------------------------------------------- */
HalfSpace::HalfSpace(double length, 
		     unsigned int nb_nodes,
		     int side_factor) :
  length(length),
  nb_nodes(nb_nodes),
  time_step(0.),
  side_factor(side_factor)
{
  std::cout << "HalfSpace constructed" << std::endl;

  this->q1 = 2*M_PI / this->length;

  this->disp_1 = new FFTableNodalField(nb_nodes);
  this->disp_2 = new FFTableNodalField(nb_nodes);
  this->nb_fft = this->disp_1->getNbFFT();

  this->U_1_r.resize(nb_fft);
  this->U_1_i.resize(nb_fft);
  this->U_2_r.resize(nb_fft);
  this->U_2_i.resize(nb_fft);

  this->H11_pi.resize(nb_fft);
  this->H12_pi.resize(nb_fft);
  this->H22_pi.resize(nb_fft);


  this->velo_1 = new NodalField(nb_nodes);
  this->velo_2 = new NodalField(nb_nodes);
  
  this->internal_1 = new FFTableNodalField(nb_nodes);
  this->internal_2 = new FFTableNodalField(nb_nodes);

  this->residual_1 = new NodalField(nb_nodes);
  this->residual_2 = new NodalField(nb_nodes);
}

/* -------------------------------------------------------------------------- */
HalfSpace::~HalfSpace() {

  delete this->disp_1;
  delete this->disp_2;

  for (unsigned int j=0; j<this->nb_fft; ++j) {
    delete this->H11_pi[j];
    delete this->H12_pi[j];
    delete this->H22_pi[j];
  }

  delete this->velo_1;
  delete this->velo_2;

  delete this->internal_1;
  delete this->internal_2;
  
  delete this->residual_1;
  delete this->residual_2;

  for (unsigned int j=0; j<this->nb_fft; ++j) {
    delete this->U_1_r[j];
    delete this->U_1_i[j];

    delete this->U_2_r[j];
    delete this->U_2_i[j];
  }
}

/* -------------------------------------------------------------------------- */
void HalfSpace::initConvolutions() {

  double q1_cs = this->q1 * this->material->getCs();

  // history for q1 is longest q = j*q1
  for (unsigned int j=1; j<this->nb_fft; ++j) {

    this->H11_pi[j] = new PreintKernel(this->H11);
    this->H12_pi[j] = new PreintKernel(this->H12);
    this->H22_pi[j] = new PreintKernel(this->H22);

    double qj_cs = j * q1_cs;
    this->H11_pi[j]->preintegrate(qj_cs, this->time_step);
    this->H12_pi[j]->preintegrate(qj_cs, this->time_step);
    this->H22_pi[j]->preintegrate(qj_cs, this->time_step);


    unsigned int nb_hist_1 = std::max(this->H11_pi[j]->getSize(),
				      this->H12_pi[j]->getSize());
    unsigned int nb_hist_2 = std::max(this->H22_pi[j]->getSize(),
				      this->H12_pi[j]->getSize());

    this->U_1_r[j] = new LimitedHistory(nb_hist_1);
    this->U_1_i[j] = new LimitedHistory(nb_hist_1);

    this->U_2_r[j] = new LimitedHistory(nb_hist_2);
    this->U_2_i[j] = new LimitedHistory(nb_hist_2);

  }
  
}

/* -------------------------------------------------------------------------- */
void HalfSpace::computeDisplacement() {
  this->computeDisplacement(this->disp_1, this->velo_1);
  this->computeDisplacement(this->disp_2, this->velo_2);
}

/* -------------------------------------------------------------------------- */
// u_i+1 = u_i + dt * v_i
void HalfSpace::computeDisplacement(FFTableNodalField * disp,
				    NodalField * velo) {

  double * disp_p = disp->storage();
  double * velo_p = velo->storage();

  for (unsigned int n=0; n<this->nb_nodes; ++n) {
    disp_p[n] += velo_p[n] * this->time_step;
  }
}
/* -------------------------------------------------------------------------- */
void HalfSpace::forwardFFT() {
  
  this->disp_1->forwardFFT();
  this->disp_2->forwardFFT();
}

/* -------------------------------------------------------------------------- */
void HalfSpace::backwardFFT() {
  
  this->internal_1->backwardFFT();
  this->internal_2->backwardFFT();

}

/* -------------------------------------------------------------------------- */
void HalfSpace::computeStressFourierCoeff() {
  
  double mu  = this->material->getShearModulus();
  double eta = this->material->getCp() / this->material->getCs();

  // access to fourier coefficients of stresses
  fftw_complex * internal_1_fd = this->internal_1->fd_storage();
  fftw_complex * internal_2_fd = this->internal_2->fd_storage();

  // set for mode zero fourier coefficients to zero
  internal_1_fd[0][0] = 0.; // real part
  internal_1_fd[0][1] = 0.; // imag part
  internal_2_fd[0][0] = 0.; // real part
  internal_2_fd[0][1] = 0.; // imag part

  // imaginary number i
  std::complex<double> imag = {0., 1.};

  //double start = omp_get_wtime();
  #pragma omp parallel for
  for (unsigned int j=1; j<this->nb_fft; ++j) {

    std::complex<double> U1 = {this->disp_1->fd(j)[0], this->disp_1->fd(j)[1]};
    std::complex<double> U2 = {this->disp_2->fd(j)[0], this->disp_2->fd(j)[1]};

    // store current displacement in history
    this->U_1_r[j]->addCurrentValue(std::real(U1));
    this->U_1_i[j]->addCurrentValue(std::imag(U1));
    this->U_2_r[j]->addCurrentValue(std::real(U2));
    this->U_2_i[j]->addCurrentValue(std::imag(U2));

    // mu * q * int(H11,U1)
    std::complex<double> F1 = - this->side_factor * mu * j * this->q1 
                            * this->H11_pi[j]->convolve(this->U_1_r[j],
							this->U_1_i[j]);

    // i *   * mu * q * U2
    F1 += mu * j*this->q1 * (2 - eta) * (imag * U2);

    // + i * mu * q * int(H12, U2)
    F1 += mu * j*this->q1 * imag * this->H12_pi[j]->convolve(this->U_2_r[j],
							     this->U_2_i[j]);

    // - mu * q * int(H22, U2)
    std::complex<double> F2 = - this->side_factor * mu * j*this->q1
                            * this->H22_pi[j]->convolve(this->U_2_r[j],
							this->U_2_i[j]);

    // - i * (2 - etq) * mu * q * U1
    F2 -= mu * j*this->q1 * (2 - eta) * (imag * U1);

    // - i * mu * q * int(H12, U1)
    F2 -= mu * j*this->q1 * imag * this->H12_pi[j]->convolve(this->U_1_r[j],
							     this->U_1_i[j]);

    // set values to internal force
    internal_1_fd[j][0] = std::real(F1); // real part
    internal_1_fd[j][1] = std::imag(F1); // imag part
    internal_2_fd[j][0] = std::real(F2); // real part
    internal_2_fd[j][1] = std::imag(F2); // imag part
  }
 // double end=omp_get_wtime();
  //time_bie += end-start;

}

/*-------------------------------*/
void HalfSpace::computeStressFourierCoeff_correction() {
    
    double mu  = this->material->getShearModulus();
    double eta = this->material->getCp() / this->material->getCs();
    
    // access to fourier coefficients of stresses
    fftw_complex * internal_1_fd = this->internal_1->fd_storage();
    fftw_complex * internal_2_fd = this->internal_2->fd_storage();
    
    // set for mode zero fourier coefficients to zero
    internal_1_fd[0][0] = 0.; // real part
    internal_1_fd[0][1] = 0.; // imag part
    internal_2_fd[0][0] = 0.; // real part
    internal_2_fd[0][1] = 0.; // imag part
    
    // imaginary number i
    std::complex<double> imag = {0., 1.};
    
    for (unsigned int j=1; j<this->nb_fft; ++j) {
        
        std::complex<double> U1 = {this->disp_1->fd(j)[0], this->disp_1->fd(j)[1]};
        std::complex<double> U2 = {this->disp_2->fd(j)[0], this->disp_2->fd(j)[1]};
        
        // store current displacement in history
//        this->U_1_r[j]->addCurrentValue_one(std::real(U1));
//        this->U_1_i[j]->addCurrentValue_one(std::imag(U1));
//        this->U_2_r[j]->addCurrentValue_one(std::real(U2));
//        this->U_2_i[j]->addCurrentValue_one(std::imag(U2));
        
        // mu * q * int(H11,U1)
        std::complex<double> F1 = - this->side_factor * mu * j * this->q1
        * this->H11_pi[j]->convolve(this->U_1_r[j],
                                    this->U_1_i[j]);
        
        // i *   * mu * q * U2
        F1 += mu * j*this->q1 * (2 - eta) * (imag * U2);
        
        // + i * mu * q * int(H12, U2)
        F1 += mu * j*this->q1 * imag * this->H12_pi[j]->convolve(this->U_2_r[j],
                                                                 this->U_2_i[j]);
        
        // - mu * q * int(H22, U2)
        std::complex<double> F2 = - this->side_factor * mu * j*this->q1
        * this->H22_pi[j]->convolve(this->U_2_r[j],
                                    this->U_2_i[j]);
        
        // - i * (2 - etq) * mu * q * U1
        F2 -= mu * j*this->q1 * (2 - eta) * (imag * U1);
        
        // - i * mu * q * int(H12, U1)
        F2 -= mu * j*this->q1 * imag * this->H12_pi[j]->convolve(this->U_1_r[j],
                                                                 this->U_1_i[j]);
        
        // set values to internal force
        internal_1_fd[j][0] = std::real(F1); // real part
        internal_1_fd[j][1] = std::imag(F1); // imag part
        internal_2_fd[j][0] = std::real(F2); // real part
        internal_2_fd[j][1] = std::imag(F2); // imag part
    }
    
}
/* -------------------------------------------------------------------------- */
// residual = (internal + external) * side_factor
void HalfSpace::computeResidual(NodalField * external_1,
				NodalField * external_2) {
  

  double * int_1_p = this->internal_1->storage();
  double * int_2_p = this->internal_2->storage();

  double * ext_1_p = external_1->storage();
  double * ext_2_p = external_2->storage();

  double * res_1_p = this->residual_1->storage();
  double * res_2_p = this->residual_2->storage();

  for (unsigned int n=0; n<this->nb_nodes; ++n) {
    res_1_p[n] = this->side_factor * (int_1_p[n] + ext_1_p[n]);
    res_2_p[n] = this->side_factor * (int_2_p[n] + ext_2_p[n]);
  }
}

/* -------------------------------------------------------------------------- */
// velocity = cs / mu       * residual (for shear components)
// velocity = cs / mu / eta * residual (for normal component)
void HalfSpace::computeVelocity() {

  double mu = this->material->getShearModulus();
  double Cs = this->material->getCs();
  double Cp = this->material->getCp();
  double eta = Cp/Cs;

  double * velo_1_p = this->velo_1->storage();
  double * velo_2_p = this->velo_2->storage();

  double * res_1_p = this->residual_1->storage();
  double * res_2_p = this->residual_2->storage();

  for (unsigned int n=0; n<this->nb_nodes; ++n) {
    velo_1_p[n] = Cs / mu * res_1_p[n];
    velo_2_p[n] = Cs / mu / eta * res_2_p[n];
  }
}

