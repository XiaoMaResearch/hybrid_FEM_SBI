/**
 * @file half_space.hh
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
#ifndef __HALF_SPACE_H__
#define __HALF_SPACE_H__
/* -------------------------------------------------------------------------- */
#include "material.hh"
#include "kernel.hh"
#include "preint_kernel.hh"
#include "nodal_field.hh"
#include "fftable_nodal_field.hh"
#include "limited_history.hh"

/* -------------------------------------------------------------------------- */
class HalfSpace {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  // side factor top=1 bot=-1
  HalfSpace(double length, unsigned int nb_nodes, int side_factor);
  virtual ~HalfSpace();
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  void computeDisplacement();
  void computeStressFourierCoeff();
  void computeStressFourierCoeff_correction();
  void computeResidual(NodalField * external_1,
		       NodalField * external_2);
  void computeVelocity();

  // init convolutions
  void initConvolutions();

  // apply fft forward on displacement
  void forwardFFT();
  // apply fft backward on internal
  void backwardFFT();

private:
  void computeDisplacement(FFTableNodalField * disp,
			   NodalField * velo);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  // set a material to half space
  void setMaterial(Material * material) { this->material = material; };
  // get material of half space
  const Material & getMaterial() const { return (*this->material); };
 
  // set Kernels to half space
  void setKernels(Kernel * H11, Kernel * H12, Kernel * H22) { 
    this->H11 = H11;
    this->H12 = H12;
    this->H22 = H22;
  };

  void setTimeStep(double time_step) { this->time_step = time_step; };

  FFTableNodalField * getDisp1() { return this->disp_1; };
  FFTableNodalField * getDisp2() { return this->disp_2; };

  NodalField * getVelo1() { return this->velo_1; };
  NodalField * getVelo2() { return this->velo_2; };

  FFTableNodalField * getInternal1() { return this->internal_1; };
  FFTableNodalField * getInternal2() { return this->internal_2; };

  NodalField * getResidual1() { return this->residual_1; };
  NodalField * getResidual2() { return this->residual_2; };

  // const accessors
  const FFTableNodalField * getDisp1() const { return this->disp_1; };
  const FFTableNodalField * getDisp2() const { return this->disp_2; };


  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  // length of domain / replication length
  double length;

  // number of nodes
  unsigned int nb_nodes;

  // time step
  double time_step;

  // used to know in which directions the tractions pull
  int side_factor;

  // some useful info that is constructed automatically
  unsigned int nb_fft;
  double q1; // fundamental mode for j=1: q1 = 2*Pi/X

  // material properties
  Material * material;
  Kernel * H11;
  Kernel * H12;
  Kernel * H22;

  // displacement 1=shear 2=normal
  FFTableNodalField * disp_1;
  FFTableNodalField * disp_2;

  // past values of displacement in frequency domain
  // each LimitedHistory is for a given wave number q 
  std::vector<LimitedHistory *> U_1_r; // U_1 real part
  std::vector<LimitedHistory *> U_1_i; // U_1 imag part
  std::vector<LimitedHistory *> U_2_r;
  std::vector<LimitedHistory *> U_2_i;

  // convolutions
  std::vector<PreintKernel *> H11_pi;
  std::vector<PreintKernel *> H12_pi;
  std::vector<PreintKernel *> H22_pi;

  // velocity 1=shear 2=normal
  NodalField * velo_1;
  NodalField * velo_2;
  // velocity 1 = shear 2= normal 2nd from the frist iteration

  // tractions due to deformation
  FFTableNodalField * internal_1;
  FFTableNodalField * internal_2;

  // all acting forces
  NodalField * residual_1;
  NodalField * residual_2;

};

//#include "half_space_impl.cc"

#endif /* __HALF_SPACE_H__ */
