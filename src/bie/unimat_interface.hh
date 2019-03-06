/**
 * @file unimat_interface.hh
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
#ifndef __UNIIMAT_INTERFACE_H__
#define __UNIIMAT_INTERFACE_H__
/* -------------------------------------------------------------------------- */
#include "half_space.hh"
#include "interface_law.hh"
#include "interface.hh"

/*
 * WARNING:
 * This unimaterial interface does only work for a central crack
 * with a normal load that is symmetric to the center of the crack
 * and a shear load that is anti-symmetric w.r.t. to the center of the crack
 */

class UnimatInterface : public Interface {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  UnimatInterface(double length, unsigned int nb_elements,
		 Material * top_material,
		 Kernel * top_H11, Kernel * top_H12, Kernel * top_H22,
		 InterfaceLaw * law);
  virtual ~UnimatInterface();
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  // get (limiting) stable time step
  virtual double getStableTimeStep();

  // initiate interface model and half spaces
  virtual void init();

  // functions used during time stepping for each half-space
  virtual void computeDisplacement();
  virtual void forwardFFT();
  virtual void computeStressFourierCoeff();
  virtual void backwardFFT();
  virtual void computeResidual();
  virtual void computeVelocity();
  
  // compute force needed to close normal gap
  virtual void closingNormalGapForce(NodalField * close_force);
  
  // compute force needed to maintain current shear gap
  virtual void maintainShearGapForce(NodalField * maintain_force);

  // compute gap in displacement
  virtual void computeGap1(NodalField * gap_1);
  virtual void computeGap2(NodalField * gap_2);

  // dumper function
  virtual void registerDumpField(const std::string & field_name);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  virtual void setTimeStep(double time_step);

  const HalfSpace & getTop() const { return this->top; };

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  // half spaces
  HalfSpace top;
};

//#include "unimat_interface_impl.cc"

#endif /* __UNIIMAT_INTERFACE_H__ */
