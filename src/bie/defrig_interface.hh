/**
 * @file defrig_interface.hh
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
#ifndef __DEFRIG_INTERFACE_H__
#define __DEFRIG_INTERFACE_H__
/* -------------------------------------------------------------------------- */
#include "half_space.hh"
#include "interface_law.hh"
#include "interface.hh"

class DefRigInterface : public Interface {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  DefRigInterface(double length, unsigned int nb_nodes,
		  Material * top_material,
		  Kernel * top_H11, Kernel * top_H12, Kernel * top_H22,
		  InterfaceLaw * law);
  virtual ~DefRigInterface();
  
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

  HalfSpace & getTop() { return this->top; };

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  // half spaces
  HalfSpace top;

  double fact_t_2;
};

//#include "defrig_interface_impl.cc"

#endif /* __DEFRIG_INTERFACE_H__ */
