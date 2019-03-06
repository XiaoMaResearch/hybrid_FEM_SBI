/**
 * @file interface.hh
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
#ifndef __INTERFACE_H__
#define __INTERFACE_H__
/* -------------------------------------------------------------------------- */
#include "interface_law.hh"
#include "dumper.hh"

class Interface : public Dumper {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  Interface(double length, unsigned int nb_nodes,
	    InterfaceLaw * law);
  virtual ~Interface();
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  // get (limiting) stable time step
  virtual double getStableTimeStep() = 0;

  // initiate interface model and half spaces
  virtual void init(bool velocity_initial_conditions = false);

  // initiate dump (overload to get coordinates)
  void initDump(const std::string & bname,
		const std::string & path);

  // iteration of advancing one time step
  virtual void advanceTimeStep();

  // functions used during time stepping for each half-space
  virtual void computeDisplacement() = 0;
  virtual void forwardFFT() = 0;
  virtual void computeStressFourierCoeff() = 0;
  virtual void backwardFFT() = 0;
  virtual void computeCohesion();
  virtual void computeResidual() = 0;
  virtual void computeVelocity() = 0;
  
  // function that combine load and cohesionw with correct signs
  void combineLoadAndCohesion(NodalField & load_and_cohesion_1,
			      NodalField & load_and_cohesion_2);

  // compute force needed to close normal gap
  virtual void closingNormalGapForce(NodalField * close_force) = 0;
  
  // compute force needed to maintain current shear gap
  virtual void maintainShearGapForce(NodalField * maintain_force) = 0;

  // compute gap in displacement
  virtual void computeGap1(NodalField * gap_1) = 0;
  virtual void computeGap2(NodalField * gap_2) = 0;

  // dumper function
  virtual void registerDumpField(const std::string & field_name);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  unsigned int getNbNodes() const { return this->nb_nodes; };

  NodalField * getShearLoad()  { return this->load_1; };
  NodalField * getNormalLoad() { return this->load_2; };

  NodalField * getCohesion1() { return this->cohesion_1; };
  NodalField * getCohesion2() { return this->cohesion_2; };

  virtual void setTimeStep(double time_step);
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
  
  // external loading
  NodalField * load_1;
  NodalField * load_2;

  // interface forces (e.g., cohesion)
  NodalField * cohesion_1;
  NodalField * cohesion_2;

  // interface law: cohesive law and contact/friction law
  InterfaceLaw * law;
};

//#include "interface_impl.cc"

#endif /* __INTERFACE_H__ */
