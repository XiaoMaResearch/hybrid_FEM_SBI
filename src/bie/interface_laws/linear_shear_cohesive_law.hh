/**
 * @file linear_shear_cohesive_law.hh
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
#ifndef __LINEAR_SHEAR_COHESIVE_LAW_H__
#define __LINEAR_SHEAR_COHESIVE_LAW_H__
/* -------------------------------------------------------------------------- */
#include "interface_law.hh"

/* 
   Linear cohesive law in shear direction only. 
   No interpenetration allowed
   but also no opening allowed.
   Thus: should only be used for pure mode II fracture
 */

class LinearShearCohesiveLaw : public InterfaceLaw {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  LinearShearCohesiveLaw(unsigned int nb_nodes, 
			 double *Gc_default,
			 double *tau_c_default,
			 double tau_r_default = 0.);
  virtual ~LinearShearCohesiveLaw() {};
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void computeCohesiveForces(NodalField * cohesion_1,
			     NodalField * cohesion_2);

  // dumper function
  virtual void registerDumpField(const std::string & field_name);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  NodalField * getGc() { return &(this->G_c); };
  NodalField * getTauc() { return &(this->tau_c); };
  NodalField * getTaur() { return &(this->tau_r); };
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  NodalField G_c;
  NodalField tau_c;
  NodalField tau_r;
};

//#include "linear_shear_cohesive_law_impl.cc"

#endif /* __LINEAR_SHEAR_COHESIVE_LAW_H__ */
