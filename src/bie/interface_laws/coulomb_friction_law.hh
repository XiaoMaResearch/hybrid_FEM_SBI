/**
 * @file coulomb_friction_law.hh
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
#ifndef __COULOMB_FRICTION_LAW_H__
#define __COULOMB_FRICTION_LAW_H__
/* -------------------------------------------------------------------------- */
#include "interface_law.hh"

/* 
   Simple Coulomb friction law, which has a constant friction coefficient
   independent of slip or slip velocity.
   Can vary in space.
 */

class CoulombFrictionLaw : public InterfaceLaw {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  CoulombFrictionLaw(unsigned int nb_nodes, 
		     double fric_coef_default);
  virtual ~CoulombFrictionLaw() {};
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void computeCohesiveForces(NodalField * cohesion_1,
			     NodalField * cohesion_2);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  NodalField * getFrictionCoefficient() { return &(this->fric_coef); };
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  NodalField fric_coef;
};

//#include "coulomb_friction_law_impl.cc"

#endif /* __COULOMB_FRICTION_LAW_H__ */
