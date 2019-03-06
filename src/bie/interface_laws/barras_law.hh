/**
 * @file barras_law.hh
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
#ifndef __BARRAS_LAW_H__
#define __BARRAS_LAW_H__
/* -------------------------------------------------------------------------- */
#include "interface_law.hh"

/* 
   This is not exactly the same law than used in Barras et al. 2014. 
   Within the cohesive zone, this law can close in normal direction, 
   whereas Barras' law will only maintain the normal gap.
   Also, here no friction behind the cohesive zone is applied.
 */

class BarrasLaw : public InterfaceLaw {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  BarrasLaw(unsigned int nb_nodes, 
	    double tau_max_default, double delta_c_default);
  virtual ~BarrasLaw() {};
  
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
  NodalField * getTauMax() { return &(this->tau_max); };
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  NodalField tau_max;
  NodalField delta_c;

};

//#include "barras_law_impl.cc"

#endif /* __BARRAS_LAW_H__ */
