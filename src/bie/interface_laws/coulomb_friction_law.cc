/**
 * @file coulomb_friction_law.cc
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
#include "coulomb_friction_law.hh"
#include "interface.hh"

#include <cmath>

/* -------------------------------------------------------------------------- */
CoulombFrictionLaw::CoulombFrictionLaw(unsigned int nb_nodes,
				       double fric_coef_default) :
  InterfaceLaw(), 
  fric_coef(nb_nodes)
{
  this->fric_coef.setAllValuesTo(fric_coef_default);
}

/* -------------------------------------------------------------------------- */
void CoulombFrictionLaw::computeCohesiveForces(NodalField * cohesion_1,
					       NodalField * cohesion_2) {

  unsigned int nb = this->interface->getNbNodes();

  // find forces needed to close normal gap
  this->interface->closingNormalGapForce(cohesion_2);
  double * coh_2_p = cohesion_2->storage();

  // find force needed to maintain shear gap
  this->interface->maintainShearGapForce(cohesion_1);
  double * coh_1_p = cohesion_1->storage();

  // interface properties
  double * fric_coef_p = this->fric_coef.storage();

  for (unsigned int n = 0; n<nb; ++n) {
    
    // avoid penetration "at any cost"
    // apply no normal cohesive force
    // coh_2_p > 0 is a adhesive force
    // coh_2_p < 0 is a contact pressure
    coh_2_p[n] = std::min(coh_2_p[n], 0.);

    double strength = fric_coef_p[n] * coh_2_p[n];

    // maximal shear cohesive force given by strength. 
    // keep orientation of shear force
    double strength_ratio = std::abs(strength / coh_1_p[n]);
    if (strength_ratio < 1)
      coh_1_p[n] = strength_ratio * coh_1_p[n];
  }
}
