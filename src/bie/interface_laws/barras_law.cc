/**
 * @file barras_law.cc
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
#include "barras_law.hh"
#include "interface.hh"

#include <cmath>

/* -------------------------------------------------------------------------- */
BarrasLaw::BarrasLaw(unsigned int nb_nodes,
		     double tau_max_default,
		     double delta_c_default) : 
  InterfaceLaw(), 
  tau_max(nb_nodes),
  delta_c(nb_nodes)
{
  this->tau_max.setAllValuesTo(tau_max_default);
  this->delta_c.setAllValuesTo(delta_c_default);
}

/* -------------------------------------------------------------------------- */
void BarrasLaw::computeCohesiveForces(NodalField * cohesion_1,
				      NodalField * cohesion_2) {

  unsigned int nb = this->interface->getNbNodes();

  // find current gap
  NodalField gap_1(nb);
  this->interface->computeGap1(&gap_1);
  double * gap_1_p = gap_1.storage();

  NodalField gap_2(nb);
  this->interface->computeGap2(&gap_2);
  double * gap_2_p = gap_2.storage();

  // find forces needed to close normal gap
  this->interface->closingNormalGapForce(cohesion_2);
  double * coh_2_p = cohesion_2->storage();

  // find force needed to maintain shear gap
  this->interface->maintainShearGapForce(cohesion_1);
  double * coh_1_p = cohesion_1->storage();

  // interface properties
  double * tau_max_p = this->tau_max.storage();
  double * delta_c_p = this->delta_c.storage();

  for (unsigned int n = 0; n<nb; ++n) {
    double relative_gap = std::sqrt(gap_1_p[n]*gap_1_p[n] + gap_2_p[n]*gap_2_p[n]) / delta_c_p[n];
    double strength = tau_max_p[n] * std::max(0., 1 - relative_gap); // which is always >= 0

    // avoid penetration "at any cost" and apply maximal allowed normal cohesive force
    // coh_2_p > 0 is a adhesive force
    // coh_2_p < 0 is a contact pressure
    coh_2_p[n] = std::min(coh_2_p[n], strength);

    // maximal shear cohesive force given by strength. Keep orientation of shear force
    double strength_ratio = std::abs(strength / coh_1_p[n]);
    if (strength_ratio < 1)
      coh_1_p[n] = strength_ratio * coh_1_p[n];
  }
}
