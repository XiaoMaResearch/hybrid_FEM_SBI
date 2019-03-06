/**
 * @file linear_shear_cohesive_law.cc
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
#include "linear_shear_cohesive_law.hh"
#include "interface.hh"

#include <cmath>

/* -------------------------------------------------------------------------- */
LinearShearCohesiveLaw::LinearShearCohesiveLaw(unsigned int nb_nodes,
					       double *Gc_default,
					       double *tau_c_default,
					       double tau_r_default) :
  InterfaceLaw(), 
  G_c(nb_nodes),
  tau_c(nb_nodes),
  tau_r(nb_nodes)
{
//  this->G_c.setAllValuesTo(Gc_default);
  this->G_c.setValuesTo(Gc_default);
//  this->tau_c.setAllValuesTo(tau_c_default);
  this->tau_c.setValuesTo(tau_c_default);
  this->tau_r.setAllValuesTo(tau_r_default);
}

/* -------------------------------------------------------------------------- */
void LinearShearCohesiveLaw::computeCohesiveForces(NodalField * cohesion_1,
						   NodalField * cohesion_2) {

  unsigned int nb = this->interface->getNbNodes();

  // find forces needed to close normal gap
  this->interface->closingNormalGapForce(cohesion_2);
  double * coh_2_p = cohesion_2->storage();

  // find force needed to maintain shear gap
  this->interface->maintainShearGapForce(cohesion_1);
  double * coh_1_p = cohesion_1->storage();

  NodalField shear_gap(nb);
  this->interface->computeGap1(&shear_gap);
  double * shear_gap_p = shear_gap.storage();

  // interface properties
  double * G_c_p = this->G_c.storage();
  double * tau_c_p = this->tau_c.storage();
  double * tau_r_p = this->tau_r.storage();

  for (unsigned int n = 0; n<nb; ++n) {

    // coh_2_p > 0 is a adhesive force
    // coh_2_p < 0 is a contact pressure
    
    // avoid penetration "at any cost"
    // avoid opening "at any cost"
    // thus: do nothing

    // avoid penetration "at any cost"
    // apply no normal cohesive force
    //coh_2_p[n] = std::min(coh_2_p[n], 0.);

    double slope = (tau_c_p[n] - tau_r_p[n]) * (tau_c_p[n] - tau_r_p[n]) / 2. / G_c_p[n];
    double strength = std::max(tau_c_p[n] - std::abs(shear_gap_p[n]) * slope, tau_r_p[n]);

    // maximal shear cohesive force given by strength. 
    // keep orientation of shear force
    double strength_ratio = std::abs(strength / coh_1_p[n]);
    if (strength_ratio < 1)
      coh_1_p[n] = strength_ratio * coh_1_p[n];
  }
}

/* -------------------------------------------------------------------------- */
void LinearShearCohesiveLaw::registerDumpField(const std::string & field_name) {

  // G_c
  if (field_name == "G_c") {
    this->interface->registerForDump(field_name,
				     &(this->G_c));
  }

  // tau_c
  else if (field_name == "tau_c") {
    this->interface->registerForDump(field_name,
				     &(this->tau_c));
  }

  // tau_r
  else if (field_name == "tau_r") {
    this->interface->registerForDump(field_name,
				     &(this->tau_r));
  }

  // do not know this field
  else {
    InterfaceLaw::registerDumpField(field_name);
  }

}
