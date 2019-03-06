/**
 * @file unimat_interface.cc
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
#include "unimat_interface.hh"

/* -------------------------------------------------------------------------- */
UnimatInterface::UnimatInterface(double length, unsigned int nb_elements,
			       Material * top_material,
			       Kernel * top_H11, Kernel * top_H12, Kernel * top_H22,
			       InterfaceLaw * law) :
  Interface(length, nb_elements, law),
  top(length, nb_elements,  1)
{
  this->top.setMaterial(top_material);
  this->top.setKernels(top_H11, top_H12, top_H22);
}

/* -------------------------------------------------------------------------- */
UnimatInterface::~UnimatInterface() { }

/* -------------------------------------------------------------------------- */
void UnimatInterface::init() {
  this->top.initConvolutions();
  Interface::init();
}

/* -------------------------------------------------------------------------- */
double UnimatInterface::getStableTimeStep() {
  double delta_x = this->length / (double)(this->nb_nodes);
  double max_cs = this->top.getMaterial().getCs();
  return delta_x / max_cs;
}

/* -------------------------------------------------------------------------- */
void UnimatInterface::setTimeStep(double time_step) { 
  Interface::setTimeStep(time_step);
  this->top.setTimeStep(time_step);
}

/* -------------------------------------------------------------------------- */
void UnimatInterface::computeDisplacement() {
  this->top.computeDisplacement();
}

/* -------------------------------------------------------------------------- */
void UnimatInterface::forwardFFT() {
  this->top.forwardFFT();
}

/* -------------------------------------------------------------------------- */
void UnimatInterface::computeStressFourierCoeff() {
  this->top.computeStressFourierCoeff();
}

/* -------------------------------------------------------------------------- */
void UnimatInterface::backwardFFT() {
  this->top.backwardFFT();
}

/* -------------------------------------------------------------------------- */
void UnimatInterface::computeResidual() {

  NodalField load_and_cohesion_1(this->nb_nodes);
  NodalField load_and_cohesion_2(this->nb_nodes);

  this->combineLoadAndCohesion(load_and_cohesion_1,
			       load_and_cohesion_2);

  this->top.computeResidual(&load_and_cohesion_1, &load_and_cohesion_2);
}

/* -------------------------------------------------------------------------- */
void UnimatInterface::computeVelocity() {
  this->top.computeVelocity();
}

/* -------------------------------------------------------------------------- */
void UnimatInterface::closingNormalGapForce(NodalField * close_force) {
  // top material information
  const Material & mat_t = this->top.getMaterial();
  double cp_t = mat_t.getCp();
  double cs_t = mat_t.getCs();
  double mu_t = mat_t.getShearModulus();
  double eta_t = cp_t / cs_t;
  double fact_t = this->time_step * cs_t / mu_t / eta_t; 

  // accessors
  double * u_2_t = this->top.getDisp2()->storage();
  double * f_2_t = this->top.getInternal2()->storage();
  double * t0_2 = this->load_2->storage();
  double * cf = close_force->storage();

  // PBC node:
  double u_2_gap = 2 * u_2_t[0];
  double f_2_gap = f_2_t[0];
  cf[0] = 0.5 * u_2_gap / fact_t + t0_2[0] + f_2_gap;
  for (unsigned int n=1; n<this->nb_nodes; ++n) {
    u_2_gap = u_2_t[n] + u_2_t[this->nb_nodes - n];
    f_2_gap = 0.5 * (f_2_t[n] + f_2_t[this->nb_nodes - n]);
    cf[n] = 0.5 * u_2_gap / fact_t + t0_2[n] + f_2_gap;
  }
}

/* -------------------------------------------------------------------------- */
void UnimatInterface::maintainShearGapForce(NodalField * maintain_force) {

  // accessors
  double * f_1_t = this->top.getInternal1()->storage();
  double * t0_1 = this->load_1->storage();
  double * mf = maintain_force->storage();

  mf[0] = t0_1[0] + f_1_t[0];
  for (unsigned int n=1; n<this->nb_nodes; ++n) {
    mf[n] = t0_1[n] + 0.5 * (f_1_t[n] + f_1_t[this->nb_nodes - n]);
  }
}

/* -------------------------------------------------------------------------- */
void UnimatInterface::computeGap1(NodalField * gap_1) {

  double * top_disp_1 = this->top.getDisp1()->storage();
  double * gap_1_p = gap_1->storage();

  gap_1_p[0] = 2*top_disp_1[0];
  for (unsigned int n=1; n<this->nb_nodes; ++n) {
    gap_1_p[n] = top_disp_1[n] + top_disp_1[this->nb_nodes-n];
  }
}

/* -------------------------------------------------------------------------- */
void UnimatInterface::computeGap2(NodalField * gap_2) {

  double * top_disp_2 = this->top.getDisp2()->storage();
  double * gap_2_p = gap_2->storage();

  gap_2_p[0] = top_disp_2[0];
  for (unsigned int n=1; n<this->nb_nodes; ++n) {
    gap_2_p[n] = top_disp_2[n] + top_disp_2[this->nb_nodes-n];
  }
}

/* -------------------------------------------------------------------------- */
void UnimatInterface::registerDumpField(const std::string & field_name) {

  // disp
  if (field_name == "top_disp_1") {
    this->registerForDump(field_name,
			  this->top.getDisp1());
  }
  else if (field_name == "top_disp_2") {
    this->registerForDump(field_name,
			  this->top.getDisp2());
  }

  // velo
  else if (field_name == "top_velo_1") {
    this->registerForDump(field_name,
			  this->top.getVelo1());
  }
  else if (field_name == "top_velo_2") {
    this->registerForDump(field_name,
			  this->top.getVelo2());
  }

  // residual
  else if (field_name == "top_residual_1") {
    this->registerForDump(field_name,
			  this->top.getResidual1());
  }
  else if (field_name == "top_residual_2") {
    this->registerForDump(field_name,
			  this->top.getResidual2());
  }

  // internal
  else if (field_name == "top_internal_1") {
    this->registerForDump(field_name,
			  this->top.getInternal1());
  }
  else if (field_name == "top_internal_2") {
    this->registerForDump(field_name,
			  this->top.getInternal2());
  }

  else {
    Interface::registerDumpField(field_name);
  }
}
