/**
 * @file defrig_interface.cc
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
#include "defrig_interface.hh"

/* -------------------------------------------------------------------------- */
DefRigInterface::DefRigInterface(double length, unsigned int nb_nodes,
				 Material * top_material,
				 Kernel * top_H11, Kernel * top_H12, Kernel * top_H22,
				 InterfaceLaw * law) :
  Interface(length, nb_nodes, law),
  top(length, nb_nodes,  1)
{
  this->top.setMaterial(top_material);
  this->top.setKernels(top_H11, top_H12, top_H22);

  // top material information
  const Material & mat_t = this->top.getMaterial();
  double cp_t = mat_t.getCp();
  double cs_t = mat_t.getCs();
  double mu_t = mat_t.getShearModulus();
  double eta_t = cp_t / cs_t;
  this->fact_t_2 = cs_t / mu_t / eta_t; 
}

/* -------------------------------------------------------------------------- */
DefRigInterface::~DefRigInterface() { }

/* -------------------------------------------------------------------------- */
void DefRigInterface::init() {
  this->top.initConvolutions();
  Interface::init();
}

/* -------------------------------------------------------------------------- */
double DefRigInterface::getStableTimeStep() {
  double delta_x = this->length / (double)(this->nb_nodes);
  double cs = this->top.getMaterial().getCs();
  return delta_x / cs;
}

/* -------------------------------------------------------------------------- */
void DefRigInterface::setTimeStep(double time_step) { 
  Interface::setTimeStep(time_step);
  this->top.setTimeStep(time_step);
}

/* -------------------------------------------------------------------------- */
void DefRigInterface::computeDisplacement() {
  this->top.computeDisplacement();
}

/* -------------------------------------------------------------------------- */
void DefRigInterface::forwardFFT() {
  this->top.forwardFFT();
}

/* -------------------------------------------------------------------------- */
void DefRigInterface::computeStressFourierCoeff() {
  this->top.computeStressFourierCoeff();
}

/* -------------------------------------------------------------------------- */
void DefRigInterface::backwardFFT() {
  this->top.backwardFFT();
}

/* -------------------------------------------------------------------------- */
void DefRigInterface::computeResidual() {

  NodalField load_and_cohesion_1(this->nb_nodes);
  NodalField load_and_cohesion_2(this->nb_nodes);

  this->combineLoadAndCohesion(load_and_cohesion_1,
			       load_and_cohesion_2);

  this->top.computeResidual(&load_and_cohesion_1, &load_and_cohesion_2);
}

/* -------------------------------------------------------------------------- */
void DefRigInterface::computeVelocity() {
  this->top.computeVelocity();
}

/* -------------------------------------------------------------------------- */
void DefRigInterface::closingNormalGapForce(NodalField * close_force) {
  // C factor of notes
  double fact_t = this->time_step * this->fact_t_2;

  // accessors
  double * f_2_t = this->top.getInternal2()->storage();
  double * t0_2 = this->load_2->storage();
  double * cf = close_force->storage();

  NodalField gap_2(this->nb_nodes);
  this->computeGap2(&gap_2);
  double * gap_2_p = gap_2.storage();

  for (unsigned int n=0; n<this->nb_nodes; ++n) {
    double u_2_gap = gap_2_p[n] / fact_t;
    double du_2_t = t0_2[n] + f_2_t[n];
    cf[n] = u_2_gap + du_2_t;
  }
}

/* -------------------------------------------------------------------------- */
void DefRigInterface::maintainShearGapForce(NodalField * maintain_force) {

  // accessors
  double * f_1_t = this->top.getInternal1()->storage();
  double * t0_1 = this->load_1->storage();
  double * mf = maintain_force->storage();

  for (unsigned int n=0; n<this->nb_nodes; ++n) {
    mf[n] = t0_1[n] + f_1_t[n];
  }
}

/* -------------------------------------------------------------------------- */
void DefRigInterface::computeGap1(NodalField * gap_1) {

  double * top_disp_1 = this->top.getDisp1()->storage();
  double * gap_1_p = gap_1->storage();

  for (unsigned int n=0; n<this->nb_nodes; ++n) {
    gap_1_p[n] = top_disp_1[n];
  }
}

/* -------------------------------------------------------------------------- */
void DefRigInterface::computeGap2(NodalField * gap_2) {

  double * top_disp_2 = this->top.getDisp2()->storage();
  double * gap_2_p = gap_2->storage();

  for (unsigned int n=0; n<this->nb_nodes; ++n) {
    gap_2_p[n] = top_disp_2[n];
  }
}

/* -------------------------------------------------------------------------- */
void DefRigInterface::registerDumpField(const std::string & field_name) {

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
