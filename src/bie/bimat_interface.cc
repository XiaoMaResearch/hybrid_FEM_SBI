/**
 * @file bimat_interface.cc
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
#include "bimat_interface.hh"

/* -------------------------------------------------------------------------- */
BimatInterface::BimatInterface(double length, unsigned int nb_nodes,
			       Material * top_material,
			       Kernel * top_H11, Kernel * top_H12, Kernel * top_H22,
			       Material * bot_material,
			       Kernel * bot_H11, Kernel * bot_H12, Kernel * bot_H22,
			       InterfaceLaw * law) :
  Interface(length, nb_nodes, law),
  top(length, nb_nodes,  1), 
  bot(length, nb_nodes, -1)
{
  this->top.setMaterial(top_material);
  this->top.setKernels(top_H11, top_H12, top_H22);
  
  this->bot.setMaterial(bot_material);
  this->bot.setKernels(bot_H11, bot_H12, bot_H22);
}

/* -------------------------------------------------------------------------- */
BimatInterface::~BimatInterface() { }

/* -------------------------------------------------------------------------- */
void BimatInterface::init() {
  this->top.initConvolutions();
  this->bot.initConvolutions();
  Interface::init();
}

/* -------------------------------------------------------------------------- */
double BimatInterface::getStableTimeStep() {
  double delta_x = this->length / (double)(this->nb_nodes);
  double max_cs = std::max(this->top.getMaterial().getCs(),
			   this->bot.getMaterial().getCs());
  return delta_x / max_cs;
}

/* -------------------------------------------------------------------------- */
void BimatInterface::setTimeStep(double time_step) { 
  Interface::setTimeStep(time_step);
  this->top.setTimeStep(time_step);
  this->bot.setTimeStep(time_step);
}

/* -------------------------------------------------------------------------- */
void BimatInterface::computeDisplacement() {
  this->top.computeDisplacement();
  this->bot.computeDisplacement();
}

/* -------------------------------------------------------------------------- */
void BimatInterface::forwardFFT() {
  this->top.forwardFFT();
  this->bot.forwardFFT();
}

/* -------------------------------------------------------------------------- */
void BimatInterface::computeStressFourierCoeff() {
  this->top.computeStressFourierCoeff();
  this->bot.computeStressFourierCoeff();
}

/* -------------------------------------------------------------------------- */
void BimatInterface::backwardFFT() {
  this->top.backwardFFT();
  this->bot.backwardFFT();
}

/* -------------------------------------------------------------------------- */
void BimatInterface::computeResidual() {

  NodalField load_and_cohesion_1(this->nb_nodes);
  NodalField load_and_cohesion_2(this->nb_nodes);

  this->combineLoadAndCohesion(load_and_cohesion_1,
			       load_and_cohesion_2);

  this->top.computeResidual(&load_and_cohesion_1, &load_and_cohesion_2);
  this->bot.computeResidual(&load_and_cohesion_1, &load_and_cohesion_2);
}

/* -------------------------------------------------------------------------- */
void BimatInterface::computeVelocity() {
  this->top.computeVelocity();
  this->bot.computeVelocity();
}

/* -------------------------------------------------------------------------- */
void BimatInterface::closingNormalGapForce(NodalField * close_force) {
  // top material information
  const Material & mat_t = this->top.getMaterial();
  double cp_t = mat_t.getCp();
  double cs_t = mat_t.getCs();
  double mu_t = mat_t.getShearModulus();
  double eta_t = cp_t / cs_t;
  double fact_t = this->time_step * cs_t / mu_t / eta_t; 

  // bot material information
  const Material & mat_b = this->bot.getMaterial();
  double cp_b = mat_b.getCp();
  double cs_b = mat_b.getCs();
  double mu_b = mat_b.getShearModulus();
  double eta_b = cp_b / cs_b;
  double fact_b = this->time_step * cs_b / mu_b / eta_b; 

  // accessors
  double * u_2_t = this->top.getDisp2()->storage();
  double * u_2_b = this->bot.getDisp2()->storage();
  double * f_2_t = this->top.getInternal2()->storage();
  double * f_2_b = this->bot.getInternal2()->storage();
  double * t0_2 = this->load_2->storage();
  double * cf = close_force->storage();

  double fact_cf = fact_t + fact_b;

  for (unsigned int n=0; n<this->nb_nodes; ++n) {
    double u_2_gap = u_2_t[n] - u_2_b[n];
    double du_2_t = fact_t * (t0_2[n] + f_2_t[n]);
    double du_2_b = fact_b * (t0_2[n] + f_2_b[n]);
    cf[n] = (u_2_gap + du_2_t + du_2_b) / fact_cf;
  }
}

/* -------------------------------------------------------------------------- */
void BimatInterface::maintainShearGapForce(NodalField * maintain_force) {
  // top material information
  const Material & mat_t = this->top.getMaterial();
  double cs_t = mat_t.getCs();
  double mu_t = mat_t.getShearModulus();
  double fact_t = cs_t / mu_t; 

  // bot material information
  const Material & mat_b = this->bot.getMaterial();
  double cs_b = mat_b.getCs();
  double mu_b = mat_b.getShearModulus();
  double fact_b = cs_b / mu_b; 

  // accessors
  double * f_1_t = this->top.getInternal1()->storage();
  double * f_1_b = this->bot.getInternal1()->storage();
  double * t0_1 = this->load_1->storage();
  double * mf = maintain_force->storage();

  double fact_mf = fact_t + fact_b;

  for (unsigned int n=0; n<this->nb_nodes; ++n) {
    mf[n] = t0_1[n] + (fact_t * f_1_t[n] + fact_b * f_1_b[n]) / fact_mf;
  }
}

/* -------------------------------------------------------------------------- */
void BimatInterface::computeGap1(NodalField * gap_1) {

  double * top_disp_1 = this->top.getDisp1()->storage();
  double * bot_disp_1 = this->bot.getDisp1()->storage();
  double * gap_1_p = gap_1->storage();

  for (unsigned int n=0; n<this->nb_nodes; ++n) {
    gap_1_p[n] = top_disp_1[n] - bot_disp_1[n];
  }
}

/* -------------------------------------------------------------------------- */
void BimatInterface::computeGap2(NodalField * gap_2) {

  double * top_disp_2 = this->top.getDisp2()->storage();
  double * bot_disp_2 = this->bot.getDisp2()->storage();
  double * gap_2_p = gap_2->storage();

  for (unsigned int n=0; n<this->nb_nodes; ++n) {
    gap_2_p[n] = top_disp_2[n] - bot_disp_2[n];
  }
}

/* -------------------------------------------------------------------------- */
void BimatInterface::registerDumpField(const std::string & field_name) {

  // disp
  if (field_name == "top_disp_1") {
    this->registerForDump(field_name,
			  this->top.getDisp1());
  }
  else if (field_name == "top_disp_2") {
    this->registerForDump(field_name,
			  this->top.getDisp2());
  }
  else if (field_name == "bot_disp_1") {
    this->registerForDump(field_name,
			  this->bot.getDisp1());
  }
  else if (field_name == "bot_disp_2") {
    this->registerForDump(field_name,
			  this->bot.getDisp2());
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
  else if (field_name == "bot_velo_1") {
    this->registerForDump(field_name,
			  this->bot.getVelo1());
  }
  else if (field_name == "bot_velo_2") {
    this->registerForDump(field_name,
			  this->bot.getVelo2());
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
  else if (field_name == "bot_residual_1") {
    this->registerForDump(field_name,
			  this->bot.getResidual1());
  }
  else if (field_name == "bot_residual_2") {
    this->registerForDump(field_name,
			  this->bot.getResidual2());
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
  else if (field_name == "bot_internal_1") {
    this->registerForDump(field_name,
			  this->bot.getInternal1());
  }
  else if (field_name == "bot_internal_2") {
    this->registerForDump(field_name,
			  this->bot.getInternal2());
  }


  else {
    Interface::registerDumpField(field_name);
  }
}
