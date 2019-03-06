/**
 * @file interface.cc
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
#include "interface.hh"

/* -------------------------------------------------------------------------- */
Interface::Interface(double length, unsigned int nb_nodes,
		     InterfaceLaw * law) :
  Dumper(),
  length(length),
  nb_nodes(nb_nodes),
  time_step(0.),
  law(law)
{
  this->load_1 = new NodalField(nb_nodes);
  this->load_2 = new NodalField(nb_nodes);

  this->cohesion_1 = new NodalField(nb_nodes);
  this->cohesion_2 = new NodalField(nb_nodes);

  this->law->setInterface(this);
}

/* -------------------------------------------------------------------------- */
Interface::~Interface() {
  delete this->load_1;
  delete this->load_2;

  delete this->cohesion_1;
  delete this->cohesion_2;
}

/* -------------------------------------------------------------------------- */
void Interface::init(bool velocity_initial_conditions) { 

  // like a typical time step (advanceTimeStep)
  // but displacement is potentially already imposed as initial condition
  this->forwardFFT();
  this->computeStressFourierCoeff();
  this->backwardFFT();
  this->computeCohesion();
  this->computeResidual();

  if (!velocity_initial_conditions) 
    this->computeVelocity();
}

/* -------------------------------------------------------------------------- */
void Interface::initDump(const std::string & bname,
			 const std::string & path) { 
  Dumper::initDump(bname, path);
  
  double dx = this->length / double(this->nb_nodes);
  double * coords = new double[this->nb_nodes];
  for (unsigned int i=0; i<this->nb_nodes; ++i)
    coords[i] = i*dx;
  Dumper::setCoords(1, this->nb_nodes, coords);
  delete[] coords;
}

/* -------------------------------------------------------------------------- */
void Interface::setTimeStep(double time_step) { 
  this->time_step = time_step; 
}

/* -------------------------------------------------------------------------- */
void Interface::advanceTimeStep() {
  // compute displacement
  this->computeDisplacement();

  // to fourier space
  this->forwardFFT();
  
  // compute convolutions
  this->computeStressFourierCoeff();
  
  // back to normal space
  this->backwardFFT();

  // compute forces due to interface law
  this->computeCohesion();
  
  // compute residual force
  this->computeResidual();

  // compute velocity
  this->computeVelocity();
}

/* -------------------------------------------------------------------------- */
void Interface::computeCohesion() {
  this->law->computeCohesiveForces(this->cohesion_1,
				   this->cohesion_2);
}

/* -------------------------------------------------------------------------- */
void Interface::combineLoadAndCohesion(NodalField & load_and_cohesion_1,
				       NodalField & load_and_cohesion_2) {

  // tau_0 - tau_coh
  for (unsigned int n=0; n<this->nb_nodes; ++n) {
    load_and_cohesion_1(n) = (*this->load_1)(n) - (*this->cohesion_1)(n);
    load_and_cohesion_2(n) = (*this->load_2)(n) - (*this->cohesion_2)(n);
  }
}

/* -------------------------------------------------------------------------- */
void Interface::registerDumpField(const std::string & field_name) {

  // load 
  if (field_name ==  "load_1") {
    this->registerForDump(field_name,
			  this->load_1);
  }
  else if (field_name == "load_2") {
    this->registerForDump(field_name,
			  this->load_2);
  }

  // cohesion 
  else if (field_name ==  "cohesion_1") {
    this->registerForDump(field_name,
			  this->cohesion_1);
  }
  else if (field_name == "cohesion_2") {
    this->registerForDump(field_name,
			  this->cohesion_2);
  }

  else {
    // used to be like this. Now give option to InterfaceLaw to registerDumpField as well.
    //std::cout << "Do not know dump field with name: " << field_name << std::endl;
    this->law->registerDumpField(field_name);
  }
}
