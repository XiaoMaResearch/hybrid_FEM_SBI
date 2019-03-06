/**
 * @file nodal_field.cc
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
 * MIT License
 * Copyright (c) 2017 David S. Kammer 
 *
 */
#include "nodal_field.hh"

/* -------------------------------------------------------------------------- */
NodalField::NodalField(unsigned int nb_nodes) :
  nb_nodes(nb_nodes) {
  
  this->field = new double[nb_nodes];
  this->zeros();
}

/* -------------------------------------------------------------------------- */
NodalField::~NodalField() {
  delete[] this->field;
}

/* -------------------------------------------------------------------------- */
void NodalField::zeros() {
  this->setAllValuesTo(0.);
}

/* -------------------------------------------------------------------------- */
void NodalField::setAllValuesTo(double value) {
  for (unsigned int n=0; n<this->nb_nodes; ++n)
      this->field[n] = value;
}
/* -------------------------------------------------------------------------- */
void NodalField::setValuesTo(double *ptr_value, int spacing) {
    for (unsigned int n=0; n<this->nb_nodes; ++n)
        this->field[n] = ptr_value[spacing*n];
}


