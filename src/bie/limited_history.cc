/**
 * @file limited_history.cc
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
#include "limited_history.hh"

/* -------------------------------------------------------------------------- */
LimitedHistory::LimitedHistory(unsigned int size) :
  nb_history_points(0),
  size(size) {

  this->index_now = 0;
  this->values = new double[size];
  
  for (unsigned int i=0; i<this->size; ++i)
    this->values[i] = 0.;
}

/* -------------------------------------------------------------------------- */
LimitedHistory::~LimitedHistory() {
  delete[] this->values;
}
