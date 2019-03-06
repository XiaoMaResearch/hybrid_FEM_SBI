/**
 * @file material.cc
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
#include <cmath>
#include "material.hh"

Material::Material(double E, double nu, double rho) :
  E(E), nu(nu), rho(rho) {

  this->computeShearModulus();
  this->computeFirstLame();
  this->computeCp();
  this->computeCs();

}

void Material::computeShearModulus() {
  this->mu = 0.5*this->E / (1+this->nu);
}

void Material::computeFirstLame() {
  // 3D or plane-strain (would be different for plane stress)
  this->lambda = this->nu*this->E / (1+this->nu) / (1-2*this->nu);
}

void Material::computeCp() {
  this->cp = std::sqrt((this->lambda + 2*this->mu)/this->rho);
}

void Material::computeCs() {
  this->cs = std::sqrt(this->mu/this->rho);
}
