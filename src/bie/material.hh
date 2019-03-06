/**
 * @file material.hh
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
#ifndef __MATERIAL_H__
#define __MATERIAL_H__
/* -------------------------------------------------------------------------- */


/* -------------------------------------------------------------------------- */
class Material {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  Material(double E, double nu, double rho);
  virtual ~Material() {};
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
private:
  void computeShearModulus();
  void computeFirstLame();
  void computeCp();
  void computeCs();

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  double getYoungsModulus() const { return this->E; };
  double getShearModulus() const { return this->mu; };
  double getPoissonRatio() const { return this->nu; };
  double getDensity() const { return this->rho; };
  double getCp() const { return this->cp; };
  double getCs() const { return this->cs; };

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  // PRIMARY
  // young's modulus
  double E;
  // poisson's ratio
  double nu;
  // density
  double rho;

  // SECONDARY
  // first lame
  double lambda;
  // shear modulus
  double mu;
  // p-wave speed
  double cp;
  // s-wave speed
  double cs;
};

//#include "material_impl.cc"

#endif /* __MATERIAL_H__ */
