/**
 * @file nodal_field.hh
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
/* -------------------------------------------------------------------------- */
#ifndef __NODAL_FIELD_H__
#define __NODAL_FIELD_H__
/* -------------------------------------------------------------------------- */
#include <vector>
/* -------------------------------------------------------------------------- */

class NodalField {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  NodalField(unsigned int nb_nodes);
  virtual ~NodalField();

private:
  // private copy constructor: NodalField cannot be copied (for now to debug)
  NodalField(NodalField & to_copy) {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void zeros();

  void setAllValuesTo(double value);
   
  void setValuesTo(double *ptr_value, int spacing = 1);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  // get number of nodes
  unsigned int getNbNodes() const { return this->nb_nodes; };

  // access the value of node n (reading and writing)
  inline double & operator()(unsigned int node);

  // access the value of node n (only reading for const nodal field)
  inline double at(unsigned int node) const ;

  // access to storage
  inline double * storage() { return this->field; };
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  // number of nodes
  unsigned int nb_nodes;

  // nodal field
  double * field;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
inline double & NodalField::operator()(unsigned int node) {
  return this->field[node];
};

inline double NodalField::at(unsigned int node) const {
  return this->field[node];
};


#endif /* __NODAL_FIELD_H__ */
