/**
 * @file interface_law.hh
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
#ifndef __INTERFACE_LAW_H__
#define __INTERFACE_LAW_H__
/* -------------------------------------------------------------------------- */
#include <sstream>

#include "nodal_field.hh"

class Interface; // <--- don't know if this works --------------------------

class InterfaceLaw {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  InterfaceLaw() {};
  virtual ~InterfaceLaw() {};
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  virtual void computeCohesiveForces(NodalField * cohesion_1,
				     NodalField * cohesion_2) = 0;

  // dumper function
  virtual void registerDumpField(const std::string & field_name);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  void setInterface(Interface * interface) { this->interface = interface; };
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  Interface * interface;
};

//#include "interface_law_impl.cc"

#endif /* __INTERFACE_LAW_H__ */
