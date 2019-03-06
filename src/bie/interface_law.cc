/**
 * @file interface_law.cc
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
#include "interface_law.hh"
#include "interface.hh"

#include <iostream>

/* -------------------------------------------------------------------------- */
void InterfaceLaw::registerDumpField(const std::string & field_name) {

  // registerDumpField went through child-interface, parent-interface, child-interface-law
  // and no one knew this field_name: No more options -> do nothing
  std::cout << "Do not know dump field with name: " << field_name << std::endl;

}
