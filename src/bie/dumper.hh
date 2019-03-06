/**
 * @file dumper.hh
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
#ifndef __DUMPER_H__
#define __DUMPER_H__
/* -------------------------------------------------------------------------- */
#include <sstream>
#include <fstream>
#include <map>

#include "nodal_field.hh"

/* -------------------------------------------------------------------------- */
class Dumper {
  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */
private:
  typedef std::map<std::ofstream *, const NodalField *> FileToFieldMap;

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  Dumper();
  virtual ~Dumper();
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  virtual void initDump(const std::string & bname,
			const std::string & path);

  virtual void registerDumpField(const std::string & field_name) = 0;

  void registerForDump(const std::string & field_name,
		       const NodalField * nodal_field);

  void dump(unsigned int step, double time);

protected:
  void setBaseName(const std::string & bname);

  void setCoords(unsigned int nb_dim,
		 unsigned int nb_points,
		 double * coords);

  void dumpField(std::ofstream * dump_file,
		 const NodalField * nodal_field);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  // base name
  std::string base_name;
  
  // path to dumped files
  std::string path;

  // has dump been initiated?
  bool initiated;

  // information based on base_name
  std::string info_file_name;
  std::string time_file_name;
  std::string coord_file_name;
  std::string field_file_name;
  std::string folder_name;

  // files corresponding to field
  FileToFieldMap files_and_fields;

  // file with time stamps
  std::ofstream * time_file;

  // file with coord
  std::ofstream * coord_file;

  // file with field infos
  std::ofstream * field_file;

  // characteristics of dumper
  std::string separator = " ";

};

//#include "dumper_impl.cc"

#endif /* __DUMPER_H__ */
