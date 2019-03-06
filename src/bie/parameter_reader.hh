/**
 * @file parameter_reader.hh
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
#ifndef __PARAMETER_READER_HH__
#define __PARAMETER_READER_HH__

/* -------------------------------------------------------------------------- */
// std
#include <set>
#include <map>

/* -------------------------------------------------------------------------- */
class ParameterReader {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  ParameterReader();
  virtual ~ParameterReader() {};
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// read input file
  void readInputFile(std::string file_name);

  /// write input file
  void writeInputFile(std::string file_name) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// 
  template<typename T>
  T get(std::string key) const;

  template<typename T>
  bool has(std::string key) const;
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// type of data available
  std::set<std::string> data_types;

  /// data
  std::map<std::string,std::string> string_data;
  std::map<std::string,int> int_data;
  std::map<std::string,unsigned int> uint_data;
  std::map<std::string,double> double_data;
  std::map<std::string,bool> bool_data;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "parameter_reader_inline_impl.cc"

#endif /* __PARAMETER_READER_HH__ */
