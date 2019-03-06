/**
 * @file parameter_reader.cc
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
// std
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
#include <algorithm>

#include "parameter_reader.hh"

/* -------------------------------------------------------------------------- */
ParameterReader::ParameterReader() : 
  data_types(),
  string_data(), 
  int_data(), 
  uint_data(), 
  double_data(), 
  bool_data()
{

  data_types.insert("string");
  data_types.insert("uint");
  data_types.insert("int");
  data_types.insert("double");
  data_types.insert("bool");
}

/* -------------------------------------------------------------------------- */
void ParameterReader::readInputFile(std::string file_name) {
  
  char comment_char = '#';
  char equal_char = '=';

  // open a file called file name 
  std::ifstream infile;
  infile.open(file_name.c_str());

  if(!infile.good()) {
    std::cerr << "Cannot open file " << file_name << "!!!" << std::endl;
    exit(EXIT_FAILURE);
  }

  std::string line;
  std::string clean_line;
  while(infile.good()) {
    getline(infile, line);
    clean_line = line;

    // take out comments
    size_t found_comment;
    found_comment = line.find_first_of(comment_char);
    if (found_comment != std::string::npos)
      clean_line = line.substr(0,found_comment);
    if (clean_line.empty())
      continue;

    std::stringstream sstr(clean_line);

    // check if data type exists
    std::string type;
    sstr >> type;
    std::transform(type.begin(),type.end(),type.begin(),::tolower);
    if (this->data_types.find(type) == this->data_types.end()) {
      std::cerr << " *** WARNING *** Data type " << type << " does not exist"
		<< " in this input data structure. Ignore line: ";
      std::cerr << clean_line << std::endl;
      continue;
    }
    
    std::string keyword;
    std::string equal;
    std::string value;

    // get keyword
    sstr >> keyword;
    size_t equal_p = keyword.find_first_of(equal_char);
    if (equal_p != std::string::npos) {
      equal = keyword.substr(equal_p,std::string::npos);
      keyword = keyword.substr(0,equal_p);
    }

    // get equal
    if (equal.empty())
      sstr >> equal;
    if (equal.length() != 1) {
      value = equal.substr(1,std::string::npos);
      equal = equal[0];
    }
    if (equal[0] != equal_char) {
      std::cerr << " *** WARNING *** Unrespected convention! Ignore line: ";
      std::cerr << clean_line << std::endl;
      continue;
    }
   
    // get value
    if (value.empty())
      sstr >> value;
    
    // no value
    if (value.empty()) {
      std::cerr << " *** WARNING *** No value given! Ignore line: ";
      std::cerr << clean_line << std::endl;
      continue;
    }
    
    // put value in map
    std::stringstream convert(value);
    if (type.compare("string") == 0) {
      this->string_data.insert(std::make_pair(keyword,value));
    }
    else if (type.compare("int") == 0) {
      int i;
      convert >> i;
      this->int_data.insert(std::make_pair(keyword,i));
    }
    else if (type.compare("uint") == 0) {
      unsigned int i;
      convert >> i;
      this->uint_data.insert(std::make_pair(keyword,i));
    }
    else if (type.compare("double") == 0) {
      double r;
      convert >> r;
      this->double_data.insert(std::make_pair(keyword,r));
    }
    else if (type.compare("bool") == 0) {
      std::transform(value.begin(),value.end(),value.begin(),::tolower);
      bool b;
      if (value.compare("true") == 0)
	b = true;
      else if (value.compare("false") == 0)
	b = false;
      else {
	std::cerr << " *** WARNING *** boolean cannot be " << value << ". Ignore line: ";
	std::cerr << clean_line << std::endl;
	continue;
      }
      this->bool_data.insert(std::make_pair(keyword,b));
    }
    else {
      std::cerr << " *** ERROR *** Could not add data to InputData for line: ";
      std::cerr << clean_line << std::endl;
      continue;
      exit(EXIT_FAILURE);
    }
  }
}

/* -------------------------------------------------------------------------- */
void ParameterReader::writeInputFile(std::string file_name) const {
  
  // open file to write input information
  std::ofstream outfile;
  outfile.open(file_name.c_str());
  
  // string
  for (std::map<std::string, std::string>::const_iterator it = string_data.begin();
       it != string_data.end(); ++it)
    outfile << "string " << it->first << " = " << it->second << std::endl;

  // Int
  for (std::map<std::string, int>::const_iterator it = int_data.begin();
       it != int_data.end(); ++it)
    outfile << "int " << it->first << " = " << it->second << std::endl;

  // UInt
  for (std::map<std::string, unsigned int>::const_iterator it = uint_data.begin();
       it != uint_data.end(); ++it)
    outfile << "uint " << it->first << " = " << it->second << std::endl;

  // Double
  for (std::map<std::string, double>::const_iterator it = double_data.begin();
       it != double_data.end(); ++it)
    outfile << "double " << it->first << " = " << it->second << std::endl;
 
  // Bool
  for (std::map<std::string, bool>::const_iterator it = bool_data.begin();
       it != bool_data.end(); ++it) {
    std::string b = "false";
    if (it->second)
      b = "true";
    outfile << "bool " << it->first << " = " << b << std::endl;
  }
}

/* -------------------------------------------------------------------------- */
template<>
unsigned int ParameterReader::get<unsigned int>(std::string key) const {
  std::map<std::string,unsigned int>::const_iterator it;
  it = this->uint_data.find(key);
  
  // if not in map
  if (it == this->uint_data.end()) {
    std::cerr << " *** ERROR *** This data was not in input file. " 
	      << "You need the following line in your input file: ";
    std::cerr << "uint " << key << " = ???" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  else 
    return it->second;
}

/* -------------------------------------------------------------------------- */
template<>
std::string ParameterReader::get<std::string>(std::string key) const {
  std::map<std::string,std::string>::const_iterator it;
  it = this->string_data.find(key);
  
  // if not in map
  if (it == this->string_data.end()) {
    std::cerr << " *** ERROR *** This data was not in input file. " 
	      << "You need the following line in your input file: ";
    std::cerr << "string " << key << " = ???" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  else 
    return it->second;
}

/* -------------------------------------------------------------------------- */
template<>
int ParameterReader::get<int>(std::string key) const {
  std::map<std::string,int>::const_iterator it;
  it = this->int_data.find(key);
  
  // if not in map
  if (it == this->int_data.end()) {
    std::cerr << " *** ERROR *** This data was not in input file. " 
	      << "You need the following line in your input file: ";
    std::cerr << "int " << key << " = ???" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  else 
    return it->second;
}

/* -------------------------------------------------------------------------- */
template<>
double ParameterReader::get<double>(std::string key) const {
  std::map<std::string,double>::const_iterator it;
  it = this->double_data.find(key);
  
  // if not in map
  if (it == this->double_data.end()) {
    std::cerr << " *** ERROR *** This data was not in input file. " 
	      << "You need the following line in your input file: ";
    std::cerr << "double " << key << " = ???" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  else 
    return it->second;
}

/* -------------------------------------------------------------------------- */
template<>
bool ParameterReader::get<bool>(std::string key) const {
  std::map<std::string,bool>::const_iterator it;
  it = this->bool_data.find(key);
  
  // if not in map
  if (it == this->bool_data.end()) {
    std::cerr << " *** ERROR *** This data was not in input file. " 
	      << "You need the following line in your input file: ";
    std::cerr << "bool " << key << " = ???" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  else 
    return it->second;
}

/* -------------------------------------------------------------------------- */
template<>
bool ParameterReader::has<bool>(std::string key) const {
  std::map<std::string,bool>::const_iterator it;
  it = this->bool_data.find(key);
  return (it != this->bool_data.end());
}
template<>
bool ParameterReader::has<std::string>(std::string key) const {
  std::map<std::string,std::string>::const_iterator it;
  it = this->string_data.find(key);
  return (it != this->string_data.end());
}
template<>
bool ParameterReader::has<int>(std::string key) const {
  std::map<std::string,int>::const_iterator it;
  it = this->int_data.find(key);
  return (it != this->int_data.end());
}
template<>
bool ParameterReader::has<unsigned int>(std::string key) const {
  std::map<std::string,unsigned int>::const_iterator it;
  it = this->uint_data.find(key);
  return (it != this->uint_data.end());
}
template<>
bool ParameterReader::has<double>(std::string key) const {
  std::map<std::string,double>::const_iterator it;
  it = this->double_data.find(key);
  return (it != this->double_data.end());
}

