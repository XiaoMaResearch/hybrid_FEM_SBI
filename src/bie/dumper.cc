/**
 * @file dumper.cc
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
#include "dumper.hh"

#include <sys/stat.h>
#include <iomanip>

/* -------------------------------------------------------------------------- */
Dumper::Dumper() {

  // default name and path
  this->setBaseName("standard-bname");
  this->path = ".";
  
  this->time_file = NULL;
  this->field_file = NULL;

  this->initiated = false;
}

/* -------------------------------------------------------------------------- */
Dumper::~Dumper() {

  FileToFieldMap::iterator it = this->files_and_fields.begin();
  FileToFieldMap::iterator end = this->files_and_fields.end();

  for (; it!=end; ++it) {
    it->first->close();
    delete it->first;
  }

  if (this->initiated) {
    this->time_file->close();
    delete this->time_file;

    this->coord_file->close();
    delete this->coord_file;

    this->field_file->close();
    delete this->field_file;
  }

}

/* -------------------------------------------------------------------------- */
void Dumper::setBaseName(const std::string & bname) {

  this->base_name = bname;

  this->info_file_name  = this->base_name + ".info";
  this->time_file_name  = this->base_name + ".time";
  this->coord_file_name = this->base_name + ".coord";
  this->field_file_name = this->base_name + ".fields";
  this->folder_name     = this->base_name + "-DataFiles";
}

/* -------------------------------------------------------------------------- */
void Dumper::initDump(const std::string & bname,
		      const std::string & path) {

  this->initiated = true;

  this->setBaseName(bname);
  this->path = path;

  // create folder for files (works only on linux)
  // read/write/search permission for owner and group
  // read/search permissions for others
  std::string full_path_to_folder = this->path + "/" + this->folder_name;
  mkdir(full_path_to_folder.c_str(), 
	S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  // info file
  std::string path_to_info_file = this->path + "/" + this->info_file_name;
  std::ofstream info_file;
  info_file.open(path_to_info_file, std::ofstream::out);
  info_file << "field_description " << this->field_file_name << std::endl;
  info_file << "time_description " << this->time_file_name << std::endl;
  info_file << "coord_description " << this->coord_file_name << std::endl;
  info_file << "folder_name " << this->folder_name << std::endl;
  info_file.close();

  // time file
  std::string path_to_time_file = this->path + "/" + this->time_file_name;

  this->time_file = new std::ofstream();
  this->time_file->open(path_to_time_file, std::ofstream::out);

  (*this->time_file) << std::scientific << std::setprecision(10);

  // coord file
  std::string path_to_coord_file = this->path + "/" + this->coord_file_name;

  this->coord_file = new std::ofstream();
  this->coord_file->open(path_to_coord_file, std::ofstream::out);

  (*this->coord_file) << std::scientific << std::setprecision(10);

  // field file
  std::string path_to_field_file = this->path + "/" + this->field_file_name;

  this->field_file = new std::ofstream();
  this->field_file->open(path_to_field_file, std::ofstream::out);
}

/* -------------------------------------------------------------------------- */
void Dumper::registerForDump(const std::string & field_name,
			     const NodalField * nodal_field) {

  // name and path
  std::string file_name = field_name + ".out";
  std::string path_to_file = this->path + "/" + this->folder_name + "/" + file_name;

  // open file
  std::ofstream * new_file = new std::ofstream();
  new_file->open(path_to_file, std::ofstream::out);

  // set parameters
  (*new_file) << std::scientific << std::setprecision(10);

  // keep reference to file
  this->files_and_fields[new_file] = nodal_field;

  // put info into field file
  (*this->field_file) << field_name << " " << file_name << std::endl;
}

/* -------------------------------------------------------------------------- */
// coords has all coordinates of all points
void Dumper::setCoords(unsigned int nb_dim,
		       unsigned int nb_points,
		       double * coords) {
  
  for (unsigned int i=0; i<nb_points; ++i) {
    for (unsigned int d=0; d<nb_dim; ++d) {
      if (d != 0)
	(*this->coord_file) << this->separator;
      (*this->coord_file) << coords[i*nb_dim + d];
    }
    (*this->coord_file) << std::endl;    
  }

}

/* -------------------------------------------------------------------------- */
void Dumper::dump(unsigned int step, double time) {

  FileToFieldMap::iterator it = this->files_and_fields.begin();
  FileToFieldMap::iterator end = this->files_and_fields.end();

  for (; it!=end; ++it) {
    this->dumpField(it->first, it->second);
  }

  (*this->time_file) << step << this->separator << time << std::endl;

}

/* -------------------------------------------------------------------------- */
void Dumper::dumpField(std::ofstream * dump_file, 
		       const NodalField * nodal_field) {

  unsigned int nb_values = nodal_field->getNbNodes();

  for (unsigned int i=0; i<nb_values; ++i) {
    if (i != 0)
      (*dump_file) << this->separator;
    (*dump_file) << nodal_field->at(i);
  }
  (*dump_file) << std::endl;

}
