/**
 * @file limited_history.hh
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
#ifndef __LIMITED_HISTORY_H__
#define __LIMITED_HISTORY_H__
/* -------------------------------------------------------------------------- */

#include <iostream>
#include <fstream>

class LimitedHistory {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  LimitedHistory(unsigned int size);
  virtual ~LimitedHistory();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  // at the current value of the history
  inline void addCurrentValue(double value);
  inline void addCurrentValue_one(double value);


  // get history value at index with index=0 : now
  inline double at(unsigned int index) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  unsigned int getSize() const { return this->size; };
  unsigned int getNbHistoryPoints() const { return std::min(this->nb_history_points,
							    this->size); };

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  // number of accumulated history points
  unsigned int nb_history_points;

  // number of history entries
  unsigned int size;
  
  // index pointing to the newest entry
  unsigned int index_now;

  // values
  double * values;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
inline void LimitedHistory::addCurrentValue(double value) {

  if (this->index_now == 0)
    this->index_now = this->size;

  this->index_now -= 1;

  this->values[this->index_now] = value;

  // increase the counter of history points
  this->nb_history_points = std::min(this->nb_history_points + 1,
				     this->size);
}
inline void LimitedHistory::addCurrentValue_one(double value) {
    
    if (this->index_now == 0)
        this->index_now = this->size;
//    
    this->index_now += 1;
//    
    this->values[this->index_now] = value;
//    
    // increase the counter of history points
    this->nb_history_points = std::min(this->nb_history_points-1,
                                       this->size);

}

/* -------------------------------------------------------------------------- */
inline double LimitedHistory::at(unsigned int index) const {
/*  if (index >= this->size) {
    std::cout << "try to access history value beyond existence" << std::endl;
    throw index;
  }
*/

  unsigned int i = (this->index_now + index) % this->size;
  return this->values[i];
}

#endif /* __LIMITED_HISTORY_H__ */
