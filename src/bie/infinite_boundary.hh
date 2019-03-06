/**
 * @file infinite_boundary.hh
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
#ifndef __INFINITE_BOUNDARY_H__
#define __INFINITE_BOUNDARY_H__
/* -------------------------------------------------------------------------- */
#include "half_space.hh"

/* -------------------------------------------------------------------------- */
class InfiniteBoundary : public HalfSpace {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  // side factor top=1 bot=-1
  InfiniteBoundary(double length, unsigned int nb_nodes, int side_factor,
		   Material * material,
		   Kernel * H11, Kernel * H12, Kernel * H22);
  virtual ~InfiniteBoundary() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  // initialize
  void init();

  // get displacement for next step
  // external_1 = in-plane shear displacement
  // external_2 = in-plane normal displacement
  void getDisplacement(NodalField * external_1,
		       NodalField * external_2,
		       NodalField * bc_disp_1,
		       NodalField * bc_disp_2,
           NodalField * bc_vel_1,
           NodalField * bc_vel_2);
  void getDisplacement_correction(NodalField * external_1,
                         NodalField * external_2,
                         NodalField * bc_disp_1,
                         NodalField * bc_disp_2,
                         NodalField * bc_vel_1,
                         NodalField * bc_vel_2);
  void setDisplacement(NodalField *disp_1_inp, NodalField * disp_2_inp);
  void setComplexU();
  void gethistorynum(int & numohistory, double & valueoffirst);


  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:


  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

};

//#include "infinite_boundary_impl.cc"

#endif /* __INFINITE_BOUNDARY_H__ */
