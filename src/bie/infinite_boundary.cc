/**
 * @file infinite_boundary.cc
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
#include "infinite_boundary.hh"
//#include <omp.h>
/* -------------------------------------------------------------------------- */
InfiniteBoundary::InfiniteBoundary(double length,
				   unsigned int nb_nodes,
				   int side_factor,
				   Material * material,
				   Kernel * H11,
				   Kernel * H12,
				   Kernel * H22):
  HalfSpace(length, nb_nodes, side_factor) {

  this->setMaterial(material);
  this->setKernels(H11,H12,H22);
}

/* -------------------------------------------------------------------------- */
void InfiniteBoundary::init() {
  this->initConvolutions();  

  // might need to do this?
  /*
  // need to choose a time step to initialize 
  double delta_x = this->length / (double)(this->nb_nodes);
  double max_cs = this->getMaterial().getCs();
  double stable_time_step = delta_x / max_cs;
  this->setTimeStep(0.5*stable_time_step);
  */
}

/* -------------------------------------------------------------------------- */
void InfiniteBoundary::getDisplacement(NodalField * external_1,
				       NodalField * external_2,
				       NodalField * bc_disp_1,
				       NodalField * bc_disp_2,
                       NodalField * bc_vel_1,
                       NodalField * bc_vel_2) {

  // previous displacement forward fft
  this->forwardFFT();
  this->computeStressFourierCoeff();
  this->backwardFFT();
  //this->computeCohesion(); // not sure what to do with this?

  // with the loads from the domain
  this->computeResidual(external_1,
			external_2);

  // update displacement
  this->computeVelocity();
  this->computeDisplacement();

    
  // copy for return
  double * disp_1_p = disp_1->storage();
  double * disp_2_p = disp_2->storage();
  double * bc_disp_1_p = bc_disp_1->storage();
  double * bc_disp_2_p = bc_disp_2->storage();
  for (unsigned int n=0; n<this->nb_nodes; ++n) {
    bc_disp_1_p[n] = disp_1_p[n];
    bc_disp_2_p[n] = disp_2_p[n];
  }
  // copy for velocity
    double * vel_1_p = velo_1->storage();
    double * vel_2_p = velo_2->storage();
    double * bc_vel_1_p = bc_vel_1->storage();
    double * bc_vel_2_p = bc_vel_2->storage();
    for (unsigned int n=0; n<this->nb_nodes; ++n) {
        bc_vel_1_p[n] = vel_1_p[n];
        bc_vel_2_p[n] = vel_2_p[n];
    }

}
/*------------------------------------------------------*/
void InfiniteBoundary::getDisplacement_correction(NodalField * external_1,
                                       NodalField * external_2,
                                       NodalField * bc_disp_1,
                                       NodalField * bc_disp_2,
                                       NodalField * bc_vel_1,
                                       NodalField * bc_vel_2) {
    
    // previous displacement forward fft
    this->forwardFFT();
    this->computeStressFourierCoeff_correction();
    this->backwardFFT();
    //this->computeCohesion(); // not sure what to do with this?
    
    // with the loads from the domain
    this->computeResidual(external_1,
                          external_2);
    
    // update displacement
    this->computeVelocity();
    this->computeDisplacement();
    
    
    // copy for return
    double * disp_1_p = disp_1->storage();
    double * disp_2_p = disp_2->storage();
    double * bc_disp_1_p = bc_disp_1->storage();
    double * bc_disp_2_p = bc_disp_2->storage();
    for (unsigned int n=0; n<this->nb_nodes; ++n) {
        bc_disp_1_p[n] = disp_1_p[n];
        bc_disp_2_p[n] = disp_2_p[n];
    }
    // copy for velocity
    double * vel_1_p = velo_1->storage();
    double * vel_2_p = velo_2->storage();
    double * bc_vel_1_p = bc_vel_1->storage();
    double * bc_vel_2_p = bc_vel_2->storage();
    for (unsigned int n=0; n<this->nb_nodes; ++n) {
        bc_vel_1_p[n] = vel_1_p[n];
        bc_vel_2_p[n] = vel_2_p[n];
    }
    
}

/* -------------------------------------------------------------------------- */
void InfiniteBoundary::setDisplacement(NodalField *disp_1_inp, NodalField *disp_2_inp){
    
    
    // copy for return
    double * disp_1_p = disp_1->storage();
    double * disp_2_p = disp_2->storage();
    double * disp_1_inp_p = disp_1_inp->storage();
    double * disp_2_inp_p = disp_2_inp->storage();
    for (unsigned int n=0; n<this->nb_nodes; ++n) {
         disp_1_p[n]=disp_1_inp_p[n];
         disp_2_p[n]=disp_2_inp_p[n];
    }
    
    
}

void InfiniteBoundary::setComplexU(){
    #pragma omp parallel for
    for (unsigned int j=1; j<this->nb_fft; ++j)
        {
  //  std::complex<double> U1 = {0.0, 0.0};
  //  std::complex<double> U2 = {0.0, 0.0};
            std::complex<double> U1 = {this->disp_1->fd(j)[0], this->disp_1->fd(j)[1]};
            std::complex<double> U2 = {this->disp_2->fd(j)[0], this->disp_2->fd(j)[1]};
////
////    // store current displacement in history
//    this->U_1_r[j]->addCurrentValue(std::real(U1));
//    this->U_1_i[j]->addCurrentValue(std::imag(U1));
//    this->U_2_r[j]->addCurrentValue(std::real(U2));
//    this->U_2_i[j]->addCurrentValue(std::imag(U2));
                this->U_1_r[j]->addCurrentValue_one(std::real(U1));
                this->U_1_i[j]->addCurrentValue_one(std::imag(U1));
                this->U_2_r[j]->addCurrentValue_one(std::real(U2));
                this->U_2_i[j]->addCurrentValue_one(std::imag(U2));
        }
    
//            numberofhistory =this->U_1_i[1]->getNbHistoryPoints();
//            valueoffirst= this->U_1_i[1]->at(1);

       // }
}

void InfiniteBoundary::gethistorynum(int & numberofhistory,double &valueoffirst){
    
    numberofhistory =this->U_1_i[1]->getNbHistoryPoints();
    valueoffirst= this->U_1_i[1]->at(1);
    
    // }
}


