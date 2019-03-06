#include <iostream>
#include <math.h>
// #include <omp.h>
// #include <vector>
// #include <complex>
// #include <algorithm>
// #include <functional>
// #include <numeric>
// #include <ctime>
//include <omp.h>
#include "Eigen/Eigen"
#include <fstream>
//include the bie header files
#include "material.hh"
#include "precomputed_kernel.hh"
#include "bimat_interface.hh"
#include "infinite_boundary.hh"
//include the fem header files
#include "Meshrectangular.hpp"
#include "mesh_Generated.hpp"
#include "bcdof.hpp"
#include "cal_BTDB.hpp"
#include "maplocal.hpp"
#include "mapglobal.hpp"
#include "BIE_correct.hpp"
#include "update_mesh_topology.hpp"

using namespace Eigen;
using namespace std;

int main() {
    int startTime = clock();
    // Domain Size
    double x_min = -5e3;
    double x_max = 5e3;
    double y_min = -1.0e3;
    double y_max = 1.0e3;
    int dim = 2.0;
    double dx = 100;
    double dy = 100;
    double nx = (x_max-x_min)/dx;
    double ny = (y_max-y_min)/dy;
    MatrixXd Node_top = MatrixXd::Zero((nx+1)*(ny+1),2);
    MatrixXd Element_top = MatrixXd::Zero(nx*ny,4);
    MatrixXd Node_bot = MatrixXd::Zero((nx+1)*(ny+1),2);
    MatrixXd Element_bot = MatrixXd::Zero(nx*ny,4);
    VectorXd fault_surf_nodes = VectorXd::Zero((nx+1),1);
    VectorXd fault_surf_nodes_new = VectorXd::Zero((nx+1),1);
    mesh_Generated(x_min,x_max,y_min,y_max,dx,dy,nx,ny,Node_top, Element_top, fault_surf_nodes, fault_surf_nodes_new);
    cout<<fault_surf_nodes<<endl;
    mesh_Generated(x_min,x_max,y_min,y_max,dx,dy,nx,ny,Node_bot, Element_bot, fault_surf_nodes, fault_surf_nodes_new);
    VectorXd x = VectorXd::LinSpaced(nx+1,x_min,x_max);
    VectorXd y = VectorXd::LinSpaced(ny+1,y_min,y_max);
    int n_nodes = Node_top.rows();
    int nel = Element_top.rows();
    // Material Property
    double density = 2670.0;
    double v_s =3.464e3;
    double v_p = 6.0e3;
    double G= pow(v_s,2)*density;
    double Lambda = pow(v_p,2)*density-2.0*G;
    double E  = G*(3.0*Lambda+2.0*G)/(Lambda+G);
    double nu = Lambda/(2.0*(Lambda+G));
    // Time step
    double alpha = 0.4;
    double dt = alpha*dx/v_p;
    // Reyleigh Damping
    double beta =0.1;
    double q = beta*dt;
    // number of time steps
    double time_run = 6.0;
    int numt = time_run/dt;
    //int numt =1;
    VectorXd time_total = dt*VectorXd::LinSpaced(numt,1,numt);
    // Number of dof per nodes
    int Ndofn = 2;
    // Number of node per elements
    int Nnel = 4;
    // Element mass
    double M=density*dx*dy*1.0;
    // External force vector
    VectorXd F_ext_global = VectorXd::Zero(n_nodes*Ndofn,1);
    // solution vector
    // + half plane
    VectorXd u_n_pos = VectorXd::Zero(n_nodes*dim,1);
    VectorXd v_n_pos = VectorXd::Zero(n_nodes*dim,1);
    VectorXd u_new_pos = VectorXd::Zero(n_nodes*dim,1);
    VectorXd v_new_pos = VectorXd::Zero(n_nodes*dim,1);
    VectorXd a_n_pos = VectorXd::Zero(n_nodes*dim,1);
    // - half plane
    VectorXd u_n_neg = VectorXd::Zero(n_nodes*dim,1);
    VectorXd v_n_neg = VectorXd::Zero(n_nodes*dim,1);
    VectorXd u_new_neg = VectorXd::Zero(n_nodes*dim,1);
    VectorXd v_new_neg = VectorXd::Zero(n_nodes*dim,1);
    VectorXd a_n_neg = VectorXd::Zero(n_nodes*dim,1);
    // slip and slip rate on the fault
    VectorXd delt_u_n_x = VectorXd::Zero(nx+1,1);
    VectorXd delt_v_n_x = VectorXd::Zero(nx+1,1);
    VectorXd delt_u_n_y = VectorXd::Zero(nx+1,1);
    VectorXd delt_v_n_y = VectorXd::Zero(nx+1,1);
    // Intial shear traction on the fault
    VectorXd Tx_0 = VectorXd::Zero(nx+1,1);
    for (int i=0 ; i<=nx; i++)
    {
        if ((x(i)<=(x_max+x_min)/2+1.5e3)&&(x(i)>=(x_max+x_min)/2-1.5e3))
        {
            Tx_0(i) = 81.6e6;
        }
        else if ((x(i)<=(x_max+x_min)/2+7.5e3+1.5e3)&&(x(i)>=(x_max+x_min)/2+7.5e3-1.5e3))
        {
            Tx_0(i) = 62.0e6;
        }
        else if ((x(i)<=(x_max+x_min)/2-7.5e3+1.5e3)&&(x(i)>=(x_max+x_min)/2-7.5e3-1.5e3))
        {
            Tx_0(i) = 78.0e6;
        }
        else
        {
            Tx_0(i) = 70.0e6;
            
        }
    }
    VectorXd sigmabar = 120.e6*VectorXd::Ones(nx+1,1);
    VectorXd Ty_0= -sigmabar;
    VectorXd T_cx = VectorXd::Zero(nx+1,1);
    VectorXd T_cy = VectorXd::Zero(nx+1,1);
    // Slip Weakening friction parameters
    double Dc = 0.4;
    double mu_d=0.525;
    VectorXd mu_s = VectorXd::Zero(nx+1,1);
    for (int i=0; i<=nx;i++)
    {
        if ((x(i)<=(x_max+x_min)/2+15e3)&&(x(i)>=(x_max+x_min)/2-15e3))
        {
            mu_s(i) = 0.677;
        }
        else
        {
            mu_s(i) = 10000.0;
        }
    }
    VectorXd tau_s = VectorXd::Zero(nx+1,1);
    VectorXd M_el_vec = M/4*VectorXd::Ones(Nnel*Ndofn,1);
    // Get the index degree of freedome for each element
    VectorXd index_el_pos = VectorXd::Zero(Ndofn*Nnel,1);
    VectorXd index_el_neg = VectorXd::Zero(Ndofn*Nnel,1);
    MatrixXd index_store_pos = MatrixXd::Zero(Ndofn*Nnel,nel);
    MatrixXd index_store_neg = MatrixXd::Zero(Ndofn*Nnel,nel);
    for (int i=0;i<nel;i++)
    {
        bcdof(Element_top.row(i),dim,index_el_pos);
        bcdof(Element_bot.row(i),dim,index_el_neg);
        index_store_pos.col(i) = index_el_pos;
        index_store_neg.col(i) = index_el_neg;
    }
    printf("finisheindexing\n");
    VectorXd top_surf_index = VectorXd::Zero(2*(nx+1),1);
    VectorXd bot_surf_index = VectorXd::Zero(2*(nx+1),1);
    top_surf_index.head(nx+1) = VectorXd::LinSpaced(nx+1,0,nx*2);
    top_surf_index.tail(nx+1)=VectorXd::LinSpaced(nx+1,1,nx*2+1);
    bot_surf_index.head(nx+1) = VectorXd::LinSpaced(nx+1,2*ny*(nx+1),2*(ny+1)*(nx+1)-2);
    bot_surf_index.tail(nx+1) = VectorXd::LinSpaced(nx+1,2*ny*(nx+1)+1,2*(ny+1)*(nx+1)-1);
    // Adding the coupling layer index
    VectorXd BIE_top_surf_index = VectorXd::Zero(2*(nx+1),1);
    VectorXd BIE_bot_surf_index = VectorXd::Zero(2*(nx+1),1);
    BIE_top_surf_index.head(nx+1) = VectorXd::LinSpaced(nx+1,2*ny*(nx+1),2*(ny+1)*(nx+1)-2);
    BIE_top_surf_index.tail(nx+1)=VectorXd::LinSpaced(nx+1,2*ny*(nx+1)+1,2*(ny+1)*(nx+1)-1);
    BIE_bot_surf_index.head(nx+1) = VectorXd::LinSpaced(nx+1,0,nx*2);
    BIE_bot_surf_index.tail(nx+1)=VectorXd::LinSpaced(nx+1,1,nx*2+1);
    // The force vector that gonna be provided to the BIE code
    VectorXd BIE_top_surf_force = VectorXd::Zero(2*(nx+1),1);
    VectorXd BIE_bot_surf_force = VectorXd::Zero(2*(nx+1),1);
    VectorXd BIE_top_surf_force_x = VectorXd::Zero(nx+1,1);
    VectorXd BIE_top_surf_force_y = VectorXd::Zero(nx+1,1);
    VectorXd BIE_bot_surf_force_x = VectorXd::Zero(nx+1,1);
    VectorXd BIE_bot_surf_force_y = VectorXd::Zero(nx+1,1);
    // Setting up the material property for the BIE code
    Material BIE_top_mat = Material(E,nu,density);
    Material BIE_bot_mat = Material(E,nu,density);
    double length = x_max-x_min;
    // infinte bc BIE call infinite_boundary.cc
    PrecomputedKernel h11("kernels/nu_.25_h11.dat");
    PrecomputedKernel h12("kernels/nu_.25_k12.dat");
    PrecomputedKernel h22("kernels/nu_.25_h22.dat");
    InfiniteBoundary BIE_inf_top(length,nx+1,1.0,&BIE_top_mat,&h11,&h12,&h22);
    InfiniteBoundary BIE_inf_bot(length,nx+1,-1.0,&BIE_bot_mat,&h11,&h12,&h22);
    // BIE setting time step
    BIE_inf_top.setTimeStep(dt);
    BIE_inf_bot.setTimeStep(dt);
    // BIE initialization
    BIE_inf_top.init();
    BIE_inf_bot.init();
    // Calculate the global Mass matrix
    VectorXd M_global=VectorXd::Zero(n_nodes*Ndofn,1);
    for (int i=0 ; i<nel;i++)
    {
        index_el_pos = index_store_pos.col(i);
        mapglobal(index_el_pos,M_global,M_el_vec);
    }
    // Calculate the stiffness matrix
    MatrixXd coord = MatrixXd::Zero(4,2);
    VectorXd Element_0= Element_top.row(0);
    coord.row(0) = Node_top.row(Element_0(0));
    coord.row(1) = Node_top.row(Element_0(1));
    coord.row(2) = Node_top.row(Element_0(2));
    coord.row(3) = Node_top.row(Element_0(3));
    MatrixXd B_T_D_B = MatrixXd::Zero(32,8);
    cal_BTDB (coord,E,nu,B_T_D_B);
    MatrixXd ke = MatrixXd::Zero(8,8);
    int k = 0;
    for (int i =0;i<4;i++)
    {
        ke = ke+ B_T_D_B.block(k,0,8,8);
        k=k+8;
    }
    printf("ready to start\n");
    // Output
    ofstream slip("slip.bin",ios::binary);
    slip.close();
    ofstream slip_rate("slip_rate.bin",ios::binary);
    slip_rate.close();
    ofstream shear("shear.bin",ios::binary);
    shear.close();
//    std::ofstream slip_x("results/slip.txt");
//    std::ofstream rate_x("results/slip_rate.txt");
//    std::ofstream shear_x("results/shear.txt");
    // Main time loop
    for (int j=0;j<numt;j++)
    {
        VectorXd fe_global_pos= VectorXd::Zero(n_nodes*Ndofn,1);
        VectorXd fe_global_neg= VectorXd::Zero(n_nodes*Ndofn,1);
#pragma omp parallel
        {
#pragma omp for
            for (int i=0;i<nel;i++)
            {
                VectorXd u_n_local_pos=VectorXd::Zero(8,1) ;
                VectorXd v_n_local_pos=VectorXd::Zero(8,1) ;
                VectorXd u_n_local_neg=VectorXd::Zero(8,1) ;
                VectorXd v_n_local_neg=VectorXd::Zero(8,1) ;
                VectorXd index_pos = index_store_pos.col(i);
                VectorXd index_neg = index_store_neg.col(i);
                maplocal(index_pos,u_n_pos,u_n_local_pos);
                maplocal(index_pos,v_n_pos,v_n_local_pos);
                maplocal(index_neg,u_n_neg,u_n_local_neg);
                maplocal(index_neg,v_n_neg,v_n_local_neg);
                VectorXd fe_int_pos = ke*(u_n_local_pos+q*v_n_local_pos);
                VectorXd fe_int_neg = ke*(u_n_local_neg+q*v_n_local_neg);
                mapglobal(index_pos,fe_global_pos,fe_int_pos);
                mapglobal(index_neg,fe_global_neg,fe_int_neg);
                
            }
        }
        //maplocal(BIE_top_surf_index,fe_global_pos,BIE_top_surf_force);
        // maplocal(BIE_bot_surf_index,fe_global_neg,BIE_bot_surf_force);
        
        VectorXd  BIE_top_surf_2nd_index = BIE_top_surf_index-2*(nx+1)*VectorXd::Ones(2*(nx+1), 1);
        VectorXd  BIE_bot_surf_2nd_index = BIE_bot_surf_index+2*(nx+1)*VectorXd::Ones(2*(nx+1), 1);
        VectorXd  BIE_top_surf_3rd_index = BIE_top_surf_index-2*2*(nx+1)*VectorXd::Ones(2*(nx+1), 1);
        VectorXd  BIE_bot_surf_3rd_index = BIE_bot_surf_index+2*2*(nx+1)*VectorXd::Ones(2*(nx+1), 1);
        
        VectorXd u_n_top=VectorXd::Zero(2*(nx+1),1);
        VectorXd u_n_bot=VectorXd::Zero(2*(nx+1),1);
        maplocal(BIE_top_surf_index,u_n_pos,u_n_top);
        maplocal(BIE_bot_surf_index,u_n_neg,u_n_bot);
        VectorXd u_n_top_2nd=VectorXd::Zero(2*(nx+1),1);
        VectorXd u_n_bot_2nd=VectorXd::Zero(2*(nx+1),1);
        maplocal(BIE_top_surf_2nd_index,u_n_pos,u_n_top_2nd);
        maplocal(BIE_bot_surf_2nd_index,u_n_neg,u_n_bot_2nd);
        VectorXd u_n_top_3rd=VectorXd::Zero(2*(nx+1),1);
        VectorXd u_n_bot_3rd=VectorXd::Zero(2*(nx+1),1);
        maplocal(BIE_top_surf_3rd_index,u_n_pos,u_n_top_3rd);
        maplocal(BIE_bot_surf_3rd_index,u_n_neg,u_n_bot_3rd);
        
        VectorXd u_n_top_x = u_n_top.head(nx+1);
        VectorXd u_n_top_y = u_n_top.tail(nx+1);
        VectorXd u_n_bot_x = u_n_bot.head(nx+1);
        VectorXd u_n_bot_y = u_n_bot.tail(nx+1);
        
        VectorXd u_n_top_2nd_x = u_n_top_2nd.head(nx+1);
        VectorXd u_n_top_2nd_y = u_n_top_2nd.tail(nx+1);
        VectorXd u_n_bot_2nd_x = u_n_bot_2nd.head(nx+1);
        VectorXd u_n_bot_2nd_y = u_n_bot_2nd.tail(nx+1);
        
        VectorXd u_n_top_3rd_x = u_n_top_3rd.head(nx+1);
        VectorXd u_n_top_3rd_y = u_n_top_3rd.tail(nx+1);
        VectorXd u_n_bot_3rd_x = u_n_bot_3rd.head(nx+1);
        VectorXd u_n_bot_3rd_y = u_n_bot_3rd.tail(nx+1);
        
        
        // top
        VectorXd u_n_top_x_leftshit = VectorXd::Zero((nx+1),1);
        u_n_top_x_leftshit.tail(nx)=u_n_top_x.head(nx);
        VectorXd u_n_top_x_leftshit_2 = VectorXd::Zero((nx+1),1);
        u_n_top_x_leftshit_2.tail(nx-1) = u_n_top_x.head(nx-1);
        // d u_x /d x
        // e_xx = d u_x / d x
        // VectorXd e_n_top_xx = 1.0/dx/2.0*(u_n_top_x_leftshit-u_n_top_x_rightshit);
        //
        // VectorXd e_n_top_xx = 1.0/dx/2.0*(3.0*u_n_top_x-4.0*u_n_top_x_leftshit+u_n_top_x_leftshit_2);
        VectorXd e_n_top_xx = 1.0/dx*(u_n_top_x-u_n_top_x_leftshit);
        // e_yy = d u_y / d y
        VectorXd e_n_top_yy = 1.0/dy/2.0*(3.0*u_n_top_y-4.0*u_n_top_2nd_y+u_n_top_3rd_y);
        // e_xy = 0.5*(d u_x /d y + d u_y / d x)
        // d u_x / d y
        VectorXd duxdy_top = 1.0/dy/2.0*(3.0*u_n_top_x-4.0*u_n_top_2nd_x+u_n_top_3rd_x);
        // d u_y / d x
        
        VectorXd u_n_top_y_leftshit = VectorXd::Zero((nx+1),1);
        u_n_top_y_leftshit.tail(nx)=u_n_top_y.head(nx);
        VectorXd u_n_top_y_leftshit_2 = VectorXd::Zero((nx+1),1);
        u_n_top_y_leftshit_2.tail(nx-1) = u_n_top_y.head(nx-1);
        
        
        //VectorXd duydx_top = 1.0/dx/2.0*(u_n_top_y_rightshit-u_n_top_y_leftshit);
        //  VectorXd duydx_top = 1.0/dx/2.0*(3.0*u_n_top_y-4.0*u_n_top_y_leftshit+u_n_top_y_leftshit_2);
        
        VectorXd duydx_top = 1.0/dx*(u_n_top_y-u_n_top_y_leftshit);
        
        // e_xy
        VectorXd e_n_top_xy = 0.5*(duxdy_top+duydx_top);
        
        // Compute stress
        // sigma_yy = Lambda*e_xx+ (2G+Lambda)*(e_yy)
        VectorXd sigma_top_yy = Lambda*(e_n_top_xx)+(2*G+Lambda)*(e_n_top_yy);
        // sigma_xy = 2G*e_xy
        VectorXd sigma_top_xy = 2*G*e_n_top_xy;
        
        // bot
        VectorXd u_n_bot_x_leftshit = VectorXd::Zero((nx+1),1);
        u_n_bot_x_leftshit.segment(1,nx)=u_n_bot_x.head(nx);
        VectorXd u_n_bot_x_rightshit = VectorXd::Zero((nx+1),1);
        // u_n_bot_x_rightshit.segment(0,nx)=u_n_bot_x.tail(nx);
        u_n_bot_x_rightshit.head(nx) = u_n_bot_x.tail(nx);
        VectorXd u_n_bot_x_rightshit_2 = VectorXd::Zero((nx+1),1);
        //  u_n_bot_x_rightshit.segment(0,nx)=u_n_bot_x.tail(nx);
        u_n_bot_x_rightshit_2.head(nx-1)=u_n_bot_x.tail(nx-1);
        
        // e_xx = d u_x / d x
        //  VectorXd e_n_bot_xx = 1.0/dx/2.0*(u_n_bot_x_leftshit-u_n_bot_x_rightshit);
        
        VectorXd e_n_bot_xx = 1.0/dx*(u_n_bot_x-u_n_bot_x_rightshit);
        // VectorXd e_n_bot_xx = 1.0/dx/2.0*(3.0*u_n_bot_x-4.0*u_n_bot_x_rightshit+u_n_bot_x_rightshit_2);
        
        // e_yy = d u_y / d y
        VectorXd e_n_bot_yy = 1.0/dy/2.0*(3.0*u_n_bot_y-4.0*u_n_bot_2nd_y+u_n_bot_3rd_y);
        // e_xy = 0.5*(d u_x /d y + d u_y / d x)
        // d u_x / d y
        VectorXd duxdy_bot = 1.0/dy/2.0*(3.0*u_n_bot_x-4.0*u_n_bot_2nd_x+u_n_bot_3rd_x);
        // d u_y / d x
        
        VectorXd u_n_bot_y_leftshit = VectorXd::Zero((nx+1),1);
        u_n_bot_y_leftshit.segment(1,nx)=u_n_bot_y.head(nx);
        VectorXd u_n_bot_y_rightshit = VectorXd::Zero((nx+1),1);
        //   u_n_bot_y_rightshit.segment(0,nx)=u_n_bot_y.tail(nx);
        u_n_bot_y_rightshit.head(nx)=u_n_bot_y.tail(nx);
        VectorXd u_n_bot_y_rightshit_2 = VectorXd::Zero((nx+1),1);
        u_n_bot_y_rightshit_2.head(nx-1)=u_n_bot_y.tail(nx-1);
        
        //  VectorXd duydx_bot = 1.0/dx/2.0*(u_n_bot_y_leftshit-u_n_bot_y_rightshit);
        
        VectorXd duydx_bot = 1.0/dx*(u_n_bot_y-u_n_bot_y_rightshit);
        // VectorXd duydx_bot = 1.0/dx/2.0*(3.0*u_n_bot_y-4.0*u_n_bot_y_rightshit+u_n_bot_y_rightshit_2);
        
        // e_xy
        VectorXd e_n_bot_xy = 0.5*(duxdy_bot+duydx_bot);
        
        // Compute stress
        // sigma_yy = Lambda*e_xx+ (2G+Lambda)*(e_yy)
        VectorXd sigma_bot_yy = Lambda*(e_n_bot_xx)+(2*G+Lambda)*(e_n_bot_yy);
        // sigma_xy = 2G*e_xy
        VectorXd sigma_bot_xy = 2*G*e_n_bot_xy;
        
        
        
        
        BIE_top_surf_force_x = -sigma_top_xy ;
        BIE_top_surf_force_y = -sigma_top_yy ;
        BIE_bot_surf_force_x = sigma_bot_xy ;
        BIE_bot_surf_force_y = sigma_bot_yy ;
        
        
        //   BIE_top_surf_force_x = 1.0/dx*(-BIE_top_surf_force.head(nx+1));
        //  BIE_top_surf_force_y = 1.0/dx*(-BIE_top_surf_force.tail(nx+1));
        // Sign need to flip??  for the bottom layer
        // BIE_bot_surf_force_x = 1.0/dx*(BIE_bot_surf_force.head(nx+1));
        //BIE_bot_surf_force_y = 1.0/dx*(BIE_bot_surf_force.tail(nx+1));
        /* -------------------------------------------------------------------------- */
        // top -x
        vector<double> std_BIE_top_surf_force_x;
        std_BIE_top_surf_force_x.resize(BIE_top_surf_force_x.size());
        VectorXd::Map(&std_BIE_top_surf_force_x[0], BIE_top_surf_force_x.size()) = BIE_top_surf_force_x;
        // bot -x
        vector<double> std_BIE_bot_surf_force_x;
        std_BIE_bot_surf_force_x.resize(BIE_bot_surf_force_x.size());
        VectorXd::Map(&std_BIE_bot_surf_force_x[0], BIE_bot_surf_force_x.size()) = BIE_bot_surf_force_x;
        // top -y
        vector<double> std_BIE_top_surf_force_y;
        std_BIE_top_surf_force_y.resize(BIE_top_surf_force_y.size());
        VectorXd::Map(&std_BIE_top_surf_force_y[0], BIE_top_surf_force_y.size()) = BIE_top_surf_force_y;
        // bot -y
        vector<double> std_BIE_bot_surf_force_y;
        std_BIE_bot_surf_force_y.resize(BIE_bot_surf_force_y.size());
        VectorXd::Map(&std_BIE_bot_surf_force_y[0], BIE_bot_surf_force_y.size()) = BIE_bot_surf_force_y;
        
        // top BIE bc
        double *x_top_ptr;
        x_top_ptr = &std_BIE_top_surf_force_x[0];
        double *y_top_ptr;
        y_top_ptr = &std_BIE_top_surf_force_y[0];
        
        NodalField *BIE_top_force_x_ptr = new NodalField(nx+1);
        BIE_top_force_x_ptr->setValuesTo(x_top_ptr);
        NodalField *BIE_top_force_y_ptr = new NodalField(nx+1);
        BIE_top_force_y_ptr->setValuesTo(y_top_ptr);
        
        NodalField *BIE_top_disp_x_ptr = new NodalField(nx+1);
        NodalField *BIE_top_disp_y_ptr = new NodalField(nx+1);
        NodalField *BIE_top_vel_x_ptr = new NodalField(nx+1);
        NodalField *BIE_top_vel_y_ptr = new NodalField(nx+1);
        
        BIE_inf_top.getDisplacement(BIE_top_force_x_ptr,BIE_top_force_y_ptr,BIE_top_disp_x_ptr,BIE_top_disp_y_ptr,BIE_top_vel_x_ptr,BIE_top_vel_y_ptr);
        
        // bot BIE bc
        double *x_bot_ptr;
        x_bot_ptr = &std_BIE_bot_surf_force_x[0];
        double *y_bot_ptr;
        y_bot_ptr = &std_BIE_bot_surf_force_y[0];
        
        NodalField *BIE_bot_force_x_ptr = new NodalField(nx+1);
        BIE_bot_force_x_ptr->setValuesTo(x_bot_ptr);
        NodalField *BIE_bot_force_y_ptr = new NodalField(nx+1);
        BIE_bot_force_y_ptr->setValuesTo(y_bot_ptr);
        
        NodalField *BIE_bot_disp_x_ptr = new NodalField(nx+1);
        NodalField *BIE_bot_disp_y_ptr = new NodalField(nx+1);
        NodalField *BIE_bot_vel_x_ptr = new NodalField(nx+1);
        NodalField *BIE_bot_vel_y_ptr = new NodalField(nx+1);
        
        BIE_inf_bot.getDisplacement(BIE_bot_force_x_ptr,BIE_bot_force_y_ptr,BIE_bot_disp_x_ptr,BIE_bot_disp_y_ptr,BIE_bot_vel_x_ptr,BIE_bot_vel_y_ptr);
        
        // get the displacement from BIE to FEM
        // top -x
        double *FEM_top_disp_x_ptr= BIE_top_disp_x_ptr->storage();
        Eigen::Map<Eigen::VectorXd> FEM_top_disp_x(FEM_top_disp_x_ptr, nx+1);
        // top -y
        double *FEM_top_disp_y_ptr= BIE_top_disp_y_ptr->storage();
        Eigen::Map<Eigen::VectorXd> FEM_top_disp_y(FEM_top_disp_y_ptr, nx+1);
        // bot -x
        double *FEM_bot_disp_x_ptr= BIE_bot_disp_x_ptr->storage();
        Eigen::Map<Eigen::VectorXd> FEM_bot_disp_x(FEM_bot_disp_x_ptr, nx+1);
        // bot -y
        double *FEM_bot_disp_y_ptr= BIE_bot_disp_y_ptr->storage();
        Eigen::Map<Eigen::VectorXd> FEM_bot_disp_y(FEM_bot_disp_y_ptr, nx+1);
        
        // get the velocity from BIE to FEM
        // top -x
        double *FEM_top_vel_x_ptr= BIE_top_vel_x_ptr->storage();
        Eigen::Map<Eigen::VectorXd> FEM_top_vel_x(FEM_top_vel_x_ptr, nx+1);
        // top -y
        double *FEM_top_vel_y_ptr= BIE_top_vel_y_ptr->storage();
        Eigen::Map<Eigen::VectorXd> FEM_top_vel_y(FEM_top_vel_y_ptr, nx+1);
        // bot -x
        double *FEM_bot_vel_x_ptr= BIE_bot_vel_x_ptr->storage();
        Eigen::Map<Eigen::VectorXd> FEM_bot_vel_x(FEM_bot_vel_x_ptr, nx+1);
        // bot -y
        double *FEM_bot_vel_y_ptr= BIE_bot_vel_y_ptr->storage();
        Eigen::Map<Eigen::VectorXd> FEM_bot_vel_y(FEM_bot_vel_y_ptr, nx+1);
        
        
        
        
        
        
        
        
        
        VectorXd M_pos = VectorXd::Zero(nx+1,1);
        M_pos(0) = M/4.0;
        M_pos.segment(1,nx)=M/2*VectorXd::Ones(nx,1);
        M_pos(nx)= M/4.0;
        VectorXd M_neg=M_pos;
        
        VectorXd fe_int_fault_pos = VectorXd::Zero(nx+1,1);
        VectorXd fe_int_fault_neg = VectorXd::Zero(nx+1,1);
        
        maplocal(top_surf_index,fe_global_pos,fe_int_fault_pos);
        maplocal(bot_surf_index,fe_global_neg,fe_int_fault_neg);
        
        
        VectorXd fe_int_pos_x=-fe_int_fault_pos.head(nx+1);
        VectorXd fe_int_neg_x = -fe_int_fault_neg.head(nx+1);
        VectorXd fe_int_pos_y = fe_int_fault_pos.tail(nx+1);
        VectorXd fe_int_neg_y = -fe_int_fault_neg.tail(nx+1);
        //   cout<<top_surf_index<<"***********\n"<<bot_surf_index;
        double a = dx*1.0;
        VectorXd Tx = (1.0/dt*M_neg.cwiseProduct(M_pos).cwiseProduct(delt_v_n_x)+\
                       (M_neg.cwiseProduct(fe_int_pos_x)-M_pos.cwiseProduct(fe_int_neg_x))).cwiseQuotient(a*(M_neg+M_pos))+Tx_0;
        VectorXd Ty = (1.0/dt*M_neg.cwiseProduct(M_pos).cwiseProduct(delt_v_n_y+1.0/dt*(delt_u_n_y))+(M_neg.cwiseProduct(fe_int_pos_y)-M_pos.cwiseProduct(fe_int_neg_y))).cwiseQuotient(a*(M_neg+M_pos))+Ty_0;
        
        
        for (int k =0;k<nx+1;k++)
        {
            if (Ty(k)<=0.0)
            {
                T_cy(k) = Ty(k);
            }
            else
            {
                T_cy(k) = 0.0;
            }
        }
        for (int k =0;k<nx+1;k++)
        {
            if (delt_u_n_x(k)<Dc)
            {
                tau_s(k) = (mu_s(k)-(mu_s(k)-mu_d)*delt_u_n_x(k)/Dc)*(-T_cy(k));
            }
            else
            {
                tau_s(k) = mu_d*(-T_cy(k));
            }
        }
        for (int k =0;k<nx+1;k++)
        {
            if (Tx(k)<=tau_s(k))
            {
                T_cx(k) = Tx(k);
            }
            else
            {
                T_cx(k) = tau_s(k);
            }
        }
        VectorXd F_total_pos = F_ext_global-fe_global_pos;
        mapglobal(top_surf_index.head(nx+1),F_total_pos,-a*(T_cx-Tx_0));
        mapglobal(top_surf_index.tail(nx+1),F_total_pos,-a*(T_cy-Ty_0));
        
        VectorXd F_total_neg = F_ext_global-fe_global_neg;
        mapglobal(bot_surf_index.head(nx+1),F_total_neg,a*(T_cx-Tx_0));
        mapglobal(bot_surf_index.tail(nx+1),F_total_neg,a*(T_cy-Ty_0));
        // Central Difference Time integration
        a_n_pos = F_total_pos.cwiseQuotient(M_global);
        v_new_pos = v_n_pos+dt*a_n_pos;
        u_new_pos = u_n_pos+dt*v_new_pos;
        a_n_neg = F_total_neg.cwiseQuotient(M_global);
        v_new_neg = v_n_neg+dt*a_n_neg;
        u_new_neg = u_n_neg+dt*v_new_neg;
        v_n_pos = v_new_pos;
        u_n_pos = u_new_pos;
        v_n_neg = v_new_neg;
        u_n_neg = u_new_neg;//
        // First BIE Prediction
        // BIE solution correct the top and bot Infitinet BC displacement
        // top -x
        BIE_correct(BIE_top_surf_index.head(nx+1),u_n_pos,FEM_top_disp_x);
        // top -y
        BIE_correct(BIE_top_surf_index.tail(nx+1),u_n_pos,FEM_top_disp_y);
        // bot -x
        BIE_correct(BIE_bot_surf_index.head(nx+1),u_n_neg,FEM_bot_disp_x);
        // bot -y
        BIE_correct(BIE_bot_surf_index.tail(nx+1),u_n_neg,FEM_bot_disp_y);
        
        // BIE solution correct the top and bot Infitinet BC velocity
        // top -x
        BIE_correct(BIE_top_surf_index.head(nx+1),v_n_pos,FEM_top_vel_x);
        // top -y
        BIE_correct(BIE_top_surf_index.tail(nx+1),v_n_pos,FEM_top_vel_y);
        // bot -x
        BIE_correct(BIE_bot_surf_index.head(nx+1),v_n_neg,FEM_bot_vel_x);
        // bot -y
        BIE_correct(BIE_bot_surf_index.tail(nx+1),v_n_neg,FEM_bot_vel_y);
        
        VectorXd u_n_fault_x_pos = VectorXd::Zero(nx+1,1);
        VectorXd v_n_fault_x_pos = VectorXd::Zero(nx+1,1);
        VectorXd u_n_fault_y_pos = VectorXd::Zero(nx+1,1);
        VectorXd v_n_fault_y_pos = VectorXd::Zero(nx+1,1);
        // neg
        VectorXd u_n_fault_x_neg = VectorXd::Zero(nx+1,1);
        VectorXd v_n_fault_x_neg = VectorXd::Zero(nx+1,1);
        VectorXd u_n_fault_y_neg = VectorXd::Zero(nx+1,1);
        VectorXd v_n_fault_y_neg = VectorXd::Zero(nx+1,1);
        // Get the local displacement on the fault
        maplocal(top_surf_index.head(nx+1),u_n_pos,u_n_fault_x_pos);
        maplocal(top_surf_index.head(nx+1),v_n_pos,v_n_fault_x_pos);
        maplocal(top_surf_index.tail(nx+1),u_n_pos,u_n_fault_y_pos);
        maplocal(top_surf_index.tail(nx+1),v_n_pos,v_n_fault_y_pos);
        //
        maplocal(bot_surf_index.head(nx+1),u_n_neg,u_n_fault_x_neg);
        maplocal(bot_surf_index.head(nx+1),v_n_neg,v_n_fault_x_neg);
        maplocal(bot_surf_index.tail(nx+1),u_n_neg,u_n_fault_y_neg);
        maplocal(bot_surf_index.tail(nx+1),v_n_neg,v_n_fault_y_neg);
        // Update Slip and slip rate
        delt_u_n_x = u_n_fault_x_pos-u_n_fault_x_neg;
        delt_u_n_y = u_n_fault_y_pos-u_n_fault_y_neg;
        delt_v_n_x = v_n_fault_x_pos-v_n_fault_x_neg;
        delt_v_n_y = v_n_fault_y_pos-v_n_fault_y_neg;
        printf("Simulation time = %f\n",time_total(j));
        //ofstream file;
        slip.open("slip.bin",ios::binary | ios::app);
        slip.write((char*)(delt_u_n_x.data()),delt_u_n_x.size()*sizeof(double));
        slip.close();
        slip_rate.open("slip_rate.bin",ios::binary | ios::app);
        slip_rate.write((char*)(delt_v_n_x.data()),delt_v_n_x.size()*sizeof(double));
        slip_rate.close();
        shear.open("shear.bin",ios::binary | ios::app);
        shear.write((char*)(T_cx.data()),T_cx.size()*sizeof(double));
        shear.close();
        }
    int endTime = clock();
    cout << "time: " << (endTime-startTime)/double(CLOCKS_PER_SEC)<< endl;
    return 0;
    
}
