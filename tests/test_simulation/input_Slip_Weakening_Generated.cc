#include <iostream>
#include <fstream>
#include <Eigen/Eigen>
//include the bie header files
#include "material.hh"
#include "precomputed_kernel.hh"
#include "bimat_interface.hh"
#include "infinite_boundary.hh"
//include the fem header files
#include "mesh_Generated.hpp"
#include "bcdof.hpp"
#include "cal_ke.hpp"
#include "cal_fe_global_const_ke.hpp"
#include "mapglobal.hpp"
#include "Slip_Weakening.hpp"
#include "cal_slip_sliprate.hpp"
#include "time_advance.hpp"
#include "BIE_correct.hpp"


using namespace Eigen;
using namespace std;

int main() {
    // Domain Size
    double x_min = -50e3;
    double x_max = 50e3;
    double y_min = -1.0e3;
    double y_max = 1.0e3;
    int dim = 2.0;
    double dx = 100;
    double dy = 100;
    int nx = (x_max-x_min)/dx;
    int ny = (y_max-y_min)/dy;
    MatrixXd Node = MatrixXd::Zero((nx+1)*(ny+1),2);
    MatrixXd Element = MatrixXd::Zero(nx*ny,4);
    VectorXd fault_surf_nodes = VectorXd::Zero((nx+1),1);
    VectorXd fault_surf_nodes_new = VectorXd::Zero((nx+1),1);
    VectorXd BIE_top_surf_nodes = VectorXd::Zero((nx+1),1);
    VectorXd BIE_bot_surf_nodes = VectorXd::Zero((nx+1),1);
    // Mesh
    mesh_Generated(x_min,x_max,y_min,y_max,dx,dy,nx,ny,Node, Element, fault_surf_nodes, fault_surf_nodes_new, BIE_top_surf_nodes, BIE_bot_surf_nodes);
    int n_nodes = Node.rows();
    int n_el = Element.rows();
    int Ndofn = 2;
    int Nnel = Element.cols();
    nx = fault_surf_nodes.size()-1;
    // Material
    double density = 2670.0;
    double v_s =3.464e3;
    double v_p = 6.0e3;
    double G= pow(v_s,2)*density;
    double Lambda = pow(v_p,2)*density-2.0*G;
    double E  = G*(3.0*Lambda+2.0*G)/(Lambda+G);
    double nu = Lambda/(2.0*(Lambda+G));
    // Time
    double alpha = 0.4;
    double dt = alpha*dx/v_p;
    // Reyleigh Damping
    double beta =0.1;
    double q = beta*dt;
    double time_run = 6.0;
    int numt = time_run/dt;
    //numt = 3;
    VectorXd time = dt*VectorXd::LinSpaced(numt,1,numt);
    // Slip weakening friction parameters
    double Dc = 0.2;
    double mu_d= 0.5;
    double mu_s = 0.6;
    // Intialization
    // disp velocity current and next time step (new)
    VectorXd u_n = VectorXd::Zero(n_nodes*Ndofn,1);
    VectorXd v_n = VectorXd::Zero(n_nodes*Ndofn,1);
    VectorXd u_new = VectorXd::Zero(n_nodes*Ndofn,1);
    VectorXd v_new = VectorXd::Zero(n_nodes*Ndofn,1);
    VectorXd a_n = VectorXd::Zero(n_nodes*Ndofn,1);
    // slip and slip-rate
    VectorXd delt_u_n = VectorXd::Zero(Ndofn*(nx+1),1);
    VectorXd delt_v_n = VectorXd::Zero(Ndofn*(nx+1),1);
    // Stress on the fault T_0 = intial stress , T = sticking force, T_c= stress critical goes into F_global
    VectorXd T_0 = VectorXd::Zero(Ndofn*(nx+1),1);
    VectorXd T = VectorXd::Zero(Ndofn*(nx+1),1);
    VectorXd T_c = VectorXd::Zero(Ndofn*(nx+1),1);
    VectorXd tau_s = VectorXd::Zero(nx+1,1);
    VectorXd F_ext_global = VectorXd::Zero(n_nodes*Ndofn,1);
    // Setting intial stress on the fault
    for (int i=0; i<T_0.size()/2; i++)
    {
        T_0(2*i+1) = -50.0e6;
    }
    VectorXd x = VectorXd::LinSpaced(nx+1,x_min,x_max);
    for (int i=0 ; i<nx+1; i++)
    {
        if ((x(i)<=(x_max+x_min)/2+0.8e3)&&(x(i)>=(x_max+x_min)/2-0.8e3))
        {
            T_0(2*i) = 31.0e6;
        }
        else
        {
            T_0(2*i) = 27.5e6;
            
        }
    }
    // Get the index degree of freedome for each element
    VectorXd index_el = VectorXd::Zero(Ndofn*Nnel,1);
    MatrixXd index_store = MatrixXd::Zero(Nnel*Ndofn,n_el);
    for (int i=0;i<n_el;i++)
    {
        bcdof(Element.row(i),dim,index_el);
        index_store.col(i) = index_el;
    }
    VectorXd top_surf_index = VectorXd::Zero(Ndofn*(nx+1),1);
    VectorXd bot_surf_index = VectorXd::Zero(Ndofn*(nx+1),1);
    VectorXd BIE_top_surf_index = VectorXd::Zero(Ndofn*(nx+1),1);
    VectorXd BIE_bot_surf_index = VectorXd::Zero(Ndofn*(nx+1),1);
    bcdof(fault_surf_nodes,dim,top_surf_index);
    bcdof(fault_surf_nodes_new,dim,bot_surf_index);
    bcdof(BIE_top_surf_nodes,dim,BIE_top_surf_index);
    bcdof(BIE_bot_surf_nodes,dim,BIE_bot_surf_index);
    // Calculating the Global Mass Vector (lumped mass)
    // Element mass
    double M=density*dx*dy*1.0;
    VectorXd M_el_vec = M/4*VectorXd::Ones(Nnel*Ndofn,1);
    VectorXd M_global_vec=VectorXd::Zero(n_nodes*Ndofn,1);
    for (int i=0 ; i<n_el;i++)
    {
        index_el = index_store.col(i);
        mapglobal(index_el,M_global_vec,M_el_vec);
    }
    // Element matrix
    MatrixXd ke = MatrixXd::Zero(8,8);
    MatrixXd coord = MatrixXd::Zero(4,2);
    VectorXd Element_0= Element.row(0);
    coord.row(0) = Node.row(Element_0(0));
    coord.row(1) = Node.row(Element_0(1));
    coord.row(2) = Node.row(Element_0(2));
    coord.row(3) = Node.row(Element_0(3));
    cal_ke (coord,E,nu,ke);
    // BIE part initiation
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

    //BIE_initiation(E, nu, density, x_max, x_min, nx, dt);
    printf("ready to start\n");
    // Output
    ofstream slip("slip.bin",ios::binary);
    slip.close();
    ofstream slip_rate("slip_rate.bin",ios::binary);
    slip_rate.close();
    ofstream shear("shear.bin",ios::binary);
    shear.close();
    std::ofstream fe("results/fe.txt");

    
    // Main time loop
    for (int j=0;j<numt;j++)
    {
        // Compute the global internal force
        VectorXd fe_global= VectorXd::Zero(n_nodes*Ndofn,1);
        cal_fe_global_const_ke(n_nodes, n_el, index_store, q, u_n, v_n, Ndofn, ke, fe_global);
        // Friction subroutine
        VectorXd F_fault = VectorXd::Zero(Ndofn*(nx+1),1);
        Slip_Weakening(M_global_vec, top_surf_index, bot_surf_index, fe_global, dt, dx, dy, nx, delt_v_n, delt_u_n, T_0, tau_s, mu_s, mu_d, Dc, Ndofn, M, F_fault, T_c);
        // Calculate the global force vector
        // Adding contribution of the fault force to the global force vector F_total
        VectorXd F_total = F_ext_global-fe_global;
        mapglobal(top_surf_index,F_total,-F_fault);
        mapglobal(bot_surf_index,F_total,F_fault);
        // Central Difference Time integration
        time_advance(u_n, v_n, F_total, M_global_vec, dt);
        // Get the slip and slip rate
        cal_slip_slip_rate(u_n, v_n, top_surf_index, bot_surf_index, Ndofn, nx, delt_u_n, delt_v_n);
        // Correct the BIE surf nodes solutions from FEM with the BIE solution
        BIE_correct(BIE_top_surf_index, BIE_bot_surf_index, fe_global, Ndofn, nx, dx, BIE_inf_top, BIE_inf_bot, u_n, v_n);
        //ofstream file;
        slip.open("slip.bin",ios::binary | ios::app);
        slip.write((char*)(delt_u_n.data()),delt_u_n.size()*sizeof(double));
        slip.close();
        slip_rate.open("slip_rate.bin",ios::binary | ios::app);
        slip_rate.write((char*)(delt_v_n.data()),delt_v_n.size()*sizeof(double));
        slip_rate.close();
        shear.open("shear.bin",ios::binary | ios::app);
        shear.write((char*)(T_c.data()),T_c.size()*sizeof(double));
        shear.close();
        printf("Simulation time = %f\n",time(j));
    }
    return 0;
}
