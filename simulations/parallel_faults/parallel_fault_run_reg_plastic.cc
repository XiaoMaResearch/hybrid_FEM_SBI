#include <iostream>
#include <fstream>
#include <Eigen/Eigen>
//include the bie header files
#include "material.hh"
#include "precomputed_kernel.hh"
#include "bimat_interface.hh"
#include "infinite_boundary.hh"
//include the fem header files
#include "mesh_Generated_multi_faults.hpp"
#include "bcdof.hpp"
#include "bcdof_ptr.hpp"
#include "cal_ke.hpp"
#include "cal_fe_global_plastic.hpp"
#include "mapglobal.hpp"
#include "Slip_Weakening_reg.hpp"
#include "Slip_Weakening.hpp"
#include "cal_slip_sliprate.hpp"
#include "time_advance.hpp"
#include "BIE_correct.hpp"
#include "igl/list_to_matrix.h"
#include "Getplastic.hpp"
#include "cal_Bmat.hpp"

//#include <omp.h>
using namespace Eigen;
using namespace std;

typedef Eigen::Matrix<int, -1, -1,RowMajor> MatrixXi_rm;

void read_matrix(std::string fileName, Eigen::MatrixXd &outputMat);

double time_fem=0;
double time_bie;
int main() {
     time_fem = 0.0;
    // Domain Size
    double x_min = -10e3;
    double x_max = 10e3;
    double y_min = -0.5e3;
    double y_max = 0.5e3;
    int dim = 2.0;
    double dx = 50;
    double dy = 50;
    int nx_el = (x_max-x_min)/dx;
    int ny = (y_max-y_min)/dy;
    MatrixXd Node = MatrixXd::Zero((nx_el+1)*(ny+1),2);
    MatrixXi_rm Element(nx_el*ny,4); Element.setZero();
    ArrayXi BIE_top_surf_nodes = ArrayXi::Zero((nx_el+1),1);
    ArrayXi BIE_bot_surf_nodes = ArrayXi::Zero((nx_el+1),1);
    int num_faults = 1;
    // Mesh
    // Position of the fault : y_pos, x_left, x_right
    MatrixXd fault_pos(num_faults,3);
    fault_pos << 0.0,  -10.0e3 ,10.0e3;

//    MatrixXd fault_pos;
//    read_matrix("fault_pos.txt", fault_pos);
//    int num_faults = fault_pos.rows();
    cout<<"test_mat="<<"\n"<<fault_pos<<"\n"<<"rows="<<fault_pos.rows()<<"\n";
    std::vector<std::vector<int>> fault_nodes(num_faults*2);
    mesh_Generated_multi_faults(x_min,x_max,y_min,y_max,dx,dy,nx_el,ny,Node, Element,BIE_top_surf_nodes, BIE_bot_surf_nodes,fault_pos, fault_nodes);
    
    int n_nodes = Node.rows();
    int n_el = Element.rows();
    int Ndofn = 2;
    int Nnel = Element.cols();
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
    double time_run = 1.0;
    int numt = time_run/dt;
    //numt = 1;
    //numt = 200;
    //numt = 3;
    //numt =4;
    VectorXd time = dt*VectorXd::LinSpaced(numt,1,numt);
    // Slip weakening friction parameters
    double Dc = 0.2;
    double mu_d= 0.1;
    double mu_s = 0.6;
    // Intialization
    // disp velocity current and next time step (new)
    VectorXd u_n = VectorXd::Zero(n_nodes*Ndofn,1);
    VectorXd v_n = VectorXd::Zero(n_nodes*Ndofn,1);
    VectorXd u_new = VectorXd::Zero(n_nodes*Ndofn,1);
    VectorXd v_new = VectorXd::Zero(n_nodes*Ndofn,1);
    VectorXd a_n = VectorXd::Zero(n_nodes*Ndofn,1);
    // slip and slip-rate
    std::vector<Eigen::Array<double, -1, 1> > delt_u_n(num_faults);
    std::vector<Eigen::Array<double, -1, 1> > delt_v_n(num_faults);
    std::vector<Eigen::Array<double, -1, 1> > T_c(num_faults);
    std::vector<Eigen::Array<double, -1, 1> > T_0(num_faults);
    std::vector<Eigen::Array<double, -1, 1> > tau_s(num_faults);


    // Vector containing number of elements on each faults (nx)
    std::vector<int> nx_faults(num_faults);
    for (int i=0;i<num_faults;i++)
    {
        nx_faults[i] =fault_nodes[2*i].size();
    }
    
    
    for (int i=0;i<num_faults;i++)
    {
        delt_u_n[i] = ArrayXd::Zero(Ndofn*(nx_faults[i]),1);
        delt_v_n[i] = ArrayXd::Zero(Ndofn*(nx_faults[i]),1);
        T_c[i] = ArrayXd::Zero(Ndofn*(nx_faults[i]),1);
        T_0[i] = ArrayXd::Zero(Ndofn*(nx_faults[i]),1);
       // tau_s[i] = ArrayXd::Zero(Ndofn*(nx_faults[i]),1);
        tau_s[i] = 0.6*ArrayXd::Ones((nx_faults[i]),1)*50.0e6;
    }
    VectorXd F_ext_global = VectorXd::Zero(n_nodes*Ndofn,1);
    // Setting intial stress on the fault
    for (int j=0; j<num_faults;j++)
    {
        for (int i=0; i<nx_faults[j]; i++)
        {
            T_0[j](2*i+1) = -50.0e6;
            T_c[j](2*i+1) = T_0[j](2*i+1);
        }
        if (j==0)
        {
            VectorXd x = VectorXd::LinSpaced(nx_faults[j], fault_pos(j,1), fault_pos(j,2));
            for (int i=0 ; i<nx_faults[j]; i++)
            {
                if ((x(i)<=(fault_pos(j,1)+fault_pos(j,2))/2+0.6e3)&&(x(i)>=(fault_pos(j,1)+fault_pos(j,2))/2-0.6e3))
                {
                    T_0[j](2*i) = 31.0e6;
                    T_c[j](2*i) = T_0[j](2*i);
                }
                else
                {
                    T_0[j](2*i) = 19.28e6;
                    T_c[j](2*i) = T_0[j](2*i);
                }
            }
        }
        else
        {
            for (int i=0 ; i<nx_faults[j]; i++)
            {
                T_0[j](2*i) = 19.28e6;
            }
        }
    }
    // Get the index degree of freedome for each element
    VectorXi index_el = VectorXi::Zero(Ndofn*Nnel,1);
    MatrixXi index_store = MatrixXi::Zero(Nnel*Ndofn,n_el);
    for (int i=0;i<n_el;i++)
    {
        bcdof(Element.row(i),dim,index_el);
        index_store.col(i) = index_el;
    }
    VectorXi BIE_top_surf_index = VectorXi::Zero(Ndofn*(BIE_top_surf_nodes.size()),1);
    VectorXi BIE_bot_surf_index = VectorXi::Zero(Ndofn*(BIE_bot_surf_nodes.size()),1);
    
    std::vector<Eigen::ArrayXi> Fault_surf_index(num_faults*2);
    for (int i=0; i <num_faults*2; i++)
    {
        Fault_surf_index[i] = ArrayXi::Zero(Ndofn*(fault_nodes[i].size()),1);
    }
  
    // Getting the faults DOF index
    for (int i=0; i<num_faults*2;i++)
    {
        bcdof_ptr(fault_nodes[i], dim, Fault_surf_index[i].data());
    }
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
    VectorXi Element_0= Element.row(0);
    coord.row(0) = Node.row(Element_0(0));
    coord.row(1) = Node.row(Element_0(1));
    coord.row(2) = Node.row(Element_0(2));
    coord.row(3) = Node.row(Element_0(3));
    cal_ke (coord,E,nu,ke);
    double detJ = 0.0;
    std::vector<Eigen::MatrixXd> B_mat(4);
    cal_Bmat(coord, E, nu, B_mat, detJ);
    // BIE part initiation
    // Setting up the material property for the BIE code
    Material BIE_top_mat = Material(E,nu,density);
    Material BIE_bot_mat = Material(E,nu,density);
    double length = x_max-x_min;
    // infinte bc BIE call infinite_boundary.cc
    PrecomputedKernel h11("kernels/nu_.25_h11.dat");
    PrecomputedKernel h12("kernels/nu_.25_k12.dat");
    PrecomputedKernel h22("kernels/nu_.25_h22.dat");
    InfiniteBoundary BIE_inf_top(length,nx_el+1,1.0,&BIE_top_mat,&h11,&h12,&h22);
    InfiniteBoundary BIE_inf_bot(length,nx_el+1,-1.0,&BIE_bot_mat,&h11,&h12,&h22);
    // BIE setting time step
    BIE_inf_top.setTimeStep(dt);
    BIE_inf_bot.setTimeStep(dt);
    // BIE initialization
    BIE_inf_top.init();
    BIE_inf_bot.init();
    printf("ready to start\n");
    // Output
    ofstream file;
    file.open("results/num_nodes_fault.bin",ios::binary);
    file.write((char*)(nx_faults.data()),nx_faults.size()*sizeof(int));
    file.close();
    
    for (int i=0;i<num_faults;i++)
    {
        std::string slip = "results/slip_"+std::to_string(i)+".bin";
        file.open(slip);
        file.close();
        std::string slip_rate = "results/slip_rate_"+std::to_string(i)+".bin";
        file.open(slip_rate);
        file.close();
        std::string shear = "results/shear_"+std::to_string(i)+".bin";
        file.open(shear);
        file.close();
    }
   // double start = omp_get_wtime();
    std::vector<MatrixXd> strain_n_store(n_el);
    std::vector<MatrixXd> stress_n_store(n_el);

    for (int i=0; i<n_el; i++)
    {
        strain_n_store[i] = MatrixXd::Zero(3, 4);
        stress_n_store[i] = MatrixXd::Zero(3, 4);

    }
    double cohes = 100.0e20;
    double blkfric = 0.6;
    double sxx_initial = -122.54e6;
    double syy_initial = -50.0e6;
    double sxy_initial = 19.285e6;
    // Main time loop
    for (int j=0;j<numt;j++)
    {
        // Compute the global internal force
        VectorXd fe_global= VectorXd::Zero(n_nodes*Ndofn,1);
      //  cal_fe_global_const_ke(n_nodes, n_el, index_store, q, u_n, v_n, Ndofn, ke, fe_global);
        cal_fe_global_plastic(n_el, index_store, ke,B_mat, detJ, E, nu, q, u_n, v_n, Ndofn, fe_global, strain_n_store,stress_n_store,
                              sxx_initial, syy_initial, sxy_initial, cohes, blkfric);
        // Friction subroutine
        VectorXd F_total = F_ext_global-fe_global;
        for (int i=0 ; i<num_faults; i++)
        {
            VectorXd F_fault = VectorXd::Zero(Ndofn*(nx_faults[i]),1);
            Slip_Weakening(M_global_vec, Fault_surf_index[2*i], Fault_surf_index[2*i+1], fe_global, dt, dx, dy, nx_faults[i]-1, delt_v_n[i], delt_u_n[i], T_0[i], tau_s[i], mu_s, mu_d, Dc, Ndofn, M, F_fault, T_c[i]);
            mapglobal(Fault_surf_index[2*i],F_total,-F_fault);
            mapglobal(Fault_surf_index[2*i+1],F_total,F_fault);
        }
        
        // Central Difference Time integration
        time_advance(u_n, v_n, F_total, M_global_vec, dt);
        // Get the slip and slip rate
        
        for (int i=0; i<num_faults; i++)
        {
            cal_slip_slip_rate(u_n, v_n,  Fault_surf_index[2*i], Fault_surf_index[2*i+1] , Ndofn, nx_faults[i]-1, delt_u_n[i], delt_v_n[i]);
        }
        // Correct the BIE surf nodes solutions from FEM with the BIE solution
        BIE_correct(BIE_top_surf_index, BIE_bot_surf_index, fe_global, Ndofn, nx_el, dx, BIE_inf_top, BIE_inf_bot, u_n, v_n);
        
//        VectorXd delt_u_n_out = delt_u_n[0];
//        VectorXd delt_v_n_out = delt_v_n[0];
//        VectorXd T_c_out = T_c[0];
//        file.open("results/slip_main.bin",ios::binary | ios::app);
//        file.write((char*)(delt_u_n[0].data()),delt_u_n[0].size()*sizeof(double));
//        file.close();
//        file.open("results/slip_rate_main.bin",ios::binary | ios::app);
//        file.write((char*)(delt_v_n[0].data()),delt_v_n[0].size()*sizeof(double));
//        file.close();
//        file.open("results/shear_main.bin",ios::binary | ios::app);
//        file.write((char*)(T_c[0].data()),T_c[0].size()*sizeof(double));
//        file.close();
        //
    //    if (j%4==3)
       // {
        for (int i=0;i<num_faults;i++)
        {
            std::string slip = "results/slip_"+std::to_string(i)+".bin";
            file.open(slip,ios::binary | ios::app);
            file.write((char*)(delt_u_n[i].data()),delt_u_n[i].size()*sizeof(double));
            file.close();
            
            std::string slip_rate = "results/slip_rate_"+std::to_string(i)+".bin";
            file.open(slip_rate,ios::binary | ios::app);
            file.write((char*)(delt_v_n[i].data()),delt_v_n[i].size()*sizeof(double));
            file.close();
            
            std::string shear = "results/shear_"+std::to_string(i)+".bin";
            file.open(shear,ios::binary | ios::app);
            file.write((char*)(T_c[i].data()),T_c[i].size()*sizeof(double));
            file.close();
        }
    //    }
        printf("Simulation time = %f\n",time(j));
      //  double end_t = omp_get_wtime();
       // std::cout<<"time_cpu_t="<<end_t-start<<std::endl;
    //   std::cout<<"time_fem="<< time_fem<<std::endl;
    //   std::cout<<"time_bie="<< time_bie<<std::endl;
    }
   // double end = omp_get_wtime();
    //std::cout<<"time_cpu="<<end-start<<std::endl;
    return 0;
}

void read_matrix(std::string fileName, Eigen::MatrixXd &outputMat) {
    fstream cin;
    cin.open(fileName.c_str());
    if (cin.fail())
    {
        std::cerr << "Failed to open file: " << fileName << std::endl;
        std::cin.get(); }
    string s;
    vector <vector <double> > matrix;
    while (getline(cin, s)) {
        stringstream input(s);
        double temp;
        vector <double> currentLine;
        while (input >> temp)
            currentLine.push_back(temp);
        matrix.push_back(currentLine);
    }
    if (!igl::list_to_matrix(matrix, outputMat))
    { std::cerr << "list tom matrix error" << std::endl; std::cin.get();
        //return false;
    }
}
