#include <iostream>
#include <fstream>
#include <Eigen/Eigen>
//include the bie header files
#include "material.hh"
#include "precomputed_kernel.hh"
#include "infinite_boundary.hh"
//include the fem header files
#include "mesh_Generated_multi_faults.hpp"
#include "bcdof.hpp"
#include "bcdof_ptr.hpp"
#include "cal_ke.hpp"
#include "cal_fe_global_plastic.hpp"
#include "cal_fe_global_const_ke.hpp"
#include "mapglobal.hpp"
#include "Slip_Weakening_reg.hpp"
#include "Slip_Weakening.hpp"
#include "cal_slip_sliprate.hpp"
#include "time_advance.hpp"
#include "BIE_correct.hpp"
//#include "igl/list_to_matrix.h"
#include "Getplastic.hpp"
#include "cal_Bmat.hpp"
#include "Slip_Weakening_lumpM_with_buffer.hpp"


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
    double x_min = -36e3;
    double x_max = 36e3;
    double y_min = -1.0e3;
    double y_max = 1.0e3;
    int dim = 2.0;
    double dx = 12.5;
    double dy = 12.5;
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
    fault_pos << 0.0,  x_min ,x_max;
    
    //    MatrixXd fault_pos;
    //    read_matrix("fault_pos.txt", fault_pos);
    //    int num_faults = fault_pos.rows();
    cout<<"test_mat="<<"\n"<<fault_pos<<"\n"<<"rows="<<fault_pos.rows()<<"\n";
    std::vector<std::vector<int>> fault_nodes(num_faults*2);
    mesh_Generated_multi_faults(x_min,x_max,y_min,y_max,dx,dy,nx_el,ny,Node, Element,BIE_top_surf_nodes, BIE_bot_surf_nodes,fault_pos, fault_nodes);
    
    std::ofstream Element_output("results/Element.txt");
    std::ofstream Node_output("results/Node.txt");
    
    std::ofstream u_n_test("results/u_n_test.txt");
    std::ofstream v_n_test("results/v_n_test.txt");
    std::ofstream u_n_test_2("results/u_n_test_2.txt");
    std::ofstream v_n_test_2("results/v_n_test_2.txt");
    Node_output<<Node;
    Element_output<<Element;
    
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
    //dt = 5e-4;
    dt = 4e-4;
    //dt = 0.0025;
    //dt = 0.000625;
    // Reyleigh Damping
    double beta =0.1;
    double q = beta*dt;
    double time_run = 0.5;
    int numt = time_run/dt;
    std::ofstream time_output("results/time.txt");
    time_output<<time_run<<"\n";
    time_output<<dt<<"\n";
    time_output<<numt<<std::endl;
    // numt = 100;
    //numt = 200;
    //numt = 3;
    // numt =8;
    VectorXd time = dt*VectorXd::LinSpaced(numt,1,numt);
    // Intialization
    // disp velocity current and next time step (new)
    VectorXd u_n = VectorXd::Zero(n_nodes*Ndofn,1);
    VectorXd v_n = VectorXd::Zero(n_nodes*Ndofn,1);
    VectorXd u_new = VectorXd::Zero(n_nodes*Ndofn,1);
    VectorXd v_new = VectorXd::Zero(n_nodes*Ndofn,1);
    VectorXd a_n = VectorXd::Zero(n_nodes*Ndofn,1);
    // Vector containing number of elements on each faults (nx)

    VectorXd F_ext_global = VectorXd::Zero(n_nodes*Ndofn,1);


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
    
    std::cout<<coord<<std::endl;
    
    cal_ke (coord,E,nu,ke);
    
    std::cout<<ke<<std::endl;
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
    //
    file.open("results/u_n.bin",ios::binary);
    file.close();
    file.open("results/v_n.bin",ios::binary);
    file.close();
    
    // output BIE boundary u_n v_n
//    VectorXd BIE_top_u_n = VectorXd::Zero(2*(nx_el+1),1);
//    VectorXd BIE_bot_u_n = VectorXd::Zero(2*(nx_el+1),1);
//    VectorXd BIE_top_v_n = VectorXd::Zero(2*(nx_el+1),1);
//    VectorXd BIE_bot_v_n = VectorXd::Zero(2*(nx_el+1),1);
    VectorXd BIE_top_u_n ;
    // Add output files BIE_top BIE_bot u_n v_n
    streampos size;
    ifstream binaryIo("results/BIE_bot_u_n.bin",ios::in|ios::binary|ios::ate);
    size = binaryIo.tellg();
    char * memblock;
    cout<<"size = "<<size <<"\n";
    memblock = new char [size];
    binaryIo.seekg(0,ios::beg);
    binaryIo.read(memblock,size);
    binaryIo.close();
    cout<<"the entire file content is in memory \n";
    double *double_u_n = (double*)memblock;  // reinterpret as dobule
    
    streampos sizev;
    ifstream binaryIov("results/BIE_bot_v_n.bin",ios::in|ios::binary|ios::ate);
    sizev = binaryIov.tellg();
    char * memblockv;
    cout<<"size = "<<sizev <<"\n";
    memblockv = new char [sizev];
    binaryIov.seekg(0,ios::beg);
    binaryIov.read(memblockv,sizev);
    binaryIov.close();
    cout<<"the entire file content is in memory \n";
    double *double_v_n = (double*)memblockv;  // reinterpret as dobule
    
    VectorXd u_n_local = VectorXd::Zero(2*(nx_el+1),1);
    VectorXd v_n_local = VectorXd::Zero(2*(nx_el+1),1);


    // Main time loop
    for (int j=0;j<numt;j++)
    {
       // mapglobal(BIE_top_surf_index, u_n, u_n_local);
       // mapglobal(BIE_top_surf_index, v_n, v_n_local);
        // Compute the global internal force
        VectorXd fe_global= VectorXd::Zero(n_nodes*Ndofn,1);
          cal_fe_global_const_ke(n_nodes, n_el, index_store, q, u_n, v_n, Ndofn, ke, fe_global);
        // Friction subroutine
        VectorXd F_total = F_ext_global-fe_global;
        // Central Difference Time integration
        time_advance(u_n, v_n, F_total, M_global_vec, dt);
        
        for (int i=0; i<v_n_local.size();i++)
        {
            double u_n_temp = double_u_n[j*(2*(nx_el+1))+i];
            double v_n_temp = double_v_n[j*(2*(nx_el+1))+i];
            u_n_local(i) = u_n_temp;
            v_n_local(i) = v_n_temp;
        }

       u_n_test<<u_n_local<<std::endl;
       v_n_test<<v_n_local<<std::endl;

        mapglobal(BIE_top_surf_index, u_n, u_n_local);
        mapglobal(BIE_top_surf_index, v_n, v_n_local);
        
        VectorXd u_n_local_2 = VectorXd::Zero(2*(nx_el+1),1);
        VectorXd v_n_local_2 = VectorXd::Zero(2*(nx_el+1),1);
        
        maplocal(BIE_top_surf_index, u_n, u_n_local_2);
      //  maplocal(BIE_top_surf_index, v_n, v_n_local_2);
        
        
        u_n_test_2<<u_n_local_2<<std::endl;
        v_n_test_2<<v_n_local_2<<std::endl;

        // Get the slip and slip rate
        if (j%10==1)
        {
            file.open("results/u_n.bin",ios::binary | ios::app);
            file.write((char*)(u_n.data()),u_n.size()*sizeof(double));
            file.close();
            file.open("results/v_n.bin",ios::binary | ios::app);
            file.write((char*)(v_n.data()),v_n.size()*sizeof(double));
            file.close();
        }
        printf("Simulation time = %f\n",time(j));
    }
    return 0;
}
