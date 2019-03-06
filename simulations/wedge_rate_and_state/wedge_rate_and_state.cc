#include <iostream>
#include <fstream>
#include <Eigen/Eigen>
//include the bie header files
#include "material.hh"
#include "precomputed_kernel.hh"
#include "bimat_interface.hh"
#include "infinite_boundary.hh"
//include the fem header files
//#include "mesh_Abaqus_multi_faults_withAB.hpp"
#include "mesh_Abaqus_multi_faults.hpp"
#include "bcdof.hpp"
#include "bcdof_ptr.hpp"
#include "cal_ke.hpp"
#include "cal_fe_global_vary_ke.hpp"
#include "cal_ab_force.hpp"
#include "mapglobal.hpp"
#include "maplocal.hpp"
#include "cal_M.hpp"
#include "cal_T_abc.hpp"
#include "cal_Bmat.hpp"
#include "Slip_Weakening_lumpM_with_buffer.hpp"
#include "Rate_and_State_Aging.hpp"
#include "cal_slip_sliprate_angle.hpp"
#include "time_advance.hpp"
#include "BIE_correct.hpp"
#include "cal_fe_global_plastic_vary.cpp"
#include "update_disp_velocity.hpp"
#include <math.h>
using namespace Eigen;
using namespace std;

typedef Eigen::Matrix<int, -1, -1,RowMajor> MatrixXi_rm;


double time_fem=0;
double time_bie;
int main() {
    
    std::string filename = "wedge_50m.inp";

    std::cout<<filename<<std::endl;
   // double fault_angle  = 0/180.0*M_PI;
    int dim = 2.0;
    double x_max = 10.0e3;
    double x_min = -10.0e3;
    double dx = 50.0;
    double dy = 50.0;
    std::vector<int> BIE_top, BIE_bot;
    //std::vector<int> AB_left, AB_right;
    MatrixXd Node;
    MatrixXi_rm Element;
    // Total number of faults icluding the main fault
    // number of secondary faults + 1
    int num_faults = 1;
    std::vector<std::vector<int>> fault_nodes(2*num_faults);
    mesh_Abaqus_multi_faults(filename, Node, Element, BIE_top, BIE_bot,fault_nodes,num_faults);
   
    int nx_BIE = BIE_top.size();
    std::cout<<"nx_BIE="<<nx_BIE;
    int n_nodes = Node.rows();
    int n_el = Element.rows();
    int Ndofn = 2;
    int Nnel = Element.cols();
    Node = Node*1e3;
    // Material
    double density = 2670.0;
    double v_s =3.464e3;
    double v_p = 6.0e3;
    double G= pow(v_s,2)*density;
    double Lambda = pow(v_p,2)*density-2.0*G;
    double E  = G*(3.0*Lambda+2.0*G)/(Lambda+G);
    double nu = Lambda/(2.0*(Lambda+G));
    //std::cout<<nu<<std::endl;
    // Slip weakening friction parameters
    double Dc = 0.2;
    double mu_d= 0.3;
    //double mu_s = 0.677;
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
    std::vector<Eigen::Array<double, -1, 1> > mu_s(num_faults);


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
        tau_s[i] = ArrayXd::Zero(Ndofn*(nx_faults[i]),1);
        mu_s[i] = ArrayXd::Zero(nx_faults[i], 1);
    }
    // Making a copy of delt_u and delt_v to calculate the frictional disspation
    std::vector<Eigen::Array<double, -1, 1> > delt_u_n_copy(num_faults);
    std::vector<Eigen::Array<double, -1, 1> > delt_v_n_copy(num_faults);
    for (int i=0;i<num_faults;i++)
    {
        delt_u_n_copy[i] = ArrayXd::Zero(Ndofn*(nx_faults[i]),1);
        delt_v_n_copy[i] = ArrayXd::Zero(Ndofn*(nx_faults[i]),1);
    }
    
    VectorXd F_ext_global = VectorXd::Zero(n_nodes*Ndofn,1);
    // Write the x coordiantes for each fault
    std::vector<double> fault_angle(num_faults);
    double x1 = 0.0;
    double x2 = 0.0;
    double y1 = 0.0;
    double y2 = 0.0;
    for (int j=0;j<num_faults; j++)
    {
        VectorXd x = VectorXd::Zero(nx_faults[j], 1);
        VectorXd y = VectorXd::Zero(nx_faults[j], 1);

        for (int i=0; i<x.size(); i++)
        {
            x(i) = Node(fault_nodes[2*j][i],0);
            y(i) = Node(fault_nodes[2*j][i],1);
        }
   
        //std::ofstream x_output("results/fault_x_coord.txt");
        std::string  x_fault= "results/x_fault_"+std::to_string(j)+".txt";
        std::ofstream x_output(x_fault);
        x_output<<x;
        
        x1 = Node(fault_nodes[2*j][0],0);
        x2 = Node(fault_nodes[2*j][nx_faults[j]-1],0);
        y1 = Node(fault_nodes[2*j][0],1);
        y2 = Node(fault_nodes[2*j][nx_faults[j]-1],1);
        fault_angle[j]=std::atan2(y2-y1,x2-x1);
        //std::cout<<x1/1e3<<" "<<x2/1e3<<" "<<y1/1e3<<" "<<y2/1e3<<std::endl;
      }
    
    
    //Adding pertubation to the main fault 0
    VectorXd x_main = VectorXd::Zero(nx_faults[0], 1);
    
    for (int i=0; i<x_main.size(); i++)
    {
        x_main(i) = Node(fault_nodes[0][i],0);
    }
    
    // Get the index degree of freedome for each element
    VectorXi index_el = VectorXi::Zero(Ndofn*Nnel,1);
    MatrixXi index_store = MatrixXi::Zero(Nnel*Ndofn,n_el);
    for (int i=0;i<n_el;i++)
    {
        bcdof(Element.row(i),dim,index_el);
        index_store.col(i) = index_el;
    }
  //  VectorXi BIE_top_surf_index = VectorXi::Zero(Ndofn*(BIE_top_surf_nodes.size()),1);
   // VectorXi BIE_bot_surf_index = VectorXi::Zero(Ndofn*(BIE_bot_surf_nodes.size()),1);
    VectorXi BIE_top_surf_index = VectorXi::Zero(Ndofn*(BIE_top.size()),1);
    VectorXi BIE_bot_surf_index = VectorXi::Zero(Ndofn*(BIE_bot.size()),1);
    std::vector<Eigen::ArrayXi> Fault_surf_index(num_faults*2);
    for (int i=0; i <num_faults*2; i++)
    {
        Fault_surf_index[i] = ArrayXi::Zero(Ndofn*(fault_nodes[i].size()),1);
    }
    
//    VectorXi AB_left_surf_index = VectorXi::Zero(Ndofn*(AB_left.size()),1);
//    VectorXi AB_right_surf_index = VectorXi::Zero(Ndofn*(AB_right.size()),1);
//
//    bcdof_ptr(AB_left,dim,AB_left_surf_index.data());
//    bcdof_ptr(AB_right,dim,AB_right_surf_index.data());

    
    // Getting the faults DOF index
    for (int i=0; i<num_faults*2;i++)
    {
        bcdof_ptr(fault_nodes[i], dim, Fault_surf_index[i].data());
    }
    bcdof_ptr(BIE_top,dim,BIE_top_surf_index.data());
    bcdof_ptr(BIE_bot,dim,BIE_bot_surf_index.data());
    // Calculating the Global Mass Vector (lumped mass)
    // Element stiffness matrix,Element mass
    std::vector<MatrixXd> ke_store(n_el);
    MatrixXd M_el = MatrixXd::Zero(Nnel*Ndofn,Nnel*Ndofn);
    VectorXd M_el_vec = VectorXd::Zero(Nnel*Ndofn,1);
    VectorXd M_global_vec=VectorXd::Zero(n_nodes*Ndofn,1);
    std::vector<double> dx_min_store;
    std::vector<double> detJ(n_el);
  //  std::vector<Eigen::MatrixXd> B_mat(4);
    std::vector<Eigen::MatrixXd> B_mat_temp(4);

    std::vector<std::vector<Eigen::MatrixXd>> B_mat;
    //B_mat[0][0] = MatrixXd::Zero(3, 8);
    for (int i=0; i<n_el;i++)
    {
        MatrixXd coord = MatrixXd::Zero(4,2);
        VectorXi Element_0= Element.row(i);
        coord.row(0) = Node.row(Element_0(0));
        coord.row(1) = Node.row(Element_0(1));
        coord.row(2) = Node.row(Element_0(2));
        coord.row(3) = Node.row(Element_0(3));
        VectorXd dx_cal = VectorXd::Zero(4,1);
        dx_cal(0)=std::sqrt(pow(std::abs(coord(0,0)-coord(1,0)),2)+pow(std::abs(coord(0,1)-coord(1,1)),2));
        dx_cal(1)=std::sqrt(pow(std::abs(coord(1,0)-coord(2,0)),2)+pow(std::abs(coord(1,1)-coord(2,1)),2));
        dx_cal(2)=std::sqrt(pow(std::abs(coord(2,0)-coord(3,0)),2)+pow(std::abs(coord(2,1)-coord(3,1)),2));
        dx_cal(3)=std::sqrt(pow(std::abs(coord(3,0)-coord(0,0)),2)+pow(std::abs(coord(3,1)-coord(0,1)),2));
        double dx_min = dx_cal.minCoeff();
        dx_min_store.push_back(dx_min);
       // std::cout<<dx_cal<<std::endl;
        cal_Bmat(coord, E, nu, B_mat_temp, detJ[i]);
        B_mat.push_back(B_mat_temp);
        cal_ke (coord,E,nu,ke_store[i]);
        cal_M(coord, density,M_el);
        M_el_vec = M_el.rowwise().sum();
        index_el = index_store.col(i);
        mapglobal(index_el,M_global_vec,M_el_vec);
        

    }
    //
    double dx_min = *std::min_element(std::begin(dx_min_store), std::end(dx_min_store));
    // Time
    double alpha = 0.1;
    //dx_min = 100;
    double dt = alpha*dx_min/v_p;
    std::cout<<"dt="<<dt<<std::endl;
    //dt = 0.001;
    // dt = 5.0e-4;
    // Reyleigh Damping
    double beta =0.1;
    double q = beta*dt;
    double time_run = 6.0;
    int numt = time_run/dt;
    // numt= 1 ;
    //
    std::ofstream time_output("results/time.txt");
    time_output<<time_run<<"\n";
    time_output<<dt<<"\n";
    time_output<<numt<<std::endl;
    
    //numt = 1;
    VectorXd time = dt*VectorXd::LinSpaced(numt,1,numt);
    // BIE part initiation
    // Setting up the material property for the BIE code
    Material BIE_top_mat = Material(E,nu,density);
    Material BIE_bot_mat = Material(E,nu,density);
    double length = x_max- x_min;
    // infinte bc BIE call infinite_boundary.cc
    PrecomputedKernel h11("kernels/nu_.25_h11.dat");
    PrecomputedKernel h12("kernels/nu_.25_k12.dat");
    PrecomputedKernel h22("kernels/nu_.25_h22.dat");
    InfiniteBoundary BIE_inf_top(length,nx_BIE,1.0,&BIE_top_mat,&h11,&h12,&h22);
    InfiniteBoundary BIE_inf_bot(length,nx_BIE,-1.0,&BIE_bot_mat,&h11,&h12,&h22);
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
    
    file.open("results/u_n.bin",ios::binary);
   // file.write((char*)(u_n.data()),u_n.size()*sizeof(double));
    file.close();
    file.open("results/v_n.bin",ios::binary);
    // file.write((char*)(u_n.data()),u_n.size()*sizeof(double));
    file.close();
    file.open("results/eq_ep_n.bin",ios::binary | ios::app);
    // file.write((char*)(u_n.data()),u_n.size()*sizeof(double));
    file.close();
    file.open("results/E_p.bin",ios::binary | ios::app);
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
    // output BIE boundary u_n v_n
    VectorXd BIE_top_u_n = VectorXd::Zero(Ndofn*(nx_faults[0]),1);
    VectorXd BIE_bot_u_n = VectorXd::Zero(Ndofn*(nx_faults[0]),1);
    VectorXd BIE_top_v_n = VectorXd::Zero(Ndofn*(nx_faults[0]),1);
    VectorXd BIE_bot_v_n = VectorXd::Zero(Ndofn*(nx_faults[0]),1);
    // Add output files BIE_top BIE_bot u_n v_n
    file.open("results/BIE_top_u_n.bin",ios::binary | ios::app);
    file.close();
    file.open("results/BIE_top_v_n.bin",ios::binary | ios::app);
    file.close();
    file.open("results/BIE_bot_u_n.bin",ios::binary | ios::app);
    file.close();
    file.open("results/BIE_bot_u_n.bin",ios::binary | ios::app);
    file.close();
    
    std::ofstream Element_output("results/Element.txt");
    std::ofstream Node_output("results/Node.txt");
    Node_output<<Node;
    Element_output<<Element;
    
    double a = 0.008;
    double b = 0.012;
    double V_0 = 1e-6;
    double f_0 = 0.6;
    double L = 0.02;
    double tau_ini = 75e6;
    double sigma_ini = 120e6;
    double V_ini  = 1e-12;
    VectorXd a_array = a*VectorXd::Zero(nx_faults[0],1);
    VectorXd b_array = b*VectorXd::Zero(nx_faults[0],1);
    VectorXd theta_ini = b*VectorXd::Zero(nx_faults[0],1);
    VectorXd theta_n = VectorXd::Zero(nx_faults[0],1);
    VectorXd theta_dot_n = VectorXd::Zero(nx_faults[0],1);
    

    double Vw_width = 5.0e3;
    for (int i =0; i<nx_faults[0]; i++){
        if (abs(x_main(i))<Vw_width/2.0)
        {
            a_array(i) = 0.008;
            b_array(i) = 0.012;
            theta_ini(i) = L/V_0*exp((a_array(i)*log(2*sinh(tau_ini/(a_array(i)*sigma_ini)))-f_0-a_array(i)*log(V_ini/V_0))/b_array(i));
        }
        else
        {
            b_array(i) = 0.012;
            a_array(i) = 0.008;
            theta_ini(i) = L/V_0*exp((a_array(i)*log(2*sinh(tau_ini/(a_array(i)*sigma_ini)))-f_0-a_array(i)*log(V_ini/V_0))/b_array(i));
        }
    }
    theta_n = theta_ini;
    double M = density*dx*dy*1.0;
   
    // Rate and state initlaized v_n delt_v_n
    for (int i=0;i<nx_faults[0];i++)
    {
        delt_v_n[0](2*i)= V_ini;
     
        v_n(Fault_surf_index[0](2*i))= V_ini/2.0;
        
        v_n(Fault_surf_index[1](2*i))= -V_ini/2.0;
    }
    
    // Setting intial stress on the fault
    for (int i=0; i<T_0[0].size()/2; i++)
    {
        T_0[0](2*i) = tau_ini;
        T_0[0](2*i+1) = -sigma_ini;
    }
    VectorXd x = VectorXd::LinSpaced(nx_faults[0],x_min,x_max);
    double d_tau_0 = 25e6;
    double R = 2.5e3;
    VectorXd F_r = VectorXd::Zero(nx_faults[0],1);
    for (int i=0; i<nx_faults[0]; i++)
    {
        if (fabs(x(i))<R)
        {
            F_r(i) = exp(pow(x(i),2)/( pow(x(i),2)-pow(R,2)));
        }
        else
        {
            F_r(i) = 0.0;
        }
    }
    VectorXd G_t = VectorXd::Zero(numt,1);
    double T_ramp = 1.0;
    for (int i=0; i<numt;i++)
    {
        if(time(i)<T_ramp)
        {
            G_t(i) = exp(pow((time(i)-T_ramp),2)/(time(i)*(time(i)-2*T_ramp)));
        }
        else
        {
            G_t(i) = 1.0;
        }
    }
    
    MatrixXd d_tau = MatrixXd::Zero(nx_faults[0],numt);
    for (int i=0;i<numt;i++)
    {
        d_tau.col(i)=75e6*VectorXd::Ones(nx_faults[0],1)+d_tau_0*F_r*(G_t(i));
        // d_tau.col(i)=75e6*VectorXd::Ones(nx_el+1,1)+3.75e6*F_r+dt*0.3125e7*i*VectorXd::Ones(nx_el+1,1);
    }
    

   // time_output<<time_run<<"\n";
   // time_output<<dt<<"\n";
    // Main time loop
    for (int j=0;j<numt;j++)
    {
        // Compute the global internal force
        VectorXd fe_global= VectorXd::Zero(n_nodes*Ndofn,1);
        cal_fe_global_vary_ke(n_nodes, n_el, index_store, q, u_n, v_n, Ndofn, ke_store, fe_global);
      //  cal_fe_global_plastic_vary(n_el, index_store, ke_store, B_mat, detJ, E, nu, q, u_n, v_n, Ndofn, fe_global, strain_n_store,stress_n_store, sxx_initial, syy_initial, sxy_initial, cohes, blkfric,eq_ep_out,E_p_out);
        // Friction subroutine
        VectorXd F_total = F_ext_global-fe_global;
        for (int i=0 ; i<num_faults; i++)
        {
            VectorXd F_fault = VectorXd::Zero(Ndofn*(nx_faults[i]),1);
            VectorXd M_pos = VectorXd::Zero(Ndofn*nx_faults[i], 1);
            VectorXd M_neg = VectorXd::Zero(Ndofn*nx_faults[i],1);
            
            maplocal(Fault_surf_index[2*i], M_global_vec, M_pos);
            maplocal(Fault_surf_index[2*i+1], M_global_vec, M_neg);
            
           // Slip_Weakening_lumpM_with_buffer(M_global_vec, Fault_surf_index[2*i], Fault_surf_index[2*i+1], fe_global, dt, dx, dy, nx_faults[i]-1, delt_v_n[i], delt_u_n[i], T_0[i], tau_s[i], mu_s[i], mu_d, Dc, Ndofn, M_pos, M_neg , F_fault, T_c[i], fault_angle[i]);
            
            for (int i=0;i<nx_faults[0];i++)
            {
                T_0[0](2*i) = d_tau(i,j);
            }
            
            Rate_and_State_Aging(M_global_vec, Fault_surf_index[0], Fault_surf_index[1], fe_global, dt, dx, dy, nx_faults[i]-1, delt_u_n[0], delt_v_n[0], T_0[i], tau_s[i], a_array, b_array , L, V_0, f_0, Ndofn, M, F_fault, T_c[i], theta_n, theta_dot_n);
            
            
            mapglobal(Fault_surf_index[2*i],F_total,-F_fault);
            mapglobal(Fault_surf_index[2*i+1],F_total,F_fault);
        }
        
        // Calculate Absorbing force
      //  cal_ab_force(T_abc, F_total, AB_left_surf_index, v_n, Ndofn);
      //  cal_ab_force(T_abc, F_total, AB_right_surf_index, v_n, Ndofn);
       // Central Difference Time integration
        time_advance(u_n, v_n, F_total, M_global_vec, dt);
        
        update_disp_velocity(u_n, v_n, Fault_surf_index[0], Fault_surf_index[1], Ndofn, nx_faults[0]-1, delt_u_n[0], delt_v_n[0]);

        delt_v_n[0](0)=0;
        delt_v_n[0](2*nx_faults[0]-2)=0;
        delt_u_n[0](0)=0;
        delt_u_n[0](2*nx_faults[0]-2)=0;
        // Correct the BIE surf nodes solutions from FEM with the BIE solution
        BIE_correct(BIE_top_surf_index, BIE_bot_surf_index, fe_global, Ndofn, nx_BIE-1, dx, BIE_inf_top, BIE_inf_bot, u_n, v_n);

    

        for (int i=0;i<1;i++)
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
//        
        if (j%100==1)
        {
            file.open("results/u_n.bin",ios::binary | ios::app);
            file.write((char*)(u_n.data()),u_n.size()*sizeof(double));
            file.close();
            file.open("results/v_n.bin",ios::binary | ios::app);
            file.write((char*)(v_n.data()),v_n.size()*sizeof(double));
            file.close();
        }
        
//  Adding virtual boundary information
        maplocal(BIE_top_surf_index,u_n,BIE_top_u_n);
        maplocal(BIE_top_surf_index,v_n,BIE_top_v_n);
        maplocal(BIE_bot_surf_index,u_n,BIE_bot_u_n);
        maplocal(BIE_bot_surf_index,v_n,BIE_bot_v_n);
        
        
        file.open("results/BIE_top_u_n.bin",ios::binary | ios::app);
        file.write((char*)(BIE_top_u_n.data()),BIE_top_u_n.size()*sizeof(double));
        file.close();
        file.open("results/BIE_top_v_n.bin",ios::binary | ios::app);
        file.write((char*)(BIE_top_v_n.data()),BIE_top_v_n.size()*sizeof(double));
        file.close();
        file.open("results/BIE_bot_u_n.bin",ios::binary | ios::app);
        file.write((char*)(BIE_bot_u_n.data()),BIE_bot_u_n.size()*sizeof(double));
        file.close();
        file.open("results/BIE_bot_v_n.bin",ios::binary | ios::app);
        file.write((char*)(BIE_bot_v_n.data()),BIE_bot_v_n.size()*sizeof(double));
        file.close();
//
        printf("Simulation time = %f\n",time(j));
    }
    return 0;
}
