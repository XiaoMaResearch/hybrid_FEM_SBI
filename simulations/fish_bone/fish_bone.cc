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
#include "cal_slip_sliprate_angle.hpp"
#include "time_advance.hpp"
#include "BIE_correct.hpp"
#include "cal_fe_global_plastic_vary.cpp"
#include <math.h>
using namespace Eigen;
using namespace std;

typedef Eigen::Matrix<int, -1, -1,RowMajor> MatrixXi_rm;


double time_fem=0;
double time_bie;
int main() {
    
    std::string filename = "fish_bone_DL30R_dx25m_Lf2R_SR_D30.inp";

    std::cout<<filename<<std::endl;
   // double fault_angle  = 0/180.0*M_PI;
    int dim = 2.0;
    double x_max = 15.0e3;
    double x_min = -15.0e3;
    double dx = 25.0;
    double dy = 25.0;
    std::vector<int> BIE_top, BIE_bot;
    //std::vector<int> AB_left, AB_right;
    MatrixXd Node;
    MatrixXi_rm Element;
    // Total number of faults icluding the main fault
    // number of secondary faults + 1
    int num_faults = 16;
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
    // Calcualte the resolving shear and nommal stress on each fault
    std::vector<double> sigma_f(num_faults);
    std::vector<double> tau_f(num_faults);
    double sigma_xx = 50.0e6;
    double sigma_yy = 50.0e6;
    double tau_xy = 20.0e6;

    for (int i=0; i<num_faults; i++)
    {
        double c = std::cos(fault_angle[i]);
        double s = std::sin(fault_angle[i]);
        sigma_f[i] = sigma_xx*s*s+sigma_yy*c*c+2*tau_xy*s*c;
        tau_f[i] = -sigma_xx*s*c+sigma_yy*s*c-tau_xy*(c*c-s*s);
      if (i!=0)
        {
            sigma_f[i] =32.0e6;
            tau_f[i] =  -18.0e6;
        }
        //double logistic =std::abs(sigma_f[i])*0.8>std::abs(tau_f[i]);
   std::cout<<"fault="<<i<<"angle="<<fault_angle[i]<<"\t sigma_f = "<<sigma_f[i]<<"\t tau_f="<<tau_f[i]<<std::endl;

    }
    
    
   // std::cout<<x<<std::endl;
    // Setting intial stress on the all fault according the the stress depcomposition.
    for (int j=0; j<num_faults;j++)
    {
        for (int i=0; i<nx_faults[j]; i++)
        {
            T_0[j](2*i+1) = -sigma_f[j];
            T_0[j](2*i) =  -tau_f[j];
            
        }
        if (j==0)
        {
            mu_s[j] =0.6*ArrayXd::Ones(nx_faults[j], 1);
        }
        else
        {
            mu_s[j] =0.6*ArrayXd::Ones(nx_faults[j], 1);
        }
    }
    //Adding pertubation to the main fault 0
    VectorXd x_main = VectorXd::Zero(nx_faults[0], 1);
    
    for (int i=0; i<x_main.size(); i++)
    {
        x_main(i) = Node(fault_nodes[0][i],0);
    }
    
    for (int k=0 ; k<nx_faults[0]; k++)
    {
        if ((x_main(k)<=(0.0+0.25e3))&&(x_main(k)>=(0.0-0.25e3)))
        {
            T_0[0](2*k) = 31.0e6;
        }
        else
        {
            T_0[0](2*k) = 20.0e6;
        }
        
        if (std::abs(x_main(k))<=5.0e3)
        {
            mu_s[0](k)= 0.6;
        }
        else
        {
            mu_s[0](k) = 0.6;
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
    double alpha = 0.4;
    //dx_min = 100;
    double dt = alpha*dx_min/v_p;
    std::cout<<"dt="<<dt<<std::endl;
    //dt = 0.001;
     dt = 5.0e-4;
    // Reyleigh Damping
    double beta =0.1;
    double q = beta*dt;
    double time_run = 4.0;
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
    
    std::ofstream M_pos_output("results/M_pos.txt");
    std::ofstream M_neg_output("results/M_neg.txt");

    //    std::ofstream shear_x("results/shear.txt");
    //double start = omp_get_wtime();
    std::vector<MatrixXd> strain_n_store(n_el);
    std::vector<MatrixXd> stress_n_store(n_el);
    
    for (int i=0; i<n_el; i++)
    {
        strain_n_store[i] = MatrixXd::Zero(3, 4);
        stress_n_store[i] = MatrixXd::Zero(3, 4);
        
    }
    double cohes = 20e6 ;
    double blkfric = 0.8;
    double sxx_initial = sigma_xx;
    double syy_initial = sigma_yy;
    double sxy_initial = -tau_xy;
    
    // Equivalent plastic strain
    std::vector<double> eq_ep_out(n_el);
    // Total Energy Dissipation
    double E_p_out;
    // Total Frictional Dissipation
    double E_fric_sec=0.0;
    double E_fric_main=0.0;
    double avg_slip = 0.0;
    double max_slip = 0.0;
    std::ofstream E_p_output("results/E_p.txt");
    std::ofstream E_fric_sec_output("results/E_fric_sec.txt");
    std::ofstream E_fric_main_output("results/E_fric_main.txt");
    std::ofstream avg_slip_main_output("results/avg_slip_main.txt");
    std::ofstream max_slip_rate_main_output("results/max_slip_rate_main.txt");

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
            
            Slip_Weakening_lumpM_with_buffer(M_global_vec, Fault_surf_index[2*i], Fault_surf_index[2*i+1], fe_global, dt, dx, dy, nx_faults[i]-1, delt_v_n[i], delt_u_n[i], T_0[i], tau_s[i], mu_s[i], mu_d, Dc, Ndofn, M_pos, M_neg , F_fault, T_c[i], fault_angle[i]);
            mapglobal(Fault_surf_index[2*i],F_total,-F_fault);
            mapglobal(Fault_surf_index[2*i+1],F_total,F_fault);
        }
        
        // Calculate Absorbing force
      //  cal_ab_force(T_abc, F_total, AB_left_surf_index, v_n, Ndofn);
      //  cal_ab_force(T_abc, F_total, AB_right_surf_index, v_n, Ndofn);
       // Central Difference Time integration
        time_advance(u_n, v_n, F_total, M_global_vec, dt);
        // Get a copy of the previous delt_u and delt_v
        for (int i=0; i<num_faults;i++)
        {
            delt_u_n_copy[i]=delt_u_n[i];
            delt_v_n_copy[i]=delt_v_n[i];
        }
        // Get the slip and slip rate
        for (int i=0; i<num_faults; i++)
        {
            
            cal_slip_slip_rate_angle(u_n, v_n,  Fault_surf_index[2*i], Fault_surf_index[2*i+1] , Ndofn, nx_faults[i]-1, delt_u_n[i], delt_v_n[i],fault_angle[i]);
        }
        delt_v_n[0](0)=0;
        delt_v_n[0](2*nx_faults[0]-2)=0;
        delt_u_n[0](0)=0;
        delt_u_n[0](2*nx_faults[0]-2)=0;
        // Correct the BIE surf nodes solutions from FEM with the BIE solution
        BIE_correct(BIE_top_surf_index, BIE_bot_surf_index, fe_global, Ndofn, nx_BIE-1, dx, BIE_inf_top, BIE_inf_bot, u_n, v_n);

        // Calculating the frictional disspation for the seondary faults and main faults
        // Secondary faults disspaltions
        for (int i=1;i<num_faults;i++)
        {
            VectorXd E_fric_sec_temp = VectorXd::Zero(nx_faults[i],1);
            for (int k=0; k<nx_faults[i];k++)
            {
                E_fric_sec_temp(k) = (delt_u_n[i](2*k)-delt_u_n_copy[i](2*k))*T_c[i](2*k);
            }
            E_fric_sec +=E_fric_sec_temp.sum();
        }
        
        E_fric_sec_output<<E_fric_sec<<std::endl;

        for (int i=0;i<1;i++)
        {
            VectorXd E_fric_main_temp = VectorXd::Zero(nx_faults[i],1);
            for (int k=0; k<nx_faults[i];k++)
            {
                E_fric_main_temp(k) = (delt_u_n[i](2*k)-delt_u_n_copy[i](2*k))*T_c[i](2*k);            }
            E_fric_main +=E_fric_main_temp.sum();
            
            VectorXd delt_u_n_temp = VectorXd::Zero(nx_faults[i],1);
            VectorXd delt_v_n_temp = VectorXd::Zero(nx_faults[i],1);

            for (int k=0; k<nx_faults[i];k++)
            {
                delt_u_n_temp(k) = delt_u_n[i](2*k);
                delt_v_n_temp(k) = delt_v_n[i](2*k);

            }
            avg_slip = delt_u_n_temp.sum()/nx_faults[i];
            max_slip = delt_v_n_temp.maxCoeff();

        }
        E_fric_main_output<<E_fric_main<<std::endl;
        avg_slip_main_output<<avg_slip<<std::endl;
        max_slip_rate_main_output<<max_slip<<std::endl;

        
        for (int i=0;i<2;i++)
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
//        if (j%10==1)
//        {
//            file.open("results/eq_ep_n.bin",ios::binary | ios::app);
//            file.write((char*)(eq_ep_out.data()),eq_ep_out.size()*sizeof(double));
//            file.close();
//        }
        printf("Simulation time = %f\n",time(j));
      //  double end_t = omp_get_wtime();
       // std::cout<<"time_cpu_t="<<end_t-start<<std::endl;
    //   std::cout<<"time_fem="<< time_fem<<std::endl;
    //   std::cout<<"time_bie="<< time_bie<<std::endl;21
        E_p_output<<E_p_out<<std::endl;
    }

    //double end = omp_get_wtime();
    //std::cout<<"time_cpu="<<end-start<<std::endl;
    return 0;
}
