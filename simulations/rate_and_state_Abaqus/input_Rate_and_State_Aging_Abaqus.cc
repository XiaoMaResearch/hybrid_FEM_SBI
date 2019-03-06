#include <iostream>
#include <fstream>
#include <Eigen/Eigen>
//include the bie header files
#include "material.hh"
#include "precomputed_kernel.hh"
#include "bimat_interface.hh"
#include "infinite_boundary.hh"
//include the fem header files
//#include "mesh_Generated.hpp"
#include "mesh_Generated_multi_faults.hpp"

#include "bcdof.hpp"
#include "bcdof_ptr.hpp"
#include "cal_ke.hpp"
#include "cal_fe_global_const_ke.hpp"
#include "mapglobal.hpp"
#include "Slip_Weakening.hpp"
#include "Rate_and_State_Aging.hpp"
#include "cal_slip_sliprate.hpp"
#include "update_disp_velocity.hpp"
#include "time_advance.hpp"
#include "BIE_correct.hpp"


using namespace Eigen;
using namespace std;

int main() {
    // Domain Size
    double x_min = -10e3;
    double x_max = 10e3;
    double y_min = -1.0e3;
    double y_max = 1.0e3;
    int dim = 2.0;
    double dx = 100;
    double dy = 100;
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
    
     double Vw_width = 20e3;
    // Finding Nodes on the diagonal
    ArrayXi left_diag_nodes =ArrayXi::Zero((ny+1),1);
    ArrayXi right_diag_nodes =ArrayXi::Zero((ny+1),1); ;
    double TOL = 1e-6;
    int m = 0;
    int n = 0;
    for(int i=0; i<Node.rows(); i++)
    {
        if ((std::abs(Node(i,0) +Vw_width/2.0)<TOL)&&(std::abs(Node(i,1))>=dx))
        {
            left_diag_nodes(m) = i ;
            m+=1;
        };
        if ((std::abs(Node(i,0) - Vw_width/2.0)<TOL)&&(std::abs(Node(i,1))>=dx))
        {
            right_diag_nodes(n) = i ;
            n+=1;
        };
    }
    VectorXi left_diag_index = VectorXi::Zero(2*(left_diag_nodes.size()),1);
    VectorXi right_diag_index = VectorXi::Zero(2*(right_diag_nodes.size()),1);

    bcdof(left_diag_nodes,dim,left_diag_index);
    bcdof(right_diag_nodes,dim,right_diag_index);
    
//    std::cout<<left_diag_nodes<<std::endl;
//    std::cout<<"********"<<std::endl;
//    std::cout<<right_diag_nodes<<std::endl;

    std::ofstream Element_output("results/Element.txt");
    std::ofstream Node_output("results/Node.txt");
    Node_output<<Node;
    Element_output<<Element;
    
    Element_output.close();
    Node_output.close();
    // Vector containing number of elements on each faults (nx)
    std::vector<int> nx_faults(num_faults);
    for (int i=0;i<num_faults;i++)
    {
        nx_faults[i] =fault_nodes[2*i].size();
    }
    
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
    }
    int n_nodes = Node.rows();
    int n_el = Element.rows();
    int Ndofn = 2;
    int Nnel = Element.cols();
   // nx_el = fault_nodes.size()-1;
    // Material
    double density = 2670.0;
    double v_s =3.464e3;
    double v_p = 6.0e3;
    double G= pow(v_s,2)*density;
    double Lambda = pow(v_p,2)*density-2.0*G;
    double E  = G*(3.0*Lambda+2.0*G)/(Lambda+G);
    double nu = Lambda/(2.0*(Lambda+G));
    // Time
    double alpha = 0.1;
    double dt = alpha*dx/v_p;
    // Reyleigh Damping 
    double beta =0.2;
    double q = beta*dt;
    double time_run = 9.0;
    int numt = time_run/dt;
    //numt =1;
    
    
    std::ofstream time_output("results/time.txt");
    time_output<<time_run<<"\n";
    time_output<<dt<<"\n";
    time_output<<numt<<std::endl;
    //numt = 3;
    VectorXd time = dt*VectorXd::LinSpaced(numt,1,numt);
    // Rate and state friction parameters
    
    double a = 0.008;
    double b = 0.012;
    
    VectorXd a_array = a*VectorXd::Zero(nx_el+1,1);
    VectorXd b_array = b*VectorXd::Zero(nx_el+1,1);

    double V_0 = 1e-6;
    double f_0 = 0.6;
    double L = 0.02;
    double tau_ini = 75e6;
    double sigma_ini = 120e6;
    double V_ini  = 1e-12;
    // Intialization
    // disp velocity current and next time step (new)
    VectorXd u_n = VectorXd::Zero(n_nodes*Ndofn,1);
    VectorXd v_n = VectorXd::Zero(n_nodes*Ndofn,1);
    VectorXd u_new = VectorXd::Zero(n_nodes*Ndofn,1);
    VectorXd v_new = VectorXd::Zero(n_nodes*Ndofn,1);
    VectorXd a_n = VectorXd::Zero(n_nodes*Ndofn,1);
    // slip and slip-rate
    VectorXd delt_u_n = VectorXd::Zero(Ndofn*(nx_el+1),1);
    VectorXd delt_v_n = VectorXd::Zero(Ndofn*(nx_el+1),1);
    std::cout<<delt_u_n.size()<<std::endl;
    std::cout<<nx_el<<std::endl;
    // Stress on the fault T_0 = intial stress , T = sticking force, T_c= stress critical goes into F_global
    VectorXd T_0 = VectorXd::Zero(Ndofn*(nx_el+1),1);
    VectorXd T = VectorXd::Zero(Ndofn*(nx_el+1),1);
    VectorXd T_c = VectorXd::Zero(Ndofn*(nx_el+1),1);
    VectorXd tau_s = VectorXd::Zero(nx_el+1,1);
    VectorXd F_ext_global = VectorXd::Zero(n_nodes*Ndofn,1);
    VectorXd theta_n =VectorXd::Zero(nx_el+1, 1);
    VectorXd theta_dot_n = VectorXd::Zero(nx_el+1, 1);
    VectorXd theta_ini =VectorXd::Zero(nx_el+1, 1);

    // Setting intial stress on the fault
    for (int i=0; i<T_0.size()/2; i++)
    {
        T_0(2*i) = tau_ini;
        T_0(2*i+1) = -sigma_ini;
    }
    VectorXd x = VectorXd::LinSpaced(nx_el+1,x_min,x_max);
    double d_tau_0 = 25e6;
    double R = 2.5e3;
    VectorXd F_r = VectorXd::Zero(nx_el+1,1);
    for (int i=0; i<nx_el+1; i++)
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
   
    for (int i =0; i<nx_el+1; i++){
        if (abs(x(i))<Vw_width/2.0)
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
    MatrixXd d_tau = MatrixXd::Zero(nx_el+1,numt);
    for (int i=0;i<numt;i++)
    {
        d_tau.col(i)=75e6*VectorXd::Ones(nx_el+1,1)+d_tau_0*F_r*(G_t(i));
       // d_tau.col(i)=75e6*VectorXd::Ones(nx_el+1,1)+3.75e6*F_r+dt*0.3125e7*i*VectorXd::Ones(nx_el+1,1);
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
    
    // Rate and state initlaized v_n delt_v_n
    for (int i=0;i<nx_el+1;i++)
    {
        delt_v_n(2*i)= V_ini;
  //      v_n(top_surf_index(2*i))= V_ini/2.0;
        
  //      v_n(bot_surf_index(2*i))= -V_ini/2.0;
        
        v_n(Fault_surf_index[0](2*i))= V_ini/2.0;
        
        v_n(Fault_surf_index[1](2*i))= -V_ini/2.0;
    }
    
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

    //BIE_initiation(E, nu, density, x_max, x_min, nx_el, dt);
    printf("ready to start\n");
    // Output
    ofstream file;
    file.open("results/num_nodes_fault.bin",ios::binary);
    file.write((char*)(nx_faults.data()),nx_faults.size()*sizeof(int));
    file.close();
    //
    file.open("results/u_n.bin",ios::binary);
    file.close();
    file.open("results/v_n.bin",ios::binary);
    file.close();
    file.open("results/eq_ep_n.bin",ios::binary);
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

    // Main time loop
    for (int j=0;j<numt;j++)
    {
        // Compute the global internal force
        VectorXd fe_global= VectorXd::Zero(n_nodes*Ndofn,1);
        cal_fe_global_const_ke(n_nodes, n_el, index_store, q, u_n, v_n, Ndofn, ke, fe_global);
        // Friction subroutine
        VectorXd F_fault = VectorXd::Zero(Ndofn*(nx_el+1),1);
      //  Slip_Weakening(M_global_vec, top_surf_index, bot_surf_index, fe_global, dt, dx, dy, nx_el, delt_v_n, delt_u_n, T_0, tau_s, mu_s, mu_d, Dc, Ndofn, M, F_fault, T_c);
        for (int i=0;i<nx_el+1;i++)
        {
            T_0(2*i) = d_tau(i,j);
        }
        Rate_and_State_Aging(M_global_vec, Fault_surf_index[0], Fault_surf_index[1], fe_global, dt, dx, dy, nx_el, delt_u_n, delt_v_n, T_0, tau_s, a_array, b_array , L, V_0, f_0, Ndofn, M, F_fault, T_c, theta_n, theta_dot_n);

        
        // Calculate the global force vector
        // Adding contribution of the fault force to the global force vector F_total
        VectorXd F_total = F_ext_global-fe_global;
        mapglobal(Fault_surf_index[0],F_total,-F_fault);
        mapglobal(Fault_surf_index[1],F_total,F_fault);
        // Central Difference Time integration
        time_advance(u_n, v_n, F_total, M_global_vec, dt);
        
        for (int i=0 ; i<left_diag_index.size();i++)
        {
            u_n(left_diag_index(i)) = 0.0;
            v_n(left_diag_index(i)) = 0.0;
            u_n(right_diag_index(i)) = 0.0;
            v_n(right_diag_index(i)) = 0.0;
        }
        
        
        // Get the slip and slip rate
        //cal_slip_slip_rate(u_n, v_n, top_surf_index, bot_surf_index, Ndofn, nx_el, delt_u_n, delt_v_n);
        update_disp_velocity(u_n, v_n, Fault_surf_index[0], Fault_surf_index[1], Ndofn, nx_el, delt_u_n, delt_v_n);

        // Correct the BIE surf nodes solutions from FEM with the BIE solution
        BIE_correct(BIE_top_surf_index, BIE_bot_surf_index, fe_global, Ndofn, nx_el, dx, BIE_inf_top, BIE_inf_bot, u_n, v_n);
        
        // Set the diagonal to be fixed
        
 
        //ofstream file;
        for (int i=0;i<1;i++)
        {
            std::string slip = "results/slip_"+std::to_string(i)+".bin";
            file.open(slip,ios::binary | ios::app);
            file.write((char*)(delt_u_n.data()),delt_u_n.size()*sizeof(double));
            file.close();
            
            std::string slip_rate = "results/slip_rate_"+std::to_string(i)+".bin";
            file.open(slip_rate,ios::binary | ios::app);
            file.write((char*)(delt_v_n.data()),delt_v_n.size()*sizeof(double));
            file.close();
            
            std::string shear = "results/shear_"+std::to_string(i)+".bin";
            file.open(shear,ios::binary | ios::app);
            file.write((char*)(T_c.data()),T_c.size()*sizeof(double));
            file.close();
        }
        
        if (j%100==1)
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
