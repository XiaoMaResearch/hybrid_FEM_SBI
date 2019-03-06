#include <iostream>
#include <fstream>
#include <Eigen/Eigen>
#include "exodusII.h"
#include "share_header.hpp"
#include "mesh_Generated_multi_faults.hpp"
#include "Meshrectangular.hpp"
using namespace Eigen;
using namespace std;

int main()
{
    // Domain Size
    double x_min = -20e3;
    double x_max = 20e3;
    double y_min = -4.0e3;
    double y_max = 4.0e3;
    int dim = 2.0;
    double dx = 100;
    double dy = 100;
    int nx = (x_max-x_min)/dx;
    int ny = (y_max-y_min)/dy;
    MatrixXd Node = MatrixXd::Zero((nx+1)*(ny+1),2);
    MatrixXi_rm Element(nx*ny,4); Element.setZero();
    ArrayXi BIE_top_surf_nodes = ArrayXi::Zero((nx+1),1);
    ArrayXi BIE_bot_surf_nodes = ArrayXi::Zero((nx+1),1);
    int num_faults = 3;
    // Mesh
    VectorXd fault_pos_y(num_faults);
    fault_pos_y << 0.0 ,1.0e3, -1.0e3;
    std::vector<std::vector<int>> fault_nodes(num_faults*2);
  //  Meshrectangular(x_min, x_max, y_min, y_max, dx, dy, nx, ny, Node, Element);
    mesh_Generated_multi_faults(x_min,x_max,y_min,y_max,dx,dy,nx,ny,Node, Element,BIE_top_surf_nodes, BIE_bot_surf_nodes,fault_pos_y, fault_nodes);


    
    int exoid, num_dim, num_nodes, num_elem, num_elem_blk;
    int num_elem_in_block[10], num_nodes_per_elem[10];
    int num_node_sets, num_side_sets, error;
    int ebids[10];
    int CPU_word_size,IO_word_size;
   // float x[100], y[100], z[100], *dummy;
    double *dummy;
    dummy = 0; /* assign this so the Cray compiler doesn’t complain */
    /* Specify compute and i/o word size */
    CPU_word_size = 8;/* float or double */
    IO_word_size = 8;/* use system default (4 bytes) */
    /* create EXODUS II file */
    exoid = ex_create ("test.exo",/* filename path */
                       EX_CLOBBER,/* create mode */
                       &CPU_word_size,/* CPU float word size in bytes */
                       &IO_word_size);/* I/O float word size in bytes */
    /* ncopts = NC_VERBOSE; */
    /* initialize file with parameters */
    num_dim = 2;
    num_nodes = Node.rows();
    num_elem = Element.rows();
    num_elem_blk = 1;
    num_node_sets = 0;
    num_side_sets = 0;
    error = ex_put_init (exoid, "This is a test", num_dim, num_nodes, num_elem,
                         num_elem_blk, num_node_sets, num_side_sets);
    /* write nodal coordinates values and names to database */
    double *z ;
    error = ex_put_coord (exoid, Node.col(0).data(), Node.col(1).data(), z);
    /* write element block parameters */
    num_elem_in_block[0] = num_elem;
    num_nodes_per_elem[0] = 4; /* elements in block #1 are 4-node quads  */
    ebids[0] = 1;
    error = ex_put_elem_block (exoid, ebids[0], "QUAD", num_elem_in_block[0],
                               num_nodes_per_elem[0], 1);
    Element.array() += 1;
//    /* write element connectivity */
    error = ex_put_elem_conn (exoid, ebids[0], Element.data());
    /* for each time step, write the analysis results;
     * the code below fills the arrays glob_var_vals,
     * nodal_var_vals, and elem_var_vals with values for debugging purposes;
     * obviously the analysis code will populate these arrays
     */
    /* write results variables parameters and names */
    char *var_names[3];
    var_names[0] = "nod_var0";
    int num_nod_vars;
    num_nod_vars = 1;
    error = ex_put_var_param (exoid, "n", num_nod_vars);
    error = ex_put_var_names (exoid, "n", num_nod_vars, var_names);
    /* for each time step, write the analysis results;
     * the code below fills the arrays glob_var_vals,
     * nodal_var_vals, and elem_var_vals with values for debugging purposes;
     * obviously the analysis code will populate these arrays
     */
    int whole_time_step, num_time_steps;
    double time_value;
    double *nodal_var_vals;
    whole_time_step = 1;
    num_time_steps = 10;
    nodal_var_vals = (double *) calloc (num_nodes, CPU_word_size);
    for (int i=0; i<num_time_steps; i++)
    {
        time_value = (double)(i+1)/100.;
        /* write time value */
        error = ex_put_time (exoid, whole_time_step, &time_value);
        /* write nodal variables */
        for (int k=1; k<=num_nod_vars; k++)
        {
            for (int j=0; j<num_nodes; j++)
            {
                nodal_var_vals[j] = (double)k + ((double)(j+1) * time_value*100);
            }
            error = ex_put_nodal_var (exoid, whole_time_step, k, num_nodes,
                                      nodal_var_vals);
        }
        whole_time_step++;
        /* update the data file; this should be done at the end of every time step
         * to ensure that no data is lost if the analysis dies
         */
        error = ex_update (exoid);
    }
    /* close the EXODUS files */
    free(nodal_var_vals);
    error = ex_close (exoid);
}
