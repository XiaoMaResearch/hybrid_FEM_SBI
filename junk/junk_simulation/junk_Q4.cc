#include <iostream>
#include <fstream>
#include <Eigen/Eigen>
#include "exodusII.h"

using namespace Eigen;
using namespace std;

int main() {
    
    int exoid, num_dim, num_nodes, num_elem, num_elem_blk;
    int num_elem_in_block[10], num_nodes_per_elem[10];
    int num_node_sets, num_side_sets, error;
    int i, j, k, m, *connect;
    int ebids[10];
    int num_glo_vars, num_nod_vars, num_ele_vars;
    int *truth_tab;
    int whole_time_step, num_time_steps;
    int CPU_word_size,IO_word_size;
    float *glob_var_vals, *nodal_var_vals, *elem_var_vals;
    float time_value;
    float x[100], y[100], z[100], *dummy;
    char *coord_names[3], *qa_record[2][4], *info[3], *var_names[3];
    char tmpstr[80];
    char *prop_names[2];
    dummy = 0; /* assign this so the Cray compiler doesn’t complain */
    /* Specify compute and i/o word size */
    CPU_word_size = 0;/* float or double */
    IO_word_size = 0;/* use system default (4 bytes) */
    /* create EXODUS II file */
    exoid = ex_create ("test.exo",/* filename path */
                       EX_CLOBBER,/* create mode */
                       &CPU_word_size,/* CPU float word size in bytes */
                       &IO_word_size);/* I/O float word size in bytes */
    /* ncopts = NC_VERBOSE; */
    /* initialize file with parameters */
    num_dim = 2;
    num_nodes = 9;
    num_elem = 4;
    num_elem_blk = 1;
    num_node_sets = 0;
    num_side_sets = 0;
    error = ex_put_init (exoid, "This is a test", num_dim, num_nodes, num_elem,
                         num_elem_blk, num_node_sets, num_side_sets);
    /* write nodal coordinates values and names to database */
    /* Quad #1 */
    x[0] = 0.0; y[0] = 0.0;
    x[1] = 1.0; y[1] = 0.0;
    x[2] = 2.0; y[2] = 0.0;
    x[3] = 0.0; y[3] = 1.0;
    x[4] = 1.0; y[4] = 1.0;
    x[5] = 2.0; y[5] = 1.0;
    x[6] = 0.0; y[6] = 2.0;
    x[7] = 1.0; y[7] = 2.0;
    x[8] = 2.0; y[8] = 2.0;


    error = ex_put_coord (exoid, x, y, z);
 //   error = ex_put_coord_names (exoid, coord_names);
    /* write element order map */
//    elem_map = (int *) calloc(num_elem, sizeof(int));
//    for (i=1; i<=num_elem; i++)
//    {
//        elem_map[i-1] = i;
//    }
//    error = ex_put_map (exoid, elem_map);
//    free (elem_map);
    /* write element block parameters */
    num_elem_in_block[0] = 4;
    num_nodes_per_elem[0] = 4; /* elements in block #1 are 4-node quads  */
    ebids[0] = 1;
    error = ex_put_elem_block (exoid, ebids[0], "QUAD", num_elem_in_block[0],
                               num_nodes_per_elem[0], 1);
    /* write element connectivity */
    connect = (int *) calloc(16, sizeof(int));
    connect[0] = 1; connect[1] = 2; connect[2] = 5; connect[3] = 4;
    connect[4] = 2; connect[5] = 3; connect[6] = 6; connect[7] = 5;
    connect[8] = 4; connect[9] = 5; connect[10] = 8; connect[11] = 7;
    connect[12] = 5; connect[13] = 6; connect[14] = 9; connect[15] = 8;

    error = ex_put_elem_conn (exoid, ebids[0], connect);
    free (connect);
    
    float *glob_var_vals, *nodal_var_vals, *elem_var_vals;

//      /* for each time step, write the analysis results;
//     * the code below fills the arrays glob_var_vals,
//     * nodal_var_vals, and elem_var_vals with values for debugging purposes;
//     * obviously the analysis code will populate these arrays
//     */
//    /* write results variables parameters and names */
//    num_glo_vars = 1;
//    var_names[0] = "glo_vars";
//    error = ex_put_var_param (exoid, "g", num_glo_vars);
//    error = ex_put_var_names (exoid, "g", num_glo_vars, var_names);
//    num_nod_vars = 2;
//    var_names[0] = "nod_var0";
//    var_names[1] = "nod_var1";
//    error = ex_put_var_param (exoid, "n", num_nod_vars);
//    error = ex_put_var_names (exoid, "n", num_nod_vars, var_names);
//    num_ele_vars = 3;
//    var_names[0] = "ele_var0";
//    var_names[1] = "ele_var1";
//    var_names[2] = "ele_var2";
//    error = ex_put_var_param (exoid, "e", num_ele_vars);
//    error = ex_put_var_names (exoid, "e", num_ele_vars, var_names);
//
//// Time variables write out
//    whole_time_step = 1;
//    num_time_steps = 10;
//    glob_var_vals = (float *) calloc (num_glo_vars, CPU_word_size);
//    nodal_var_vals = (float *) calloc (num_nodes, CPU_word_size);
//    elem_var_vals = (float *) calloc (4, CPU_word_size);
//    for (i=0; i<num_time_steps; i++)
//    {
//        time_value = (float)(i+1)/100.;
//        /* write time value */
//        error = ex_put_time (exoid, whole_time_step, &time_value);
//        /* write global variables */
//        for (j=0; j<num_glo_vars; j++)
//        {
//            glob_var_vals[j] = (float)(j+2) * time_value;
//        }
//        error = ex_put_glob_vars (exoid, whole_time_step, num_glo_vars,
//                                  glob_var_vals);
//        /* write nodal variables */
//        for (k=1; k<=num_nod_vars; k++)
//        {
//            for (j=0; j<num_nodes; j++)
//            {
//                nodal_var_vals[j] = (float)k + ((float)(j+1) * time_value);
//            }
//            error = ex_put_nodal_var (exoid, whole_time_step, k, num_nodes,
//                                      nodal_var_vals);
//        }
//        /* write element variables */
//        for (k=1; k<=num_ele_vars; k++)
//        {
//            for (j=0; j<num_elem_blk; j++)
//            {
//                for (m=0; m<num_elem_in_block[j]; m++)
//                {
//                    elem_var_vals[m] = (float)(k+1) + (float)(j+2) +
//                    ((float)(m+1)*time_value);
//                }
//                error = ex_put_elem_var (exoid, whole_time_step, k, ebids[j],num_elem_in_block[j], elem_var_vals);
//            }
//            whole_time_step++;
//        }
//                /* update the data file; this should be done at the end of every time step
//                 * to ensure that no data is lost if the analysis dies
//                 */
//                error = ex_update (exoid);
//            }
//            free(glob_var_vals);
//            free(nodal_var_vals);
//            free(elem_var_vals);
//
//            /* close the EXODUS files
//             */
            error = ex_close (exoid);
        }

