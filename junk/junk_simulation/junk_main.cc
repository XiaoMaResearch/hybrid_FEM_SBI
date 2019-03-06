#include <iostream>
#include <fstream>
#include <Eigen/Eigen>
#include "exodusII.h"

using namespace Eigen;
using namespace std;

int main() {
    
    int exoid, num_dim, num_nodes, num_elem, num_elem_blk;
    int num_elem_in_block[10], num_nodes_per_elem[10];
    int num_node_sets, num_sides, num_side_sets, error;
    int i, j, k, m, *elem_map, *connect;
    int node_list[100],elem_list[100],side_list[100];
    int ebids[10], ids[10];
    int num_sides_per_set[10], num_nodes_per_set[10], num_elem_per_set[10];
    int num_df_per_set[10];
    int df_ind[10], node_ind[10], elem_ind[10], side_ind[10];
    int  num_qa_rec, num_info;
    int num_glo_vars, num_nod_vars, num_ele_vars;
    int *truth_tab;
    int whole_time_step, num_time_steps;
    int ndims, nvars, ngatts, recdim;
    int CPU_word_size,IO_word_size;
    int prop_array[2];
    float *glob_var_vals, *nodal_var_vals, *elem_var_vals;
    float time_value;
    float x[100], y[100], z[100], *dummy;
    float attrib[1], dist_fact[100];
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
    num_dim = 3;
    num_nodes = 26;
    num_elem = 5;
    num_elem_blk = 5;
    num_node_sets = 2;
    num_side_sets = 5;
    error = ex_put_init (exoid, "This is a test", num_dim, num_nodes, num_elem,
                         num_elem_blk, num_node_sets, num_side_sets);
    /* write nodal coordinates values and names to database */
    /* Quad #1 */
    x[0] = 0.0; y[0] = 0.0; z[0] = 0.0;
    x[1] = 1.0; y[1] = 0.0; z[1] = 0.0;
    x[2] = 1.0; y[2] = 1.0; z[2] = 0.0;
    x[3] = 0.0; y[3] = 1.0; z[3] = 0.0;
    /* Quad #2 */
    x[4]  =  1.0; y[4]  =  0.0; z[4]  =  0.0;
    x[5]  =  2.0; y[5]  =  0.0; z[5]  =  0.0;
    x[6]  =  2.0; y[6]  =  1.0; z[6]  =  0.0;
    x[7]  =  1.0; y[7]  =  1.0; z[7]  =  0.0;
    /* Hex #1 */
    x[8]  =  0.0; y[8]  =  0.0; z[8]  =  0.0;
    x[9]  = 10.0; y[9]  =  0.0; z[9]  =  0.0;
    x[10] = 10.0; y[10] =  0.0; z[10] =-10.0;
    x[11] =  1.0; y[11] =  0.0; z[11] =-10.0;
    x[12] =  1.0; y[12] = 10.0; z[12] =  0.0;
    x[13] = 10.0; y[13] = 10.0; z[13] =  0.0;
    x[14] = 10.0; y[14] = 10.0; z[14] =-10.0;
    x[15] =  1.0; y[15] = 10.0; z[15] =-10.0;
    /* Tetra #1 */
    x[16] =  0.0; y[16] =  0.0; z[16] =  0.0;
    x[17] =  1.0; y[17] =  0.0; z[17] =  5.0;
    x[18] = 10.0; y[18] =  0.0; z[18] =  2.0;
    x[19] =  7.0; y[19] =  5.0; z[19] =  3.0;
    /* Wedge #1 */
    x[20] =  3.0; y[20] =  0.0; z[20] =  6.0;
    x[21] =  6.0; y[21] =  0.0; z[21] =  0.0;
    x[22] =  0.0; y[22] =  0.0; z[22] =  0.0;
    x[23] =  3.0; y[23] =  2.0; z[23] =  6.0;
    x[24] =  6.0; y[24] =  2.0; z[24] =  2.0;
    x[25] =  0.0; y[25] =  2.0; z[25] =  0.0;
    error = ex_put_coord (exoid, x, y, z);
    coord_names[0] = "xcoor";
    coord_names[1] = "ycoor";
    coord_names[2] = "zcoor";
    error = ex_put_coord_names (exoid, coord_names);
    /* write element order map */
    elem_map = (int *) calloc(num_elem, sizeof(int));
    for (i=1; i<=num_elem; i++)
    {
        elem_map[i-1] = i;
    }
    error = ex_put_map (exoid, elem_map);
    free (elem_map);
    /* write element block parameters */
    num_elem_in_block[0] = 1;
    num_elem_in_block[1] = 1;
    num_elem_in_block[2] = 1;
    num_elem_in_block[3] = 1;
    num_elem_in_block[4] = 1;
    num_nodes_per_elem[0] = 4; /* elements in block #1 are 4-node quads  */
    num_nodes_per_elem[1] = 4; /* elements in block #2 are 4-node quads  */
    num_nodes_per_elem[2] = 8; /* elements in block #3 are 8-node hexes  */
    num_nodes_per_elem[3] = 4; /* elements in block #3 are 4-node tetras */
    num_nodes_per_elem[4] = 6; /* elements in block #3 are 6-node wedges */
    ebids[0] = 10;
    ebids[1] = 11;
    ebids[2] = 12;
    ebids[3] = 13;
    ebids[4] = 14;
    error = ex_put_elem_block (exoid, ebids[0], "QUAD", num_elem_in_block[0],
                               num_nodes_per_elem[0], 1);
    error = ex_put_elem_block (exoid, ebids[1], "QUAD", num_elem_in_block[1],
                               num_nodes_per_elem[1], 1);
    error = ex_put_elem_block (exoid, ebids[2], "HEX", num_elem_in_block[2],
                               num_nodes_per_elem[2], 1);
    error = ex_put_elem_block (exoid, ebids[3], "TETRA", num_elem_in_block[3],
                               num_nodes_per_elem[3], 1);
    error = ex_put_elem_block (exoid, ebids[4], "WEDGE", num_elem_in_block[4],
                               num_nodes_per_elem[4], 1);
    /* write element block properties */
    prop_names[0] = "TOP";
    prop_names[1] = "RIGHT";
    error = ex_put_prop_names(exoid,EX_ELEM_BLOCK,2,prop_names);
    error = ex_put_prop(exoid, EX_ELEM_BLOCK, ebids[0], "TOP", 1);
    error = ex_put_prop(exoid, EX_ELEM_BLOCK, ebids[1], "TOP", 1);
    error = ex_put_prop(exoid, EX_ELEM_BLOCK, ebids[2], "RIGHT", 1);
    error = ex_put_prop(exoid, EX_ELEM_BLOCK, ebids[3], "RIGHT", 1);
    error = ex_put_prop(exoid, EX_ELEM_BLOCK, ebids[4], "RIGHT", 1);
    /* write element connectivity */
    connect = (int *) calloc(8, sizeof(int));
    connect[0] = 1; connect[1] = 2; connect[2] = 3; connect[3] = 4;
    error = ex_put_elem_conn (exoid, ebids[0], connect);
    connect[0] = 5; connect[1] = 6; connect[2] = 7; connect[3] = 8;
    error = ex_put_elem_conn (exoid, ebids[1], connect);
    connect[0] = 9; connect[1] = 10; connect[2] = 11; connect[3] = 12;
    connect[4] = 13; connect[5] = 14; connect[6] = 15; connect[7] = 16;
    error = ex_put_elem_conn (exoid, ebids[2], connect);
    connect[0] = 17; connect[1] = 18; connect[2] = 19; connect[3] = 20;
    error = ex_put_elem_conn (exoid, ebids[3], connect);
    connect[0] = 21; connect[1] = 22; connect[2] = 23;
    connect[3] = 24; connect[4] = 25; connect[5] = 26;
    error = ex_put_elem_conn (exoid, ebids[4], connect);
    free (connect);
    /* write element block attributes */
    attrib[0] = 3.14159;
    error = ex_put_elem_attr (exoid, ebids[0], attrib);
    attrib[0] = 6.14159;
    error = ex_put_elem_attr (exoid, ebids[1], attrib);
    error = ex_put_elem_attr (exoid, ebids[2], attrib);
    error = ex_put_elem_attr (exoid, ebids[3], attrib);
    error = ex_put_elem_attr (exoid, ebids[4], attrib);
    /* write individual node sets */
    error = ex_put_node_set_param (exoid, 20, 5, 5);
    node_list[0] = 100; node_list[1] = 101; node_list[2] = 102;
    node_list[3] = 103; node_list[4] = 104;
    dist_fact[0] = 1.0; dist_fact[1] = 2.0; dist_fact[2] = 3.0;
    dist_fact[3] = 4.0; dist_fact[4] = 5.0;
    error = ex_put_node_set (exoid, 20, node_list);
    error = ex_put_node_set_dist_fact (exoid, 20, dist_fact);
    error = ex_put_node_set_param (exoid, 21, 3, 3);
    node_list[0] = 200; node_list[1] = 201; node_list[2] = 202;
    dist_fact[0] = 1.1; dist_fact[1] = 2.1; dist_fact[2] = 3.1;
    error = ex_put_node_set (exoid, 21, node_list);
    error = ex_put_node_set_dist_fact (exoid, 21, dist_fact);
    error = ex_put_prop(exoid, EX_NODE_SET, 20, "FACE", 4);
    error = ex_put_prop(exoid, EX_NODE_SET, 21, "FACE", 5);
    prop_array[0] = 1000;
    prop_array[1] = 2000;
    error = ex_put_prop_array(exoid, EX_NODE_SET, "VELOCITY", prop_array);
    /* write concatenated node sets; this produces the same information as
     * the above code which writes individual node sets
     */
    /* THIS SECTION IS COMMENTED OUT
     ids[0] = 20; ids[1] = 21;
     num_nodes_per_set[0] = 5; num_nodes_per_set[1] = 3;
     node_ind[0] = 0; node_ind[1] = 5;
     node_list[0] = 100; node_list[1] = 101; node_list[2] = 102;
     node_list[3] = 103; node_list[4] = 104;
     node_list[5] = 200; node_list[6] = 201; node_list[7] = 202;
     num_df_per_set[0] = 5; num_df_per_set[1] = 3;
     df_ind[0] = 0; df_ind[1] = 5;
     dist_fact[0] = 1.0; dist_fact[1] = 2.0; dist_fact[2] = 3.0;
     dist_fact[3] = 4.0; dist_fact[4] = 5.0;
     dist_fact[5] = 1.1; dist_fact[6] = 2.1; dist_fact[7] = 3.1;
     error = ex_put_concat_node_sets (exoid, ids, num_nodes_per_set,
     num_df_per_set, node_ind,
     df_ind, node_list, dist_fact);
     error = ex_put_prop(exoid, EX_NODE_SET, 20, “FACE”, 4);
     error = ex_put_prop(exoid, EX_NODE_SET, 21, “FACE”, 5);
     D-4
     prop_array[0] = 1000;
     prop_array[1] = 2000;
     error = ex_put_prop_array(exoid, EX_NODE_SET, “VELOCITY”, prop_array);
     END COMMENTED OUT SECTION */
    /* write individual side sets */
    /* side set #1  - quad */
    error = ex_put_side_set_param (exoid, 30, 2, 4);
    elem_list[0] = 2; elem_list[1] = 2;
    side_list[0] = 4; side_list[1] = 2;
    dist_fact[0] = 30.0; dist_fact[1] = 30.1; dist_fact[2] = 30.2;
    dist_fact[3] = 30.3;
    error = ex_put_side_set (exoid, 30, elem_list, side_list);
    error = ex_put_side_set_dist_fact (exoid, 30, dist_fact);
    /* side set #2  - quad, spanning 2 elements  */
    error = ex_put_side_set_param (exoid, 31, 2, 4);
    elem_list[0] = 1; elem_list[1] = 2;
    side_list[0] = 2; side_list[1] = 3;
    dist_fact[0] = 31.0; dist_fact[1] = 31.1; dist_fact[2] = 31.2;
    dist_fact[3] = 31.3;
    error = ex_put_side_set (exoid, 31, elem_list, side_list);
    error = ex_put_side_set_dist_fact (exoid, 31, dist_fact);
    /* side set #3  - hex */
    error = ex_put_side_set_param (exoid, 32, 7, 0);
    elem_list[0] = 3; elem_list[1] = 3;
    elem_list[2] = 3; elem_list[3] = 3;
    elem_list[4] = 3; elem_list[5] = 3;
    elem_list[6] = 3;
    side_list[0] = 5; side_list[1] = 3;
    side_list[2] = 3; side_list[3] = 2;
    side_list[4] = 4; side_list[5] = 1;
    side_list[6] = 6;
    error = ex_put_side_set (exoid, 32, elem_list, side_list);
    /* side set #4  - tetras */
    error = ex_put_side_set_param (exoid, 33, 4, 0);
    elem_list[0] = 4; elem_list[1] = 4;
    elem_list[2] = 4; elem_list[3] = 4;
    side_list[0] = 1; side_list[1] = 2;
    side_list[2] = 3; side_list[3] = 4;
    error = ex_put_side_set (exoid, 33, elem_list, side_list);
    /* side set #5  - wedges */
    error = ex_put_side_set_param (exoid, 34, 5, 0);
    elem_list[0] = 5; elem_list[1] = 5;
    elem_list[2] = 5; elem_list[3] = 5;
    elem_list[4] = 5;
    side_list[0] = 1; side_list[1] = 2;
    side_list[2] = 3; side_list[3] = 4;
    side_list[4] = 5;
    error = ex_put_side_set (exoid, 34, elem_list, side_list);
    /* write concatenated side sets; side set node lists (which is how side sets
     * were described in EXODUS I) are converted to side set side lists and then
     * written out; this produces the same information as the above code which
     * writes individual side sets
     */
    /* THIS SECTION IS COMMENTED OUT
     ids[0] = 30;
     ids[1] = 31;
     ids[2] = 32;
     ids[3] = 33;
     ids[4] = 34;
     node_list[0] = 8; node_list[1] = 5;
     node_list[2] = 6; node_list[3] = 7;
     node_list[4] = 2; node_list[5] = 3;
     node_list[6] = 7; node_list[7] = 8;
     node_list[8] = 9; node_list[9] = 12;
     node_list[10] = 11; node_list[11] = 10;
     node_list[12] = 11; node_list[13] = 12;
     node_list[14] = 16; node_list[15] = 15;
     node_list[16] = 16; node_list[17] = 15;
     node_list[18] = 11; node_list[19] = 12;
     node_list[20] = 10; node_list[21] = 11;
     node_list[22] = 15; node_list[23] = 14;
     node_list[24] = 13; node_list[25] = 16;
     node_list[26] = 12; node_list[27] =  9;
     node_list[28] = 14; node_list[29] = 13;
     node_list[30] =  9; node_list[31] = 10;
     node_list[32] = 16; node_list[33] = 13;
     node_list[34] = 14; node_list[35] = 15;
     node_list[36] = 17; node_list[37] = 18;
     node_list[38] = 20;
     node_list[39] = 18; node_list[40] = 19;
     node_list[41] = 20;
     node_list[42] = 20; node_list[43] = 19;
     node_list[44] = 17;
     node_list[45] = 19; node_list[46] = 18;
     node_list[47] = 17;
     node_list[48] = 25; node_list[49] = 24;
     node_list[50] = 21; node_list[51] = 22;
     node_list[52] = 26; node_list[53] = 25;
     node_list[54] = 22; node_list[55] = 23;
     node_list[56] = 26; node_list[57] = 23;
     node_list[58] = 21; node_list[59] = 24;
     D-6
     node_list[60] = 23; node_list[61] = 22;
     node_list[62] = 21;
     node_list[63] = 24; node_list[64] = 25;
     node_list[65] = 26;
     node_ind[0] = 0;
     node_ind[1] = 4;
     node_ind[2] = 8;
     node_ind[3] = 36;
     node_ind[4] = 47;
     num_elem_per_set[0] = 2;
     num_elem_per_set[1] = 2;
     num_elem_per_set[2] = 7;
     num_elem_per_set[3] = 4;
     num_elem_per_set[4] = 5;
     num_nodes_per_set[0] = 4;
     num_nodes_per_set[1] = 4;
     num_nodes_per_set[2] = 28;
     num_nodes_per_set[3] = 12;
     num_nodes_per_set[4] = 18;
     elem_ind[0] = 0;
     elem_ind[1] = 2;
     elem_ind[2] = 4;
     elem_ind[3] = 11;
     elem_ind[4] = 15;
     elem_list[0] = 2; elem_list[1] = 2;
     elem_list[2] = 1; elem_list[3] = 2;
     elem_list[4] = 3; elem_list[5] = 3;
     elem_list[6] = 3; elem_list[7] = 3;
     elem_list[8] = 3; elem_list[9] = 3;
     elem_list[10] = 3; elem_list[11] = 4;
     elem_list[12] = 4; elem_list[13] = 4;
     elem_list[14] = 4; elem_list[15] = 5;
     elem_list[16] = 5; elem_list[17] = 5;
     elem_list[18] = 5; elem_list[19] = 5;
     error = ex_cvt_nodes_to_sides(exoid,
     num_elem_per_set,
     num_nodes_per_set,
     elem_ind,
     node_ind,
     elem_list,
     node_list,
     side_list);
     num_df_per_set[0] = 4;
     num_df_per_set[1] = 4;
     num_df_per_set[2] = 0;
     num_df_per_set[3] = 0;
     num_df_per_set[4] = 0;
     df_ind[0] = 0;
     df_ind[1] = 4;
     dist_fact[0] = 30.0; dist_fact[1] = 30.1;
     dist_fact[2] = 30.2; dist_fact[3] = 30.3;
     dist_fact[4] = 31.0; dist_fact[5] = 31.1;
     dist_fact[6] = 31.2; dist_fact[7] = 31.3;
     error = ex_put_concat_side_sets (exoid, ids, num_elem_per_set,
     num_df_per_set, elem_ind, df_ind,
     elem_list, side_list, dist_fact);
     END COMMENTED OUT SECTION */
    error = ex_put_prop(exoid, EX_SIDE_SET, 30, "COLOR", 100);
    error = ex_put_prop(exoid, EX_SIDE_SET, 31, "COLOR", 101);
    /* write QA records */
    num_qa_rec = 2;
    qa_record[0][0] = "TESTWT";
    qa_record[0][1] = "testwt";
    qa_record[0][2] = "07/07/93";
    qa_record[0][3] = "15:41:33";
    qa_record[1][0] = "FASTQ";
    qa_record[1][1] = "fastq";
    qa_record[1][2] = "07/07/93";
    qa_record[1][3] = "16:41:33";
    error = ex_put_qa (exoid, num_qa_rec, qa_record);
    /* write information records */
    num_info = 3;
    info[0] = "This is the first information record.";
    info[1] = "This is the second information record.";
    info[2] = "This is the third information record.";
    error = ex_put_info (exoid, num_info, info);
    /* write results variables parameters and names */
    num_glo_vars = 1;
    var_names[0] = "glo_vars";
    error = ex_put_var_param (exoid, "g", num_glo_vars);
    error = ex_put_var_names (exoid, "g", num_glo_vars, var_names);
    num_nod_vars = 2;
    var_names[0] = "nod_var0";
    var_names[1] = "nod_var1";
    error = ex_put_var_param (exoid, "n", num_nod_vars);
    error = ex_put_var_names (exoid, "n", num_nod_vars, var_names);
    num_ele_vars = 3;
    var_names[0] = "ele_var0";
    var_names[1] = "ele_var1";
    var_names[2] = "ele_var2";
    error = ex_put_var_param (exoid, "e", num_ele_vars);
    error = ex_put_var_names (exoid, "e", num_ele_vars, var_names);
    /* write element variable truth table */
    truth_tab = (int *) calloc ((num_elem_blk*num_ele_vars), sizeof(int));
    k = 0;
    for (i=0; i<num_elem_blk; i++)
    {
        for (j=0; j<num_ele_vars; j++)
        {
            truth_tab[k++] = 1;
        }
    }
    error = ex_put_elem_var_tab (exoid, num_elem_blk, num_ele_vars, truth_tab);
    free (truth_tab);
    /* for each time step, write the analysis results;
     * the code below fills the arrays glob_var_vals,
     D-8
     * nodal_var_vals, and elem_var_vals with values for debugging purposes;
     * obviously the analysis code will populate these arrays
     */
    whole_time_step = 1;
    num_time_steps = 10;
    glob_var_vals = (float *) calloc (num_glo_vars, CPU_word_size);
    nodal_var_vals = (float *) calloc (num_nodes, CPU_word_size);
    elem_var_vals = (float *) calloc (4, CPU_word_size);
    for (i=0; i<num_time_steps; i++)
    {
        time_value = (float)(i+1)/100.;
        /* write time value */
        error = ex_put_time (exoid, whole_time_step, &time_value);
        /* write global variables */
        for (j=0; j<num_glo_vars; j++)
        {
            glob_var_vals[j] = (float)(j+2) * time_value;
        }
        error = ex_put_glob_vars (exoid, whole_time_step, num_glo_vars,
                                  glob_var_vals);
        /* write nodal variables */
        for (k=1; k<=num_nod_vars; k++)
        {
            for (j=0; j<num_nodes; j++)
            {
                nodal_var_vals[j] = (float)k + ((float)(j+1) * time_value);
            }
            error = ex_put_nodal_var (exoid, whole_time_step, k, num_nodes,
                                      nodal_var_vals);
        }
        /* write element variables */
        for (k=1; k<=num_ele_vars; k++)
        {
            for (j=0; j<num_elem_blk; j++)
            {
                for (m=0; m<num_elem_in_block[j]; m++)
                {
                    elem_var_vals[m] = (float)(k+1) + (float)(j+2) +
                    ((float)(m+1)*time_value);
                }
                error = ex_put_elem_var (exoid, whole_time_step, k, ebids[j],num_elem_in_block[j], elem_var_vals);
            }
            whole_time_step++;
        }
                /* update the data file; this should be done at the end of every time step
                 * to ensure that no data is lost if the analysis dies
                 */
                error = ex_update (exoid);
            }
            free(glob_var_vals);
            free(nodal_var_vals);
            free(elem_var_vals);
            /* close the EXODUS files
             */
            error = ex_close (exoid);
        }

