/* -----------------------------------------------------------------------------
 
 Unified STO projection / AOM overlap code
 
 Alpha version: 0.1
 1-Oct-2019
 
 Orestis George Ziogos, UCL
 o.ziogos@ucl.ac.uk
 
 For more information, examine the README file in the parent directory.
 
-------------------------------------------------------------------------------- */
#define cmax_length 1000
#define sub_length 200
#define AngToBohr 1.8897259886
#define equal_diff 1e-6
// cube related
#define offset 5.0
#define print_thres 1.0e-20
#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<math.h>
#include<string.h>

struct simplex_obj_f_input {
    int MOLog_channel;                  //
    int simplex_entries;                //
    char **simplex_name;                //
    char **simplex_MO;                  //
    char **simplex_cube;                //
    int *pos_storage_L;                 //
    int *pos_storage_R;                 //
    int *CGTOs_storage_L;               //
    int *CGTOs_storage_R;               //
    int *PCGTOs_storage_L;              //
    int *PCGTOs_storage_R;              //
    int *MO_storage_L;                  //
    int *MO_storage_R;                  //
    int *atoms_aggr;                    //
    double *x_aggr;                     //
    double *y_aggr;                     //
    double *z_aggr;                     //
    char **species_aggr;                //
    int *CGTOs_agg;                     //
    int *PCGTOs_agg;                    //
    int *bfnPerAtom_agg;                //
    int *GTO_depth_agg;                 //
    double *palpha_agg;                 //
    double *pcoeff_agg;                 //
    int **pqn_agg;                      //
    double *MOcoeffs_aggr_alpha;        //
    double *MOcoeffs_aggr_beta;         //
    char current_folder[cmax_length];   //
    char STOproj_cube_grid[cmax_length];//
    int *simplex_id;                    //
    int HOMOs;                          //
    int LUMOs;                          //
    
    int *AOM_type_map;
    int *AOM_frag_atoms;
    int **AOM_indices;
    double *AOM_pi_values;
    double *HAB;
    
    int mode;
};


struct simplex_params {
    int steps;
    double acc;
    double rho;
    double xi;
    double gamma;
    double sigma;
    double tau;
    double tau_prime;
    double p_mu_min[16];
    double p_mu_max[16];
    
    double p_mu_sample_step[16];
    double p_mu_sample_max[16];
    
};
