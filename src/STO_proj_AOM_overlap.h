/* -----------------------------------------------------------------------------
 
 Unified STO projection / AOM overlap code
 
 Beta version: 1.0
 1-Mar-2019
 
 Orestis George Ziogos, UCL
 o.ziogos@ucl.ac.uk
 
 For more information, examine the README file in the parent directory.
 
-------------------------------------------------------------------------------- */

void read_config_file(char *current_folder,char *config_file,char *mode,char *STOproj_name,char *STOproj_molecule,char *STOproj_basis,char *STOproj_MO_file,char *STOproj_MO,char *STOproj_basis_file,char *STOproj_cube_grid,double *smu_per_species,double *pmu_per_species,double *AOMsmu_per_species,double *AOMpmu_per_species,int *verb,int *cube,char *STOproj_AOM_include, int *frag1atoms,int *MOLog_channel,int *phase1,int *phase2,
    int *simplex_entries,char ***simplex_name,char ***simplex_molecule,char ***simplex_MO_file,char ***simplex_MO,char ***simplex_cube,char ***simplex_type,struct simplex_params *simplex,
    char ***simplex_dimer,char ***simplex_AOM_include,int **simplex_frag1_atoms,double **simplex_HAB);

void read_xyz_convert_Ang_to_Bohr(char *STOproj_molecule,int *atoms,double **x,double **y,double **z,char ***species,int verb);

void initialize_STO(int atoms,char **species,double *smu_per_species,double *pmu_per_species,int *STOs,int **STO_id_array,int **STO_type_array,double **STO_mu_array,int verb);

void calculate_overlap_S_matrix(int STOs,double *x,double *y,double *z,int *STO_id_array,int *STO_type_array,double *STO_mu_array,double ***Smatrix,int verb);

void resolve_unique_species(int atoms,char **species,int *atom_types,char ***unique_species,int verb);

void read_CP2K_GTOs(char *current_folder,int atoms,int atom_types,char **species,char **unique_species,char *STOproj_basis_file,char *STOproj_basis,int *CGTOs,int **bfnPerAtom,int **GTO_depth,int *PCGTOs,double **palpha,double **pcoeff,int ***pqn,int verb);

void locate_MO(char *STOproj_MO_file,char *STOproj_MO,int CGTOs,double ***MOcoeffs,int verb);

double GTO_overlap(int atoms,double *x,double *y,double *z,int *bfnPerAtom,int *GTO_depth,double **MOcoeffs,double *pcoeff,double *palpha,int **pqn,int MOLog_channel);

void STO_GTO_projection(int atoms,double *x,double *y,double *z,int STOs,int *STO_type_array,int *STO_id_array,double *STO_mu_array,double **Smatrix,int *bfnPerAtom,int *GTO_depth,double **MOcoeffs,double *pcoeff,double *palpha,int **pqn,double **res,int verb,int MOLog_channel,double **V_array,double **s_array);

void create_cube_file(char *current_folder,char *STOproj_cube_grid,char *STOproj_name,char *STOproj_MO,int atoms,char **species,double *x,double *y,double *z,double **STO_matrix,double *smu_per_species,double *pmu_per_species,int verb);

void resolve_pvecs(int atoms,double *x,double *y,double *z,char **species,double **px,double **py,double **pz);

int resolve_atomic_Z(char *species);

double AOM_overlap_calculation(int istart,int istop,int jstart,int jstop,double *x,double *y,double *z,int *STO_id_array,int *STO_type_array,double *STO_mu_array,double **STO_matrix);

void resolve_pvecs_per_fragment(int atoms,double *x,double *y,double *z,char **species,double **px,double **py,double **pz,int frag1atoms);

void AOM_frag_cubes(char *current_folder,int atoms,int frag1atoms,double *x,double *y,double *z,char **species,double **STO_matrix,double *smu_per_species,double *pmu_per_species,int verb);

double simplex_obj_f(double *smu_per_species,double *pmu_per_species,struct simplex_obj_f_input current_input,double **compl_array);

void run_simplex(double *smu_per_species,double *pmu_per_species,struct simplex_obj_f_input current_input,struct simplex_params simplex);

double AOM_simplex_obj_f(double *smu_per_species,double *pmu_per_species,struct simplex_obj_f_input current_input,double **res);

void write_state(char *current_folder,char *STOproj_name,char *STOproj_MO,
                 int atoms,int STOs,double *px,double *py,double *pz,double **Smatrix,double **STO_matrix,
                 double orb_compl,double *V_array,double *s_array);
