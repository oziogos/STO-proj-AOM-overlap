/* -----------------------------------------------------------------------------
 
 Unified STO projection / AOM overlap code
 
 Alpha version: 0.1
 1-Oct-2019
 
 Orestis George Ziogos, UCL
 o.ziogos@ucl.ac.uk
 
 For more information, examine the README file in the parent directory.
 
-------------------------------------------------------------------------------- */
// -llapack -lblas -lgfortran
#include"general.h"
#include"STO_proj_AOM_overlap.h"
int main(int argc,char **argv)
{
    FILE *fp;
    char current_folder[cmax_length],file_path[cmax_length],**species,**unique_species;
    int atoms=0,i,j,atom_types;
    double S,*res,*px,*py,*pz,*V_array,*s_array;
    char mode[cmax_length],STOproj_name[cmax_length],STOproj_molecule[cmax_length],STOproj_basis[cmax_length];
    char STOproj_MO_file[cmax_length],STOproj_MO[cmax_length],STOproj_basis_file[cmax_length],STOproj_cube_grid[cmax_length];
    int verb,cube;
    double smu_per_species[16],pmu_per_species[16],AOMsmu_per_species[16],AOMpmu_per_species[16];
    int STOs,*STO_id_array,*STO_type_array;
    double orb_compl,*x,*y,*z,*STO_mu_array,**Smatrix,**STO_matrix;
    int CGTOs,PCGTOs,*GTO_depth,*bfnPerAtom,**pqn;
    double *palpha,*pcoeff;
    double **MOcoeffs;
    char STOproj_AOM_include[cmax_length],word[cmax_length],buffer[cmax_length];
    int frag1atoms,frag1STOs,MOLog_channel;
    double *AOM_pi;
    int phase1=1,phase2=1;
    
    int simplex_entries,k;
    char **simplex_name,**simplex_molecule,**simplex_MO_file,**simplex_cube,**simplex_MO,**simplex_type;
    
    struct simplex_params simplex;
    
    char **simplex_dimer,**simplex_AOM_include;
    int *simplex_frag1_atoms;
    double *simplex_HAB;
    
    // store working directory
    getcwd(current_folder,cmax_length);
    // read config file
    read_config_file(current_folder,argv[1],mode,
                     STOproj_name,STOproj_molecule,STOproj_basis,STOproj_MO_file,STOproj_MO,STOproj_basis_file,STOproj_cube_grid,smu_per_species,pmu_per_species,
                     AOMsmu_per_species,AOMpmu_per_species,&verb,&cube,STOproj_AOM_include,&frag1atoms,&MOLog_channel,&phase1,&phase2,
                     &simplex_entries,&simplex_name,&simplex_molecule,&simplex_MO_file,&simplex_MO,&simplex_cube,&simplex_type,&simplex,
                     &simplex_dimer,&simplex_AOM_include,&simplex_frag1_atoms,&simplex_HAB);
    // single molecule mode
    if(strcmp(mode,"molecule")==0)
    {
        // read single molecule xyz file (distance: Ang)
        read_xyz_convert_Ang_to_Bohr(STOproj_molecule,&atoms,&x,&y,&z,&species,verb);
        // resolve number of STO basis functions according to species information and populate STO-related arrays
        initialize_STO(atoms,species,smu_per_species,pmu_per_species,&STOs,&STO_id_array,&STO_type_array,&STO_mu_array,verb);
        // calculate overlap S matrix
        calculate_overlap_S_matrix(STOs,x,y,z,STO_id_array,STO_type_array,STO_mu_array,&Smatrix,verb);
        // deduce unique species
        resolve_unique_species(atoms,species,&atom_types,&unique_species,verb);
        // read and process CP2K internal basis ASCII file
        read_CP2K_GTOs(current_folder,atoms,atom_types,species,unique_species,STOproj_basis_file,STOproj_basis,&CGTOs,&bfnPerAtom,&GTO_depth,&PCGTOs,&palpha,&pcoeff,&pqn,verb);
        // read MOLog and store target MO coefficients
        locate_MO(STOproj_MO_file,STOproj_MO,CGTOs,&MOcoeffs,verb);
        // calculate GTO overlap - console output
        if(verb==-1){
            S=GTO_overlap(atoms,x,y,z,bfnPerAtom,GTO_depth,MOcoeffs,pcoeff,palpha,pqn,MOLog_channel);
            printf("Single molecule GTO overlap based on CP2K MO %s: %lf\n",STOproj_MO,S);
            printf("\n--------------------------------------------------------------------------\n\n");
        }
        // STO-GTO projection
        STO_GTO_projection(atoms,x,y,z,STOs,STO_type_array,STO_id_array,STO_mu_array,Smatrix,bfnPerAtom,GTO_depth,MOcoeffs,pcoeff,palpha,pqn,&res,verb,MOLog_channel,&V_array,&s_array);
        // calculate orbital completeness
        orb_compl=0.0;for(i=0;i<STOs;++i)for(j=0;j<STOs;++j)orb_compl=orb_compl+res[i]*res[j]*Smatrix[i][j];
        // console output
        if(verb==-1||verb==1){
            printf("Projection completeness: %lf\n",orb_compl);
            printf("\n--------------------------------------------------------------------------\n\n");
        }
        // normalize
        for(i=0;i<STOs;++i)res[i]=res[i]/sqrt(fabs(orb_compl));
        // STO orbital decomposition
        STO_matrix=(double**)malloc(atoms*sizeof(double*));for(i=0;i<atoms;++i)STO_matrix[i]=(double*)malloc(4*sizeof(double));
        // console output
        if(verb==-1){printf("STO orbital decomposition (input for getorbsgen); normalized:\n");}
        for(i=0;i<STOs;++i)
        {
            if(STO_type_array[i]==1){
                STO_matrix[STO_id_array[i]-1][0]=res[i];
                STO_matrix[STO_id_array[i]-1][1]=0.0;
                STO_matrix[STO_id_array[i]-1][2]=0.0;
                STO_matrix[STO_id_array[i]-1][3]=0.0;
            }
            else{
                STO_matrix[STO_id_array[i]-1][(STO_type_array[i]-2)%4]=res[i];
            }
        }
        // console output
        if(verb==-1){
            for(i=0;i<atoms;++i)printf("%s\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",species[i],x[i],y[i],z[i],STO_matrix[i][0],STO_matrix[i][1],STO_matrix[i][2],STO_matrix[i][3]);
            printf("\n--------------------------------------------------------------------------\n\n");
        }
        // STO cube file
        if(cube==1)create_cube_file(current_folder,STOproj_cube_grid,STOproj_name,STOproj_MO,atoms,species,x,y,z,STO_matrix,smu_per_species,pmu_per_species,verb);
        // pvecs
        resolve_pvecs(atoms,x,y,z,species,&px,&py,&pz);
        // console output
        if(verb==-1){
            printf("pvecs | AOM pi contribution:\n");
            for(i=0;i<atoms;++i)printf("[%d]\t%s\t%lf\t%lf\t%lf\t|\t%lf\n",i+1,species[i],px[i],py[i],pz[i],px[i]*STO_matrix[i][1]+py[i]*STO_matrix[i][2]+pz[i]*STO_matrix[i][3]);
            printf("\n--------------------------------------------------------------------------\n\n");
        }
        // write AOM include file
        sprintf(file_path,"%s/AOM_COEFF_single_%s_MO_%s.include",current_folder,STOproj_name,STOproj_MO);
        fp=fopen(file_path,"w+");
        for(i=0;i<atoms;++i)fprintf(fp,"%s %d %d %.1lf %lf\n",species[i],resolve_atomic_Z(species[i]),1,0.0,px[i]*STO_matrix[i][1]+py[i]*STO_matrix[i][2]+pz[i]*STO_matrix[i][3]);
        fclose(fp);
        // console output
        if(verb==-1||verb==1){
            printf("Generated single molecule AOM file: %s\n",file_path);
            printf("\n--------------------------------------------------------------------------\n\n");
        }
        // write state
        write_state(current_folder,STOproj_name,STOproj_MO,atoms,STOs,px,py,pz,Smatrix,STO_matrix,orb_compl,V_array,s_array);
        // free memory
        for(i=0;i<atom_types;++i)free(unique_species[i]);free(unique_species);
        for(i=0;i<STOs;++i)free(Smatrix[i]);free(Smatrix);
        free(bfnPerAtom);free(GTO_depth);free(palpha);free(pcoeff);for(i=0;i<PCGTOs;++i)free(pqn[i]);free(pqn);
        for(i=0;i<CGTOs;++i)free(MOcoeffs[i]);free(MOcoeffs);
        free(res);
        free(V_array);free(s_array);
        // free common memory
        for(i=0;i<atoms;++i)free(species[i]);free(species);free(x);free(y);free(z);
        free(STO_type_array);free(STO_id_array);free(STO_mu_array);
        free(px);free(py);free(pz);
        for(i=0;i<atoms;++i)free(STO_matrix[i]);free(STO_matrix);
    }
    else if(strcmp(mode,"dimer")==0)
    {
        read_xyz_convert_Ang_to_Bohr(STOproj_molecule,&atoms,&x,&y,&z,&species,verb);
        initialize_STO(atoms,species,AOMsmu_per_species,AOMpmu_per_species,&STOs,&STO_id_array,&STO_type_array,&STO_mu_array,verb);
        // read AOM include file
        sprintf(file_path,"%s",STOproj_AOM_include);fp=fopen(file_path,"r");if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
        AOM_pi=(double*)malloc(atoms*sizeof(double));for(i=0;i<atoms;++i){fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%s\t%s\t%s\t%lf",word,word,word,word,&AOM_pi[i]);}
        for(i=0;i<frag1atoms;++i)AOM_pi[i]=AOM_pi[i]*phase1;
        for(i=frag1atoms;i<atoms;++i)AOM_pi[i]=AOM_pi[i]*phase2;
        fclose(fp);
        // pvecs: calculates pvecs using topology derivation of the whole dimer
        //resolve_pvecs(atoms,x,y,z,species,&px,&py,&pz);
        // pvecs: treat fragments
        resolve_pvecs_per_fragment(atoms,x,y,z,species,&px,&py,&pz,frag1atoms);
        // console output
        if(verb==-1){
            printf("pvecs:\n");
            for(i=0;i<atoms;++i)printf("[%d]\t%s\t%lf\t%lf\t%lf\n",i+1,species[i],px[i],py[i],pz[i]);
            printf("\n--------------------------------------------------------------------------\n\n");
        }
        // STO orbital decomposition
        STO_matrix=(double**)malloc(atoms*sizeof(double*));for(i=0;i<atoms;++i)STO_matrix[i]=(double*)malloc(4*sizeof(double));
        // console output
        if(verb==-1){printf("STO orbital decomposition:\n");}
        for(i=0;i<STOs;++i)
        {
            if(STO_type_array[i]==1){
                STO_matrix[STO_id_array[i]-1][0]=0.0;
                STO_matrix[STO_id_array[i]-1][1]=0.0;
                STO_matrix[STO_id_array[i]-1][2]=0.0;
                STO_matrix[STO_id_array[i]-1][3]=0.0;
            }
            else{
                j=(STO_type_array[i]-2)%4;
                if(j==0)STO_matrix[STO_id_array[i]-1][j]=0.0;
                if(j==1)STO_matrix[STO_id_array[i]-1][j]=px[STO_id_array[i]-1]*AOM_pi[STO_id_array[i]-1];
                if(j==2)STO_matrix[STO_id_array[i]-1][j]=py[STO_id_array[i]-1]*AOM_pi[STO_id_array[i]-1];
                if(j==3)STO_matrix[STO_id_array[i]-1][j]=pz[STO_id_array[i]-1]*AOM_pi[STO_id_array[i]-1];
            }
        }
        // console output
        if(verb==-1){
            for(i=0;i<atoms;++i)printf("%s\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",species[i],x[i],y[i],z[i],STO_matrix[i][0],STO_matrix[i][1],STO_matrix[i][2],STO_matrix[i][3]);
            printf("\n--------------------------------------------------------------------------\n\n");
        }
        // resolve number of STO orbitals belonging to fragment 1
        frag1STOs=0;for(i=0;i<STOs;++i)if(STO_id_array[i]<=frag1atoms)++frag1STOs;
        // calculate overlap for frag1
        S=AOM_overlap_calculation(0,frag1STOs,0,frag1STOs,x,y,z,STO_id_array,STO_type_array,STO_mu_array,STO_matrix);
        if(verb==-1||verb==1){printf("Fragment 1 AOM overlap prior to normalization: %lf\n",S);}
        // normalize
        S=fabs(sqrt(S));for(i=0;i<frag1atoms;++i)for(j=1;j<=3;++j)STO_matrix[i][j]=STO_matrix[i][j]/S;
        // check normalization
        if(verb==-1||verb==1){
            S=AOM_overlap_calculation(0,frag1STOs,0,frag1STOs,x,y,z,STO_id_array,STO_type_array,STO_mu_array,STO_matrix);
            printf("Fragment 1 AOM overlap after normalization: %lf\n",S);printf("\n--------------------------------------------------------------------------\n\n");}
        // calculate overlap for frag2
        S=AOM_overlap_calculation(frag1STOs,STOs,frag1STOs,STOs,x,y,z,STO_id_array,STO_type_array,STO_mu_array,STO_matrix);
        if(verb==-1||verb==1){printf("Fragment 2 AOM overlap prior to normalization: %lf\n",S);}
        // normalize
        S=fabs(sqrt(S));for(i=frag1atoms;i<atoms;++i)for(j=1;j<=3;++j)STO_matrix[i][j]=STO_matrix[i][j]/S;
        // check normalization
        if(verb==-1||verb==1){
            S=AOM_overlap_calculation(frag1STOs,STOs,frag1STOs,STOs,x,y,z,STO_id_array,STO_type_array,STO_mu_array,STO_matrix);
            printf("Fragment 2 AOM overlap after normalization: %lf\n",S);printf("\n--------------------------------------------------------------------------\n\n");}
        // calculate AOM overlap between frag1 and frag2
        S=AOM_overlap_calculation(0,frag1STOs,frag1STOs,STOs,x,y,z,STO_id_array,STO_type_array,STO_mu_array,STO_matrix);
        // console output
        if(verb==-1||verb==1){printf("Dimer AOM overlap: %lf\n",S);printf("\n--------------------------------------------------------------------------\n\n");}
        if(verb==0)printf("%lf\n",S);
        // frag AOM cuve files
        if(verb==-1)AOM_frag_cubes(current_folder,atoms,frag1atoms,x,y,z,species,STO_matrix,AOMsmu_per_species,AOMpmu_per_species,verb);
        // free memory
        free(AOM_pi);
        // free common memory
        for(i=0;i<atoms;++i)free(species[i]);free(species);free(x);free(y);free(z);
        free(STO_type_array);free(STO_id_array);free(STO_mu_array);
        free(px);free(py);free(pz);
        for(i=0;i<atoms;++i)free(STO_matrix[i]);free(STO_matrix);
    }
    else if(strcmp(mode,"molecule_simplex")==0)
    {
        printf("Caching simplex obj function memory:\n");
        
        struct simplex_obj_f_input current_input;
                        
        current_input.MOLog_channel=MOLog_channel;
        
        current_input.simplex_entries=simplex_entries;
        current_input.simplex_name=(char**)malloc(simplex_entries*sizeof(char*));
        for(i=0;i<simplex_entries;++i)current_input.simplex_name[i]=(char*)malloc(cmax_length*sizeof(char));
        for(i=0;i<simplex_entries;++i)strcpy(current_input.simplex_name[i],simplex_name[i]);
        
        current_input.simplex_MO=(char**)malloc(simplex_entries*sizeof(char*));
        for(i=0;i<simplex_entries;++i)current_input.simplex_MO[i]=(char*)malloc(cmax_length*sizeof(char));
        for(i=0;i<simplex_entries;++i)strcpy(current_input.simplex_MO[i],simplex_MO[i]);
        
        current_input.simplex_cube=(char**)malloc(simplex_entries*sizeof(char*));
        for(i=0;i<simplex_entries;++i)current_input.simplex_cube[i]=(char*)malloc(cmax_length*sizeof(char));
        for(i=0;i<simplex_entries;++i)strcpy(current_input.simplex_cube[i],simplex_cube[i]);
        
        strcpy(current_input.current_folder,current_folder);
        strcpy(current_input.STOproj_cube_grid,STOproj_cube_grid);
        
        current_input.simplex_id=(int*)malloc(simplex_entries*sizeof(int));
        current_input.HOMOs=0;
        for(i=0;i<simplex_entries;++i)
            if(strcmp(simplex_type[i],"HOMO")==0)
            {
                current_input.simplex_id[i]=2;
                current_input.HOMOs=current_input.HOMOs+1;
            }
            else
            {
                current_input.simplex_id[i]=0;
            }
        current_input.LUMOs=simplex_entries-current_input.HOMOs;
        
        current_input.MO_storage_L=(int*)malloc(simplex_entries*sizeof(int));
        current_input.MO_storage_R=(int*)malloc(simplex_entries*sizeof(int));
        current_input.pos_storage_L=(int*)malloc(simplex_entries*sizeof(int));
        current_input.pos_storage_R=(int*)malloc(simplex_entries*sizeof(int));
        current_input.CGTOs_storage_L=(int*)malloc(simplex_entries*sizeof(int));
        current_input.CGTOs_storage_R=(int*)malloc(simplex_entries*sizeof(int));
        current_input.PCGTOs_storage_L=(int*)malloc(simplex_entries*sizeof(int));
        current_input.PCGTOs_storage_R=(int*)malloc(simplex_entries*sizeof(int));
        
        current_input.atoms_aggr=(int*)malloc(simplex_entries*sizeof(int));

        current_input.CGTOs_agg=(int*)malloc(simplex_entries*sizeof(int));
        current_input.PCGTOs_agg=(int*)malloc(simplex_entries*sizeof(int));

        
        for(i=0;i<simplex_entries;++i)
        {
            read_xyz_convert_Ang_to_Bohr(simplex_molecule[i],&atoms,&x,&y,&z,&species,0);
            current_input.atoms_aggr[i]=atoms;
            if(i==0)
            {
                current_input.pos_storage_L[i]=0;current_input.pos_storage_R[i]=atoms-1;
                current_input.x_aggr=(double*)malloc(atoms*sizeof(double));
                current_input.y_aggr=(double*)malloc(atoms*sizeof(double));
                current_input.z_aggr=(double*)malloc(atoms*sizeof(double));
                current_input.species_aggr=(char**)malloc(atoms*sizeof(char*));
                for(j=0;j<atoms;++j)current_input.species_aggr[j]=(char*)malloc(cmax_length*sizeof(char));
                for(j=0;j<atoms;++j)
                {
                    current_input.x_aggr[j]=x[j];
                    current_input.y_aggr[j]=y[j];
                    current_input.z_aggr[j]=z[j];
                    sprintf(current_input.species_aggr[j],"%s",species[j]);
                }
            }
            else
            {
                current_input.pos_storage_L[i]=current_input.pos_storage_R[i-1]+1;
                current_input.pos_storage_R[i]=current_input.pos_storage_L[i]+atoms-1;
                current_input.x_aggr=(double*)realloc(current_input.x_aggr,(current_input.pos_storage_R[i]+1)*sizeof(double));
                current_input.y_aggr=(double*)realloc(current_input.y_aggr,(current_input.pos_storage_R[i]+1)*sizeof(double));
                current_input.z_aggr=(double*)realloc(current_input.z_aggr,(current_input.pos_storage_R[i]+1)*sizeof(double));
                current_input.species_aggr=(char**)realloc(current_input.species_aggr,(current_input.pos_storage_R[i]+1)*sizeof(char*));
                k=-1;
                for(j=current_input.pos_storage_L[i];j<=current_input.pos_storage_R[i];++j)
                {
                    current_input.species_aggr[j]=(char*)malloc(cmax_length*sizeof(char));
                    k=k+1;
                    current_input.x_aggr[j]=x[k];
                    current_input.y_aggr[j]=y[k];
                    current_input.z_aggr[j]=z[k];
                    sprintf(current_input.species_aggr[j],"%s",species[k]);
                }
            }
            
            initialize_STO(atoms,species,smu_per_species,pmu_per_species,&STOs,&STO_id_array,&STO_type_array,&STO_mu_array,0);
            calculate_overlap_S_matrix(STOs,x,y,z,STO_id_array,STO_type_array,STO_mu_array,&Smatrix,0);
            resolve_unique_species(atoms,species,&atom_types,&unique_species,0);
            
            read_CP2K_GTOs(current_folder,atoms,atom_types,species,unique_species,STOproj_basis_file,STOproj_basis,&CGTOs,&bfnPerAtom,&GTO_depth,&PCGTOs,&palpha,&pcoeff,&pqn,0);
            current_input.CGTOs_agg[i]=CGTOs;
            current_input.PCGTOs_agg[i]=PCGTOs;
            if(i==0)
            {
                current_input.CGTOs_storage_L[i]=0;current_input.CGTOs_storage_R[i]=CGTOs-1;
                current_input.PCGTOs_storage_L[i]=0;current_input.PCGTOs_storage_R[i]=PCGTOs-1;
                current_input.bfnPerAtom_agg=(int*)malloc(atoms*sizeof(int));
                current_input.GTO_depth_agg=(int*)malloc(CGTOs*sizeof(int));
                current_input.palpha_agg=(double*)malloc(PCGTOs*sizeof(double));
                current_input.pcoeff_agg=(double*)malloc(PCGTOs*sizeof(double));
                current_input.pqn_agg=(int**)malloc(PCGTOs*sizeof(int*));
                for(j=0;j<PCGTOs;++j)current_input.pqn_agg[j]=(int*)malloc(3*sizeof(int));
                for(j=0;j<atoms;++j)
                {
                    current_input.bfnPerAtom_agg[j]=bfnPerAtom[j];
                }
                for(j=0;j<CGTOs;++j)
                {
                    current_input.GTO_depth_agg[j]=GTO_depth[j];
                }
                for(j=0;j<PCGTOs;++j)
                {
                    current_input.palpha_agg[j]=palpha[j];
                    current_input.pcoeff_agg[j]=pcoeff[j];
                    current_input.pqn_agg[j][0]=pqn[j][0];
                    current_input.pqn_agg[j][1]=pqn[j][1];
                    current_input.pqn_agg[j][2]=pqn[j][2];
                }
            }
            else
            {
                current_input.CGTOs_storage_L[i]=current_input.CGTOs_storage_R[i-1]+1;
                current_input.CGTOs_storage_R[i]=current_input.CGTOs_storage_L[i]+CGTOs-1;
                current_input.PCGTOs_storage_L[i]=current_input.PCGTOs_storage_R[i-1]+1;
                current_input.PCGTOs_storage_R[i]=current_input.PCGTOs_storage_L[i]+PCGTOs-1;
                current_input.bfnPerAtom_agg=(int*)realloc(current_input.bfnPerAtom_agg,(current_input.pos_storage_R[i]+1)*sizeof(int));
                current_input.GTO_depth_agg=(int*)realloc(current_input.GTO_depth_agg,(current_input.CGTOs_storage_R[i]+1)*sizeof(int));
                current_input.palpha_agg=(double*)realloc(current_input.palpha_agg,(current_input.PCGTOs_storage_R[i]+1)*sizeof(double));
                current_input.pcoeff_agg=(double*)realloc(current_input.pcoeff_agg,(current_input.PCGTOs_storage_R[i]+1)*sizeof(double));
                current_input.pqn_agg=(int**)realloc(current_input.pqn_agg,(current_input.PCGTOs_storage_R[i]+1)*sizeof(int*));
                k=-1;
                for(j=current_input.pos_storage_L[i];j<=current_input.pos_storage_R[i];++j)
                {
                    k=k+1;
                    current_input.bfnPerAtom_agg[j]=bfnPerAtom[k];
                }
                k=-1;
                for(j=current_input.CGTOs_storage_L[i];j<=current_input.CGTOs_storage_R[i];++j)
                {
                    k=k+1;
                    current_input.GTO_depth_agg[j]=GTO_depth[k];
                }
                k=-1;
                for(j=current_input.PCGTOs_storage_L[i];j<=current_input.PCGTOs_storage_R[i];++j)
                {
                    current_input.pqn_agg[j]=(int*)malloc(3*sizeof(int));
                    k=k+1;
                    current_input.palpha_agg[j]=palpha[k];
                    current_input.pcoeff_agg[j]=pcoeff[k];
                    current_input.pqn_agg[j][0]=pqn[k][0];
                    current_input.pqn_agg[j][1]=pqn[k][1];
                    current_input.pqn_agg[j][2]=pqn[k][2];
                }
            }
            
            locate_MO(simplex_MO_file[i],simplex_MO[i],CGTOs,&MOcoeffs,0);
            if(i==0)
            {
                current_input.MO_storage_L[i]=0;current_input.MO_storage_R[i]=CGTOs-1;
                current_input.MOcoeffs_aggr_alpha=(double*)malloc(CGTOs*sizeof(double));
                current_input.MOcoeffs_aggr_beta=(double*)malloc(CGTOs*sizeof(double));
                for(j=0;j<CGTOs;++j){current_input.MOcoeffs_aggr_alpha[j]=MOcoeffs[j][0];current_input.MOcoeffs_aggr_beta[j]=MOcoeffs[j][1];}
            }
            else
            {
                current_input.MO_storage_L[i]=current_input.MO_storage_R[i-1]+1;
                current_input.MO_storage_R[i]=current_input.MO_storage_L[i]+CGTOs-1;
                current_input.MOcoeffs_aggr_alpha=(double*)realloc(current_input.MOcoeffs_aggr_alpha,(current_input.MO_storage_R[i]+1)*sizeof(double));
                current_input.MOcoeffs_aggr_beta=(double*)realloc(current_input.MOcoeffs_aggr_beta,(current_input.MO_storage_R[i]+1)*sizeof(double));
                k=-1;
                for(j=current_input.MO_storage_L[i];j<=current_input.MO_storage_R[i];++j)
                {
                    k=k+1;
                    current_input.MOcoeffs_aggr_alpha[j]=MOcoeffs[k][0];
                    current_input.MOcoeffs_aggr_beta[j]=MOcoeffs[k][1];
                }
            }
            printf("%d/%d\t%s\n",i+1,simplex_entries,current_input.simplex_name[i]);
            
            for(j=0;j<atoms;++j)free(species[j]);free(species);free(x);free(y);free(z);
            free(STO_type_array);free(STO_id_array);free(STO_mu_array);
            for(j=0;j<STOs;++j)free(Smatrix[j]);free(Smatrix);
            for(j=0;j<atom_types;++j)free(unique_species[j]);free(unique_species);
            free(bfnPerAtom);free(GTO_depth);free(palpha);free(pcoeff);for(j=0;j<PCGTOs;++j)free(pqn[j]);free(pqn);
            for(j=0;j<CGTOs;++j)free(MOcoeffs[j]);free(MOcoeffs);
        }
        current_input.mode=0;
        printf("Done!\n");
        
        run_simplex(smu_per_species,pmu_per_species,current_input,simplex);
        
        // free simplex variables
        for(i=0;i<simplex_entries;++i)
        {
            free(simplex_name[i]);
            free(simplex_molecule[i]);
            free(simplex_MO_file[i]);
            free(simplex_cube[i]);
            free(simplex_MO[i]);
            free(simplex_type[i]);
        }
        free(simplex_name);
        free(simplex_molecule);
        free(simplex_MO_file);
        free(simplex_cube);
        free(simplex_MO);
        free(simplex_type);
        
        for(i=0;i<simplex_entries;++i)
        {
            free(current_input.simplex_name[i]);
            free(current_input.simplex_MO[i]);
            free(current_input.simplex_cube[i]);
        }
        free(current_input.simplex_name);
        free(current_input.simplex_MO);
        free(current_input.simplex_cube);
        free(current_input.atoms_aggr);
        free(current_input.CGTOs_agg);
        free(current_input.PCGTOs_agg);
        free(current_input.simplex_id);
        
        free(current_input.MOcoeffs_aggr_alpha);free(current_input.MOcoeffs_aggr_beta);
        for(i=0;i<=current_input.pos_storage_R[simplex_entries-1];++i)free(current_input.species_aggr[i]);free(current_input.species_aggr);
        free(current_input.x_aggr);free(current_input.y_aggr);free(current_input.z_aggr);
        free(current_input.bfnPerAtom_agg);free(current_input.GTO_depth_agg);free(current_input.palpha_agg);free(current_input.pcoeff_agg);
        for(i=0;i<=current_input.PCGTOs_storage_R[simplex_entries-1];++i)free(current_input.pqn_agg[i]);free(current_input.pqn_agg);
        
        free(current_input.MO_storage_L);free(current_input.MO_storage_R);
        free(current_input.pos_storage_L);free(current_input.pos_storage_R);
        free(current_input.CGTOs_storage_L);free(current_input.CGTOs_storage_R);
        free(current_input.PCGTOs_storage_L);free(current_input.PCGTOs_storage_R);
        
    }
    else if(strcmp(mode,"dimer_simplex")==0)
    {
        
        int molecule_types,type_map[simplex_entries];
        char **unique_molecules;
        
        struct simplex_obj_f_input current_input;
        
        // use the resolve_unique_species() function for unique dimer identification
        resolve_unique_species(simplex_entries,simplex_name,&molecule_types,&unique_molecules,0);
        //for(i=0;i<molecule_types;++i)printf("%s\n",unique_molecules[i]);
        
        int atoms_per_dimer_type[molecule_types];
        int frag1_atoms_per_dimer_type[molecule_types];
        char AOM_src[molecule_types][cmax_length];
        
        for(i=0;i<simplex_entries;++i)
        {
            for(j=0;j<molecule_types;++j)
            {
                if(strcmp(simplex_name[i],unique_molecules[j])==0)
                {
                    type_map[i]=j;
                    sprintf(AOM_src[j],"%s",simplex_AOM_include[i]);
                    frag1_atoms_per_dimer_type[j]=simplex_frag1_atoms[i];
                    break;
                }
            }
        }
        
        current_input.simplex_entries=simplex_entries;
        
        current_input.simplex_name=(char**)malloc(simplex_entries*sizeof(char*));
        for(i=0;i<simplex_entries;++i)current_input.simplex_name[i]=(char*)malloc(cmax_length*sizeof(char));
        for(i=0;i<simplex_entries;++i)strcpy(current_input.simplex_name[i],simplex_name[i]);
        
        current_input.pos_storage_L=(int*)malloc(simplex_entries*sizeof(int));
        current_input.pos_storage_R=(int*)malloc(simplex_entries*sizeof(int));
        
        current_input.atoms_aggr=(int*)malloc(simplex_entries*sizeof(int));
        
        current_input.HAB=(double*)malloc(simplex_entries*sizeof(double));
        
        printf("Caching simplex obj function memory:\n");
        
        for(i=0;i<simplex_entries;++i)
        {
            
            current_input.HAB[i]=simplex_HAB[i];
            
            read_xyz_convert_Ang_to_Bohr(simplex_dimer[i],&atoms,&x,&y,&z,&species,0);
            current_input.atoms_aggr[i]=atoms;
            
            atoms_per_dimer_type[type_map[i]]=atoms;
            
            if(i==0)
            {
                current_input.pos_storage_L[i]=0;current_input.pos_storage_R[i]=atoms-1;
                current_input.x_aggr=(double*)malloc(atoms*sizeof(double));
                current_input.y_aggr=(double*)malloc(atoms*sizeof(double));
                current_input.z_aggr=(double*)malloc(atoms*sizeof(double));
                current_input.species_aggr=(char**)malloc(atoms*sizeof(char*));
                for(j=0;j<atoms;++j)current_input.species_aggr[j]=(char*)malloc(cmax_length*sizeof(char));
                for(j=0;j<atoms;++j)
                {
                    current_input.x_aggr[j]=x[j];
                    current_input.y_aggr[j]=y[j];
                    current_input.z_aggr[j]=z[j];
                    sprintf(current_input.species_aggr[j],"%s",species[j]);
                }
            }
            else
            {
                current_input.pos_storage_L[i]=current_input.pos_storage_R[i-1]+1;
                current_input.pos_storage_R[i]=current_input.pos_storage_L[i]+atoms-1;
                current_input.x_aggr=(double*)realloc(current_input.x_aggr,(current_input.pos_storage_R[i]+1)*sizeof(double));
                current_input.y_aggr=(double*)realloc(current_input.y_aggr,(current_input.pos_storage_R[i]+1)*sizeof(double));
                current_input.z_aggr=(double*)realloc(current_input.z_aggr,(current_input.pos_storage_R[i]+1)*sizeof(double));
                current_input.species_aggr=(char**)realloc(current_input.species_aggr,(current_input.pos_storage_R[i]+1)*sizeof(char*));
                k=-1;
                for(j=current_input.pos_storage_L[i];j<=current_input.pos_storage_R[i];++j)
                {
                    current_input.species_aggr[j]=(char*)malloc(cmax_length*sizeof(char));
                    k=k+1;
                    current_input.x_aggr[j]=x[k];
                    current_input.y_aggr[j]=y[k];
                    current_input.z_aggr[j]=z[k];
                    sprintf(current_input.species_aggr[j],"%s",species[k]);
                }
            }
            for(j=0;j<atoms;++j)free(species[j]);free(species);free(x);free(y);free(z);
            
            printf("%d/%d\t%s\n",i+1,simplex_entries,simplex_dimer[i]);
            
        }
        
        int AOM_indices[molecule_types][2];
        
        AOM_indices[0][0]=0;AOM_indices[0][1]=atoms_per_dimer_type[0]-1;
        for(i=1;i<molecule_types;++i)
        {
            AOM_indices[i][0]=AOM_indices[i-1][1]+1;
            AOM_indices[i][1]=AOM_indices[i][0]+atoms_per_dimer_type[i]-1;
        }
        
        int AOM_length;
        
        AOM_length=AOM_indices[molecule_types-1][1]+1;
        
        double AOM_pi_values[AOM_length];
        
        j=-1;
        for(i=0;i<molecule_types;++i)
        {
            fp=fopen(AOM_src[i],"r");
            while(fgets(buffer,cmax_length,fp)!=NULL)
            {
                j=j+1;
                sscanf(buffer,"%s\t%s\t%s\t%s\t%lf",word,word,word,word,&AOM_pi_values[j]);
            }
            
            fclose(fp);
        }
        
        current_input.AOM_type_map=(int*)malloc(simplex_entries*sizeof(int));
        current_input.AOM_frag_atoms=(int*)malloc(molecule_types*sizeof(int));
        current_input.AOM_pi_values=(double*)malloc(AOM_length*sizeof(double));
        current_input.AOM_indices=(int**)malloc(molecule_types*sizeof(int*));
        for(i=0;i<molecule_types;++i)current_input.AOM_indices[i]=(int*)malloc(2*sizeof(int));
        
        for(i=0;i<simplex_entries;++i)current_input.AOM_type_map[i]=type_map[i];
        for(i=0;i<molecule_types;++i)
        {
            current_input.AOM_frag_atoms[i]=frag1_atoms_per_dimer_type[i];
            current_input.AOM_indices[i][0]=AOM_indices[i][0];
            current_input.AOM_indices[i][1]=AOM_indices[i][1];
        }
        for(i=0;i<AOM_length;++i)current_input.AOM_pi_values[i]=AOM_pi_values[i];
        current_input.mode=1;
        
        printf("Done!\n");
        
        /*
        for(i=0;i<simplex_entries;++i)
        {
            printf("xyz bounds: %d\t%d\n",current_input.pos_storage_L[i],current_input.pos_storage_R[i]);
            printf("atoms: %d\n",current_input.atoms_aggr[i]);
            printf("frag1_atoms: %d\n",frag1_atoms_per_dimer_type[type_map[i]]);
            for(j=current_input.pos_storage_L[i];j<=current_input.pos_storage_R[i];++j)
                printf("%s\t%lf\t%lf\t%lf\n",current_input.species_aggr[j],current_input.x_aggr[j],current_input.y_aggr[j],current_input.z_aggr[j]);
            for(j=AOM_indices[type_map[i]][0];j<=AOM_indices[type_map[i]][1];++j)
                printf("%lf\n",AOM_pi_values[j]);
            printf("-----------------------------------------------------------\n");
        }
        */
        
        run_simplex(AOMsmu_per_species,AOMpmu_per_species,current_input,simplex);
        
        // free memory
        
        free(simplex_HAB);
        
        free(current_input.AOM_type_map);free(current_input.AOM_frag_atoms);free(current_input.AOM_pi_values);
        for(i=0;i<molecule_types;++i)free(current_input.AOM_indices[i]);free(current_input.AOM_indices);
        
        for(i=0;i<simplex_entries;++i)free(current_input.simplex_name[i]);
        free(current_input.simplex_name);
        
        free(current_input.atoms_aggr);
        
        free(current_input.HAB);
        
        for(i=0;i<=current_input.pos_storage_R[simplex_entries-1];++i)free(current_input.species_aggr[i]);free(current_input.species_aggr);
        free(current_input.x_aggr);free(current_input.y_aggr);free(current_input.z_aggr);
        
        free(current_input.pos_storage_L);free(current_input.pos_storage_R);
        
        for(i=0;i<molecule_types;++i)free(unique_molecules[i]);free(unique_molecules);
        
        // free simplex variables
        for(i=0;i<simplex_entries;++i)
        {
            free(simplex_name[i]);
            free(simplex_dimer[i]);
            free(simplex_AOM_include[i]);
        }
        free(simplex_name);
        free(simplex_dimer);
        free(simplex_AOM_include);
        free(simplex_frag1_atoms);

    }
    return 0;
}
