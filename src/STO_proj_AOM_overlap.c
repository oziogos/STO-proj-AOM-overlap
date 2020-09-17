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
    double S,*res,*px,*py,*pz;
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
    // store working directory
    getcwd(current_folder,cmax_length);
    // read config file
    read_config_file(current_folder,argv[1],mode,
                     STOproj_name,STOproj_molecule,STOproj_basis,STOproj_MO_file,STOproj_MO,STOproj_basis_file,STOproj_cube_grid,smu_per_species,pmu_per_species,
                     AOMsmu_per_species,AOMpmu_per_species,&verb,&cube,STOproj_AOM_include,&frag1atoms,&MOLog_channel,&phase1,&phase2);
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
        STO_GTO_projection(atoms,x,y,z,STOs,STO_type_array,STO_id_array,STO_mu_array,Smatrix,bfnPerAtom,GTO_depth,MOcoeffs,pcoeff,palpha,pqn,&res,verb,MOLog_channel);
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
        // free memory
        for(i=0;i<atom_types;++i)free(unique_species[i]);free(unique_species);
        for(i=0;i<STOs;++i)free(Smatrix[i]);free(Smatrix);
        free(bfnPerAtom);free(GTO_depth);free(palpha);free(pcoeff);for(i=0;i<PCGTOs;++i)free(pqn[i]);free(pqn);
        for(i=0;i<CGTOs;++i)free(MOcoeffs[i]);free(MOcoeffs);
        free(res);
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
    }
    // free common memory
    for(i=0;i<atoms;++i)free(species[i]);free(species);free(x);free(y);free(z);
    free(STO_type_array);free(STO_id_array);free(STO_mu_array);
    free(px);free(py);free(pz);
    for(i=0;i<atoms;++i)free(STO_matrix[i]);free(STO_matrix);
    return 0;
}
