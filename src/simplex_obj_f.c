
#include"general.h"

void initialize_STO(int atoms,char **species,double *smu_per_species,double *pmu_per_species,int *STOs,int **STO_id_array,int **STO_type_array,double **STO_mu_array,int verb);
void calculate_overlap_S_matrix(int STOs,double *x,double *y,double *z,int *STO_id_array,int *STO_type_array,double *STO_mu_array,double ***Smatrix,int verb);
void resolve_unique_species(int atoms,char **species,int *atom_types,char ***unique_species,int verb);
void STO_GTO_projection(int atoms,double *x,double *y,double *z,int STOs,int *STO_type_array,int *STO_id_array,double *STO_mu_array,double **Smatrix,int *bfnPerAtom,int *GTO_depth,double **MOcoeffs,double *pcoeff,double *palpha,int **pqn,double **res,int verb,int MOLog_channel,double **V_array,double **s_array);
void create_cube_file(char *current_folder,char *STOproj_cube_grid,char *STOproj_name,char *STOproj_MO,int atoms,char **species,double *x,double *y,double *z,double **STO_matrix,double *smu_per_species,double *pmu_per_species,int verb);

void resolve_pvecs(int atoms,double *x,double *y,double *z,char **species,double **px,double **py,double **pz);

int resolve_atomic_Z(char *species);

void write_state(char *current_folder,char *STOproj_name,char *STOproj_MO,
                 int atoms,int STOs,double *px,double *py,double *pz,double **Smatrix,double **STO_matrix,
                 double orb_compl,double *V_array,double *s_array);

double simplex_obj_f(double *smu_per_species,double *pmu_per_species,struct simplex_obj_f_input current_input,double **compl_array);
double simplex_obj_f(double *smu_per_species,double *pmu_per_species,struct simplex_obj_f_input current_input,double **compl_array)
{
    FILE *fp;
    char file_path[cmax_length];
    int i,j,k;
    int atoms;
    double *x,*y,*z;
    char **species;
    int STOs;
    int *STO_id_array,*STO_type_array;
    double *STO_mu_array;
    double **Smatrix;
    int atom_types;
    char **unique_species;
    int CGTOs,PCGTOs,*bfnPerAtom,*GTO_depth,**pqn;
    double *palpha,*pcoeff,**MOcoeffs,*res,orb_compl,**STO_matrix,*V_array,*s_array;
    
    //double compl_array[current_input.simplex_entries];
    double compl_array_HOMO[current_input.HOMOs];
    double compl_array_LUMO[current_input.LUMOs];
    int Hcounter=-1,Lcounter=-1;
    
    double HOMO_mean=0.0,LUMO_mean=0.0;
    double HOMO_var,LUMO_var;
    
    double metric;
    
    double *px,*py,*pz;
    
    for(i=0;i<current_input.simplex_entries;++i)
    {
        //read_xyz_convert_Ang_to_Bohr(simplex_molecule[i],&atoms,&x,&y,&z,&species,0);
        atoms=current_input.atoms_aggr[i];
        x=(double*)malloc(atoms*sizeof(double));
        y=(double*)malloc(atoms*sizeof(double));
        z=(double*)malloc(atoms*sizeof(double));
        species=(char**)malloc(atoms*sizeof(char*));
        for(j=0;j<atoms;++j)species[j]=(char*)malloc(cmax_length*sizeof(char));
        k=-1;
        for(j=current_input.pos_storage_L[i];j<=current_input.pos_storage_R[i];++j)
        {
            k=k+1;
            x[k]=current_input.x_aggr[j];
            y[k]=current_input.y_aggr[j];
            z[k]=current_input.z_aggr[j];
            sprintf(species[k],"%s",current_input.species_aggr[j]);
        }
        
        initialize_STO(atoms,species,smu_per_species,pmu_per_species,&STOs,&STO_id_array,&STO_type_array,&STO_mu_array,0);
        
        calculate_overlap_S_matrix(STOs,x,y,z,STO_id_array,STO_type_array,STO_mu_array,&Smatrix,0);
        
        resolve_unique_species(atoms,species,&atom_types,&unique_species,0);
        
        //read_CP2K_GTOs(current_folder,atoms,atom_types,species,unique_species,STOproj_basis_file,STOproj_basis,&CGTOs,&bfnPerAtom,&GTO_depth,&PCGTOs,&palpha,&pcoeff,&pqn,0);
        CGTOs=current_input.CGTOs_agg[i];
        PCGTOs=current_input.PCGTOs_agg[i];
        bfnPerAtom=(int*)malloc(atoms*sizeof(int));
        GTO_depth=(int*)malloc(CGTOs*sizeof(int));
        palpha=(double*)malloc(PCGTOs*sizeof(double));
        pcoeff=(double*)malloc(PCGTOs*sizeof(double));
        pqn=(int**)malloc(PCGTOs*sizeof(int*));
        for(j=0;j<PCGTOs;++j)pqn[j]=(int*)malloc(3*sizeof(int));
        k=-1;
        for(j=current_input.pos_storage_L[i];j<=current_input.pos_storage_R[i];++j)
        {
            k=k+1;
            bfnPerAtom[k]=current_input.bfnPerAtom_agg[j];
        }
        k=-1;
        for(j=current_input.CGTOs_storage_L[i];j<=current_input.CGTOs_storage_R[i];++j)
        {
            k=k+1;
            GTO_depth[k]=current_input.GTO_depth_agg[j];
        }
        k=-1;
        for(j=current_input.PCGTOs_storage_L[i];j<=current_input.PCGTOs_storage_R[i];++j)
        {
            k=k+1;
            palpha[k]=current_input.palpha_agg[j];
            pcoeff[k]=current_input.pcoeff_agg[j];
            pqn[k][0]=current_input.pqn_agg[j][0];
            pqn[k][1]=current_input.pqn_agg[j][1];
            pqn[k][2]=current_input.pqn_agg[j][2];
        }
        
        //locate_MO(simplex_MO_file[i],simplex_MO[i],CGTOs,&MOcoeffs,0);
        MOcoeffs=(double**)malloc(CGTOs*sizeof(double*));
        for(j=0;j<CGTOs;++j)MOcoeffs[j]=(double*)malloc(2*sizeof(double));
        k=-1;
        for(j=current_input.MO_storage_L[i];j<=current_input.MO_storage_R[i];++j)
        {
            k=k+1;
            MOcoeffs[k][0]=current_input.MOcoeffs_aggr_alpha[j];
            MOcoeffs[k][1]=current_input.MOcoeffs_aggr_beta[j];
        }
        
        STO_GTO_projection(atoms,x,y,z,STOs,STO_type_array,STO_id_array,STO_mu_array,Smatrix,bfnPerAtom,GTO_depth,MOcoeffs,pcoeff,palpha,pqn,&res,0,current_input.MOLog_channel,&V_array,&s_array);
        orb_compl=0.0;for(j=0;j<STOs;++j)for(k=0;k<STOs;++k)orb_compl=orb_compl+res[j]*res[k]*Smatrix[j][k];
        
        //
        //printf("%d\t%s\t%lf\n",i+1,simplex_name[i],orb_compl);
        //
        
        (*compl_array)[i]=orb_compl;
        if(current_input.simplex_id[i]==2)
        {
            Hcounter=Hcounter+1;
            compl_array_HOMO[Hcounter]=orb_compl;
            HOMO_mean=HOMO_mean+orb_compl;
        }
        else
        {
            Lcounter=Lcounter+1;
            compl_array_LUMO[Lcounter]=orb_compl;
            LUMO_mean=LUMO_mean+orb_compl;
        }
        
        
        if(strcmp("yes",current_input.simplex_cube[i])==0)
        {
            for(j=0;j<STOs;++j)res[j]=res[j]/sqrt(fabs(orb_compl));
            STO_matrix=(double**)malloc(atoms*sizeof(double*));for(j=0;j<atoms;++j)STO_matrix[j]=(double*)malloc(4*sizeof(double));
            for(j=0;j<STOs;++j)
            {
                if(STO_type_array[j]==1){
                    STO_matrix[STO_id_array[j]-1][0]=res[j];
                    STO_matrix[STO_id_array[j]-1][1]=0.0;
                    STO_matrix[STO_id_array[j]-1][2]=0.0;
                    STO_matrix[STO_id_array[j]-1][3]=0.0;
                }
                else{
                    STO_matrix[STO_id_array[j]-1][(STO_type_array[j]-2)%4]=res[j];
                }
            }
            /*
            sprintf(file_path,"%s/STO_matrix_%s.dat",current_input.current_folder,current_input.simplex_name[i]);
            fp=fopen(file_path,"w+");
            for(j=0;j<atoms;++j)
            {
                fprintf(fp,"%lf\t%lf\t%lf\t%lf\n",STO_matrix[j][0],STO_matrix[j][1],STO_matrix[j][2],STO_matrix[j][3]);
            }
            fclose(fp);
            */
            create_cube_file(current_input.current_folder,current_input.STOproj_cube_grid,current_input.simplex_name[i],current_input.simplex_MO[i],atoms,species,x,y,z,STO_matrix,smu_per_species,pmu_per_species,0);
            
            resolve_pvecs(atoms,x,y,z,species,&px,&py,&pz);
            
            sprintf(file_path,"%s/AOM_COEFF_single_%s_MO_%s.include",current_input.current_folder,current_input.simplex_name[i],current_input.simplex_MO[i]);
            fp=fopen(file_path,"w+");
            for(j=0;j<atoms;++j)fprintf(fp,"%s %d %d %.1lf %lf\n",species[j],resolve_atomic_Z(species[j]),1,0.0,px[j]*STO_matrix[j][1]+py[j]*STO_matrix[j][2]+pz[j]*STO_matrix[j][3]);
            fclose(fp);
            
            // write state
            write_state(current_input.current_folder,current_input.simplex_name[i],current_input.simplex_MO[i],atoms,STOs,px,py,pz,Smatrix,STO_matrix,orb_compl,V_array,s_array);
            
            for(j=0;j<atoms;++j)free(STO_matrix[j]);free(STO_matrix);
            free(px);free(py);free(pz);
        }
        
        for(j=0;j<atoms;++j)free(species[j]);free(species);free(x);free(y);free(z);
        free(STO_type_array);free(STO_id_array);free(STO_mu_array);
        for(j=0;j<STOs;++j)free(Smatrix[j]);free(Smatrix);
        for(j=0;j<atom_types;++j)free(unique_species[j]);free(unique_species);
        free(bfnPerAtom);free(GTO_depth);free(palpha);free(pcoeff);for(j=0;j<PCGTOs;++j)free(pqn[j]);free(pqn);
        for(j=0;j<CGTOs;++j)free(MOcoeffs[j]);free(MOcoeffs);
        free(res);
        free(V_array);free(s_array);
    }
    HOMO_mean=HOMO_mean/current_input.HOMOs;
    LUMO_mean=LUMO_mean/current_input.LUMOs;
    HOMO_var=0.0;
    for(i=0;i<current_input.HOMOs;++i)HOMO_var=HOMO_var+(compl_array_HOMO[i]-HOMO_mean)*(compl_array_HOMO[i]-HOMO_mean);HOMO_var=HOMO_var/current_input.HOMOs;
    LUMO_var=0.0;
    for(i=0;i<current_input.LUMOs;++i)LUMO_var=LUMO_var+(compl_array_LUMO[i]-LUMO_mean)*(compl_array_LUMO[i]-LUMO_mean);LUMO_var=LUMO_var/current_input.LUMOs;
    //printf("%lf\t%lf\t%lf\t%lf\n",HOMO_mean,HOMO_var,LUMO_mean,LUMO_var);
    metric=sqrt((HOMO_mean-1.0)*(HOMO_mean-1.0)+(HOMO_var-0.0)*(HOMO_var-0.0)+(LUMO_mean-1.0)*(LUMO_mean-1.0)+(LUMO_var-0.0)*(LUMO_var-0.0));
    //printf("%lf\n",metric);
    return metric;
}
