/* -----------------------------------------------------------------------------
 
 Unified STO projection / AOM overlap code
 
 Beta version: 1.0
 1-Mar-2019
 
 Orestis George Ziogos, UCL
 o.ziogos@ucl.ac.uk
 
 For more information, examine the README file in the parent directory.
 
 -------------------------------------------------------------------------------- */
#include"general.h"
int resolve_atomic_Z(char *species);
void create_cube_file(char *current_folder,char *STOproj_cube_grid,char *STOproj_name,char *STOproj_MO,int atoms,char **species,double *x,double *y,double *z,double **STO_matrix,double *smu_per_species,double *pmu_per_species,int verb);
void AOM_frag_cubes(char *current_folder,int atoms,int frag1atoms,double *x,double *y,double *z,char **species,double **STO_matrix,double *smu_per_species,double *pmu_per_species,int verb);
void AOM_frag_cubes(char *current_folder,int atoms,int frag1atoms,double *x,double *y,double *z,char **species,double **STO_matrix,double *smu_per_species,double *pmu_per_species,int verb)
{
    FILE *fp,*fpw;
    char buffer[cmax_length];
    int i,j,localatoms;
    double d1,d2,d3;
    double *localx,*localy,*localz,**localSTO_matrix;
    char **localspecies;
    // frag 1
    localatoms=frag1atoms;
    localx=(double*)malloc(localatoms*sizeof(double));
    localy=(double*)malloc(localatoms*sizeof(double));
    localz=(double*)malloc(localatoms*sizeof(double));
    localSTO_matrix=(double**)malloc(localatoms*sizeof(double*));
    for(i=0;i<localatoms;++i)localSTO_matrix[i]=(double*)malloc(4*sizeof(double));
    localspecies=(char**)malloc(localatoms*sizeof(char*));for(i=0;i<localatoms;++i)localspecies[i]=(char*)malloc(sub_length*sizeof(char));
    for(i=0;i<localatoms;++i)
    {
        localx[i]=x[i];
        localy[i]=y[i];
        localz[i]=z[i];
        sprintf(localspecies[i],"%s",species[i]);
        for(j=0;j<4;++j)localSTO_matrix[i][j]=STO_matrix[i][j];
        //printf("$$ [%d]\t%s\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",i+1,species[i],x[i],y[i],z[i],STO_matrix[i][0],STO_matrix[i][1],STO_matrix[i][2],STO_matrix[i][3]);
    }
    create_cube_file(current_folder,"0.25","single_frag","1",localatoms,localspecies,localx,localy,localz,localSTO_matrix,smu_per_species,pmu_per_species,verb);
    fp=fopen("single_frag_STO_MO_1.cube","r");if(fp==NULL){printf("Could not locate %s\n","single_frag_STO_MO_1.cube");exit(-1);}
    fpw=fopen("frag_STO_MO_1.cube","w+");
    for(i=0;i<2;++i){fgets(buffer,cmax_length,fp);fprintf(fpw,"%s",buffer);}
    fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%lf\t%lf\t%lf",&j,&d1,&d2,&d3);fprintf(fpw,"%d\t%lf\t%lf\t%lf\n",atoms,d1,d2,d3);
    for(i=0;i<3;++i){fgets(buffer,cmax_length,fp);fprintf(fpw,"%s",buffer);}
    for(i=0;i<localatoms;++i)fgets(buffer,cmax_length,fp);
    for(i=0;i<atoms;++i){fprintf(fpw,"%d\t%lf\t%lf\t%lf\t%lf\n",resolve_atomic_Z(species[i]),(double)resolve_atomic_Z(species[i]),x[i],y[i],z[i]);}
    while(fgets(buffer,cmax_length,fp)!=NULL)fprintf(fpw,"%s",buffer);
    fclose(fpw);
    fclose(fp);
    free(localx);free(localy);free(localz);for(i=0;i<localatoms;++i){free(localSTO_matrix[i]);free(localspecies[i]);}free(localSTO_matrix);free(localspecies);
    // frag 2
    localatoms=atoms-frag1atoms;
    localx=(double*)malloc(localatoms*sizeof(double));
    localy=(double*)malloc(localatoms*sizeof(double));
    localz=(double*)malloc(localatoms*sizeof(double));
    localSTO_matrix=(double**)malloc(localatoms*sizeof(double*));
    for(i=0;i<localatoms;++i)localSTO_matrix[i]=(double*)malloc(4*sizeof(double));
    localspecies=(char**)malloc(localatoms*sizeof(char*));for(i=0;i<localatoms;++i)localspecies[i]=(char*)malloc(sub_length*sizeof(char));
    for(i=0;i<localatoms;++i)
    {
        localx[i]=x[i+frag1atoms];
        localy[i]=y[i+frag1atoms];
        localz[i]=z[i+frag1atoms];
        sprintf(localspecies[i],"%s",species[i+frag1atoms]);
        for(j=0;j<4;++j)localSTO_matrix[i][j]=STO_matrix[i+frag1atoms][j];
        //printf("$$ [%d]\t%s\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",i+1,species[i+frag1atoms],x[i+frag1atoms],y[i+frag1atoms],z[i+frag1atoms],STO_matrix[i+frag1atoms][0],STO_matrix[i+frag1atoms][1],STO_matrix[i+frag1atoms][2],STO_matrix[i+frag1atoms][3]);
    }
    create_cube_file(current_folder,"0.25","single_frag","2",localatoms,localspecies,localx,localy,localz,localSTO_matrix,smu_per_species,pmu_per_species,verb);
    fp=fopen("single_frag_STO_MO_2.cube","r");if(fp==NULL){printf("Could not locate %s\n","single_frag_STO_MO_2.cube");exit(-1);}
    fpw=fopen("frag_STO_MO_2.cube","w+");
    for(i=0;i<2;++i){fgets(buffer,cmax_length,fp);fprintf(fpw,"%s",buffer);}
    fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%lf\t%lf\t%lf",&j,&d1,&d2,&d3);fprintf(fpw,"%d\t%lf\t%lf\t%lf\n",atoms,d1,d2,d3);
    for(i=0;i<3;++i){fgets(buffer,cmax_length,fp);fprintf(fpw,"%s",buffer);}
    for(i=0;i<localatoms;++i)fgets(buffer,cmax_length,fp);
    for(i=0;i<atoms;++i){fprintf(fpw,"%d\t%lf\t%lf\t%lf\t%lf\n",resolve_atomic_Z(species[i]),(double)resolve_atomic_Z(species[i]),x[i],y[i],z[i]);}
    while(fgets(buffer,cmax_length,fp)!=NULL)fprintf(fpw,"%s",buffer);
    fclose(fpw);
    fclose(fp);
    free(localx);free(localy);free(localz);for(i=0;i<localatoms;++i){free(localSTO_matrix[i]);free(localspecies[i]);}free(localSTO_matrix);free(localspecies);
}
