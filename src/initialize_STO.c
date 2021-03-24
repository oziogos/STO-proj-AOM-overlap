/* -----------------------------------------------------------------------------
 
 Unified STO projection / AOM overlap code
 
 Alpha version: 0.1
 1-Oct-2019
 
 Orestis George Ziogos, UCL
 o.ziogos@ucl.ac.uk
 
 For more information, examine the README file in the parent directory.
 
-------------------------------------------------------------------------------- */
#include"general.h"
int resolveSTOnum(char *species);
double resolveSTOmu(char *species,int type,double *smu_per_species,double *pmu_per_species);
int resolve_atomic_Z(char *species);
void initialize_STO(int atoms,char **species,double *smu_per_species,double *pmu_per_species,int *STOs,int **STO_id_array,int **STO_type_array,double **STO_mu_array,int verb);
void initialize_STO(int atoms,char **species,double *smu_per_species,double *pmu_per_species,int *STOs,int **STO_id_array,int **STO_type_array,double **STO_mu_array,int verb)
{
    int i,j;
    
    *STOs=0;for(i=0;i<atoms;++i)*STOs=*STOs+resolveSTOnum(species[i]);
    *STO_id_array=(int*)malloc(*STOs*sizeof(int));*STO_type_array=(int*)malloc(*STOs*sizeof(int));*STO_mu_array=(double*)malloc(*STOs*sizeof(double));
    j=-1;
    for(i=0;i<atoms;++i)
    {
        if(strcmp(species[i],"H")==0){
            ++j;(*STO_id_array)[j]=i+1;(*STO_type_array)[j]=1;
        }
        if(strcmp(species[i],"C")==0||
           strcmp(species[i],"N")==0||
           strcmp(species[i],"O")==0||
           strcmp(species[i],"F")==0){
            ++j;(*STO_id_array)[j]=i+1;(*STO_type_array)[j]=2;
            ++j;(*STO_id_array)[j]=i+1;(*STO_type_array)[j]=3;
            ++j;(*STO_id_array)[j]=i+1;(*STO_type_array)[j]=4;
            ++j;(*STO_id_array)[j]=i+1;(*STO_type_array)[j]=5;
        }
        if(strcmp(species[i],"S")==0){
            ++j;(*STO_id_array)[j]=i+1;(*STO_type_array)[j]=6;
            ++j;(*STO_id_array)[j]=i+1;(*STO_type_array)[j]=7;
            ++j;(*STO_id_array)[j]=i+1;(*STO_type_array)[j]=8;
            ++j;(*STO_id_array)[j]=i+1;(*STO_type_array)[j]=9;
        }
    }
    for(i=0;i<*STOs;++i)(*STO_mu_array)[i]=resolveSTOmu(species[(*STO_id_array)[i]-1],(*STO_type_array)[i],smu_per_species,pmu_per_species);
    if(verb==-1){
        // console output
        printf("Number of total STO basis functions: %d\n",*STOs);
        printf("Atomic decomposition:\n");
        printf("[id]\tSTO\torbital\tmu\ttype\n");
        for(i=0;i<*STOs;++i){printf("[%d]\t%d\t",(*STO_id_array)[i],i+1);
            if((*STO_type_array)[i]==1)printf("1s\t");
            if((*STO_type_array)[i]==2)printf("2s\t");
            if((*STO_type_array)[i]==3)printf("2px\t");
            if((*STO_type_array)[i]==4)printf("2py\t");
            if((*STO_type_array)[i]==5)printf("2pz\t");
            if((*STO_type_array)[i]==6)printf("3s\t");
            if((*STO_type_array)[i]==7)printf("3px\t");
            if((*STO_type_array)[i]==8)printf("3py\t");
            if((*STO_type_array)[i]==9)printf("3pz\t");
            printf("%lf\t%d\n",(*STO_mu_array)[i],(*STO_type_array)[i]);
        }
        printf("\n--------------------------------------------------------------------------\n\n");
    }
}
int resolveSTOnum(char *species);
int resolveSTOnum(char *species)
{
    int STOs=0;
    if(strcmp(species,"H")==0)STOs=1;
    if(strcmp(species,"C")==0)STOs=1+3;
    if(strcmp(species,"N")==0)STOs=1+3;
    if(strcmp(species,"O")==0)STOs=1+3;
    if(strcmp(species,"F")==0)STOs=1+3;
    if(strcmp(species,"S")==0)STOs=1+3;
    return STOs;
}
double resolveSTOmu(char *species,int type,double *smu_per_species,double *pmu_per_species);
double resolveSTOmu(char *species,int type,double *smu_per_species,double *pmu_per_species)
{
    double res=0.0;
    if(strcmp(species,"H")==0&&type==1)res=smu_per_species[resolve_atomic_Z(species)-1];
    if(strcmp(species,"C")==0&&type==2)res=smu_per_species[resolve_atomic_Z(species)-1];
    if(strcmp(species,"C")==0&&type==3)res=pmu_per_species[resolve_atomic_Z(species)-1];
    if(strcmp(species,"C")==0&&type==4)res=pmu_per_species[resolve_atomic_Z(species)-1];
    if(strcmp(species,"C")==0&&type==5)res=pmu_per_species[resolve_atomic_Z(species)-1];
    if(strcmp(species,"O")==0&&type==2)res=smu_per_species[resolve_atomic_Z(species)-1];
    if(strcmp(species,"O")==0&&type==3)res=pmu_per_species[resolve_atomic_Z(species)-1];
    if(strcmp(species,"O")==0&&type==4)res=pmu_per_species[resolve_atomic_Z(species)-1];
    if(strcmp(species,"O")==0&&type==5)res=pmu_per_species[resolve_atomic_Z(species)-1];
    if(strcmp(species,"N")==0&&type==2)res=smu_per_species[resolve_atomic_Z(species)-1];
    if(strcmp(species,"N")==0&&type==3)res=pmu_per_species[resolve_atomic_Z(species)-1];
    if(strcmp(species,"N")==0&&type==4)res=pmu_per_species[resolve_atomic_Z(species)-1];
    if(strcmp(species,"N")==0&&type==5)res=pmu_per_species[resolve_atomic_Z(species)-1];
    if(strcmp(species,"F")==0&&type==2)res=smu_per_species[resolve_atomic_Z(species)-1];
    if(strcmp(species,"F")==0&&type==3)res=pmu_per_species[resolve_atomic_Z(species)-1];
    if(strcmp(species,"F")==0&&type==4)res=pmu_per_species[resolve_atomic_Z(species)-1];
    if(strcmp(species,"F")==0&&type==5)res=pmu_per_species[resolve_atomic_Z(species)-1];
    if(strcmp(species,"S")==0&&type==6)res=smu_per_species[resolve_atomic_Z(species)-1];
    if(strcmp(species,"S")==0&&type==7)res=pmu_per_species[resolve_atomic_Z(species)-1];
    if(strcmp(species,"S")==0&&type==8)res=pmu_per_species[resolve_atomic_Z(species)-1];
    if(strcmp(species,"S")==0&&type==9)res=pmu_per_species[resolve_atomic_Z(species)-1];
    return res;
}
