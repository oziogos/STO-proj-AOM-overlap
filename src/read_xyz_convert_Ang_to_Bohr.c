/* -----------------------------------------------------------------------------
 
 Unified STO projection / AOM overlap code
 
 Alpha version: 0.1
 1-Oct-2019
 
 Orestis George Ziogos, UCL
 o.ziogos@ucl.ac.uk
 
 For more information, examine the README file in the parent directory.
 
-------------------------------------------------------------------------------- */
#include"general.h"
void read_xyz_convert_Ang_to_Bohr(char *STOproj_molecule,int *atoms,double **x,double **y,double **z,char ***species,int verb);
void read_xyz_convert_Ang_to_Bohr(char *STOproj_molecule,int *atoms,double **x,double **y,double **z,char ***species,int verb)
{
    FILE *fp;
    char file_path[cmax_length],buffer[cmax_length];
    int i;
    
    sprintf(file_path,"%s",STOproj_molecule);
    fp=fopen(file_path,"r");if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
    fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&(*atoms));
    *x=(double*)malloc((*atoms)*sizeof(double));
    *y=(double*)malloc((*atoms)*sizeof(double));
    *z=(double*)malloc((*atoms)*sizeof(double));
    *species=(char**)malloc((*atoms)*sizeof(char*));for(i=0;i<(*atoms);++i)(*species)[i]=(char*)malloc(sub_length*sizeof(char));
    fgets(buffer,cmax_length,fp);
    for(i=0;i<(*atoms);++i){fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%lf\t%lf\t%lf",(*species)[i],&(*x)[i],&(*y)[i],&(*z)[i]);}
    fclose(fp);
    if(verb==-1){
        // console output
        printf("Number of atoms: %d\n",*atoms);
        printf("[Atomic id], species and coordinates (Ang|Bohr):\n");for(i=0;i<*atoms;++i)printf("[%d]\t%s\t%lf\t%lf\t%lf\t|\t%lf\t%lf\t%lf\n",i+1,(*species)[i],(*x)[i],(*y)[i],(*z)[i],(*x)[i]*AngToBohr,(*y)[i]*AngToBohr,(*z)[i]*AngToBohr);
        printf("\n--------------------------------------------------------------------------\n\n");
    }
    // convert Ang to Bohr
    for(i=0;i<*atoms;++i){(*x)[i]=(*x)[i]*AngToBohr;(*y)[i]=(*y)[i]*AngToBohr;(*z)[i]=(*z)[i]*AngToBohr;}
}
