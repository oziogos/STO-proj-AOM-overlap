/* -----------------------------------------------------------------------------
 
 Unified STO projection / AOM overlap code
 
 Alpha version: 0.1
 1-Oct-2019
 
 Orestis George Ziogos, UCL
 o.ziogos@ucl.ac.uk
 
 For more information, examine the README file in the parent directory.
 
-------------------------------------------------------------------------------- */
#include"general.h"
void resolve_pvecs(int atoms,double *x,double *y,double *z,char **species,double **px,double **py,double **pz);
void resolve_pvecs_per_fragment(int atoms,double *x,double *y,double *z,char **species,double **px,double **py,double **pz,int frag1atoms);
void resolve_pvecs_per_fragment(int atoms,double *x,double *y,double *z,char **species,double **px,double **py,double **pz,int frag1atoms)
{
    int i,j,frag2atoms,localatoms;
    double *xfrag,*yfrag,*zfrag,*pxfrag,*pyfrag,*pzfrag;
    char **speciesfrag;
    
    frag2atoms=atoms-frag1atoms;
    *px=(double*)malloc(atoms*sizeof(double));
    *py=(double*)malloc(atoms*sizeof(double));
    *pz=(double*)malloc(atoms*sizeof(double));
    
    localatoms=frag1atoms;
    speciesfrag=(char**)malloc(localatoms*sizeof(char*));for(i=0;i<localatoms;++i)speciesfrag[i]=(char*)malloc(sub_length*sizeof(char));
    xfrag=(double*)malloc(localatoms*sizeof(double));
    yfrag=(double*)malloc(localatoms*sizeof(double));
    zfrag=(double*)malloc(localatoms*sizeof(double));
    j=-1;for(i=0;i<frag1atoms;++i){++j;xfrag[j]=x[i];yfrag[j]=y[i];zfrag[j]=z[i];sprintf(speciesfrag[j],"%s",species[i]);}
    resolve_pvecs(localatoms,xfrag,yfrag,zfrag,speciesfrag,&pxfrag,&pyfrag,&pzfrag);
    j=-1;for(i=0;i<frag1atoms;++i){++j;(*px)[i]=pxfrag[j];(*py)[i]=pyfrag[j];(*pz)[i]=pzfrag[j];}
    for(i=0;i<localatoms;++i)free(speciesfrag[i]);free(speciesfrag);
    free(xfrag);free(yfrag);free(zfrag);
    free(pxfrag);free(pyfrag);free(pzfrag);
    
    localatoms=frag2atoms;
    speciesfrag=(char**)malloc(localatoms*sizeof(char*));for(i=0;i<localatoms;++i)speciesfrag[i]=(char*)malloc(sub_length*sizeof(char));
    xfrag=(double*)malloc(localatoms*sizeof(double));
    yfrag=(double*)malloc(localatoms*sizeof(double));
    zfrag=(double*)malloc(localatoms*sizeof(double));
    j=-1;for(i=frag1atoms;i<atoms;++i){++j;xfrag[j]=x[i];yfrag[j]=y[i];zfrag[j]=z[i];sprintf(speciesfrag[j],"%s",species[i]);}
    resolve_pvecs(localatoms,xfrag,yfrag,zfrag,speciesfrag,&pxfrag,&pyfrag,&pzfrag);
    j=-1;for(i=frag1atoms;i<atoms;++i){++j;(*px)[i]=pxfrag[j];(*py)[i]=pyfrag[j];(*pz)[i]=pzfrag[j];}
    for(i=0;i<localatoms;++i)free(speciesfrag[i]);free(speciesfrag);
    free(xfrag);free(yfrag);free(zfrag);
    free(pxfrag);free(pyfrag);free(pzfrag);
}
