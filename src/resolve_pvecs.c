/* -----------------------------------------------------------------------------
 
 Unified STO projection / AOM overlap code
 
 Alpha version: 0.1
 1-Oct-2019
 
 Orestis George Ziogos, UCL
 o.ziogos@ucl.ac.uk
 
 For more information, examine the README file in the parent directory.
 
-------------------------------------------------------------------------------- */
#include"general.h"
int resolve_atomic_Z(char *species);
void resolve_pvecs(int atoms,double *x,double *y,double *z,char **species,double **px,double **py,double **pz);
void resolve_pvecs(int atoms,double *x,double *y,double *z,char **species,double **px,double **py,double **pz)
{
    int i,j,k,aoinum,*atomlist,*active;
    double neighbourlist[4][4],dist;
    double veca_x,veca_y,veca_z,vecb_x,vecb_y,vecb_z,ux,uy,uz;
    
    // deduce active atoms: right now just consider non-H atoms...
    aoinum=0;
    atomlist=(int*)malloc(atoms*sizeof(int));
    active=(int*)malloc(atoms*sizeof(int));
    for(i=0;i<atoms;++i){active[i]=0;if(resolve_atomic_Z(species[i])!=1){active[i]=1;++aoinum;atomlist[aoinum-1]=i;}}
    
    // calculate p vectors
    *px=(double*)malloc(atoms*sizeof(double));
    *py=(double*)malloc(atoms*sizeof(double));
    *pz=(double*)malloc(atoms*sizeof(double));
    for(i=0;i<atoms;++i){(*px)[i]=0.0;(*py)[i]=0.0;(*pz)[i]=0.0;}
    for(i=0;i<aoinum;++i)
    {
        j=0;
        for(k=0;k<atoms;++k)
        {
            dist=(x[atomlist[i]]-x[k])*(x[atomlist[i]]-x[k])+(y[atomlist[i]]-y[k])*(y[atomlist[i]]-y[k])+(z[atomlist[i]]-z[k])*(z[atomlist[i]]-z[k]);dist=sqrt(dist);
            if(atomlist[i]!=k && dist<3.5)
            {
                neighbourlist[1][j]=x[k];
                neighbourlist[2][j]=y[k];
                neighbourlist[3][j]=z[k];
                ++j;
            }
        }
        if(j==2)
        {
            neighbourlist[1][2]=x[atomlist[i]];
            neighbourlist[2][2]=y[atomlist[i]];
            neighbourlist[3][2]=z[atomlist[i]];
        }
        veca_x=neighbourlist[1][0]-neighbourlist[1][2];
        veca_y=neighbourlist[2][0]-neighbourlist[2][2];
        veca_z=neighbourlist[3][0]-neighbourlist[3][2];
        vecb_x=neighbourlist[1][1]-neighbourlist[1][2];
        vecb_y=neighbourlist[2][1]-neighbourlist[2][2];
        vecb_z=neighbourlist[3][1]-neighbourlist[3][2];
        ux=-veca_z*vecb_y + veca_y*vecb_z;
        uy= veca_z*vecb_x - veca_x*vecb_z;
        uz=-veca_y*vecb_x + veca_x*vecb_y;
        dist=ux*ux+uy*uy+uz*uz;
        dist=sqrt(dist);
        ux=ux/dist;uy=uy/dist;uz=uz/dist;
        (*px)[atomlist[i]]=ux;
        (*py)[atomlist[i]]=uy;
        (*pz)[atomlist[i]]=uz;
    }
    free(atomlist);free(active);
}
