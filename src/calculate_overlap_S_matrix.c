/* -----------------------------------------------------------------------------
 
 Unified STO projection / AOM overlap code
 
 Alpha version: 0.1
 1-Oct-2019
 
 Orestis George Ziogos, UCL
 o.ziogos@ucl.ac.uk
 
 For more information, examine the README file in the parent directory.
 
-------------------------------------------------------------------------------- */
#include"general.h"
double overlap(double X1,double Y1,double Z1,double X2,double Y2,double Z2,double mu1,double mu2,int type1,int type2);
void calculate_overlap_S_matrix(int STOs,double *x,double *y,double *z,int *STO_id_array,int *STO_type_array,double *STO_mu_array,double ***Smatrix,int verb);
void calculate_overlap_S_matrix(int STOs,double *x,double *y,double *z,int *STO_id_array,int *STO_type_array,double *STO_mu_array,double ***Smatrix,int verb)
{
    int i,j;
    // NOTE for beta version: the matrix is symmetric; use overlap() for the upper triangular part and then mirror...
    *Smatrix=(double**)malloc(STOs*sizeof(double*));for(i=0;i<STOs;++i)(*Smatrix)[i]=(double*)malloc(STOs*sizeof(double));
    for(i=0;i<STOs;++i)for(j=0;j<STOs;++j)if(i!=j){(*Smatrix)[i][j]=0.0;}else{(*Smatrix)[i][j]=1.0;}
    for(i=0;i<STOs;++i)
    {
        for(j=0;j<STOs;++j)
        {
            if(STO_id_array[i]!=STO_id_array[j])
            {
                (*Smatrix)[i][j]=overlap(
                                         x[STO_id_array[i]-1],y[STO_id_array[i]-1],z[STO_id_array[i]-1],
                                         x[STO_id_array[j]-1],y[STO_id_array[j]-1],z[STO_id_array[j]-1],
                                         STO_mu_array[i],
                                         STO_mu_array[j],
                                         STO_type_array[i],
                                         STO_type_array[j]);
            }
        }
    }
    if(verb==-1){
        // console output
        printf("Single molecule STO overlap matrix (Smatrix):\n");
        for(i=0;i<STOs;++i){for(j=0;j<STOs;++j){printf("%.4lf\t",(*Smatrix)[i][j]);}printf("\n");}
        printf("\n--------------------------------------------------------------------------\n\n");
    }
}
