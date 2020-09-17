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
double AOM_overlap_calculation(int istart,int istop,int jstart,int jstop,double *x,double *y,double *z,int *STO_id_array,int *STO_type_array,double *STO_mu_array,double **STO_matrix);
double AOM_overlap_calculation(int istart,int istop,int jstart,int jstop,double *x,double *y,double *z,int *STO_id_array,int *STO_type_array,double *STO_mu_array,double **STO_matrix)
{
    int locali,localj,i,j;
    double S;
    S=0.0;
    for(i=istart;i<istop;++i)
    {
        locali=(STO_type_array[i]-2)%4;
        for(j=jstart;j<jstop;++j)
        {
            localj=(STO_type_array[j]-2)%4;
            if(locali>0&&localj>0)
            {
                if(STO_id_array[i]!=STO_id_array[j])
                {
                    S=S+STO_matrix[STO_id_array[i]-1][locali]*STO_matrix[STO_id_array[j]-1][localj]*overlap(x[STO_id_array[i]-1],y[STO_id_array[i]-1],z[STO_id_array[i]-1],
                                                                                                            x[STO_id_array[j]-1],y[STO_id_array[j]-1],z[STO_id_array[j]-1],
                                                                                                            STO_mu_array[i],
                                                                                                            STO_mu_array[j],
                                                                                                            STO_type_array[i],
                                                                                                            STO_type_array[j]);
                }
                else
                {
                    if(STO_type_array[i]==STO_type_array[j])S=S+STO_matrix[STO_id_array[i]-1][locali]*STO_matrix[STO_id_array[j]-1][localj];
                }
            }
        }
    }
    return S;
}
