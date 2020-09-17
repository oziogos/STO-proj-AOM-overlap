/* -----------------------------------------------------------------------------
 
 Unified STO projection / AOM overlap code
 
 Alpha version: 0.1
 1-Oct-2019
 
 Orestis George Ziogos, UCL
 o.ziogos@ucl.ac.uk
 
 For more information, examine the README file in the parent directory.
 
-------------------------------------------------------------------------------- */
#include"general.h"
double anIres(double X1,double Y1,double Z1,double a1,int lx1,int ly1,int lz1,double X2,double Y2,double Z2,double a2,int lx2,int ly2,int lz2);
double GTO_overlap(int atoms,double *x,double *y,double *z,int *bfnPerAtom,int *GTO_depth,double **MOcoeffs,double *pcoeff,double *palpha,int **pqn,int MOLog_channel);
double GTO_overlap(int atoms,double *x,double *y,double *z,int *bfnPerAtom,int *GTO_depth,double **MOcoeffs,double *pcoeff,double *palpha,int **pqn,int MOLog_channel)
{
    int i,j,k,l,m,ii,jj,kk,ll,mm;
    double S;
    S=0.0;
    kk=-1;mm=-1;
    for(ii=0;ii<atoms;++ii)
    {
        for(jj=0;jj<bfnPerAtom[ii];++jj)
        {
            ++kk;
            for(ll=0;ll<GTO_depth[kk];++ll)
            {
                ++mm;
                k=-1;m=-1;
                for(i=0;i<atoms;++i)
                {
                    for(j=0;j<bfnPerAtom[i];++j)
                    {
                        ++k;
                        for(l=0;l<GTO_depth[k];++l)
                        {
                            ++m;
                            S=S+MOcoeffs[kk][MOLog_channel]*pcoeff[mm]*MOcoeffs[k][MOLog_channel]*pcoeff[m]*anIres(x[ii],y[ii],z[ii],palpha[mm],pqn[mm][0],pqn[mm][1],pqn[mm][2],x[i],y[i],z[i],palpha[m],pqn[m][0],pqn[m][1],pqn[m][2]);
                        }
                    }
                }
            }
        }
    }
    return S;
}
