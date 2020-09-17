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

extern void dgesv_( int* n, int* nrhs, double* a, int* lda, int* ipiv,
                  double* b, int* ldb, int* info );
extern void store_matrix( int m, int n, double* a, int lda , double *res );
void store_matrix( int m, int n, double* a, int lda , double *res ) {
    int i, j, k=-1;
    for( i = 0; i < m; i++ ) {
        for( j = 0; j < n; j++ ) { ++k; res[k]=a[i+j*lda]; }
    }
}
void STO_GTO_projection(int atoms,double *x,double *y,double *z,int STOs,int *STO_type_array,int *STO_id_array,double *STO_mu_array,double **Smatrix,int *bfnPerAtom,int *GTO_depth,double **MOcoeffs,double *pcoeff,double *palpha,int **pqn,double **res,int verb,int MOLog_channel);
void STO_GTO_projection(int atoms,double *x,double *y,double *z,int STOs,int *STO_type_array,int *STO_id_array,double *STO_mu_array,double **Smatrix,int *bfnPerAtom,int *GTO_depth,double **MOcoeffs,double *pcoeff,double *palpha,int **pqn,double **res,int verb,int MOLog_channel)
{
    
    // O-ohata GTO decomposition arrays
    int specialdepth[9] = {10, 10, 6, 6, 6, 8, 8, 8, 8};
    double oohatac1s[10] = {0.0154199, 0.203216, 0.401002, 0.313305, 0.144647, 0.049609, 0.0139344, 0.00318156, 0.000622761, 0.0000886225};
    double oohataa1s[10] = {0.0372356, 0.0848076, 0.191073, 0.458596, 1.18834, 3.34954, 10.514, 37.9276, 156.411, 1188.35};
    double oohatac2s[10] = {0.00578359, 0.225661, 0.527807, 0.265263, 0.0457229, -0.0488578, -0.0223909, -0.00600006, -0.00120626, -0.000163319};
    double oohataa2s[10] = {0.0214081, 0.0497851, 0.0991734, 0.2010660, 0.263514, 1.3073, 3.98351, 14.1623, 60.6176, 432.099};
    double oohatac2p[6] = {0.162114, 0.490328, 0.362480, 0.118961, 0.0237179, 0.00306629};
    double oohataa2p[6] = {0.0554971, 0.127302, 0.31, 0.834323, 2.56, 10.3};
    double oohatac3s[8] = {0.00621804, 0.434003, 0.678264, 0.00507343, -0.163049, -0.049145, -0.00574937, -0.000273434};
    double oohataa3s[8] = {0.014179, 0.0408735, 0.0803231, 0.201715, 0.376225, 0.995934, 3.2689, 15.8};
    double oohatac3p[8] = {0.0196080, 0.324519, 0.521925, 0.207248, 0.00504818, -0.0139732, -0.00394290, -0.000523609};
    double oohataa3p[8] = {0.0247242, 0.0510310, 0.100108, 0.207164, 0.424917, 1.24038, 3.81650, 15.8360};
    
    double oohatacLOCAL[10],oohataaLOCAL[10];
    double *V_array,S;
    int isto,Lx,Ly,Lz,ii,jj,mm,ll,i,j,k,l,m;
    
    // calculate STO-GTO overlap: V array
    V_array=(double*)malloc(STOs*sizeof(double));
    for(isto=0;isto<STOs;++isto)
    {
        if(STO_type_array[isto]==1){for(i=0;i<specialdepth[STO_type_array[isto]-1];++i){oohatacLOCAL[i]=oohatac1s[i];oohataaLOCAL[i]=oohataa1s[i];Lx=0;Ly=0;Lz=0;}}
        if(STO_type_array[isto]==2){for(i=0;i<specialdepth[STO_type_array[isto]-1];++i){oohatacLOCAL[i]=oohatac2s[i];oohataaLOCAL[i]=oohataa2s[i];Lx=0;Ly=0;Lz=0;}}
        if(STO_type_array[isto]==3){for(i=0;i<specialdepth[STO_type_array[isto]-1];++i){oohatacLOCAL[i]=oohatac2p[i];oohataaLOCAL[i]=oohataa2p[i];Lx=1;Ly=0;Lz=0;}}
        if(STO_type_array[isto]==4){for(i=0;i<specialdepth[STO_type_array[isto]-1];++i){oohatacLOCAL[i]=oohatac2p[i];oohataaLOCAL[i]=oohataa2p[i];Lx=0;Ly=1;Lz=0;}}
        if(STO_type_array[isto]==5){for(i=0;i<specialdepth[STO_type_array[isto]-1];++i){oohatacLOCAL[i]=oohatac2p[i];oohataaLOCAL[i]=oohataa2p[i];Lx=0;Ly=0;Lz=1;}}
        if(STO_type_array[isto]==6){for(i=0;i<specialdepth[STO_type_array[isto]-1];++i){oohatacLOCAL[i]=oohatac3s[i];oohataaLOCAL[i]=oohataa3s[i];Lx=0;Ly=0;Lz=0;}}
        if(STO_type_array[isto]==7){for(i=0;i<specialdepth[STO_type_array[isto]-1];++i){oohatacLOCAL[i]=oohatac3p[i];oohataaLOCAL[i]=oohataa3p[i];Lx=1;Ly=0;Lz=0;}}
        if(STO_type_array[isto]==8){for(i=0;i<specialdepth[STO_type_array[isto]-1];++i){oohatacLOCAL[i]=oohatac3p[i];oohataaLOCAL[i]=oohataa3p[i];Lx=0;Ly=1;Lz=0;}}
        if(STO_type_array[isto]==9){for(i=0;i<specialdepth[STO_type_array[isto]-1];++i){oohatacLOCAL[i]=oohatac3p[i];oohataaLOCAL[i]=oohataa3p[i];Lx=0;Ly=0;Lz=1;}}
        S=0.0;
        mm=-1;
        for(ll=0;ll<specialdepth[STO_type_array[isto]-1];++ll)
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
                        S=S+oohatacLOCAL[mm]*MOcoeffs[k][MOLog_channel]*pcoeff[m]*anIres(x[STO_id_array[isto]-1],y[STO_id_array[isto]-1],z[STO_id_array[isto]-1],STO_mu_array[STO_type_array[isto]-1]*STO_mu_array[STO_type_array[isto]-1]*oohataaLOCAL[mm],Lx,Ly,Lz,x[i],y[i],z[i],palpha[m],pqn[m][0],pqn[m][1],pqn[m][2]);
                    }
                }
            }
        }
        V_array[isto]=S;
    }
    if(verb==-1){
        // console output
        printf("STO-GTO overlap constant array V:\n");for(i=0;i<STOs;++i)printf("%d\t%.16lf\n",i+1,V_array[i]);
        printf("\n--------------------------------------------------------------------------\n\n");
    }
    
    int N,NRHS,LDA,LDB;
    double *a,*b;
    N=STOs;
    NRHS=1;
    LDA=N;
    LDB=N;
    int n = N, nrhs = NRHS, lda = LDA, ldb = LDB, info;
    int ipiv[N];
    a=(double*)malloc((N*N)*sizeof(double));
    b=(double*)malloc(N*sizeof(double));
    k=-1;
    for(jj=0;jj<N;++jj)
    {
        for(ii=0;ii<N;++ii)
        {
            k=k+1;
            a[k]=Smatrix[ii][jj];
            
        }
    }
    for(jj=0;jj<N;++jj)b[jj]=V_array[jj];
    dgesv_( &n, &nrhs, a, &lda, ipiv, b, &ldb, &info );
    *res=(double*)malloc(STOs*sizeof(double));
    store_matrix( n, nrhs, b, ldb, *res );
    if(verb==-1){
        // console output
        printf("Solution to Smatrix * X = V:\n");
        for(i=0;i<STOs;++i)printf("%.16lf\n",(*res)[i]);
        printf("\n--------------------------------------------------------------------------\n\n");
    }
    free(a);free(b);
    
    free(V_array);
}
