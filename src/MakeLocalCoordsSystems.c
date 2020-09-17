/* -----------------------------------------------------------------------------
 
 Unified STO projection / AOM overlap code
 
 Alpha version: 0.1
 1-Oct-2019
 
 Orestis George Ziogos, UCL
 o.ziogos@ucl.ac.uk
 
 For more information, examine the README file in the parent directory.
 
-------------------------------------------------------------------------------- */
#include"general.h"
int equal(double A,double B);
void cross_product(double ax,double ay,double az,double bx,double by,double bz,double *crossx,double *crossy,double *crossz);
void MakeLocalCoordsSystems(double X1,double Y1,double Z1,double X2,double Y2,double Z2,double *exb_X,double *exb_Y,double *exb_Z,double *eyb_X,double *eyb_Y,double *eyb_Z,double *ezb_X,double *ezb_Y,double *ezb_Z);
void MakeLocalCoordsSystems(double X1,double Y1,double Z1,double X2,double Y2,double Z2,double *exb_X,double *exb_Y,double *exb_Z,double *eyb_X,double *eyb_Y,double *eyb_Z,double *ezb_X,double *ezb_Y,double *ezb_Z)
{
    int align=0,maxpos,i;
    double ez[3],norm,ey[3],absmax,t[3],ex[3];
    ez[0]=X2-X1;ez[1]=Y2-Y1;ez[2]=Z2-Z1;
    norm=ez[0]*ez[0]+ez[1]*ez[1]+ez[2]*ez[2];norm=sqrt(norm);ez[0]=ez[0]/norm;ez[1]=ez[1]/norm;ez[2]=ez[2]/norm;
    if(equal(ez[1],0.0)==1 && equal(ez[2],0.0)==1){ey[0]=0.0;ey[1]=0.0;ey[2]=1.0;align=1;}
    if(equal(ez[0],0.0)==1 && equal(ez[2],0.0)==1){ey[0]=1.0;ey[1]=0.0;ey[2]=0.0;align=1;}
    if(equal(ez[0],0.0)==1 && equal(ez[1],0.0)==1){ey[0]=0.0;ey[1]=1.0;ey[2]=0.0;align=1;}
    if(align==0)
    {
        maxpos=0;absmax=fabs(ez[0]);for(i=1;i<3;++i){if(fabs(ez[i])>absmax){maxpos=i;absmax=fabs(ez[i]);}}
        for(i=0;i<3;++i)
        {
            if(i!=maxpos){t[i]=1.0;}else{t[i]=-1.0*(ez[(i+1)%3+1-1]+ez[(i+1+1)%3+1-1])/ez[i];}
        }
        ey[0]=t[0];ey[1]=t[1];ey[2]=t[2];norm=ey[0]*ey[0]+ey[1]*ey[1]+ey[2]*ey[2];norm=sqrt(norm);ey[0]=ey[0]/norm;ey[1]=ey[1]/norm;ey[2]=ey[2]/norm;
    }
    cross_product(ey[0],ey[1],ey[2],ez[0],ez[1],ez[2],&ex[0],&ex[1],&ex[2]);
    *exb_X=ex[0];*exb_Y=ex[1];*exb_Z=ex[2];*eyb_X=ey[0];*eyb_Y=ey[1];*eyb_Z=ey[2];*ezb_X=-ez[0];*ezb_Y=-ez[1];*ezb_Z=-ez[2];
}
void cross_product(double ax,double ay,double az,double bx,double by,double bz,double *crossx,double *crossy,double *crossz);
void cross_product(double ax,double ay,double az,double bx,double by,double bz,double *crossx,double *crossy,double *crossz)
{
    *crossx=ay*bz-az*by;
    *crossy=az*bx-ax*bz;
    *crossz=ax*by-ay*bx;
}
