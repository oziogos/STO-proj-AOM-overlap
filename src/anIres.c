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
double anIres(double X1,double Y1,double Z1,double a1,int lx1,int ly1,int lz1,double X2,double Y2,double Z2,double a2,int lx2,int ly2,int lz2)
{
    double res;
    double E=2.718281828459;
    //0 0 0 - 0 0 0
    if(lx1==0&&ly1==0&&lz1==0&&lx2==0&&ly2==0&&lz2==0)res=(2*sqrt (2)*pow (a1,0.75)*pow (a2,0.75))/(pow (a1 + a2,1.5)*pow (E,(a1*a2*(pow (X1 - X2,2) + pow (Y1 - Y2,2) + pow (Z1 - Z2,2)))/(a1 + a2)));
    //0 0 0 - 0 0 1
    if(lx1==0&&ly1==0&&lz1==0&&lx2==0&&ly2==0&&lz2==1)res=(4*sqrt (2)*pow (a1,1.75)*pow (a2,1.25)*(Z1 - Z2))/(pow (a1 + a2,2.5)*pow (E,(a1*a2*(pow (X1 - X2,2) + pow (Y1 - Y2,2) + pow (Z1 - Z2,2)))/(a1 + a2)));
    //0 0 0 - 0 0 2
    if(lx1==0&&ly1==0&&lz1==0&&lx2==0&&ly2==0&&lz2==2)res=(4*sqrt (0.6666666666666666)*pow (a1,0.75)*pow (a2,0.75)*sqrt (pow (a2,2))*(a1 + a2 + 2*pow (a1,2)*pow (Z1 - Z2,2)))/(pow (a1 + a2,3.5)*pow (E,(a1*a2*(pow (X1 - X2,2) + pow (Y1 - Y2,2) + pow (Z1 - Z2,2)))/(a1 + a2)));
    //0 0 0 - 0 1 0
    if(lx1==0&&ly1==0&&lz1==0&&lx2==0&&ly2==1&&lz2==0)res=(4*sqrt (2)*pow (a1,1.75)*pow (a2,1.25)*(Y1 - Y2))/(pow (a1 + a2,2.5)*pow (E,(a1*a2*(pow (X1 - X2,2) + pow (Y1 - Y2,2) + pow (Z1 - Z2,2)))/(a1 + a2)));
    //0 0 0 - 0 1 1
    if(lx1==0&&ly1==0&&lz1==0&&lx2==0&&ly2==1&&lz2==1)res=(8*sqrt (2)*pow (a1,2.75)*pow (a2,0.75)*sqrt (pow (a2,2))*(Y1 - Y2)*(Z1 - Z2))/(pow (a1 + a2,3.5)*pow (E,(a1*a2*(pow (X1 - X2,2) + pow (Y1 - Y2,2) + pow (Z1 - Z2,2)))/(a1 + a2)));
    //0 0 0 - 0 2 0
    if(lx1==0&&ly1==0&&lz1==0&&lx2==0&&ly2==2&&lz2==0)res=(4*sqrt (0.6666666666666666)*pow (a1,0.75)*pow (a2,0.75)*sqrt (pow (a2,2))*(a1 + a2 + 2*pow (a1,2)*pow (Y1 - Y2,2)))/(pow (a1 + a2,3.5)*pow (E,(a1*a2*(pow (X1 - X2,2) + pow (Y1 - Y2,2) + pow (Z1 - Z2,2)))/(a1 + a2)));
    //0 0 0 - 1 0 0
    if(lx1==0&&ly1==0&&lz1==0&&lx2==1&&ly2==0&&lz2==0)res=(4*sqrt (2)*pow (a1,1.75)*pow (a2,1.25)*(X1 - X2))/(pow (a1 + a2,2.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //0 0 0 - 1 0 1
    if(lx1==0&&ly1==0&&lz1==0&&lx2==1&&ly2==0&&lz2==1)res=(8*sqrt (2)*pow (a1,2.75)*pow (a2,0.75)*sqrt (pow (a2,2))*(X1 - X2)*(Z1 - Z2))/(pow (a1 + a2,3.5)*pow (E,(a1*a2*(pow (X1 - X2,2) + pow (Y1 - Y2,2) + pow (Z1 - Z2,2)))/(a1 + a2)));
    //0 0 0 - 1 1 0
    if(lx1==0&&ly1==0&&lz1==0&&lx2==1&&ly2==1&&lz2==0)res=(8*sqrt (2)*pow (a1,2.75)*pow (a2,0.75)*sqrt (pow (a2,2))*(X1 - X2)*(Y1 - Y2))/(pow (a1 + a2,3.5)*pow (E,(a1*a2*(pow (X1 - X2,2) + pow (Y1 - Y2,2) + pow (Z1 - Z2,2)))/(a1 + a2)));
    //0 0 0 - 2 0 0
    if(lx1==0&&ly1==0&&lz1==0&&lx2==2&&ly2==0&&lz2==0)res=(4*sqrt (0.6666666666666666)*pow (a1,0.75)*pow (a2,0.75)*sqrt (pow (a2,2))*(a1 + a2 + 2*pow (a1,2)*pow (X1 - X2,2)))/(pow (a1 + a2,3.5)*pow (E,(a1*a2*(pow (X1 - X2,2) + pow (Y1 - Y2,2) + pow (Z1 - Z2,2)))/(a1 + a2)));
    //0 0 1 - 0 0 0
    if(lx1==0&&ly1==0&&lz1==1&&lx2==0&&ly2==0&&lz2==0)res=(4*sqrt (2)*pow (a1,1.25)*pow (a2,1.75)*(-Z1 + Z2))/(pow (a1 + a2,2.5)*pow (E,(a1*a2*(pow (X1 - X2,2) + pow (Y1 - Y2,2) + pow (Z1 - Z2,2)))/(a1 + a2)));
    //0 0 1 - 0 0 1
    if(lx1==0&&ly1==0&&lz1==1&&lx2==0&&ly2==0&&lz2==1)res=(4*sqrt (2)*pow (a1,1.25)*pow (a2,1.25)*(a2 + a1*(1 - 2*a2*pow (Z1 - Z2,2))))/(pow (a1 + a2,3.5)*pow (E,(a1*a2*(pow (X1 - X2,2) + pow (Y1 - Y2,2) + pow (Z1 - Z2,2)))/(a1 + a2)));
    //0 0 1 - 0 0 2
    if(lx1==0&&ly1==0&&lz1==1&&lx2==0&&ly2==0&&lz2==2)res=(-8*sqrt (0.6666666666666666)*pow (a1,1.25)*pow (a2,0.75)*sqrt (pow (a2,2))*(-(a1*a2) + pow (a2,2) + 2*pow (a1,2)*(-1 + a2*pow (Z1 - Z2,2)))*(Z1 - Z2))/(pow (a1 + a2,4.5)*pow (E,(a1*a2*(pow (X1 - X2,2) + pow (Y1 - Y2,2) + pow (Z1 - Z2,2)))/(a1 + a2)));
    //0 0 1 - 0 1 0
    if(lx1==0&&ly1==0&&lz1==1&&lx2==0&&ly2==1&&lz2==0)res=(-8*sqrt (2)*pow (a1,2.25)*pow (a2,2.25)*(Y1 - Y2)*(Z1 - Z2))/(pow (a1 + a2,3.5)*pow (E,(a1*a2*(pow (X1 - X2,2) + pow (Y1 - Y2,2) + pow (Z1 - Z2,2)))/(a1 + a2)));
    //0 0 1 - 0 1 1
    if(lx1==0&&ly1==0&&lz1==1&&lx2==0&&ly2==1&&lz2==1)res=(8*sqrt (2)*pow (a1,2.25)*pow (a2,0.75)*sqrt (pow (a2,2))*(Y1 - Y2)*(a2 + a1*(1 - 2*a2*pow (Z1 - Z2,2))))/(pow (a1 + a2,4.5)*pow (E,(a1*a2*(pow (X1 - X2,2) + pow (Y1 - Y2,2) + pow (Z1 - Z2,2)))/(a1 + a2)));
    //0 0 1 - 0 2 0
    if(lx1==0&&ly1==0&&lz1==1&&lx2==0&&ly2==2&&lz2==0)res=(-8*sqrt (0.6666666666666666)*pow (a1,1.25)*pow (a2,1.75)*sqrt (pow (a2,2))*(a1 + a2 + 2*pow (a1,2)*pow (Y1 - Y2,2))*(Z1 - Z2))/(pow (a1 + a2,4.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //0 0 1 - 1 0 0
    if(lx1==0&&ly1==0&&lz1==1&&lx2==1&&ly2==0&&lz2==0)res=(-8*sqrt (2)*pow (a1,2.25)*pow (a2,2.25)*(X1 - X2)*(Z1 - Z2))/(pow (a1 + a2,3.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //0 0 1 - 1 0 1
    if(lx1==0&&ly1==0&&lz1==1&&lx2==1&&ly2==0&&lz2==1)res=(8*sqrt (2)*pow (a1,2.25)*pow (a2,0.75)*sqrt (pow (a2,2))*(X1 - X2)*(a2 + a1*(1 - 2*a2*pow (Z1 - Z2,2))))/(pow (a1 + a2,4.5)*pow (E,(a1*a2*(pow (X1 - X2,2) + pow (Y1 - Y2,2) + pow (Z1 - Z2,2)))/(a1 + a2)));
    //0 0 1 - 1 1 0
    if(lx1==0&&ly1==0&&lz1==1&&lx2==1&&ly2==1&&lz2==0)res=(-16*sqrt (2)*pow (a1,3.25)*pow (a2,1.75)*sqrt (pow (a2,2))*(X1 - X2)*(Y1 - Y2)*(Z1 - Z2))/(pow (a1 + a2,4.5)*pow (E,(a1*a2*(pow (X1 - X2,2) + pow (Y1 - Y2,2) + pow (Z1 - Z2,2)))/(a1 + a2)));
    //0 0 1 - 2 0 0
    if(lx1==0&&ly1==0&&lz1==1&&lx2==2&&ly2==0&&lz2==0)res=(-8*sqrt (0.6666666666666666)*pow (a1,1.25)*pow (a2,1.75)*sqrt (pow (a2,2))*(a1 + a2 + 2*pow (a1,2)*pow (X1 - X2,2))*(Z1 - Z2))/(pow (a1 + a2,4.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //0 0 2 - 0 0 0
    if(lx1==0&&ly1==0&&lz1==2&&lx2==0&&ly2==0&&lz2==0)res=(4*sqrt (0.6666666666666666)*pow (a1,0.75)*sqrt (pow (a1,2))*pow (a2,0.75)*(a1 + a2*(1 + 2*a2*pow (Z1 - Z2,2))))/(pow (a1 + a2,3.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //0 0 2 - 0 0 1
    if(lx1==0&&ly1==0&&lz1==2&&lx2==0&&ly2==0&&lz2==1)res=(8*sqrt (0.6666666666666666)*pow (a1,0.75)*sqrt (pow (a1,2))*pow (a2,1.25)*(pow (a1,2) - 2*pow (a2,2) + a1*a2*(-1 + 2*a2*pow (Z1 - Z2,2)))*(Z1 - Z2))/(pow (a1 + a2,4.5)*pow (E,(a1*a2*(pow (X1 - X2,2) + pow (Y1 - Y2,2) + pow (Z1 - Z2,2)))/(a1 + a2)));
    //0 0 2 - 0 0 2
    if(lx1==0&&ly1==0&&lz1==2&&lx2==0&&ly2==0&&lz2==2)res=(8*sqrt (2)*pow (a1,0.75)*sqrt (pow (a1,2))*pow (a2,0.75)*sqrt (pow (a2,2))*(-6*a1*a2*(-1 + a2*pow (Z1 - Z2,2)) + pow (a2,2)*(3 + 2*a2*pow (Z1 - Z2,2)) + pow (a1,2)*(3 - 6*a2*pow (Z1 - Z2,2) + 4*pow (a2,2)*pow (Z1 - Z2,4)) + 2*pow (a1,3)*pow (Z1 - Z2,2)))/(3.*pow (a1 + a2,5.5)*pow (E,(a1*a2*(pow (X1 - X2,2) + pow (Y1 - Y2,2) + pow (Z1 - Z2,2)))/(a1 + a2)));
    //0 0 2 - 0 1 0
    if(lx1==0&&ly1==0&&lz1==2&&lx2==0&&ly2==1&&lz2==0)res=(8*sqrt (0.6666666666666666)*pow (a1,1.75)*sqrt (pow (a1,2))*pow (a2,1.25)*(Y1 - Y2)*(a1 + a2*(1 + 2*a2*pow (Z1 - Z2,2))))/(pow (a1 + a2,4.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //0 0 2 - 0 1 1
    if(lx1==0&&ly1==0&&lz1==2&&lx2==0&&ly2==1&&lz2==1)res=(16*sqrt (0.6666666666666666)*pow (a1,1.75)*sqrt (pow (a1,2))*pow (a2,0.75)*sqrt (pow (a2,2))*(Y1 - Y2)*(pow (a1,2) - 2*pow (a2,2) + a1*a2*(-1 + 2*a2*pow (Z1 - Z2,2)))*(Z1 - Z2))/(pow (a1 + a2,5.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //0 0 2 - 0 2 0
    if(lx1==0&&ly1==0&&lz1==2&&lx2==0&&ly2==2&&lz2==0)res=(8*sqrt (2)*pow (a1,0.75)*sqrt (pow (a1,2))*pow (a2,0.75)*sqrt (pow (a2,2))*(a1 + a2 + 2*pow (a1,2)*pow (Y1 - Y2,2))*(a1 + a2*(1 + 2*a2*pow (Z1 - Z2,2))))/(3.*pow (a1 + a2,5.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //0 0 2 - 1 0 0
    if(lx1==0&&ly1==0&&lz1==2&&lx2==1&&ly2==0&&lz2==0)res=(8*sqrt (0.6666666666666666)*pow (a1,1.75)*sqrt (pow (a1,2))*pow (a2,1.25)*(X1 - X2)*(a1 + a2*(1 + 2*a2*pow (Z1 - Z2,2))))/(pow (a1 + a2,4.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //0 0 2 - 1 0 1
    if(lx1==0&&ly1==0&&lz1==2&&lx2==1&&ly2==0&&lz2==1)res=(16*sqrt (0.6666666666666666)*pow (a1,1.75)*sqrt (pow (a1,2))*pow (a2,0.75)*sqrt (pow (a2,2))*(X1 - X2)*(pow (a1,2) - 2*pow (a2,2) + a1*a2*(-1 + 2*a2*pow (Z1 - Z2,2)))*(Z1 - Z2))/(pow (a1 + a2,5.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //0 0 2 - 1 1 0
    if(lx1==0&&ly1==0&&lz1==2&&lx2==1&&ly2==1&&lz2==0)res=(16*sqrt (0.6666666666666666)*pow (a1,2.75)*sqrt (pow (a1,2))*pow (a2,0.75)*sqrt (pow (a2,2))*(X1 - X2)*(Y1 - Y2)*(a1 + a2*(1 + 2*a2*pow (Z1 - Z2,2))))/(pow (a1 + a2,5.5)*pow (E,(a1*a2*(pow (X1 - X2,2) + pow (Y1 - Y2,2) + pow (Z1 - Z2,2)))/(a1 + a2)));
    //0 0 2 - 2 0 0
    if(lx1==0&&ly1==0&&lz1==2&&lx2==2&&ly2==0&&lz2==0)res=(8*sqrt (2)*pow (a1,0.75)*sqrt (pow (a1,2))*pow (a2,0.75)*sqrt (pow (a2,2))*(a1 + a2 + 2*pow (a1,2)*pow (X1 - X2,2))*(a1 + a2*(1 + 2*a2*pow (Z1 - Z2,2))))/(3.*pow (a1 + a2,5.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //0 1 0 - 0 0 0
    if(lx1==0&&ly1==1&&lz1==0&&lx2==0&&ly2==0&&lz2==0)res=(4*sqrt (2)*pow (a1,1.25)*pow (a2,1.75)*(-Y1 + Y2))/(pow (a1 + a2,2.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //0 1 0 - 0 0 1
    if(lx1==0&&ly1==1&&lz1==0&&lx2==0&&ly2==0&&lz2==1)res=(-8*sqrt (2)*pow (a1,2.25)*pow (a2,2.25)*(Y1 - Y2)*(Z1 - Z2))/(pow (a1 + a2,3.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //0 1 0 - 0 0 2
    if(lx1==0&&ly1==1&&lz1==0&&lx2==0&&ly2==0&&lz2==2)res=(8*sqrt (0.6666666666666666)*pow (a1,1.25)*pow (a2,1.75)*sqrt (pow (a2,2))*(-Y1 + Y2)*(a1 + a2 + 2*pow (a1,2)*pow (Z1 - Z2,2)))/(pow (a1 + a2,4.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //0 1 0 - 0 1 0
    if(lx1==0&&ly1==1&&lz1==0&&lx2==0&&ly2==1&&lz2==0)res=(4*sqrt (2)*pow (a1,1.25)*pow (a2,1.25)*(a2 + a1*(1 - 2*a2*pow (Y1 - Y2,2))))/(pow (a1 + a2,3.5)*pow (E,(a1*a2*(pow (X1 - X2,2) + pow (Y1 - Y2,2) + pow (Z1 - Z2,2)))/(a1 + a2)));
    //0 1 0 - 0 1 1
    if(lx1==0&&ly1==1&&lz1==0&&lx2==0&&ly2==1&&lz2==1)res=(8*sqrt (2)*pow (a1,2.25)*pow (a2,0.75)*sqrt (pow (a2,2))*(a2 + a1*(1 - 2*a2*pow (Y1 - Y2,2)))*(Z1 - Z2))/(pow (a1 + a2,4.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //0 1 0 - 0 2 0
    if(lx1==0&&ly1==1&&lz1==0&&lx2==0&&ly2==2&&lz2==0)res=(-8*sqrt (0.6666666666666666)*pow (a1,1.25)*pow (a2,0.75)*sqrt (pow (a2,2))*(-(a1*a2) + pow (a2,2) + 2*pow (a1,2)*(-1 + a2*pow (Y1 - Y2,2)))*(Y1 - Y2))/(pow (a1 + a2,4.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //0 1 0 - 1 0 0
    if(lx1==0&&ly1==1&&lz1==0&&lx2==1&&ly2==0&&lz2==0)res=(-8*sqrt (2)*pow (a1,2.25)*pow (a2,2.25)*(X1 - X2)*(Y1 - Y2))/(pow (a1 + a2,3.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //0 1 0 - 1 0 1
    if(lx1==0&&ly1==1&&lz1==0&&lx2==1&&ly2==0&&lz2==1)res=(16*sqrt (2)*pow (a1,3.25)*pow (a2,1.75)*sqrt (pow (a2,2))*(X1 - X2)*(-Y1 + Y2)*(Z1 - Z2))/(pow (a1 + a2,4.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //0 1 0 - 1 1 0
    if(lx1==0&&ly1==1&&lz1==0&&lx2==1&&ly2==1&&lz2==0)res=(8*sqrt (2)*pow (a1,2.25)*pow (a2,0.75)*sqrt (pow (a2,2))*(X1 - X2)*(a2 + a1*(1 - 2*a2*pow (Y1 - Y2,2))))/(pow (a1 + a2,4.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //0 1 0 - 2 0 0
    if(lx1==0&&ly1==1&&lz1==0&&lx2==2&&ly2==0&&lz2==0)res=(8*sqrt (0.6666666666666666)*pow (a1,1.25)*pow (a2,1.75)*sqrt (pow (a2,2))*(a1 + a2 + 2*pow (a1,2)*pow (X1 - X2,2))*(-Y1 + Y2))/(pow (a1 + a2,4.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //0 1 1 - 0 0 0
    if(lx1==0&&ly1==1&&lz1==1&&lx2==0&&ly2==0&&lz2==0)res=(8*sqrt (2)*pow (a1,0.75)*sqrt (pow (a1,2))*pow (a2,2.75)*(Y1 - Y2)*(Z1 - Z2))/(pow (a1 + a2,3.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //0 1 1 - 0 0 1
    if(lx1==0&&ly1==1&&lz1==1&&lx2==0&&ly2==0&&lz2==1)res=(8*sqrt (2)*pow (a1,0.75)*sqrt (pow (a1,2))*pow (a2,2.25)*(-Y1 + Y2)*(a2 + a1*(1 - 2*a2*pow (Z1 - Z2,2))))/(pow (a1 + a2,4.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //0 1 1 - 0 0 2
    if(lx1==0&&ly1==1&&lz1==1&&lx2==0&&ly2==0&&lz2==2)res=(16*sqrt (0.6666666666666666)*pow (a1,0.75)*sqrt (pow (a1,2))*pow (a2,1.75)*sqrt (pow (a2,2))*(Y1 - Y2)*(-(a1*a2) + pow (a2,2) + 2*pow (a1,2)*(-1 + a2*pow (Z1 - Z2,2)))*(Z1 - Z2))/(pow (a1 + a2,5.5)*pow (E,(a1*a2*(pow (X1 - X2,2) + pow (Y1 - Y2,2) + pow (Z1 - Z2,2)))/(a1 + a2)));
    //0 1 1 - 0 1 0
    if(lx1==0&&ly1==1&&lz1==1&&lx2==0&&ly2==1&&lz2==0)res=(-8*sqrt (2)*pow (a1,0.75)*sqrt (pow (a1,2))*pow (a2,2.25)*(a2 + a1*(1 - 2*a2*pow (Y1 - Y2,2)))*(Z1 - Z2))/(pow (a1 + a2,4.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //0 1 1 - 0 1 1
    if(lx1==0&&ly1==1&&lz1==1&&lx2==0&&ly2==1&&lz2==1)res=(8*sqrt (2)*pow (a1,0.75)*sqrt (pow (a1,2))*pow (a2,0.75)*sqrt (pow (a2,2))*(a2 + a1*(1 - 2*a2*pow (Y1 - Y2,2)))*(a2 + a1*(1 - 2*a2*pow (Z1 - Z2,2))))/(pow (a1 + a2,5.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //0 1 1 - 0 2 0
    if(lx1==0&&ly1==1&&lz1==1&&lx2==0&&ly2==2&&lz2==0)res=(16*sqrt (0.6666666666666666)*pow (a1,0.75)*sqrt (pow (a1,2))*pow (a2,1.75)*sqrt (pow (a2,2))*(-(a1*a2) + pow (a2,2) + 2*pow (a1,2)*(-1 + a2*pow (Y1 - Y2,2)))*(Y1 - Y2)*(Z1 - Z2))/(pow (a1 + a2,5.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //0 1 1 - 1 0 0
    if(lx1==0&&ly1==1&&lz1==1&&lx2==1&&ly2==0&&lz2==0)res=(16*sqrt (2)*pow (a1,1.75)*sqrt (pow (a1,2))*pow (a2,3.25)*(X1 - X2)*(Y1 - Y2)*(Z1 - Z2))/(pow (a1 + a2,4.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //0 1 1 - 1 0 1
    if(lx1==0&&ly1==1&&lz1==1&&lx2==1&&ly2==0&&lz2==1)res=(16*sqrt (2)*pow (a1,1.75)*sqrt (pow (a1,2))*pow (a2,1.75)*sqrt (pow (a2,2))*(X1 - X2)*(-Y1 + Y2)*(a2 + a1*(1 - 2*a2*pow (Z1 - Z2,2))))/(pow (a1 + a2,5.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //0 1 1 - 1 1 0
    if(lx1==0&&ly1==1&&lz1==1&&lx2==1&&ly2==1&&lz2==0)res=(-16*sqrt (2)*pow (a1,1.75)*sqrt (pow (a1,2))*pow (a2,1.75)*sqrt (pow (a2,2))*(X1 - X2)*(a2 + a1*(1 - 2*a2*pow (Y1 - Y2,2)))*(Z1 - Z2))/(pow (a1 + a2,5.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //0 1 1 - 2 0 0
    if(lx1==0&&ly1==1&&lz1==1&&lx2==2&&ly2==0&&lz2==0)res=(16*sqrt (0.6666666666666666)*pow (a1,0.75)*sqrt (pow (a1,2))*pow (a2,2.75)*sqrt (pow (a2,2))*(a1 + a2 + 2*pow (a1,2)*pow (X1 - X2,2))*(Y1 - Y2)*(Z1 - Z2))/(pow (a1 + a2,5.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //0 2 0 - 0 0 0
    if(lx1==0&&ly1==2&&lz1==0&&lx2==0&&ly2==0&&lz2==0)res=(4*sqrt (0.6666666666666666)*pow (a1,0.75)*sqrt (pow (a1,2))*pow (a2,0.75)*(a1 + a2*(1 + 2*a2*pow (Y1 - Y2,2))))/(pow (a1 + a2,3.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //0 2 0 - 0 0 1
    if(lx1==0&&ly1==2&&lz1==0&&lx2==0&&ly2==0&&lz2==1)res=(8*sqrt (0.6666666666666666)*pow (a1,1.75)*sqrt (pow (a1,2))*pow (a2,1.25)*(a1 + a2*(1 + 2*a2*pow (Y1 - Y2,2)))*(Z1 - Z2))/(pow (a1 + a2,4.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //0 2 0 - 0 0 2
    if(lx1==0&&ly1==2&&lz1==0&&lx2==0&&ly2==0&&lz2==2)res=(8*sqrt (2)*pow (a1,0.75)*sqrt (pow (a1,2))*pow (a2,0.75)*sqrt (pow (a2,2))*(a1 + a2*(1 + 2*a2*pow (Y1 - Y2,2)))*(a1 + a2 + 2*pow (a1,2)*pow (Z1 - Z2,2)))/(3.*pow (a1 + a2,5.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //0 2 0 - 0 1 0
    if(lx1==0&&ly1==2&&lz1==0&&lx2==0&&ly2==1&&lz2==0)res=(8*sqrt (0.6666666666666666)*pow (a1,0.75)*sqrt (pow (a1,2))*pow (a2,1.25)*(pow (a1,2) - 2*pow (a2,2) + a1*a2*(-1 + 2*a2*pow (Y1 - Y2,2)))*(Y1 - Y2))/(pow (a1 + a2,4.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //0 2 0 - 0 1 1
    if(lx1==0&&ly1==2&&lz1==0&&lx2==0&&ly2==1&&lz2==1)res=(16*sqrt (0.6666666666666666)*pow (a1,1.75)*sqrt (pow (a1,2))*pow (a2,0.75)*sqrt (pow (a2,2))*(pow (a1,2) - 2*pow (a2,2) + a1*a2*(-1 + 2*a2*pow (Y1 - Y2,2)))*(Y1 - Y2)*(Z1 - Z2))/(pow (a1 + a2,5.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //0 2 0 - 0 2 0
    if(lx1==0&&ly1==2&&lz1==0&&lx2==0&&ly2==2&&lz2==0)res=(8*sqrt (2)*pow (a1,0.75)*sqrt (pow (a1,2))*pow (a2,0.75)*sqrt (pow (a2,2))*(-6*a1*a2*(-1 + a2*pow (Y1 - Y2,2)) + pow (a2,2)*(3 + 2*a2*pow (Y1 - Y2,2)) + pow (a1,2)*(3 - 6*a2*pow (Y1 - Y2,2) + 4*pow (a2,2)*pow (Y1 - Y2,4)) + 2*pow (a1,3)*pow (Y1 - Y2,2)))/(3.*pow (a1 + a2,5.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //0 2 0 - 1 0 0
    if(lx1==0&&ly1==2&&lz1==0&&lx2==1&&ly2==0&&lz2==0)res=(8*sqrt (0.6666666666666666)*pow (a1,1.75)*sqrt (pow (a1,2))*pow (a2,1.25)*(X1 - X2)*(a1 + a2*(1 + 2*a2*pow (Y1 - Y2,2))))/(pow (a1 + a2,4.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //0 2 0 - 1 0 1
    if(lx1==0&&ly1==2&&lz1==0&&lx2==1&&ly2==0&&lz2==1)res=(16*sqrt (0.6666666666666666)*pow (a1,2.75)*sqrt (pow (a1,2))*pow (a2,0.75)*sqrt (pow (a2,2))*(X1 - X2)*(a1 + a2*(1 + 2*a2*pow (Y1 - Y2,2)))*(Z1 - Z2))/(pow (a1 + a2,5.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //0 2 0 - 1 1 0
    if(lx1==0&&ly1==2&&lz1==0&&lx2==1&&ly2==1&&lz2==0)res=(16*sqrt (0.6666666666666666)*pow (a1,1.75)*sqrt (pow (a1,2))*pow (a2,0.75)*sqrt (pow (a2,2))*(X1 - X2)*(pow (a1,2) - 2*pow (a2,2) + a1*a2*(-1 + 2*a2*pow (Y1 - Y2,2)))*(Y1 - Y2))/(pow (a1 + a2,5.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //0 2 0 - 2 0 0
    if(lx1==0&&ly1==2&&lz1==0&&lx2==2&&ly2==0&&lz2==0)res=(8*sqrt (2)*pow (a1,0.75)*sqrt (pow (a1,2))*pow (a2,0.75)*sqrt (pow (a2,2))*(a1 + a2 + 2*pow (a1,2)*pow (X1 - X2,2))*(a1 + a2*(1 + 2*a2*pow (Y1 - Y2,2))))/(3.*pow (a1 + a2,5.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //1 0 0 - 0 0 0
    if(lx1==1&&ly1==0&&lz1==0&&lx2==0&&ly2==0&&lz2==0)res=(-4*sqrt (2)*pow (a1,1.25)*pow (a2,1.75)*(X1 - X2))/(pow (a1 + a2,2.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //1 0 0 - 0 0 1
    if(lx1==1&&ly1==0&&lz1==0&&lx2==0&&ly2==0&&lz2==1)res=(-8*sqrt (2)*pow (a1,2.25)*pow (a2,2.25)*(X1 - X2)*(Z1 - Z2))/(pow (a1 + a2,3.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //1 0 0 - 0 0 2
    if(lx1==1&&ly1==0&&lz1==0&&lx2==0&&ly2==0&&lz2==2)res=(8*sqrt (0.6666666666666666)*pow (a1,1.25)*pow (a2,1.75)*sqrt (pow (a2,2))*(-X1 + X2)*(a1 + a2 + 2*pow (a1,2)*pow (Z1 - Z2,2)))/(pow (a1 + a2,4.5)*pow (E,(a1*a2*(pow (X1 - X2,2) + pow (Y1 - Y2,2) + pow (Z1 - Z2,2)))/(a1 + a2)));
    //1 0 0 - 0 1 0
    if(lx1==1&&ly1==0&&lz1==0&&lx2==0&&ly2==1&&lz2==0)res=(8*sqrt (2)*pow (a1,2.25)*pow (a2,2.25)*(-X1 + X2)*(Y1 - Y2))/(pow (a1 + a2,3.5)*pow (E,(a1*a2*(pow (X1 - X2,2) + pow (Y1 - Y2,2) + pow (Z1 - Z2,2)))/(a1 + a2)));
    //1 0 0 - 0 1 1
    if(lx1==1&&ly1==0&&lz1==0&&lx2==0&&ly2==1&&lz2==1)res=(-16*sqrt (2)*pow (a1,3.25)*pow (a2,1.75)*sqrt (pow (a2,2))*(X1 - X2)*(Y1 - Y2)*(Z1 - Z2))/(pow (a1 + a2,4.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //1 0 0 - 0 2 0
    if(lx1==1&&ly1==0&&lz1==0&&lx2==0&&ly2==2&&lz2==0)res=(-8*sqrt (0.6666666666666666)*pow (a1,1.25)*pow (a2,1.75)*sqrt (pow (a2,2))*(X1 - X2)*(a1 + a2 + 2*pow (a1,2)*pow (Y1 - Y2,2)))/(pow (a1 + a2,4.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //1 0 0 - 1 0 0
    if(lx1==1&&ly1==0&&lz1==0&&lx2==1&&ly2==0&&lz2==0)res=(4*sqrt (2)*pow (a1,1.25)*pow (a2,1.25)*(a2 + a1*(1 - 2*a2*pow (X1 - X2,2))))/(pow (a1 + a2,3.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //1 0 0 - 1 0 1
    if(lx1==1&&ly1==0&&lz1==0&&lx2==1&&ly2==0&&lz2==1)res=(8*sqrt (2)*pow (a1,2.25)*pow (a2,0.75)*sqrt (pow (a2,2))*(a2 + a1*(1 - 2*a2*pow (X1 - X2,2)))*(Z1 - Z2))/(pow (a1 + a2,4.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //1 0 0 - 1 1 0
    if(lx1==1&&ly1==0&&lz1==0&&lx2==1&&ly2==1&&lz2==0)res=(-8*sqrt (2)*pow (a1,2.25)*pow (a2,0.75)*sqrt (pow (a2,2))*(-a2 + a1*(-1 + 2*a2*pow (X1 - X2,2)))*(Y1 - Y2))/(pow (a1 + a2,4.5)*pow (E,(a1*a2*(pow (X1 - X2,2) + pow (Y1 - Y2,2) + pow (Z1 - Z2,2)))/(a1 + a2)));
    //1 0 0 - 2 0 0
    if(lx1==1&&ly1==0&&lz1==0&&lx2==2&&ly2==0&&lz2==0)res=(-8*sqrt (0.6666666666666666)*pow (a1,1.25)*pow (a2,0.75)*sqrt (pow (a2,2))*(-(a1*a2) + pow (a2,2) + 2*pow (a1,2)*(-1 + a2*pow (X1 - X2,2)))*(X1 - X2))/(pow (a1 + a2,4.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //1 0 1 - 0 0 0
    if(lx1==1&&ly1==0&&lz1==1&&lx2==0&&ly2==0&&lz2==0)res=(8*sqrt (2)*pow (a1,0.75)*sqrt (pow (a1,2))*pow (a2,2.75)*(X1 - X2)*(Z1 - Z2))/(pow (a1 + a2,3.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //1 0 1 - 0 0 1
    if(lx1==1&&ly1==0&&lz1==1&&lx2==0&&ly2==0&&lz2==1)res=(-8*sqrt (2)*pow (a1,0.75)*sqrt (pow (a1,2))*pow (a2,2.25)*(X1 - X2)*(a2 + a1*(1 - 2*a2*pow (Z1 - Z2,2))))/(pow (a1 + a2,4.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //1 0 1 - 0 0 2
    if(lx1==1&&ly1==0&&lz1==1&&lx2==0&&ly2==0&&lz2==2)res=(16*sqrt (0.6666666666666666)*pow (a1,0.75)*sqrt (pow (a1,2))*pow (a2,1.75)*sqrt (pow (a2,2))*(X1 - X2)*(-(a1*a2) + pow (a2,2) + 2*pow (a1,2)*(-1 + a2*pow (Z1 - Z2,2)))*(Z1 - Z2))/(pow (a1 + a2,5.5)*pow (E,(a1*a2*(pow (X1 - X2,2) + pow (Y1 - Y2,2) + pow (Z1 - Z2,2)))/(a1 + a2)));
    //1 0 1 - 0 1 0
    if(lx1==1&&ly1==0&&lz1==1&&lx2==0&&ly2==1&&lz2==0)res=(16*sqrt (2)*pow (a1,1.75)*sqrt (pow (a1,2))*pow (a2,3.25)*(X1 - X2)*(Y1 - Y2)*(Z1 - Z2))/(pow (a1 + a2,4.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //1 0 1 - 0 1 1
    if(lx1==1&&ly1==0&&lz1==1&&lx2==0&&ly2==1&&lz2==1)res=(-16*sqrt (2)*pow (a1,1.75)*sqrt (pow (a1,2))*pow (a2,1.75)*sqrt (pow (a2,2))*(X1 - X2)*(Y1 - Y2)*(a2 + a1*(1 - 2*a2*pow (Z1 - Z2,2))))/(pow (a1 + a2,5.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //1 0 1 - 0 2 0
    if(lx1==1&&ly1==0&&lz1==1&&lx2==0&&ly2==2&&lz2==0)res=(16*sqrt (0.6666666666666666)*pow (a1,0.75)*sqrt (pow (a1,2))*pow (a2,2.75)*sqrt (pow (a2,2))*(X1 - X2)*(a1 + a2 + 2*pow (a1,2)*pow (Y1 - Y2,2))*(Z1 - Z2))/(pow (a1 + a2,5.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //1 0 1 - 1 0 0
    if(lx1==1&&ly1==0&&lz1==1&&lx2==1&&ly2==0&&lz2==0)res=(-8*sqrt (2)*pow (a1,0.75)*sqrt (pow (a1,2))*pow (a2,2.25)*(a2 + a1*(1 - 2*a2*pow (X1 - X2,2)))*(Z1 - Z2))/(pow (a1 + a2,4.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //1 0 1 - 1 0 1
    if(lx1==1&&ly1==0&&lz1==1&&lx2==1&&ly2==0&&lz2==1)res=(8*sqrt (2)*pow (a1,0.75)*sqrt (pow (a1,2))*pow (a2,0.75)*sqrt (pow (a2,2))*(a2 + a1*(1 - 2*a2*pow (X1 - X2,2)))*(a2 + a1*(1 - 2*a2*pow (Z1 - Z2,2))))/(pow (a1 + a2,5.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //1 0 1 - 1 1 0
    if(lx1==1&&ly1==0&&lz1==1&&lx2==1&&ly2==1&&lz2==0)res=(-16*sqrt (2)*pow (a1,1.75)*sqrt (pow (a1,2))*pow (a2,1.75)*sqrt (pow (a2,2))*(a2 + a1*(1 - 2*a2*pow (X1 - X2,2)))*(Y1 - Y2)*(Z1 - Z2))/(pow (a1 + a2,5.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //1 0 1 - 2 0 0
    if(lx1==1&&ly1==0&&lz1==1&&lx2==2&&ly2==0&&lz2==0)res=(16*sqrt (0.6666666666666666)*pow (a1,0.75)*sqrt (pow (a1,2))*pow (a2,1.75)*sqrt (pow (a2,2))*(-(a1*a2) + pow (a2,2) + 2*pow (a1,2)*(-1 + a2*pow (X1 - X2,2)))*(X1 - X2)*(Z1 - Z2))/(pow (a1 + a2,5.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //1 1 0 - 0 0 0
    if(lx1==1&&ly1==1&&lz1==0&&lx2==0&&ly2==0&&lz2==0)res=(8*sqrt (2)*pow (a1,0.75)*sqrt (pow (a1,2))*pow (a2,2.75)*(X1 - X2)*(Y1 - Y2))/(pow (a1 + a2,3.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //1 1 0 - 0 0 1
    if(lx1==1&&ly1==1&&lz1==0&&lx2==0&&ly2==0&&lz2==1)res=(16*sqrt (2)*pow (a1,1.75)*sqrt (pow (a1,2))*pow (a2,3.25)*(X1 - X2)*(Y1 - Y2)*(Z1 - Z2))/(pow (a1 + a2,4.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //1 1 0 - 0 0 2
    if(lx1==1&&ly1==1&&lz1==0&&lx2==0&&ly2==0&&lz2==2)res=(16*sqrt (0.6666666666666666)*pow (a1,0.75)*sqrt (pow (a1,2))*pow (a2,2.75)*sqrt (pow (a2,2))*(X1 - X2)*(Y1 - Y2)*(a1 + a2 + 2*pow (a1,2)*pow (Z1 - Z2,2)))/(pow (a1 + a2,5.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //1 1 0 - 0 1 0
    if(lx1==1&&ly1==1&&lz1==0&&lx2==0&&ly2==1&&lz2==0)res=(-8*sqrt (2)*pow (a1,0.75)*sqrt (pow (a1,2))*pow (a2,2.25)*(X1 - X2)*(a2 + a1*(1 - 2*a2*pow (Y1 - Y2,2))))/(pow (a1 + a2,4.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //1 1 0 - 0 1 1
    if(lx1==1&&ly1==1&&lz1==0&&lx2==0&&ly2==1&&lz2==1)res=(-16*sqrt (2)*pow (a1,1.75)*sqrt (pow (a1,2))*pow (a2,1.75)*sqrt (pow (a2,2))*(X1 - X2)*(a2 + a1*(1 - 2*a2*pow (Y1 - Y2,2)))*(Z1 - Z2))/(pow (a1 + a2,5.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //1 1 0 - 0 2 0
    if(lx1==1&&ly1==1&&lz1==0&&lx2==0&&ly2==2&&lz2==0)res=(16*sqrt (0.6666666666666666)*pow (a1,0.75)*sqrt (pow (a1,2))*pow (a2,1.75)*sqrt (pow (a2,2))*(X1 - X2)*(-(a1*a2) + pow (a2,2) + 2*pow (a1,2)*(-1 + a2*pow (Y1 - Y2,2)))*(Y1 - Y2))/(pow (a1 + a2,5.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //1 1 0 - 1 0 0
    if(lx1==1&&ly1==1&&lz1==0&&lx2==1&&ly2==0&&lz2==0)res=(8*sqrt (2)*pow (a1,0.75)*sqrt (pow (a1,2))*pow (a2,2.25)*(a2 + a1*(1 - 2*a2*pow (X1 - X2,2)))*(-Y1 + Y2))/(pow (a1 + a2,4.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //1 1 0 - 1 0 1
    if(lx1==1&&ly1==1&&lz1==0&&lx2==1&&ly2==0&&lz2==1)res=(16*sqrt (2)*pow (a1,1.75)*sqrt (pow (a1,2))*pow (a2,1.75)*sqrt (pow (a2,2))*(a2 + a1*(1 - 2*a2*pow (X1 - X2,2)))*(-Y1 + Y2)*(Z1 - Z2))/(pow (a1 + a2,5.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //1 1 0 - 1 1 0
    if(lx1==1&&ly1==1&&lz1==0&&lx2==1&&ly2==1&&lz2==0)res=(8*sqrt (2)*pow (a1,0.75)*sqrt (pow (a1,2))*pow (a2,0.75)*sqrt (pow (a2,2))*(a2 + a1*(1 - 2*a2*pow (X1 - X2,2)))*(a2 + a1*(1 - 2*a2*pow (Y1 - Y2,2))))/(pow (a1 + a2,5.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //1 1 0 - 2 0 0
    if(lx1==1&&ly1==1&&lz1==0&&lx2==2&&ly2==0&&lz2==0)res=(16*sqrt (0.6666666666666666)*pow (a1,0.75)*sqrt (pow (a1,2))*pow (a2,1.75)*sqrt (pow (a2,2))*(-(a1*a2) + pow (a2,2) + 2*pow (a1,2)*(-1 + a2*pow (X1 - X2,2)))*(X1 - X2)*(Y1 - Y2))/(pow (a1 + a2,5.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //2 0 0 - 0 0 0
    if(lx1==2&&ly1==0&&lz1==0&&lx2==0&&ly2==0&&lz2==0)res=(4*sqrt (0.6666666666666666)*pow (a1,0.75)*sqrt (pow (a1,2))*pow (a2,0.75)*(a1 + a2*(1 + 2*a2*pow (X1 - X2,2))))/(pow (a1 + a2,3.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //2 0 0 - 0 0 1
    if(lx1==2&&ly1==0&&lz1==0&&lx2==0&&ly2==0&&lz2==1)res=(8*sqrt (0.6666666666666666)*pow (a1,1.75)*sqrt (pow (a1,2))*pow (a2,1.25)*(a1 + a2*(1 + 2*a2*pow (X1 - X2,2)))*(Z1 - Z2))/(pow (a1 + a2,4.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //2 0 0 - 0 0 2
    if(lx1==2&&ly1==0&&lz1==0&&lx2==0&&ly2==0&&lz2==2)res=(8*sqrt (2)*pow (a1,0.75)*sqrt (pow (a1,2))*pow (a2,0.75)*sqrt (pow (a2,2))*(a1 + a2*(1 + 2*a2*pow (X1 - X2,2)))*(a1 + a2 + 2*pow (a1,2)*pow (Z1 - Z2,2)))/(3.*pow (a1 + a2,5.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //2 0 0 - 0 1 0
    if(lx1==2&&ly1==0&&lz1==0&&lx2==0&&ly2==1&&lz2==0)res=(8*sqrt (0.6666666666666666)*pow (a1,1.75)*sqrt (pow (a1,2))*pow (a2,1.25)*(a1 + a2*(1 + 2*a2*pow (X1 - X2,2)))*(Y1 - Y2))/(pow (a1 + a2,4.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //2 0 0 - 0 1 1
    if(lx1==2&&ly1==0&&lz1==0&&lx2==0&&ly2==1&&lz2==1)res=(16*sqrt (0.6666666666666666)*pow (a1,2.75)*sqrt (pow (a1,2))*pow (a2,0.75)*sqrt (pow (a2,2))*(a1 + a2*(1 + 2*a2*pow (X1 - X2,2)))*(Y1 - Y2)*(Z1 - Z2))/(pow (a1 + a2,5.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //2 0 0 - 0 2 0
    if(lx1==2&&ly1==0&&lz1==0&&lx2==0&&ly2==2&&lz2==0)res=(8*sqrt (2)*pow (a1,0.75)*sqrt (pow (a1,2))*pow (a2,0.75)*sqrt (pow (a2,2))*(a1 + a2*(1 + 2*a2*pow (X1 - X2,2)))*(a1 + a2 + 2*pow (a1,2)*pow (Y1 - Y2,2)))/(3.*pow (a1 + a2,5.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //2 0 0 - 1 0 0
    if(lx1==2&&ly1==0&&lz1==0&&lx2==1&&ly2==0&&lz2==0)res=(8*sqrt (0.6666666666666666)*pow (a1,0.75)*sqrt (pow (a1,2))*pow (a2,1.25)*(pow (a1,2) - 2*pow (a2,2) + a1*a2*(-1 + 2*a2*pow (X1 - X2,2)))*(X1 - X2))/(pow (a1 + a2,4.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //2 0 0 - 1 0 1
    if(lx1==2&&ly1==0&&lz1==0&&lx2==1&&ly2==0&&lz2==1)res=(16*sqrt (0.6666666666666666)*pow (a1,1.75)*sqrt (pow (a1,2))*pow (a2,0.75)*sqrt (pow (a2,2))*(pow (a1,2) - 2*pow (a2,2) + a1*a2*(-1 + 2*a2*pow (X1 - X2,2)))*(X1 - X2)*(Z1 - Z2))/(pow (a1 + a2,5.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //2 0 0 - 1 1 0
    if(lx1==2&&ly1==0&&lz1==0&&lx2==1&&ly2==1&&lz2==0)res=(16*sqrt (0.6666666666666666)*pow (a1,1.75)*sqrt (pow (a1,2))*pow (a2,0.75)*sqrt (pow (a2,2))*(pow (a1,2) - 2*pow (a2,2) + a1*a2*(-1 + 2*a2*pow (X1 - X2,2)))*(X1 - X2)*(Y1 - Y2))/(pow (a1 + a2,5.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    //2 0 0 - 2 0 0
    if(lx1==2&&ly1==0&&lz1==0&&lx2==2&&ly2==0&&lz2==0)res=(8*sqrt (2)*pow (a1,0.75)*sqrt (pow (a1,2))*pow (a2,0.75)*sqrt (pow (a2,2))*(-6*a1*a2*(-1 + a2*pow (X1 - X2,2)) + pow (a2,2)*(3 + 2*a2*pow (X1 - X2,2)) + pow (a1,2)*(3 - 6*a2*pow (X1 - X2,2) + 4*pow (a2,2)*pow (X1 - X2,4)) + 2*pow (a1,3)*pow (X1 - X2,2)))/(3.*pow (a1 + a2,5.5)*pow (E,(a1*a2*(pow (X1,2) - 2*X1*X2 + pow (X2,2) + pow (Y1,2) - 2*Y1*Y2 + pow (Y2,2) + pow (Z1,2) - 2*Z1*Z2 + pow (Z2,2)))/(a1 + a2)));
    return res;
}
