/* -----------------------------------------------------------------------------
 
 Unified STO projection / AOM overlap code
 
 Beta version: 1.0
 1-Mar-2019
 
 Orestis George Ziogos, UCL
 o.ziogos@ucl.ac.uk
 
 For more information, examine the README file in the parent directory.
 
-------------------------------------------------------------------------------- */
#include"general.h"
long int factorial(int n);
void MakeLocalCoordsSystems(double X1,double Y1,double Z1,double X2,double Y2,double Z2,double *exb_X,double *exb_Y,double *exb_Z,double *eyb_X,double *eyb_Y,double *eyb_Z,double *ezb_X,double *ezb_Y,double *ezb_Z);
int equal(double A,double B);
double Afunction(int k,double p);
double Afunction(int k,double p)
{
    double res;
    int mu,n;
    double prefac,sum,nom,denom;
    prefac=exp(-p);
    sum=0.0;
    for(mu=1;mu<=k+1;++mu)
    {
        nom=(double)factorial(k);
        denom=pow(p,mu);
        n=k-mu+1;
        denom=denom*factorial(n);
        sum=sum+nom/denom;
    }
    res=prefac*sum;
    return res;
}
double Bfunction(int k,double p,double t);
double Bfunction(int k,double p,double t)
{
    double res;
    int mu,n;
    double prefac1,prefac2,sum1,sum2,nom,denom,u;
    prefac1=-exp(-p*t);
    prefac2=-exp(p*t);
    sum1=0.0;sum2=0.0;
    for(mu=1;mu<=k+1;++mu)
    {
        nom=(double)factorial(k);
        u=p*t;
        denom=pow(u,mu);
        n=k-mu+1;
        denom=denom*factorial(n);
        u=nom/denom;
        sum1=sum1+u;
        n=k-mu;
        sum2=sum2+u*pow(-1,n);
    }
    res=prefac1*sum1+prefac2*sum2;
    return res;
}
double pvalue(double X1,double Y1,double Z1,double X2,double Y2,double Z2,double mu1,double mu2);
double pvalue(double X1,double Y1,double Z1,double X2,double Y2,double Z2,double mu1,double mu2)
{
    double res,norm;
    norm=(X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2)+(Z1-Z2)*(Z1-Z2);norm=sqrt(norm);
    res=0.5*norm*(mu1+mu2);
    return res;
}
double tvalue(double mu1,double mu2);
double tvalue(double mu1,double mu2)
{
    double res;
    res=(mu1-mu2);
    res=res/(mu1+mu2);
    return res;
}
double Smulliken_s_s(double p,double t,int type1,int type2);
double Smulliken_s_s(double p,double t,int type1,int type2)
{
    double res=0.0;
    if(type1==1&&type2==1) //1s-1s
    {
        if(equal(t,0.0)==1)
        {res=(1.0/6.0)*(p*p*p)*(3.0*Afunction(2,p)-Afunction(0,p));}
        else
        {res=(1.0/4.0)*(p*p*p)*pow((1.0-t*t),1.5)*(Afunction(2,p)*Bfunction(0,p,t)-Afunction(0,p)*Bfunction(2,p,t));}
    }
    if(type1==2&&type2==2) //2s-2s
    {
        if(equal(t,0.0)==1)
        {res=(1.0/360.0)*(p*p*p*p*p)*(15.0*Afunction(4,p)-10.0*Afunction(2,p)+3.0*Afunction(0,p));}
        else
        {res=(1.0/48.0)*(p*p*p*p*p)*pow((1.0-t*t),2.5)*(Afunction(4,p)*Bfunction(0,p,t)-2.0*Afunction(2,p)*Bfunction(2,p,t)+Afunction(0,p)*Bfunction(4,p,t));}
    }
    if((type1==1&&type2==2)||(type1==2&&type2==1)) //1s-2s
    {
        if(equal(t,0.0)==1)
        {res=(1.0/12.0)*pow(3.0,-0.5)*(p*p*p*p)*(3.0*Afunction(3,p)-Afunction(1,p));}
        else
        {res=(1.0/8.0)*pow(3.0,-0.5)*(p*p*p*p)*pow((1.0+t),1.5)*pow((1.0-t),2.5)*(Afunction(3,p)*Bfunction(0,p,t)-Afunction(2,p)*Bfunction(1,p,t)-Afunction(1,p)*Bfunction(2,p,t)+Afunction(0,p)*Bfunction(3,p,t));}
    }
    if((type1==1&&type2==6)||(type1==6&&type2==1)) //1s-3s
    {
        if(equal(t,0.0)==1)
        {res=(1.0/60.0)*(1.0/sqrt(10.0))*(p*p*p*p*p)*(5.0*Afunction(4,p)-Afunction(0,p));}
        else
        {res=(1.0/24.0)*(1.0/sqrt(10.0))*(p*p*p*p*p)*pow((1.0+t),1.5)*pow((1.0-t),3.5)*(Afunction(4,p)*Bfunction(0,p,t)-2.0*Afunction(3,p)*Bfunction(1,p,t)+2.0*Afunction(1,p)*Bfunction(3,p,t)-Afunction(0,p)*Bfunction(4,p,t));}
    }
    if((type1==2&&type2==6)||(type1==6&&type2==2)) //2s-3s
    {
        if(equal(t,0.0)==1)
        {res=(1.0/360.0)*(1.0/sqrt(30.0))*(p*p*p*p*p*p)*(15.0*Afunction(5,p)-10.0*Afunction(3,p)+3.0*Afunction(1,p));}
        else
        {res=(1.0/48.0)*(1.0/sqrt(30.0))*(p*p*p*p*p*p)*pow((1.0+t),2.5)*pow((1.0-t),3.5)*(Afunction(5,p)*Bfunction(0,p,t)-Afunction(4,p)*Bfunction(1,p,t)-2.0*Afunction(3,p)*Bfunction(2,p,t)+2.0*Afunction(2,p)*Bfunction(3,p,t)+Afunction(1,p)*Bfunction(4,p,t)-Afunction(0,p)*Bfunction(5,p,t));}
    }
    if(type1==6&&type2==6) //3s-3s
    {
        if(equal(t,0.0)==1)
        {res=(1.0/25200.0)*(p*p*p*p*p*p*p)*(35.0*Afunction(6,p)-35.0*Afunction(4,p)+21.0*Afunction(2,p)-5.0*Afunction(0,p));}
        else
        {res=(1.0/1440.0)*(p*p*p*p*p*p*p)*pow((1.0-t*t),3.5)*(Afunction(6,p)*Bfunction(0,p,t)-3.0*Afunction(4,p)*Bfunction(2,p,t)+3.0*Afunction(2,p)*Bfunction(4,p,t)-Afunction(0,p)*Bfunction(6,p,t));}
    }
    return res;
}
double Smulliken_s_psigma(double p,double t,int type1,int type2);
double Smulliken_s_psigma(double p,double t,int type1,int type2)
{
    double res=0.0;
    if((type1==1&&(type2>=3&&type2<=5))||(type2==1&&(type1>=3&&type1<=5))) //1s-2psigma
    {
        if(equal(t,0.0)==1)
        {res=(1.0/12.0)*(p*p*p*p)*(3.0*Afunction(2,p)-Afunction(0,p));}
        else
        {res=(1.0/8.0)*(p*p*p*p)*pow((1.0+t),1.5)*pow((1.0-t),2.5)*(-1.0*Afunction(3,p)*Bfunction(1,p,t)+Afunction(2,p)*Bfunction(0,p,t)+Afunction(1,p)*Bfunction(3,p,t)-Afunction(0,p)*Bfunction(2,p,t));}
    }
    if((type1==2&&(type2>=3&&type2<=5))||(type2==2&&(type1>=3&&type1<=5))) //2s-2psigma
    {
        if(equal(t,0.0)==1)
        {res=(1.0/60.0)*(1.0/sqrt(3.0))*(p*p*p*p*p)*(5.0*Afunction(3,p)-Afunction(1,p));}
        else
        {res=(1.0/16.0)*(1.0/sqrt(3.0))*(p*p*p*p*p)*pow((1.0-t*t),2.5)*(Afunction(3,p)*(Bfunction(0,p,t)-Bfunction(2,p,t))+Afunction(1,p)*(Bfunction(4,p,t)-Bfunction(2,p,t))+Bfunction(1,p,t)*(Afunction(2,p)-Afunction(4,p))+Bfunction(3,p,t)*(Afunction(2,p)-Afunction(0,p)));}
    }
    if((type1==1&&(type2>=7&&type2<=9))||(type2==1&&(type1>=7&&type1<=9))) //1s-3psigma
    {
        if(equal(t,0.0)==1)
        {res=(1.0/15.0)*(1.0/sqrt(30.0))*(p*p*p*p*p)*(5.0*Afunction(3,p)-2.0*Afunction(1,p));}
        else
        {res=(1.0/8.0)*(1.0/sqrt(30.0))*(p*p*p*p*p)*pow((1.0+t),1.5)*pow((1-t),3.5)*(Afunction(3,p)*(Bfunction(0,p,t)+Bfunction(2,p,t))-Afunction(1,p)*(Bfunction(2,p,t)+Bfunction(4,p,t))-Bfunction(1,p,t)*(Afunction(2,p)+Afunction(4,p))+Bfunction(3,p,t)*(Afunction(0,p)+Afunction(2,p)));}
    }
    if((type1==2&&(type2>=7&&type2<=9))||(type2==2&&(type1>=7&&type1<=9))) //2s-3psigma
    {
        if(equal(t,0.0)==1)
        {res=(1.0/360.0)*(1.0/sqrt(10.0))*(p*p*p*p*p*p)*(15.0*Afunction(4,p)-10.0*Afunction(2,p)+3.0*Afunction(0,p));}
        else
        {res=(1.0/48.0)*(1.0/sqrt(10.0))*(p*p*p*p*p*p)*pow((1.0+t),2.5)*pow((1-t),3.5)*(-1.0*Afunction(5,p)*Bfunction(1,p,t)+Afunction(4,p)*Bfunction(0,p,t)+2.0*Afunction(3,p)*Bfunction(3,p,t)-2.0*Afunction(2,p)*Bfunction(2,p,t)-Afunction(1,p)*Bfunction(5,p,t)+Afunction(0,p)*Bfunction(4,p,t));}
    }
    if((type1==6&&(type2>=3&&type2<=5))||(type2==6&&(type1>=3&&type1<=5))) //3s-2psigma
    {
        if(equal(t,0.0)==1)
        {res=(1.0/360.0)*(1.0/sqrt(10.0))*(p*p*p*p*p*p)*(5.0*Afunction(4,p)+6.0*Afunction(2,p)-3.0*Afunction(0,p));}
        else
        {res=(1.0/48.0)*(1.0/sqrt(10.0))*(p*p*p*p*p*p)*pow((1.0+t),2.5)*pow((1.0-t),3.5)*(Afunction(4,p)*(Bfunction(0,p,t)-2.0*Bfunction(2,p,t))+Afunction(1,p)*(2.0*Bfunction(3,p,t)-Bfunction(5,p,t))+Bfunction(1,p,t)*(Afunction(5,p)-2.0*Afunction(3,p))+Bfunction(4,p,t)*(2.0*Afunction(2,p)-Afunction(0,p)));}
    }
    if((type1==6&&(type2>=7&&type2<=9))||(type2==6&&(type1>=7&&type1<=9))) //3s-3psigma
    {
        if(equal(t,0.0)==1)
        {res=(1.0/12600.0)*(1.0/sqrt(3.0))*(p*p*p*p*p*p*p)*(35.0*Afunction(5,p)-14.0*Afunction(3,p)+3.0*Afunction(1,p));}
        else
        {res=(1.0/480.0)*(1.0/sqrt(3.0))*(p*p*p*p*p*p*p)*pow((1.0-t*t),3.5)*(-1.0*Afunction(6,p)*Bfunction(1,p,t)+Afunction(5,p)*(Bfunction(0,p,t)-Bfunction(2,p,t))+Afunction(4,p)*(Bfunction(1,p,t)+2.0*Bfunction(3,p,t))+2.0*Afunction(3,p)*(Bfunction(4,p,t)-Bfunction(2,p,t))-Afunction(2,p)*(2.0*Bfunction(3,p,t)+Bfunction(5,p,t))+Afunction(1,p)*(Bfunction(4,p,t)-Bfunction(6,p,t))+Afunction(0,p)*Bfunction(5,p,t));}
    }
    return res;
}
double Smulliken_psigma_psigma(double p,double t,int type1,int type2);
double Smulliken_psigma_psigma(double p,double t,int type1,int type2)
{
    double res=0.0;
    if((type1>=3&&type1<=5)&&(type2>=3&&type2<=5)) //2psigma-2psigma
    {
        if(equal(t,0.0)==1)
        {res=(1.0/120.0)*(p*p*p*p*p)*(5.0*Afunction(4,p)-18.0*Afunction(2,p)+5.0*Afunction(0,p));}
        else
        {res=(1.0/16.0)*(p*p*p*p*p)*pow((1.0-t*t),2.5)*(Bfunction(2,p,t)*(Afunction(0,p)+Afunction(4,p))-Afunction(2,p)*(Bfunction(0,p,t)+Bfunction(4,p,t)));}
    }
    if((type1>=3&&type1<=5)&&(type2>=7&&type2<=9)) //2psigma-3psigma
    {
        if(equal(t,0.0)==1)
        {res=(1.0/120.0)*(1.0/sqrt(30.0))*(p*p*p*p*p*p)*(5.0*Afunction(5,p)-18.0*Afunction(3,p)+5.0*Afunction(1,p));}
        else
        {res=(1.0/16.0)*(1.0/sqrt(30.0))*(p*p*p*p*p*p)*pow((1.0+t),2.5)*pow((1.0-t),3.5)*(Afunction(2,p)*(Bfunction(1,p,t)+Bfunction(5,p,t))-Afunction(3,p)*(Bfunction(0,p,t)+Bfunction(4,p,t))-Bfunction(3,p,t)*(Afunction(0,p)+Afunction(4,p))+Bfunction(2,p,t)*(Afunction(1,p)+Afunction(5,p)));}
    }
    if((type1>=7&&type1<=9)&&(type2>=7&&type2<=9)) //3psigma-3psigma
    {
        if(equal(t,0.0)==1)
        {res=(1.0/25200.0)*(p*p*p*p*p*p*p)*(35.0*Afunction(6,p)-147.0*Afunction(4,p)+85.0*Afunction(2,p)-21.0*Afunction(0,p));}
        else
        {res=(1.0/480.0)*(p*p*p*p*p*p*p)*pow((1.0-t*t),3.5)*(Afunction(6,p)*Bfunction(2,p,t)-Afunction(4,p)*(Bfunction(0,p,t)+2.0*Bfunction(4,p,t))+Afunction(2,p)*(Bfunction(6,p,t)+2.0*Bfunction(2,p,t))-Afunction(0,p)*Bfunction(4,p,t));}
    }
    return res;
}
double Smulliken_ppi_ppi(double p,double t,int type1,int type2);
double Smulliken_ppi_ppi(double p,double t,int type1,int type2)
{
    double res=0.0;
    if((type1>=3&&type1<=5)&&(type2>=3&&type2<=5)) //2ppi-2ppi
    {
        if(equal(t,0.0)==1)
        {res=(1.0/120.0)*(p*p*p*p*p)*(5.0*Afunction(4,p)-6.0*Afunction(2,p)+Afunction(0,p));}
        else
        {res=(1.0/32.0)*(p*p*p*p*p)*pow((1.0-t*t),2.5)*(Afunction(4,p)*(Bfunction(0,p,t)-Bfunction(2,p,t))+Afunction(2,p)*(Bfunction(4,p,t)-Bfunction(0,p,t))+Afunction(0,p)*(Bfunction(2,p,t)-Bfunction(4,p,t)));}
    }
    if((type1>=3&&type1<=5)&&(type2>=7&&type2<=9)) //2ppi-3ppi
    {
        if(equal(t,0.0)==1)
        {res=(1.0/120.0)*(1.0/sqrt(30.0))*(p*p*p*p*p*p)*(5.0*Afunction(5,p)-6.0*Afunction(3,p)+Afunction(1,p));}
        else
        {res=(1.0/32.0)*(1.0/sqrt(30.0))*(p*p*p*p*p*p)*pow((1.0+t),2.5)*pow((1.0-t),3.5)*(Afunction(5,p)*(Bfunction(0,p,t)-Bfunction(2,p,t))+Afunction(4,p)*(Bfunction(3,p,t)-Bfunction(1,p,t))+Afunction(3,p)*(Bfunction(4,p,t)-Bfunction(0,p,t))+Afunction(2,p)*(Bfunction(1,p,t)-Bfunction(5,p,t))+Afunction(1,p)*(Bfunction(2,p,t)-Bfunction(4,p,t))+Afunction(0,p)*(Bfunction(5,p,t)-Bfunction(3,p,t)));}
    }
    if((type1>=7&&type1<=9)&&(type2>=7&&type2<=9)) //3ppi-3ppi
    {
        if(equal(t,0.0)==1)
        {res=(1.0/25200.0)*(p*p*p*p*p*p*p)*(35.0*Afunction(6,p)-49.0*Afunction(4,p)+17.0*Afunction(2,p)-3.0*Afunction(0,p));}
        else
        {res=(1.0/960.0)*(p*p*p*p*p*p*p)*pow((1.0-t*t),3.5)*(Afunction(6,p)*(Bfunction(0,p,t)-Bfunction(2,p,t))+Afunction(4,p)*(2.0*Bfunction(4,p,t)-Bfunction(0,p,t)-Bfunction(2,p,t))+Afunction(2,p)*(2.0*Bfunction(2,p,t)-Bfunction(4,p,t)-Bfunction(6,p,t))+Afunction(0,p)*(Bfunction(6,p,t)-Bfunction(4,p,t)));}
    }
    return res;
}
double overlap(double X1,double Y1,double Z1,double X2,double Y2,double Z2,double mu1,double mu2,int type1,int type2);
double overlap(double X1,double Y1,double Z1,double X2,double Y2,double Z2,double mu1,double mu2,int type1,int type2)
{
    double exb[3],eyb[3],ezb[3];
    double res=0.0,S,Spi,r1,r2,r3;
    double vector[3][3]={{1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};
    // s-s
    if(type1==1&&type2==1){res=Smulliken_s_s(pvalue(X1,Y1,Z1,X2,Y2,Z2,mu1,mu2),tvalue(mu1,mu2),type1,type2);}
    if(type1==2&&type2==2){res=Smulliken_s_s(pvalue(X1,Y1,Z1,X2,Y2,Z2,mu1,mu2),tvalue(mu1,mu2),type1,type2);}
    if(type1==6&&type2==6){res=Smulliken_s_s(pvalue(X1,Y1,Z1,X2,Y2,Z2,mu1,mu2),tvalue(mu1,mu2),type1,type2);}
    if(type1==1&&type2==2){res=Smulliken_s_s(pvalue(X1,Y1,Z1,X2,Y2,Z2,mu1,mu2),tvalue(mu1,mu2),type1,type2);}
    if(type1==2&&type2==1){res=Smulliken_s_s(pvalue(X2,Y2,Z2,X1,Y1,Z1,mu2,mu1),tvalue(mu2,mu1),type2,type1);}
    if(type1==1&&type2==6){res=Smulliken_s_s(pvalue(X1,Y1,Z1,X2,Y2,Z2,mu1,mu2),tvalue(mu1,mu2),type1,type2);}
    if(type1==6&&type2==1){res=Smulliken_s_s(pvalue(X2,Y2,Z2,X1,Y1,Z1,mu2,mu1),tvalue(mu2,mu1),type2,type1);}
    if(type1==2&&type2==6){res=Smulliken_s_s(pvalue(X1,Y1,Z1,X2,Y2,Z2,mu1,mu2),tvalue(mu1,mu2),type1,type2);}
    if(type1==6&&type2==2){res=Smulliken_s_s(pvalue(X2,Y2,Z2,X1,Y1,Z1,mu2,mu1),tvalue(mu2,mu1),type2,type1);}
    // s-p
    if(type1==1&&(type2>=3&&type2<=5)){
        S=Smulliken_s_psigma(pvalue(X1,Y1,Z1,X2,Y2,Z2,mu1,mu2),tvalue(mu1,mu2),type1,type2);
        MakeLocalCoordsSystems(X1,Y1,Z1,X2,Y2,Z2,&exb[0],&exb[1],&exb[2],&eyb[0],&eyb[1],&eyb[2],&ezb[0],&ezb[1],&ezb[2]);
        res=vector[type2-2-1][0]*ezb[0]+vector[type2-2-1][1]*ezb[1]+vector[type2-2-1][2]*ezb[2];
        res=res*S;
    }
    if(type2==1&&(type1>=3&&type1<=5)){
        S=Smulliken_s_psigma(pvalue(X2,Y2,Z2,X1,Y1,Z1,mu2,mu1),tvalue(mu2,mu1),type2,type1);
        MakeLocalCoordsSystems(X2,Y2,Z2,X1,Y1,Z1,&exb[0],&exb[1],&exb[2],&eyb[0],&eyb[1],&eyb[2],&ezb[0],&ezb[1],&ezb[2]);
        res=vector[type1-2-1][0]*ezb[0]+vector[type1-2-1][1]*ezb[1]+vector[type1-2-1][2]*ezb[2];
        res=res*S;
    }
    if(type1==2&&(type2>=3&&type2<=5)){
        S=Smulliken_s_psigma(pvalue(X1,Y1,Z1,X2,Y2,Z2,mu1,mu2),tvalue(mu1,mu2),type1,type2);
        MakeLocalCoordsSystems(X1,Y1,Z1,X2,Y2,Z2,&exb[0],&exb[1],&exb[2],&eyb[0],&eyb[1],&eyb[2],&ezb[0],&ezb[1],&ezb[2]);
        res=vector[type2-2-1][0]*ezb[0]+vector[type2-2-1][1]*ezb[1]+vector[type2-2-1][2]*ezb[2];
        res=res*S;
    }
    if(type2==2&&(type1>=3&&type1<=5)){
        S=Smulliken_s_psigma(pvalue(X2,Y2,Z2,X1,Y1,Z1,mu2,mu1),tvalue(mu2,mu1),type2,type1);
        MakeLocalCoordsSystems(X2,Y2,Z2,X1,Y1,Z1,&exb[0],&exb[1],&exb[2],&eyb[0],&eyb[1],&eyb[2],&ezb[0],&ezb[1],&ezb[2]);
        res=vector[type1-2-1][0]*ezb[0]+vector[type1-2-1][1]*ezb[1]+vector[type1-2-1][2]*ezb[2];
        res=res*S;
    }
    if(type1==1&&(type2>=7&&type2<=9)){
        S=Smulliken_s_psigma(pvalue(X1,Y1,Z1,X2,Y2,Z2,mu1,mu2),tvalue(mu1,mu2),type1,type2);
        MakeLocalCoordsSystems(X1,Y1,Z1,X2,Y2,Z2,&exb[0],&exb[1],&exb[2],&eyb[0],&eyb[1],&eyb[2],&ezb[0],&ezb[1],&ezb[2]);
        res=vector[type2-6-1][0]*ezb[0]+vector[type2-6-1][1]*ezb[1]+vector[type2-6-1][2]*ezb[2];
        res=res*S;
    }
    if(type2==1&&(type1>=7&&type1<=9)){
        S=Smulliken_s_psigma(pvalue(X2,Y2,Z2,X1,Y1,Z1,mu2,mu1),tvalue(mu2,mu1),type2,type1);
        MakeLocalCoordsSystems(X2,Y2,Z2,X1,Y1,Z1,&exb[0],&exb[1],&exb[2],&eyb[0],&eyb[1],&eyb[2],&ezb[0],&ezb[1],&ezb[2]);
        res=vector[type1-6-1][0]*ezb[0]+vector[type1-6-1][1]*ezb[1]+vector[type1-6-1][2]*ezb[2];
        res=res*S;
    }
    if(type1==2&&(type2>=7&&type2<=9)){
        S=Smulliken_s_psigma(pvalue(X1,Y1,Z1,X2,Y2,Z2,mu1,mu2),tvalue(mu1,mu2),type1,type2);
        MakeLocalCoordsSystems(X1,Y1,Z1,X2,Y2,Z2,&exb[0],&exb[1],&exb[2],&eyb[0],&eyb[1],&eyb[2],&ezb[0],&ezb[1],&ezb[2]);
        res=vector[type2-6-1][0]*ezb[0]+vector[type2-6-1][1]*ezb[1]+vector[type2-6-1][2]*ezb[2];
        res=res*S;
    }
    if(type2==2&&(type1>=7&&type1<=9)){
        S=Smulliken_s_psigma(pvalue(X2,Y2,Z2,X1,Y1,Z1,mu2,mu1),tvalue(mu2,mu1),type2,type1);
        MakeLocalCoordsSystems(X2,Y2,Z2,X1,Y1,Z1,&exb[0],&exb[1],&exb[2],&eyb[0],&eyb[1],&eyb[2],&ezb[0],&ezb[1],&ezb[2]);
        res=vector[type1-6-1][0]*ezb[0]+vector[type1-6-1][1]*ezb[1]+vector[type1-6-1][2]*ezb[2];
        res=res*S;
    }
    if((type1>=3&&type1<=5)&&type2==6){
        S=Smulliken_s_psigma(pvalue(X1,Y1,Z1,X2,Y2,Z2,mu1,mu2),tvalue(mu1,mu2),type1,type2);
        MakeLocalCoordsSystems(X1,Y1,Z1,X2,Y2,Z2,&exb[0],&exb[1],&exb[2],&eyb[0],&eyb[1],&eyb[2],&ezb[0],&ezb[1],&ezb[2]);
        res=vector[type1-2-1][0]*ezb[0]+vector[type1-2-1][1]*ezb[1]+vector[type1-2-1][2]*ezb[2];
        res=res*S;
    }
    if((type2>=3&&type2<=5)&&type1==6){
        S=Smulliken_s_psigma(pvalue(X2,Y2,Z2,X1,Y1,Z1,mu2,mu1),tvalue(mu2,mu1),type2,type1);
        MakeLocalCoordsSystems(X2,Y2,Z2,X1,Y1,Z1,&exb[0],&exb[1],&exb[2],&eyb[0],&eyb[1],&eyb[2],&ezb[0],&ezb[1],&ezb[2]);
        res=vector[type2-2-1][0]*ezb[0]+vector[type2-2-1][1]*ezb[1]+vector[type2-2-1][2]*ezb[2];
        res=res*S;
    }
    if(type1==6&&(type2>=7&&type2<=9)){
        S=Smulliken_s_psigma(pvalue(X1,Y1,Z1,X2,Y2,Z2,mu1,mu2),tvalue(mu1,mu2),type1,type2);
        MakeLocalCoordsSystems(X1,Y1,Z1,X2,Y2,Z2,&exb[0],&exb[1],&exb[2],&eyb[0],&eyb[1],&eyb[2],&ezb[0],&ezb[1],&ezb[2]);
        res=vector[type2-6-1][0]*ezb[0]+vector[type2-6-1][1]*ezb[1]+vector[type2-6-1][2]*ezb[2];
        res=res*S;
    }
    if(type2==6&&(type1>=7&&type1<=9)){
        S=Smulliken_s_psigma(pvalue(X2,Y2,Z2,X1,Y1,Z1,mu2,mu1),tvalue(mu2,mu1),type2,type1);
        MakeLocalCoordsSystems(X2,Y2,Z2,X1,Y1,Z1,&exb[0],&exb[1],&exb[2],&eyb[0],&eyb[1],&eyb[2],&ezb[0],&ezb[1],&ezb[2]);
        res=vector[type1-6-1][0]*ezb[0]+vector[type1-6-1][1]*ezb[1]+vector[type1-6-1][2]*ezb[2];
        res=res*S;
    }
    if((type1>=3&&type1<=5)&&(type2>=3&&type2<=5)){
        S=Smulliken_psigma_psigma(pvalue(X1,Y1,Z1,X2,Y2,Z2,mu1,mu2),tvalue(mu1,mu2),type1,type2);
        Spi=Smulliken_ppi_ppi(pvalue(X1,Y1,Z1,X2,Y2,Z2,mu1,mu2),tvalue(mu1,mu2),type1,type2);
        MakeLocalCoordsSystems(X1,Y1,Z1,X2,Y2,Z2,&exb[0],&exb[1],&exb[2],&eyb[0],&eyb[1],&eyb[2],&ezb[0],&ezb[1],&ezb[2]);
        r1=(vector[type1-2-1][0]*ezb[0]+vector[type1-2-1][1]*ezb[1]+vector[type1-2-1][2]*ezb[2]);
        r1=r1*(vector[type2-2-1][0]*ezb[0]+vector[type2-2-1][1]*ezb[1]+vector[type2-2-1][2]*ezb[2]);
        r2=(vector[type1-2-1][0]*eyb[0]+vector[type1-2-1][1]*eyb[1]+vector[type1-2-1][2]*eyb[2]);
        r2=r2*(vector[type2-2-1][0]*eyb[0]+vector[type2-2-1][1]*eyb[1]+vector[type2-2-1][2]*eyb[2]);
        r3=(vector[type1-2-1][0]*exb[0]+vector[type1-2-1][1]*exb[1]+vector[type1-2-1][2]*exb[2]);
        r3=r3*(vector[type2-2-1][0]*exb[0]+vector[type2-2-1][1]*exb[1]+vector[type2-2-1][2]*exb[2]);
        res=r1*S+(r2+r3)*Spi;
    }
    if((type1>=3&&type1<=5)&&(type2>=7&&type2<=9)){
        S=Smulliken_psigma_psigma(pvalue(X1,Y1,Z1,X2,Y2,Z2,mu1,mu2),tvalue(mu1,mu2),type1,type2);
        Spi=Smulliken_ppi_ppi(pvalue(X1,Y1,Z1,X2,Y2,Z2,mu1,mu2),tvalue(mu1,mu2),type1,type2);
        MakeLocalCoordsSystems(X1,Y1,Z1,X2,Y2,Z2,&exb[0],&exb[1],&exb[2],&eyb[0],&eyb[1],&eyb[2],&ezb[0],&ezb[1],&ezb[2]);
        r1=   (vector[type1-2-1][0]*ezb[0]+vector[type1-2-1][1]*ezb[1]+vector[type1-2-1][2]*ezb[2]);
        r1=r1*(vector[type2-6-1][0]*ezb[0]+vector[type2-6-1][1]*ezb[1]+vector[type2-6-1][2]*ezb[2]);
        r2=   (vector[type1-2-1][0]*eyb[0]+vector[type1-2-1][1]*eyb[1]+vector[type1-2-1][2]*eyb[2]);
        r2=r2*(vector[type2-6-1][0]*eyb[0]+vector[type2-6-1][1]*eyb[1]+vector[type2-6-1][2]*eyb[2]);
        r3=   (vector[type1-2-1][0]*exb[0]+vector[type1-2-1][1]*exb[1]+vector[type1-2-1][2]*exb[2]);
        r3=r3*(vector[type2-6-1][0]*exb[0]+vector[type2-6-1][1]*exb[1]+vector[type2-6-1][2]*exb[2]);
        res=r1*S+(r2+r3)*Spi;
    }
    if((type2>=3&&type2<=5)&&(type1>=7&&type1<=9)){
        S=Smulliken_psigma_psigma(pvalue(X2,Y2,Z2,X1,Y1,Z1,mu2,mu1),tvalue(mu2,mu1),type2,type1);
        Spi=Smulliken_ppi_ppi(pvalue(X2,Y2,Z2,X1,Y1,Z1,mu2,mu1),tvalue(mu2,mu1),type2,type1);
        MakeLocalCoordsSystems(X2,Y2,Z2,X1,Y1,Z1,&exb[0],&exb[1],&exb[2],&eyb[0],&eyb[1],&eyb[2],&ezb[0],&ezb[1],&ezb[2]);
        r1=   (vector[type2-2-1][0]*ezb[0]+vector[type2-2-1][1]*ezb[1]+vector[type2-2-1][2]*ezb[2]);
        r1=r1*(vector[type1-6-1][0]*ezb[0]+vector[type1-6-1][1]*ezb[1]+vector[type1-6-1][2]*ezb[2]);
        r2=   (vector[type2-2-1][0]*eyb[0]+vector[type2-2-1][1]*eyb[1]+vector[type2-2-1][2]*eyb[2]);
        r2=r2*(vector[type1-6-1][0]*eyb[0]+vector[type1-6-1][1]*eyb[1]+vector[type1-6-1][2]*eyb[2]);
        r3=   (vector[type2-2-1][0]*exb[0]+vector[type2-2-1][1]*exb[1]+vector[type2-2-1][2]*exb[2]);
        r3=r3*(vector[type1-6-1][0]*exb[0]+vector[type1-6-1][1]*exb[1]+vector[type1-6-1][2]*exb[2]);
        res=r1*S+(r2+r3)*Spi;
    }
    if((type1>=7&&type1<=9)&&(type2>=7&&type2<=9)){
        S=Smulliken_psigma_psigma(pvalue(X1,Y1,Z1,X2,Y2,Z2,mu1,mu2),tvalue(mu1,mu2),type1,type2);
        Spi=Smulliken_ppi_ppi(pvalue(X1,Y1,Z1,X2,Y2,Z2,mu1,mu2),tvalue(mu1,mu2),type1,type2);
        MakeLocalCoordsSystems(X1,Y1,Z1,X2,Y2,Z2,&exb[0],&exb[1],&exb[2],&eyb[0],&eyb[1],&eyb[2],&ezb[0],&ezb[1],&ezb[2]);
        r1=   (vector[type1-6-1][0]*ezb[0]+vector[type1-6-1][1]*ezb[1]+vector[type1-6-1][2]*ezb[2]);
        r1=r1*(vector[type2-6-1][0]*ezb[0]+vector[type2-6-1][1]*ezb[1]+vector[type2-6-1][2]*ezb[2]);
        r2=   (vector[type1-6-1][0]*eyb[0]+vector[type1-6-1][1]*eyb[1]+vector[type1-6-1][2]*eyb[2]);
        r2=r2*(vector[type2-6-1][0]*eyb[0]+vector[type2-6-1][1]*eyb[1]+vector[type2-6-1][2]*eyb[2]);
        r3=   (vector[type1-6-1][0]*exb[0]+vector[type1-6-1][1]*exb[1]+vector[type1-6-1][2]*exb[2]);
        r3=r3*(vector[type2-6-1][0]*exb[0]+vector[type2-6-1][1]*exb[1]+vector[type2-6-1][2]*exb[2]);
        res=r1*S+(r2+r3)*Spi;
    }
    return res;
}
long int factorial(int n);
long int factorial(int n)
{
    if (n >= 1)
        return n*factorial(n-1);
    else
        return 1;
}
