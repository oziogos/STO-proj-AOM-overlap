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
double atomic_contrib(double x,double y,double z,double X,double Y,double Z,double c1s,double c2s,double c2px,double c2py,double c2pz,double c3s,double c3px,double c3py,double c3pz,double mu1s,double mu2s,double mu2p,double mu3s,double mu3p);
void create_cube_file(char *current_folder,char *STOproj_cube_grid,char *STOproj_name,char *STOproj_MO,int atoms,char **species,double *x,double *y,double *z,double **STO_matrix,double *smu_per_species,double *pmu_per_species,int verb);
void create_cube_file(char *current_folder,char *STOproj_cube_grid,char *STOproj_name,char *STOproj_MO,int atoms,char **species,double *x,double *y,double *z,double **STO_matrix,double *smu_per_species,double *pmu_per_species,int verb)
{
    FILE *fp;
    char file_path[cmax_length];
    int i,j;
    
    double xmin,xmax,ymin,ymax,zmin,zmax;
    int resx,resy,resz;
    double length;
    int ix,iy,iz,NX,NY,NZ;
    double sum;
    double X,Y,Z;
    int atomicZ;
    double c1s,c2s,c2px,c2py,c2pz,c3s,c3px,c3py,c3pz,mu1s,mu2s,mu2p,mu3s,mu3p;
    
    length=atof(STOproj_cube_grid);
    
    xmin=x[0];xmax=x[0];ymin=y[0];ymax=y[0];zmin=z[0];zmax=z[0];
    for(j=0;j<atoms;++j)
    {
        if(x[j]>xmax)xmax=x[j];
        if(x[j]<xmin)xmin=x[j];
        if(y[j]>ymax)ymax=y[j];
        if(y[j]<ymin)ymin=y[j];
        if(z[j]>zmax)zmax=z[j];
        if(z[j]<zmin)zmin=z[j];
    }
    
    xmin=xmin-offset;xmax=xmax+offset;
    ymin=ymin-offset;ymax=ymax+offset;
    zmin=zmin-offset;zmax=zmax+offset;
    resx=(xmax-xmin)/length+1;
    resy=(ymax-ymin)/length+1;
    resz=(zmax-zmin)/length+1;
    
    sprintf(file_path,"%s/%s_STO_MO_%s.cube",current_folder,STOproj_name,STOproj_MO);
    fp=fopen(file_path,"w+");
    
    fprintf(fp,"My CUBE\n\n");
    fprintf(fp,"%d\t%lf\t%lf\t%lf\n",atoms,xmin,ymin,zmin);
    fprintf(fp,"%d\t%lf\t%lf\t%lf\n",resx,length,0.0,0.0);
    fprintf(fp,"%d\t%lf\t%lf\t%lf\n",resy,0.0,length,0.0);
    fprintf(fp,"%d\t%lf\t%lf\t%lf\n",resz,0.0,0.0,length);
    for(i=0;i<atoms;++i){fprintf(fp,"%d\t%lf\t%lf\t%lf\t%lf\n",resolve_atomic_Z(species[i]),(double)resolve_atomic_Z(species[i]),x[i],y[i],z[i]);}
    NX=resx;NY=resy;NZ=resz;
    for (ix=0;ix<NX;ix++) {
        for (iy=0;iy<NY;iy++) {
            for (iz=0;iz<NZ;iz++) {
                
                X=xmin+ix*length;
                Y=ymin+iy*length;
                Z=zmin+iz*length;
                sum=0.0;
                for(i=0;i<atoms;++i)
                {
                    atomicZ=resolve_atomic_Z(species[i]);
                    if(atomicZ==1){
                        c1s=STO_matrix[i][0];
                        c2s=0.0;
                        c2px=0.0;
                        c2py=0.0;
                        c2pz=0.0;
                        c3s=0.0;
                        c3px=0.0;
                        c3py=0.0;
                        c3pz=0.0;
                        mu1s=smu_per_species[atomicZ-1];
                        mu2s=0.0;
                        mu2p=0.0;
                        mu3s=0.0;
                        mu3p=0.0;
                    }
                    if(atomicZ==6 || atomicZ==7 || atomicZ==8 || atomicZ==9){
                        c1s=0.0;
                        c2s=STO_matrix[i][0];
                        c2px=STO_matrix[i][1];
                        c2py=STO_matrix[i][2];
                        c2pz=STO_matrix[i][3];
                        c3s=0.0;
                        c3px=0.0;
                        c3py=0.0;
                        c3pz=0.0;
                        mu1s=0.0;
                        mu2s=smu_per_species[atomicZ-1];
                        mu2p=pmu_per_species[atomicZ-1];
                        mu3s=0.0;
                        mu3p=0.0;
                    }
                    if(atomicZ==16){
                        c1s=0.0;
                        c2s=0.0;
                        c2px=0.0;
                        c2py=0.0;
                        c2pz=0.0;
                        c3s=STO_matrix[i][0];
                        c3px=STO_matrix[i][1];
                        c3py=STO_matrix[i][2];
                        c3pz=STO_matrix[i][3];
                        mu1s=0.0;
                        mu2s=0.0;
                        mu2p=0.0;
                        mu3s=smu_per_species[atomicZ-1];
                        mu3p=pmu_per_species[atomicZ-1];
                    }
                    
                    sum=sum+atomic_contrib(
                                           X,Y,Z,x[i],y[i],z[i],
                                           c1s,c2s,c2px,c2py,c2pz,c3s,c3px,c3py,c3pz,
                                           mu1s,mu2s,mu2p,mu3s,mu3p);
                }
                if(sum*sum>print_thres)
                    if(sum<0.0)
                        fprintf(fp,"%.3g ",-fabs(sum));
                    else
                        fprintf(fp,"%.3g ",fabs(sum));
                    else
                        fprintf(fp,"%.3g ",0.0);
                if (iz % 6 == 5)
                    fprintf(fp,"\n");
            }
            fprintf(fp,"\n");
        }
    }
    
    fclose(fp);
    if(verb==1||verb==-1){
        // console output
        printf("Generated STO cube file %s for preview...\n",file_path);
        printf("\n--------------------------------------------------------------------------\n\n");
    }
}
// x,y,z are the GENERIC coordinates
// X,Y,Z are the ATOMIC coordinates
double atomic_contrib(double x,double y,double z,double X,double Y,double Z,double c1s,double c2s,double c2px,double c2py,double c2pz,double c3s,double c3px,double c3py,double c3pz,double mu1s,double mu2s,double mu2p,double mu3s,double mu3p);
double atomic_contrib(double x,double y,double z,double X,double Y,double Z,double c1s,double c2s,double c2px,double c2py,double c2pz,double c3s,double c3px,double c3py,double c3pz,double mu1s,double mu2s,double mu2p,double mu3s,double mu3p)
{
    double res,R;
    R=sqrt((x-X)*(x-X)+(y-Y)*(y-Y)+(z-Z)*(z-Z));
    res=(0.5641895835477563*c1s*pow(mu1s,1.5))/exp(mu1s*R) + (0.32573500793527993*c2s*pow(mu2s,2.5)*R)/exp(mu2s*R) + (0.11894160774351806*c3s*pow(mu3s,3.5)*R*R)/exp(mu3s*R) + (0.5641895835477563*pow(mu2p,2.5)*(c2px*(x - X) + c2py*(y - Y) + c2pz*(z - Z)))/exp(mu2p*R) + (0.20601290774570113*pow(mu3p,3.5)*R*(c3px*(x - X) + c3py*(y - Y) + c3pz*(z - Z)))/exp(mu3p*R);
    return res;
}

