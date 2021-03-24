#include"general.h"

#define simplex_equal_diff 1.0e-16
#define pr_acc ".6"

int simplex_equal(double A,double B);

double simplex_obj_f(double *smu_per_species,double *pmu_per_species,struct simplex_obj_f_input current_input,double **compl_array);

double AOM_simplex_obj_f(double *smu_per_species,double *pmu_per_species,struct simplex_obj_f_input current_input,double **res);

void run_simplex(double *smu_per_species,double *pmu_per_species,struct simplex_obj_f_input current_input,struct simplex_params simplex);
void run_simplex(double *smu_per_species,double *pmu_per_species,struct simplex_obj_f_input current_input,struct simplex_params simplex)

{

    char mask[cmax_length],action[cmax_length];
    int n,i,j,flag,shrink;
    double rho=1.0,xi=2.0,gamma=0.5,sigma=0.5;
    double d_buffer,fr,fe,foc,fic;
    double **x,*f;
    double *xr,*xC,*xe,*xoc,*xic;
    double *array;
    int break_sum;
    double tau=0.05,tau_prime=0.00025;
    double *xmin,*xmax;
    int k,kmax;
    int outside;
    double **bc;
    
    int *index;
    double break_acc;
    
    double *compl_array,metric;
    
    double *res;
    
    if(current_input.mode==0)
        compl_array=(double*)malloc(current_input.simplex_entries*sizeof(double));
    else
        res=(double*)malloc(current_input.simplex_entries*sizeof(double));
    
    //
    kmax=simplex.steps;
    break_acc=simplex.acc;
    rho=simplex.rho;
    xi=simplex.xi;
    gamma=simplex.gamma;
    sigma=simplex.sigma;
    tau=simplex.tau;
    tau_prime=simplex.tau_prime;
    //
    
    // resolve active STO species
    n=0;
    for(i=0;i<16;++i)
        if(pmu_per_species[i]>1.0e-6)
            ++n;
    index=(int*)malloc(n*sizeof(int));
    j=-1;
    for(i=0;i<16;++i)
        if(pmu_per_species[i]>1.0e-6)
        {
            ++j;
            index[j]=i;
        }
    
    // preallocations
    f=(double*)malloc((n+1)*sizeof(double));
    x=(double**)malloc((n+1)*sizeof(double*));
    for(i=0;i<n+1;++i)x[i]=(double*)malloc(n*sizeof(double));
    xr=(double*)malloc(n*sizeof(double));
    xC=(double*)malloc(n*sizeof(double));
    xe=(double*)malloc(n*sizeof(double));
    xoc=(double*)malloc(n*sizeof(double));
    xic=(double*)malloc(n*sizeof(double));
    array=(double*)malloc(n*sizeof(double));
    bc=(double**)malloc(5*sizeof(double*));
    for(i=0;i<5;++i)bc[i]=(double*)malloc(n*sizeof(double));
    
    // read starting point and set constraints
    xmin=(double*)malloc(n*sizeof(double));
    xmax=(double*)malloc(n*sizeof(double));
    j=-1;
    for(i=0;i<16;++i)
        if(pmu_per_species[i]>1.0e-6)
        {
            ++j;
            x[0][j]=pmu_per_species[i];
            xmin[j]=simplex.p_mu_min[i];
            xmax[j]=simplex.p_mu_max[i];
        }
    
    //--------------------------------------------------------------------------
    
    // initialize
    
    // initialize the simplex
    for(i=1;i<n+1;++i)
    {
        for(j=0;j<n;++j)
        {
            x[i][j]=x[0][j];
            if(i-1==j)
            {
                if(simplex_equal(x[i][j],0.0)==1)
                    x[i][j]=tau_prime;
                else
                    x[i][j]=(1.0+tau)*x[0][j];
            }
        }
    }
    
    if(kmax>0)
    {
        // console out
        printf("\n*** Initialization...\n\n");
        printf("##\tInitial coordinates:\n");
        for(i=0;i<n+1;++i){printf("##\t");for(j=0;j<n;++j){sprintf(mask,"%%%slf\t",pr_acc);printf(mask,x[i][j]);}printf("\n");}
        printf("##\tCoordinates constraints:\n");
        for(i=0;i<n;++i){sprintf(mask,"##\txmin_%%d = %%%slf\txmax_%%d = %%%slf\n",pr_acc,pr_acc);printf(mask,i+1,xmin[i],i+1,xmax[i]);}
        printf("##\tTolerance = %e\n",break_acc);
        printf("\n");
        printf("*** Starting simplex...\n\n");
        
        // order
        for(i=0;i<n+1;++i)  // evaluate function and store in array f
        {
            for(j=0;j<n;++j)array[j]=x[i][j];
            
            for(j=0;j<n;++j)pmu_per_species[index[j]]=array[j];
            if(current_input.mode==0)
                f[i]=simplex_obj_f( smu_per_species,pmu_per_species,current_input,&compl_array);
            else
                f[i]=AOM_simplex_obj_f( smu_per_species,pmu_per_species,current_input,&res);
            
        }
        flag=1;             // bubble sort
        while(flag==1)
        {
            flag=0;
            for(i=1;i<n+1;++i)
            {
                if(f[i]<f[i-1])
                {
                    d_buffer=f[i-1];f[i-1]=f[i];f[i]=d_buffer;
                    for(j=0;j<n;++j)
                    {
                        d_buffer=x[i-1][j];x[i-1][j]=x[i][j];x[i][j]=d_buffer;
                    }
                    flag=1;
                }
            }
        }
        
        printf("%d\t",1);
        for(i=0;i<n;++i)printf("%e\t",x[0][i]);
        printf("f = %e\n",f[0]);
    }
    else if (kmax==0)
    {
        printf("\n*** Single evaluation...\n\n");
        if(current_input.mode==0)
        {
            metric=simplex_obj_f( smu_per_species,pmu_per_species,current_input,&compl_array);
            for(i=0;i<current_input.simplex_entries;++i)printf("%.6lf;",compl_array[i]);printf("\n");
        }
        else
        {
            metric=AOM_simplex_obj_f( smu_per_species,pmu_per_species,current_input,&res);
            for(i=0;i<current_input.simplex_entries;++i)printf("%.6lf;",res[i]);printf("\n");
            printf("%lf\n",metric);
        }
    }
    else
    {
        printf("\n*** Sampling mode...\n\n");
        //for(i=0;i<16;++i)printf("%lf\t%lf\t%lf\n",pmu_per_species[i],simplex.p_mu_sample_step[i],simplex.p_mu_sample_max[i]);
        
        int sampling_species=0;
        for(i=0;i<16;++i)
            if(simplex.p_mu_sample_step[i]>0.0)
                ++sampling_species;
        int sampling_Z[sampling_species];
        int sampling_N[sampling_species];
        double sampling_min[sampling_species],sampling_max[sampling_species],sampling_step[sampling_species];
        
        double pvalues[16];
        
        for(i=0;i<16;++i)pvalues[i]=pmu_per_species[i];
        
        j=-1;
        for(i=0;i<16;++i)
        {
            if(simplex.p_mu_sample_step[i]>0.0)
            {
                ++j;
                sampling_Z[j]=i+1;
                sampling_min[j]=pmu_per_species[i];
                sampling_max[j]=simplex.p_mu_sample_max[i];
                sampling_step[j]=simplex.p_mu_sample_step[i];
                sampling_N[j]=(int)floor((sampling_max[j]-sampling_min[j])/sampling_step[j])+1;
            }
        }
        
        for(i=0;i<sampling_species;++i)
        {
            printf("*** Atomic number %d: m%lf\tM%lf\ts%lf\tN%d\n",sampling_Z[i],sampling_min[i],sampling_max[i],sampling_step[i],sampling_N[i]);
            for(j=0;j<sampling_N[i];++j)printf("%.2lf,",sampling_min[i]+j*sampling_step[i]);printf("\n");
        }
        
        if(sampling_species==2)
        {
            printf("*** 2D mode\n");
            
            for(i=0;i<sampling_N[0];++i)
            {
                for(j=0;j<sampling_N[1];++j)
                {
                    pvalues[sampling_Z[0]-1]=sampling_min[0]+i*sampling_step[0];
                    pvalues[sampling_Z[1]-1]=sampling_min[1]+j*sampling_step[1];
                    if(current_input.mode==0)
                    {
                        metric=simplex_obj_f( smu_per_species,pvalues,current_input,&compl_array);
                        printf("$ %lf\t%lf\t%lf\t",sampling_min[0]+i*sampling_step[0],sampling_min[1]+j*sampling_step[1],metric);
                        for(int l=0;l<current_input.simplex_entries;++l)printf("%.6lf;",compl_array[l]);printf("\n");
                    }
                    else
                    {
                        metric=AOM_simplex_obj_f( smu_per_species,pvalues,current_input,&res);
                        printf("$ %lf\t%lf\t%lf\n",sampling_min[0]+i*sampling_step[0],sampling_min[1]+j*sampling_step[1],metric);
                    }
                }
            }
            
        }
        
    }
    // k-loop
    for(k=0;k<kmax;++k)
    {
        // reflect
        
        // calculate centroid excluding the n+1 point
        for(i=0;i<n;++i)    // loop on variables
        {
            xC[i]=0.0;
            for(j=0;j<n;++j)    // loop on points excluding n+1
            {
                xC[i]=xC[i]+x[j][i];
            }
            xC[i]=(double)xC[i]/n;  // normalize
        }
        // apply reflection
        for(i=0;i<n;++i)xr[i]=(1.0+rho)*xC[i]-rho*x[n][i];
        
        // check boundary conditions
        outside=0;
        for(i=0;i<n;++i)
            if(xr[i]<xmin[i]||xr[i]>xmax[i])
            {outside=1;break;}
        if(outside==1)
            fr=1.0e12;
        else
        {
            for(j=0;j<n;++j)pmu_per_species[index[j]]=xr[j];
            if(current_input.mode==0)
                fr=simplex_obj_f( smu_per_species,pmu_per_species,current_input,&compl_array);
            else
                fr=AOM_simplex_obj_f( smu_per_species,pmu_per_species,current_input,&res);
        }
        // evaluate function at xr
        
        // check the value of fr with respect to f[0] and f[n-1]
        shrink=0;   // shrink flag
        if((fr<f[n-1]&&fr>f[0])||(fr<f[n-1]&&simplex_equal(fr,f[0])==1))    // if .true., the iteration is completed
        {
            // overwrite point with xr
            for(i=0;i<n;++i)x[n][i]=xr[i];
            sprintf(action,"reflect");
        }
        else
        {
            if(fr<f[0]) // expand
            {
                // expansion
                for(i=0;i<n;++i)xe[i]=xC[i]+xi*(xr[i]-xC[i]);
                
                // check boundary conditions
                outside=0;
                for(i=0;i<n;++i)
                    if(xe[i]<xmin[i]||xe[i]>xmax[i])
                    {outside=1;break;}
                if(outside==1)
                    fe=1.0e12;
                else
                {
                    for(j=0;j<n;++j)pmu_per_species[index[j]]=xe[j];
                    if(current_input.mode==0)
                        fe=simplex_obj_f( smu_per_species,pmu_per_species,current_input,&compl_array);
                    else
                        fe=AOM_simplex_obj_f( smu_per_species,pmu_per_species,current_input,&res);
                }
                // evaluate function at xe
                
                if(fe<fr)
                {
                    // overwrite point with xe
                    for(i=0;i<n;++i)x[n][i]=xe[i];
                    sprintf(action,"expand");
                }
                else
                {
                    // overwrite point with xr
                    for(i=0;i<n;++i)x[n][i]=xr[i];
                    sprintf(action,"reflect");
                }
            }
            else
            {
                if((fr<f[n]&&fr>f[n-1])||(fr<f[n]&&simplex_equal(fr,f[n-1])==1))
                {
                    // outside contraction
                    for(i=0;i<n;++i)xoc[i]=xC[i]+gamma*(xr[i]-xC[i]);
                    
                    // check boundary conditions
                    outside=0;
                    for(i=0;i<n;++i)
                        if(xoc[i]<xmin[i]||xoc[i]>xmax[i])
                        {outside=1;break;}
                    if(outside==1)
                        foc=1.0e12;
                    else
                    {
                        for(j=0;j<n;++j)pmu_per_species[index[j]]=xoc[j];
                        if(current_input.mode==0)
                            foc=simplex_obj_f( smu_per_species,pmu_per_species,current_input,&compl_array);
                        else
                            foc=AOM_simplex_obj_f( smu_per_species,pmu_per_species,current_input,&res);
                    }
                    // evaluate function at xoc
                    
                    if(foc<fr||simplex_equal(foc,fr)==1)
                    {
                        // overwrite point with xoc
                        for(i=0;i<n;++i)x[n][i]=xoc[i];
                        sprintf(action,"contract outside");
                    }
                    else
                    {
                        shrink=1;
                    }
                }
                else
                {
                    // inside contraction
                    for(i=0;i<n;++i)xic[i]=xC[i]-gamma*(xr[i]-xC[i]);
                    
                    // check boundary conditions
                    outside=0;
                    for(i=0;i<n;++i)
                        if(xic[i]<xmin[i]||xic[i]>xmax[i])
                        {outside=1;break;}
                    if(outside==1)
                        fic=1.0e12;
                    else
                    {
                        for(j=0;j<n;++j)pmu_per_species[index[j]]=xic[j];
                        if(current_input.mode==0)
                            fic=simplex_obj_f( smu_per_species,pmu_per_species,current_input,&compl_array);
                        else
                            fic=AOM_simplex_obj_f( smu_per_species,pmu_per_species,current_input,&res);
                    }
                    // evaluate function at xic
                    
                    if(fic<f[n])
                    {
                        // overwrite point with xic
                        for(i=0;i<n;++i)x[n][i]=xic[i];
                        sprintf(action,"contract inside");
                    }
                    else
                    {
                        shrink=1;
                    }
                }
                // shrink
                if(shrink==1)
                {
                    for(i=1;i<n+1;++i)
                    {
                        for(j=0;j<n;++j)
                        {
                            x[i][j]=x[0][j]+sigma*(x[i][j]-x[0][j]);
                        }
                    }
                    sprintf(action,"shrink");
                }
            }
        }
        
        // order
        for(i=0;i<n+1;++i)  // evaluate function and store in array f
        {
            for(j=0;j<n;++j)array[j]=x[i][j];
                                                              // n+1 evaluations
            for(j=0;j<n;++j)pmu_per_species[index[j]]=array[j];
            if(current_input.mode==0)
                f[i]=simplex_obj_f( smu_per_species,pmu_per_species,current_input,&compl_array);
            else
                f[i]=AOM_simplex_obj_f( smu_per_species,pmu_per_species,current_input,&res);
        }
        flag=1;             // bubble sort
        while(flag==1)
        {
            flag=0;
            for(i=1;i<n+1;++i)
            {
                if(f[i]<f[i-1])
                {
                    d_buffer=f[i-1];f[i-1]=f[i];f[i]=d_buffer;
                    for(j=0;j<n;++j)
                    {
                        d_buffer=x[i-1][j];x[i-1][j]=x[i][j];x[i][j]=d_buffer;
                    }
                    flag=1;
                }
            }
        }
        
        // console out
        printf("%d\t",k+2);
        for(i=0;i<n;++i)printf("%e\t",x[0][i]);
        printf("f = %e\t",f[0]);
        for(i=0;i<n;++i)bc[k%5][i]=f[i];
        printf("%s\n",action);
        
        // breaking condition
        if(k>4)
        {
            break_sum=0;
            for(j=0;j<n;++j)
            {
                for(i=1;i<5;++i)
                {
                    if(fabs(bc[i][j]-bc[0][j])<break_acc)break_sum=break_sum+1;
                }
            }
            if(break_sum==n*4)break;
        }
        if(current_input.mode==0)
        {
            for(i=0;i<current_input.simplex_entries;++i)printf("%.6lf;",compl_array[i]);printf("\n");
        }
        
    }
    
    
    //--------------------------------------------------------------------------
    for(i=0;i<n+1;++i)free(x[i]);free(x);
    free(f);free(xr);free(xC);free(xe);free(xoc);free(xic);
    free(array);
    free(xmin);free(xmax);
    for(i=0;i<5;++i)free(bc[i]);free(bc);
    
    free(index);
    if(current_input.mode==0)
        free(compl_array);
    else
        free(res);
}

int simplex_equal(double A,double B);
int simplex_equal(double A,double B)
{
    int value;
    if(fabs(A-B)<simplex_equal_diff){value=1;}else{value=0;}
    return value;
}
