/* -----------------------------------------------------------------------------
 
 Unified STO projection / AOM overlap code
 
 Beta version: 1.0
 1-Mar-2019
 
 Orestis George Ziogos, UCL
 o.ziogos@ucl.ac.uk
 
 For more information, examine the README file in the parent directory.
 
-------------------------------------------------------------------------------- */
#include"general.h"

void initialize_STO(int atoms,char **species,double *smu_per_species,double *pmu_per_species,int *STOs,int **STO_id_array,int **STO_type_array,double **STO_mu_array,int verb);
void resolve_pvecs_per_fragment(int atoms,double *x,double *y,double *z,char **species,double **px,double **py,double **pz,int frag1atoms);
double AOM_overlap_calculation(int istart,int istop,int jstart,int jstop,double *x,double *y,double *z,int *STO_id_array,int *STO_type_array,double *STO_mu_array,double **STO_matrix);


double AOM_simplex_obj_f(double *smu_per_species,double *pmu_per_species,struct simplex_obj_f_input current_input,double **res);
double AOM_simplex_obj_f(double *smu_per_species,double *pmu_per_species,struct simplex_obj_f_input current_input,double **res)
{
    
    int i,j,k;
    
    int atoms;
    double *x,*y,*z;
    char **species;
    
    int STOs,*STO_id_array,*STO_type_array,frag1atoms,frag1STOs;
    double *STO_mu_array,*AOM_pi,*px,*py,*pz,**STO_matrix,S;
    
    int n;
    double A0,A1,A2;//,B0,B1,B2;
    double beta2,chisq,beta1,chisq_log;
    
    n=current_input.simplex_entries;
    
    for(i=0;i<current_input.simplex_entries;++i)
    {
        /*
        printf("xyz bounds: %d\t%d\n",current_input.pos_storage_L[i],current_input.pos_storage_R[i]);
        printf("atoms: %d\n",current_input.atoms_aggr[i]);
        printf("frag1_atoms: %d\n",current_input.AOM_frag_atoms[current_input.AOM_type_map[i]]);
        for(j=current_input.pos_storage_L[i];j<=current_input.pos_storage_R[i];++j)
            printf("%s\t%lf\t%lf\t%lf\n",current_input.species_aggr[j],current_input.x_aggr[j],current_input.y_aggr[j],current_input.z_aggr[j]);
        for(j=current_input.AOM_indices[current_input.AOM_type_map[i]][0];j<=current_input.AOM_indices[current_input.AOM_type_map[i]][1];++j)
            printf("%lf\n",current_input.AOM_pi_values[j]);
        printf("-----------------------------------------------------------\n");
        */
        
        atoms=current_input.atoms_aggr[i];
        x=(double*)malloc(atoms*sizeof(double));
        y=(double*)malloc(atoms*sizeof(double));
        z=(double*)malloc(atoms*sizeof(double));
        species=(char**)malloc(atoms*sizeof(char*));
        for(j=0;j<atoms;++j)species[j]=(char*)malloc(cmax_length*sizeof(char));
        
        k=-1;
        for(j=current_input.pos_storage_L[i];j<=current_input.pos_storage_R[i];++j)
        {
            k=k+1;
            x[k]=current_input.x_aggr[j];
            y[k]=current_input.y_aggr[j];
            z[k]=current_input.z_aggr[j];
            strcpy(species[k],current_input.species_aggr[j]);
        }
        
        initialize_STO(atoms,species,smu_per_species,pmu_per_species,&STOs,&STO_id_array,&STO_type_array,&STO_mu_array,0);
        
        AOM_pi=(double*)malloc(atoms*sizeof(double));
        k=-1;
        for(j=current_input.AOM_indices[current_input.AOM_type_map[i]][0];j<=current_input.AOM_indices[current_input.AOM_type_map[i]][1];++j)
        {
            k=k+1;
            AOM_pi[k]=current_input.AOM_pi_values[j];
        }
        
        frag1atoms=current_input.AOM_frag_atoms[current_input.AOM_type_map[i]];
        
        resolve_pvecs_per_fragment(atoms,x,y,z,species,&px,&py,&pz,frag1atoms);
        
        STO_matrix=(double**)malloc(atoms*sizeof(double*));for(j=0;j<atoms;++j)STO_matrix[j]=(double*)malloc(4*sizeof(double));

        for(j=0;j<STOs;++j)
        {
            if(STO_type_array[j]==1){
                STO_matrix[STO_id_array[j]-1][0]=0.0;
                STO_matrix[STO_id_array[j]-1][1]=0.0;
                STO_matrix[STO_id_array[j]-1][2]=0.0;
                STO_matrix[STO_id_array[j]-1][3]=0.0;
            }
            else{
                k=(STO_type_array[j]-2)%4;
                if(k==0)STO_matrix[STO_id_array[j]-1][k]=0.0;
                if(k==1)STO_matrix[STO_id_array[j]-1][k]=px[STO_id_array[j]-1]*AOM_pi[STO_id_array[j]-1];
                if(k==2)STO_matrix[STO_id_array[j]-1][k]=py[STO_id_array[j]-1]*AOM_pi[STO_id_array[j]-1];
                if(k==3)STO_matrix[STO_id_array[j]-1][k]=pz[STO_id_array[j]-1]*AOM_pi[STO_id_array[j]-1];
            }
        }

        frag1STOs=0;for(j=0;j<STOs;++j)if(STO_id_array[j]<=frag1atoms)++frag1STOs;

        S=AOM_overlap_calculation(0,frag1STOs,0,frag1STOs,x,y,z,STO_id_array,STO_type_array,STO_mu_array,STO_matrix);
        
        S=sqrt(fabs(S));for(j=0;j<frag1atoms;++j)for(k=1;k<=3;++k)STO_matrix[j][k]=STO_matrix[j][k]/S;
        
        S=AOM_overlap_calculation(frag1STOs,STOs,frag1STOs,STOs,x,y,z,STO_id_array,STO_type_array,STO_mu_array,STO_matrix);
        
        S=sqrt(fabs(S));for(j=frag1atoms;j<atoms;++j)for(k=1;k<=3;++k)STO_matrix[j][k]=STO_matrix[j][k]/S;
        
        S=AOM_overlap_calculation(0,frag1STOs,frag1STOs,STOs,x,y,z,STO_id_array,STO_type_array,STO_mu_array,STO_matrix);
        
        (*res)[i]=S;
        
        free(x);free(y);free(z);for(j=0;j<atoms;++j)free(species[j]);free(species);
        free(STO_id_array);free(STO_type_array);free(STO_mu_array);
        free(AOM_pi);free(px);free(py);free(pz);
        for(j=0;j<atoms;++j)free(STO_matrix[j]);free(STO_matrix);
    }
    /*
    A0=0.0;for(i=0;i<n;++i)A0=A0+fabs(current_input.HAB[i]);A0=-2.0*A0;
    A1=2.0*n;
    A2=0.0;for(i=0;i<n;++i)A2=A2+fabs((*res)[i]);A2=2.0*A2;
    B0=0.0;for(i=0;i<n;++i)B0=B0+fabs((*res)[i])*fabs(current_input.HAB[i]);B0=-2.0*B0;
    B1=A2;
    B2=0.0;for(i=0;i<n;++i)B2=B2+fabs((*res)[i])*fabs((*res)[i]);B2=2.0*B2;
    beta1=((A2*B0)/(A1*B2)-A0/A1)/(1.0-(A2*B1)/(A1*B2));
    beta2=-B0/B2-(B1*beta1)/B2;
    chisq=0.0;for(i=0;i<n;++i)chisq=chisq+(fabs(current_input.HAB[i])-(beta1+beta2*fabs((*res)[i])))*(fabs(current_input.HAB[i])-(beta1+beta2*fabs((*res)[i])));
    */
    /*
    A0=0.0;for(i=0;i<n;++i)A0=A0+fabs((*res)[i])*fabs(current_input.HAB[i]);
    A1=0.0;for(i=0;i<n;++i)A1=A1+fabs((*res)[i])*fabs((*res)[i]);
    beta2=A0/A1;
    chisq=0.0;for(i=0;i<n;++i)chisq=chisq+(fabs(current_input.HAB[i])-(beta2*fabs((*res)[i])))*(fabs(current_input.HAB[i])-(beta2*fabs((*res)[i])));
    */
    
    if(current_input.eval_metric!=0)
    {
        int nlog=0;
        double xlog[n],ylog[n];
        for(i=0;i<n;++i)
        {
            if(fabs((*res)[i])>0.0&&fabs(current_input.HAB[i])>0.0)
            {
                xlog[nlog]=log10(fabs((*res)[i]));
                ylog[nlog]=log10(fabs(current_input.HAB[i]));
                nlog=nlog+1;
            }
        }
        
        A1=2.0*nlog;
        A0=0.0;A2=0.0;for(i=0;i<nlog;++i){A0=A0+ylog[i];A2=A2+xlog[i];}A0=-2.0*A0;A2=2.0*A2;
        beta1=-(A0+A2)/A1;
        chisq_log=0.0;for(i=0;i<nlog;++i)chisq_log=chisq_log+(ylog[i]-(beta1+xlog[i]))*(ylog[i]-(beta1+xlog[i]));
        chisq_log=chisq_log/nlog;
    }
    else
    {
        chisq_log=-1.0;
    }
    return chisq_log;
}
