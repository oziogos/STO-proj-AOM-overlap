/* -----------------------------------------------------------------------------
 
 Unified STO projection / AOM overlap code
 
 Beta version: 1.0
 1-Mar-2019
 
 Orestis George Ziogos, UCL
 o.ziogos@ucl.ac.uk
 
 For more information, examine the README file in the parent directory.
 
-------------------------------------------------------------------------------- */
#include"general.h"
int resolve_atomic_Z(char *species);
void read_config_file(char *current_folder,char *config_file,char *mode,char *STOproj_name,char *STOproj_molecule,char *STOproj_basis,char *STOproj_MO_file,char *STOproj_MO,char *STOproj_basis_file,char *STOproj_cube_grid,double *smu_per_species,double *pmu_per_species,double *AOMsmu_per_species,double *AOMpmu_per_species,int *verb,int *cube,char *STOproj_AOM_include, int *frag1atoms,int *MOLog_channel,int *phase1,int *phase2,
    int *simplex_entries,char ***simplex_name,char ***simplex_molecule,char ***simplex_MO_file,char ***simplex_MO,char ***simplex_cube,char ***simplex_type,struct simplex_params *simplex,
    char ***simplex_dimer,char ***simplex_AOM_include,int **simplex_frag1_atoms,double **simplex_HAB);
void read_config_file(char *current_folder,char *config_file,char *mode,char *STOproj_name,char *STOproj_molecule,char *STOproj_basis,char *STOproj_MO_file,char *STOproj_MO,char *STOproj_basis_file,char *STOproj_cube_grid,double *smu_per_species,double *pmu_per_species,double *AOMsmu_per_species,double *AOMpmu_per_species,int *verb,int *cube,char *STOproj_AOM_include, int *frag1atoms,int *MOLog_channel,int *phase1,int *phase2,
    int *simplex_entries,char ***simplex_name,char ***simplex_molecule,char ***simplex_MO_file,char ***simplex_MO,char ***simplex_cube,char ***simplex_type,struct simplex_params *simplex,
    char ***simplex_dimer,char ***simplex_AOM_include,int **simplex_frag1_atoms,double **simplex_HAB)
{
    FILE *fp;
    char file_path[cmax_length],buffer[cmax_length],word[cmax_length],local_species[sub_length],verbosity[cmax_length],word2[cmax_length];
    int i;
    
    sprintf(file_path,"%s/%s",current_folder,config_file);
    fp=fopen(file_path,"r");if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
    
    while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s"," ");sscanf(buffer,"%s",word);if(strcmp(word,"verb")==0){sscanf(buffer,"%s\t%s",word,verbosity);break;}}rewind(fp);
    *verb=1;
    if(strcmp(verbosity,"none")==0)*verb=0;
    if(strcmp(verbosity,"default")==0)*verb=1;
    if(strcmp(verbosity,"debug")==0)*verb=-1;
    if(*verb==-1||*verb==1)printf("%s",buffer);
    
    while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s"," ");sscanf(buffer,"%s",word);if(strcmp(word,"mode")==0){sscanf(buffer,"%s\t%s",word,mode);if(*verb==-1||*verb==1)printf("%s",buffer);break;}}rewind(fp);
    
    // resolve STO mu per element for projection
    if(strcmp(mode,"molecule")==0)
    {
    
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s"," ");sscanf(buffer,"%s",word);if(strcmp(word,"name")==0){sscanf(buffer,"%s\t%s",word,STOproj_name);if(*verb==-1||*verb==1)printf("%s",buffer);break;}}rewind(fp);
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s"," ");sscanf(buffer,"%s",word);if(strcmp(word,"molecule")==0){sscanf(buffer,"%s\t%s",word,STOproj_molecule);if(*verb==-1||*verb==1)printf("%s",buffer);break;}}rewind(fp);
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s"," ");sscanf(buffer,"%s",word);if(strcmp(word,"basis")==0){sscanf(buffer,"%s\t%s",word,STOproj_basis);if(*verb==-1||*verb==1)printf("%s",buffer);break;}}rewind(fp);
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s"," ");sscanf(buffer,"%s",word);if(strcmp(word,"MO_file")==0){sscanf(buffer,"%s\t%s",word,STOproj_MO_file);if(*verb==-1||*verb==1)printf("%s",buffer);break;}}rewind(fp);
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s"," ");sscanf(buffer,"%s",word);if(strcmp(word,"MO")==0){sscanf(buffer,"%s\t%s",word,STOproj_MO);if(*verb==-1||*verb==1)printf("%s",buffer);break;}}rewind(fp);
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s"," ");sscanf(buffer,"%s",word);if(strcmp(word,"basis_file")==0){sscanf(buffer,"%s\t%s",word,STOproj_basis_file);if(*verb==-1||*verb==1)printf("%s",buffer);break;}}rewind(fp);
        sprintf(STOproj_cube_grid,"%s","0.25");
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s"," ");sscanf(buffer,"%s",word);if(strcmp(word,"cube_grid")==0){sscanf(buffer,"%s\t%s",word,STOproj_cube_grid);if(*verb==-1||*verb==1)printf("%s",buffer);break;}}rewind(fp);
        
        *cube=0;
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s"," ");sscanf(buffer,"%s",word);if(strcmp(word,"cube")==0){sscanf(buffer,"%s\t%s",word,word2);if(strcmp(word2,"yes")==0)*cube=1;if(*verb==-1||*verb==1)printf("%s",buffer);break;}}rewind(fp);
        
        *MOLog_channel=0;
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s"," ");sscanf(buffer,"%s",word);if(strcmp(word,"MO_channel")==0){sscanf(buffer,"%s\t%d",word,&(*MOLog_channel));*MOLog_channel=*MOLog_channel-1;if(*verb==-1||*verb==1)printf("%s",buffer);break;}}rewind(fp);

        for(i=0;i<16;++i){smu_per_species[i]=0.0;pmu_per_species[i]=0.0;}
        while(fgets(buffer,cmax_length,fp)!=NULL)
        {
            sprintf(word,"%s"," ");sscanf(buffer,"%s",word);
            if(strcmp(word,"proj_mu")==0)
            {
                sscanf(buffer,"%s\t%s",word,local_species);
                if(*verb==-1||*verb==1)printf("%s",buffer);
                i=resolve_atomic_Z(local_species);
                if(i==1){sscanf(buffer,"%s\t%s\t%lf",word,local_species,&smu_per_species[i-1]);}
                else{sscanf(buffer,"%s\t%s\t%lf\t%lf",word,local_species,&smu_per_species[i-1],&pmu_per_species[i-1]);}
            }
        }
        rewind(fp);
    }
    
    // resolve STO mu per element for AOM overlap calculations
    else if(strcmp(mode,"dimer")==0)
    {
        
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s"," ");sscanf(buffer,"%s",word);if(strcmp(word,"dimer")==0){sscanf(buffer,"%s\t%s",word,STOproj_molecule);if(*verb==-1||*verb==1)printf("%s",buffer);break;}}rewind(fp);
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s"," ");sscanf(buffer,"%s",word);if(strcmp(word,"AOM_include")==0){sscanf(buffer,"%s\t%s",word,STOproj_AOM_include);if(*verb==-1||*verb==1)printf("%s",buffer);break;}}rewind(fp);
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s"," ");sscanf(buffer,"%s",word);if(strcmp(word,"atoms_frag1")==0){sscanf(buffer,"%s\t%d",word,&(*frag1atoms));if(*verb==-1||*verb==1)printf("%s",buffer);break;}}rewind(fp);
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s"," ");sscanf(buffer,"%s",word);if(strcmp(word,"phase")==0){sscanf(buffer,"%s\t%d\t%d",word,&(*phase1),&(*phase2));if(*verb==-1||*verb==1)printf("%s",buffer);break;}}rewind(fp);

        for(i=0;i<16;++i){AOMsmu_per_species[i]=0.0;AOMpmu_per_species[i]=0.0;}
        while(fgets(buffer,cmax_length,fp)!=NULL)
        {
            sprintf(word,"%s"," ");sscanf(buffer,"%s",word);
            if(strcmp(word,"AOM_mu")==0)
            {
                sscanf(buffer,"%s\t%s",word,local_species);
                if(*verb==-1||*verb==1)printf("%s",buffer);
                i=resolve_atomic_Z(local_species);
                if(i==1){sscanf(buffer,"%s\t%s\t%lf",word,local_species,&AOMsmu_per_species[i-1]);}
                else{sscanf(buffer,"%s\t%s\t%lf\t%lf",word,local_species,&AOMsmu_per_species[i-1],&AOMpmu_per_species[i-1]);}
            }
        }
    }
    
    // simplex: single mol
    else if(strcmp(mode,"molecule_simplex")==0)
    {
        
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s"," ");sscanf(buffer,"%s",word);if(strcmp(word,"basis")==0){sscanf(buffer,"%s\t%s",word,STOproj_basis);if(*verb==-1||*verb==1)printf("%s",buffer);break;}}rewind(fp);
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s"," ");sscanf(buffer,"%s",word);if(strcmp(word,"basis_file")==0){sscanf(buffer,"%s\t%s",word,STOproj_basis_file);if(*verb==-1||*verb==1)printf("%s",buffer);break;}}rewind(fp);
        sprintf(STOproj_cube_grid,"%s","0.25");
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s"," ");sscanf(buffer,"%s",word);if(strcmp(word,"cube_grid")==0){sscanf(buffer,"%s\t%s",word,STOproj_cube_grid);if(*verb==-1||*verb==1)printf("%s",buffer);break;}}rewind(fp);
        *MOLog_channel=0;
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s"," ");sscanf(buffer,"%s",word);if(strcmp(word,"MO_channel")==0){sscanf(buffer,"%s\t%d",word,&(*MOLog_channel));*MOLog_channel=*MOLog_channel-1;if(*verb==-1||*verb==1)printf("%s",buffer);break;}}rewind(fp);
        
        
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s"," ");sscanf(buffer,"%s",word);if(strcmp(word,"entries")==0){sscanf(buffer,"%s\t%d",word,&(*simplex_entries));if(*verb==-1||*verb==1)printf("%s",buffer);break;}}
        
        
        *simplex_name=(char**)malloc(*simplex_entries*sizeof(char**));
        for(i=0;i<*simplex_entries;++i)(*simplex_name)[i]=(char*)malloc(cmax_length*sizeof(char*));
        *simplex_molecule=(char**)malloc(*simplex_entries*sizeof(char**));
        for(i=0;i<*simplex_entries;++i)(*simplex_molecule)[i]=(char*)malloc(cmax_length*sizeof(char*));
        *simplex_MO_file=(char**)malloc(*simplex_entries*sizeof(char**));
        for(i=0;i<*simplex_entries;++i)(*simplex_MO_file)[i]=(char*)malloc(cmax_length*sizeof(char*));
        *simplex_MO=(char**)malloc(*simplex_entries*sizeof(char**));
        for(i=0;i<*simplex_entries;++i)(*simplex_MO)[i]=(char*)malloc(cmax_length*sizeof(char*));
        *simplex_cube=(char**)malloc(*simplex_entries*sizeof(char**));
        for(i=0;i<*simplex_entries;++i)(*simplex_cube)[i]=(char*)malloc(cmax_length*sizeof(char*));
        *simplex_type=(char**)malloc(*simplex_entries*sizeof(char**));
        for(i=0;i<*simplex_entries;++i)(*simplex_type)[i]=(char*)malloc(cmax_length*sizeof(char*));
        
        for(i=0;i<*simplex_entries;++i)
        {
            fgets(buffer,cmax_length,fp);
            printf("%s",buffer);
            sscanf(buffer,"%s\t%s\t%s\t%s\t%s\t%s",(*simplex_name)[i],(*simplex_molecule)[i],(*simplex_MO_file)[i],(*simplex_MO)[i],(*simplex_cube)[i],(*simplex_type)[i]);
        }
        
        rewind(fp);
        
        // load STO mu coeffs
        for(i=0;i<16;++i){smu_per_species[i]=0.0;pmu_per_species[i]=0.0;}
        while(fgets(buffer,cmax_length,fp)!=NULL)
        {
            sprintf(word,"%s"," ");sscanf(buffer,"%s",word);
            if(strcmp(word,"proj_mu")==0)
            {
                sscanf(buffer,"%s\t%s",word,local_species);
                if(*verb==-1||*verb==1)printf("%s",buffer);
                i=resolve_atomic_Z(local_species);
                if(i==1){sscanf(buffer,"%s\t%s\t%lf",word,local_species,&smu_per_species[i-1]);}
                else{sscanf(buffer,"%s\t%s\t%lf\t%lf",word,local_species,&smu_per_species[i-1],&pmu_per_species[i-1]);}
            }
        }
        rewind(fp);
        
        // load simplex bounds
        double min,max;
        for(i=0;i<16;++i){simplex->p_mu_min[i]=0.0;simplex->p_mu_max[i]=0.0;}
        while(fgets(buffer,cmax_length,fp)!=NULL)
        {
            sprintf(word,"%s"," ");sscanf(buffer,"%s",word);
            if(strcmp(word,"bounds")==0)
            {
                sscanf(buffer,"%s\t%s\t%s\t%lf\t%lf",word,word2,local_species,&min,&max);
                i=resolve_atomic_Z(local_species);
                if(i!=1 && strcmp(word2,"p_mu")==0)
                {
                    if(*verb==-1||*verb==1)printf("%s",buffer);
                    simplex->p_mu_min[i-1]=min;
                    simplex->p_mu_max[i-1]=max;
                }
            }
        }
        rewind(fp);
        
        // load simplex parameters
        /*
        simplex_steps	1000
        simplex_acc	1.0e-16
        simplex_rho	1.0
        simplex_xi	2.0
        simplex_gamma	0.5
        simplex_sigma	0.5
        simplex_tau	0.05
        simplex_tau_prime	0.00025
        */
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s"," ");sscanf(buffer,"%s",word);if(strcmp(word,"simplex_steps")==0){sscanf(buffer,"%s\t%d",word,&(simplex->steps));if(*verb==-1||*verb==1)printf("%s",buffer);break;}}rewind(fp);
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s"," ");sscanf(buffer,"%s",word);if(strcmp(word,"simplex_acc")==0){sscanf(buffer,"%s\t%lf",word,&(simplex->acc));if(*verb==-1||*verb==1)printf("%s",buffer);break;}}rewind(fp);
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s"," ");sscanf(buffer,"%s",word);if(strcmp(word,"simplex_rho")==0){sscanf(buffer,"%s\t%lf",word,&(simplex->rho));if(*verb==-1||*verb==1)printf("%s",buffer);break;}}rewind(fp);
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s"," ");sscanf(buffer,"%s",word);if(strcmp(word,"simplex_xi")==0){sscanf(buffer,"%s\t%lf",word,&(simplex->xi));if(*verb==-1||*verb==1)printf("%s",buffer);break;}}rewind(fp);
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s"," ");sscanf(buffer,"%s",word);if(strcmp(word,"simplex_gamma")==0){sscanf(buffer,"%s\t%lf",word,&(simplex->gamma));if(*verb==-1||*verb==1)printf("%s",buffer);break;}}rewind(fp);
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s"," ");sscanf(buffer,"%s",word);if(strcmp(word,"simplex_sigma")==0){sscanf(buffer,"%s\t%lf",word,&(simplex->sigma));if(*verb==-1||*verb==1)printf("%s",buffer);break;}}rewind(fp);
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s"," ");sscanf(buffer,"%s",word);if(strcmp(word,"simplex_tau")==0){sscanf(buffer,"%s\t%lf",word,&(simplex->tau));if(*verb==-1||*verb==1)printf("%s",buffer);break;}}rewind(fp);
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s"," ");sscanf(buffer,"%s",word);if(strcmp(word,"simplex_tau_prime")==0){sscanf(buffer,"%s\t%lf",word,&(simplex->tau_prime));if(*verb==-1||*verb==1)printf("%s",buffer);break;}}rewind(fp);
    
        if(simplex->steps<0)
        {
            for(i=0;i<16;++i){simplex->p_mu_sample_step[i]=0.0;simplex->p_mu_sample_max[i]=0.0;}
            while(fgets(buffer,cmax_length,fp)!=NULL)
            {
                sprintf(word,"%s"," ");sscanf(buffer,"%s",word);
                if(strcmp(word,"sample")==0)
                {
                    sscanf(buffer,"%s\t%s\t%s\t%lf\t%lf",word,word2,local_species,&min,&max);
                    i=resolve_atomic_Z(local_species);
                    if(i!=1 && strcmp(word2,"p_mu")==0)
                    {
                        if(*verb==-1||*verb==1)printf("%s",buffer);
                        simplex->p_mu_sample_step[i-1]=min;
                        simplex->p_mu_sample_max[i-1]=max;
                    }
                }
            }
            rewind(fp);
        }
    
    }
    
    else if(strcmp(mode,"dimer_simplex")==0)
    {
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s"," ");sscanf(buffer,"%s",word);if(strcmp(word,"entries")==0){sscanf(buffer,"%s\t%d",word,&(*simplex_entries));if(*verb==-1||*verb==1)printf("%s",buffer);break;}}
        
        *simplex_name=(char**)malloc(*simplex_entries*sizeof(char**));
        for(i=0;i<*simplex_entries;++i)(*simplex_name)[i]=(char*)malloc(cmax_length*sizeof(char*));
        
        *simplex_dimer=(char**)malloc(*simplex_entries*sizeof(char**));
        for(i=0;i<*simplex_entries;++i)(*simplex_dimer)[i]=(char*)malloc(cmax_length*sizeof(char*));
        
        *simplex_AOM_include=(char**)malloc(*simplex_entries*sizeof(char**));
        for(i=0;i<*simplex_entries;++i)(*simplex_AOM_include)[i]=(char*)malloc(cmax_length*sizeof(char*));
        
        *simplex_frag1_atoms=(int*)malloc(*simplex_entries*sizeof(int));
        
        *simplex_HAB=(double*)malloc(*simplex_entries*sizeof(double));
        
        simplex->eval_metric=1;
        
        for(i=0;i<*simplex_entries;++i)
        {
            fgets(buffer,cmax_length,fp);
            //sscanf(buffer,"%s\t%s\t%s\t%d\t%lf",(*simplex_name)[i],(*simplex_dimer)[i],(*simplex_AOM_include)[i],&(*simplex_frag1_atoms)[i],&(*simplex_HAB)[i]);
            sscanf(buffer,"%s\t%s\t%s\t%d\t%s",(*simplex_name)[i],(*simplex_dimer)[i],(*simplex_AOM_include)[i],&(*simplex_frag1_atoms)[i],word);
            if(strcmp(word,"N/A")==0)
            {
                simplex->eval_metric=0;
                (*simplex_HAB)[i]=0.0;
            }
            else
            {
                printf("%s",buffer);
                sscanf(word,"%lf",&(*simplex_HAB)[i]);
            }
        }
        
        rewind(fp);
        
        for(i=0;i<16;++i){AOMsmu_per_species[i]=0.0;AOMpmu_per_species[i]=0.0;}
        while(fgets(buffer,cmax_length,fp)!=NULL)
        {
            sprintf(word,"%s"," ");sscanf(buffer,"%s",word);
            if(strcmp(word,"AOM_mu")==0)
            {
                sscanf(buffer,"%s\t%s",word,local_species);
                if(*verb==-1||*verb==1)printf("%s",buffer);
                i=resolve_atomic_Z(local_species);
                if(i==1){sscanf(buffer,"%s\t%s\t%lf",word,local_species,&AOMsmu_per_species[i-1]);}
                else{sscanf(buffer,"%s\t%s\t%lf\t%lf",word,local_species,&AOMsmu_per_species[i-1],&AOMpmu_per_species[i-1]);}
            }
        }
        
        rewind(fp);

        // load simplex bounds
        double min,max;
        for(i=0;i<16;++i){simplex->p_mu_min[i]=0.0;simplex->p_mu_max[i]=0.0;}
        while(fgets(buffer,cmax_length,fp)!=NULL)
        {
            sprintf(word,"%s"," ");sscanf(buffer,"%s",word);
            if(strcmp(word,"bounds")==0)
            {
                sscanf(buffer,"%s\t%s\t%s\t%lf\t%lf",word,word2,local_species,&min,&max);
                i=resolve_atomic_Z(local_species);
                if(i!=1 && strcmp(word2,"p_mu")==0)
                {
                    if(*verb==-1||*verb==1)printf("%s",buffer);
                    simplex->p_mu_min[i-1]=min;
                    simplex->p_mu_max[i-1]=max;
                }
            }
        }
        rewind(fp);
        
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s"," ");sscanf(buffer,"%s",word);if(strcmp(word,"simplex_steps")==0){sscanf(buffer,"%s\t%d",word,&(simplex->steps));if(*verb==-1||*verb==1)printf("%s",buffer);break;}}rewind(fp);
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s"," ");sscanf(buffer,"%s",word);if(strcmp(word,"simplex_acc")==0){sscanf(buffer,"%s\t%lf",word,&(simplex->acc));if(*verb==-1||*verb==1)printf("%s",buffer);break;}}rewind(fp);
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s"," ");sscanf(buffer,"%s",word);if(strcmp(word,"simplex_rho")==0){sscanf(buffer,"%s\t%lf",word,&(simplex->rho));if(*verb==-1||*verb==1)printf("%s",buffer);break;}}rewind(fp);
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s"," ");sscanf(buffer,"%s",word);if(strcmp(word,"simplex_xi")==0){sscanf(buffer,"%s\t%lf",word,&(simplex->xi));if(*verb==-1||*verb==1)printf("%s",buffer);break;}}rewind(fp);
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s"," ");sscanf(buffer,"%s",word);if(strcmp(word,"simplex_gamma")==0){sscanf(buffer,"%s\t%lf",word,&(simplex->gamma));if(*verb==-1||*verb==1)printf("%s",buffer);break;}}rewind(fp);
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s"," ");sscanf(buffer,"%s",word);if(strcmp(word,"simplex_sigma")==0){sscanf(buffer,"%s\t%lf",word,&(simplex->sigma));if(*verb==-1||*verb==1)printf("%s",buffer);break;}}rewind(fp);
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s"," ");sscanf(buffer,"%s",word);if(strcmp(word,"simplex_tau")==0){sscanf(buffer,"%s\t%lf",word,&(simplex->tau));if(*verb==-1||*verb==1)printf("%s",buffer);break;}}rewind(fp);
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s"," ");sscanf(buffer,"%s",word);if(strcmp(word,"simplex_tau_prime")==0){sscanf(buffer,"%s\t%lf",word,&(simplex->tau_prime));if(*verb==-1||*verb==1)printf("%s",buffer);break;}}rewind(fp);
        
        if(simplex->steps<0)
        {
            for(i=0;i<16;++i){simplex->p_mu_sample_step[i]=0.0;simplex->p_mu_sample_max[i]=0.0;}
            while(fgets(buffer,cmax_length,fp)!=NULL)
            {
                sprintf(word,"%s"," ");sscanf(buffer,"%s",word);
                if(strcmp(word,"sample")==0)
                {
                    sscanf(buffer,"%s\t%s\t%s\t%lf\t%lf",word,word2,local_species,&min,&max);
                    i=resolve_atomic_Z(local_species);
                    if(i!=1 && strcmp(word2,"p_mu")==0)
                    {
                        if(*verb==-1||*verb==1)printf("%s",buffer);
                        simplex->p_mu_sample_step[i-1]=min;
                        simplex->p_mu_sample_max[i-1]=max;
                    }
                }
            }
            rewind(fp);
        }
        
    }
    
    fclose(fp);
    
    if(*verb==-1||*verb==1){
        // console output
        printf("\n--------------------------------------------------------------------------\n\n");
    }
}
