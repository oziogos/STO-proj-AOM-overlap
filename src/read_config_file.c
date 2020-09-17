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
void read_config_file(char *current_folder,char *config_file,char *mode,char *STOproj_name,char *STOproj_molecule,char *STOproj_basis,char *STOproj_MO_file,char *STOproj_MO,char *STOproj_basis_file,char *STOproj_cube_grid,double *smu_per_species,double *pmu_per_species,double *AOMsmu_per_species,double *AOMpmu_per_species,int *verb,int *cube,char *STOproj_AOM_include, int *frag1atoms,int *MOLog_channel,int *phase1,int *phase2);
void read_config_file(char *current_folder,char *config_file,char *mode,char *STOproj_name,char *STOproj_molecule,char *STOproj_basis,char *STOproj_MO_file,char *STOproj_MO,char *STOproj_basis_file,char *STOproj_cube_grid,double *smu_per_species,double *pmu_per_species,double *AOMsmu_per_species,double *AOMpmu_per_species,int *verb,int *cube,char *STOproj_AOM_include, int *frag1atoms,int *MOLog_channel,int *phase1,int *phase2)
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
    
    fclose(fp);
    
    if(*verb==-1||*verb==1){
        // console output
        printf("\n--------------------------------------------------------------------------\n\n");
    }
}
