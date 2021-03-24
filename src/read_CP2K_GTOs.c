/* -----------------------------------------------------------------------------
 
 Unified STO projection / AOM overlap code
 
 Beta version: 1.0
 1-Mar-2019
 
 Orestis George Ziogos, UCL
 o.ziogos@ucl.ac.uk
 
 For more information, examine the README file in the parent directory.
 
-------------------------------------------------------------------------------- */
#include"general.h"
void tokenize2(char *buffer,int *elements_out,int *depth,char ***res);
void clear_tok(int elements,char ***res);
void read_basis(char *input_cp2k_basis_file,char *input_species,char *input_basis,int verb);
void read_CP2K_GTOs(char *current_folder,int atoms,int atom_types,char **species,char **unique_species,char *STOproj_basis_file,char *STOproj_basis,int *CGTOs,int **bfnPerAtom,int **GTO_depth,int *PCGTOs,double **palpha,double **pcoeff,int ***pqn,int verb);
void read_CP2K_GTOs(char *current_folder,int atoms,int atom_types,char **species,char **unique_species,char *STOproj_basis_file,char *STOproj_basis,int *CGTOs,int **bfnPerAtom,int **GTO_depth,int *PCGTOs,double **palpha,double **pcoeff,int ***pqn,int verb)
{
    FILE *fp;
    char file_path[cmax_length],buffer[cmax_length],word[cmax_length];
    int GTOs,i,j,cgtos_counter,pcgtos_counter;
    char **GTO_data;
    int elements_out,depth;
    char **values;
    
    if(verb==-1){
        // console output
        printf("Processing %s according to unique species information...\n",STOproj_basis_file);
        printf("Utilized Gaussian basis set: %s\n\n",STOproj_basis);
    }
    for(i=0;i<atom_types;++i){
        // console output
        if(verb==-1)printf("*** Species: %s\n*** Basis set breakdown:\n\n",unique_species[i]);
        read_basis(STOproj_basis_file,unique_species[i],STOproj_basis,verb);
        if(verb==-1)printf("*** Generated log file for %s_%s\n",unique_species[i],STOproj_basis);
        if(verb==-1)printf("\n--------------------------------------------------------------------------\n\n");
    }
    
    // generate GTO orbital information
    (*CGTOs)=0;GTOs=0;
    (*bfnPerAtom)=(int*)malloc(atoms*sizeof(int));for(i=0;i<atoms;++i)(*bfnPerAtom)[i]=0;
    for(i=0;i<atoms;++i)
    {
        sprintf(file_path,"%s/log_%s_%s.dat",current_folder,species[i],STOproj_basis);
        fp=fopen(file_path,"r");if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
        while(fgets(buffer,cmax_length,fp)!=NULL){++(*bfnPerAtom)[i];++GTOs;sscanf(buffer,"%s\t%s\t%d",word,word,&j);if(j==-2)++(*bfnPerAtom)[i];}
        fclose(fp);
        (*CGTOs)=(*CGTOs)+(*bfnPerAtom)[i];
    }
    GTO_data=(char**)malloc(GTOs*sizeof(char*));for(i=0;i<GTOs;++i)GTO_data[i]=(char*)malloc(cmax_length*sizeof(char));
    j=-1;
    for(i=0;i<atoms;++i)
    {
        sprintf(file_path,"%s/log_%s_%s.dat",current_folder,species[i],STOproj_basis);
        fp=fopen(file_path,"r");if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
        while(fgets(buffer,cmax_length,fp)!=NULL){++j;sprintf(GTO_data[j],"%s",buffer);}
        fclose(fp);
    }
    (*GTO_depth)=(int*)malloc((*CGTOs)*sizeof(int));
    cgtos_counter=-1;
    for(i=0;i<GTOs;++i)
    {
        sprintf(buffer,"%s",GTO_data[i]);
        tokenize2(buffer,&elements_out,&depth,&values);
        if(atoi(values[2])==0)
        {
            //s-type
            if(atoi(values[1])==0){++cgtos_counter;(*GTO_depth)[cgtos_counter]=atoi(values[3]);}
            //p-type
            if(atoi(values[1])==1)for(j=0;j<3;++j){++cgtos_counter;(*GTO_depth)[cgtos_counter]=atoi(values[3]);}
            //d-type
            if(atoi(values[1])==2)for(j=0;j<6;++j){++cgtos_counter;(*GTO_depth)[cgtos_counter]=atoi(values[3]);}
        }
        clear_tok(elements_out,&values);
    }
    (*PCGTOs)=0;for(i=0;i<(*CGTOs);++i)(*PCGTOs)=(*PCGTOs)+(*GTO_depth)[i];
    (*palpha)=(double*)malloc((*PCGTOs)*sizeof(double));
    (*pcoeff)=(double*)malloc((*PCGTOs)*sizeof(double));
    (*pqn)=(int**)malloc((*PCGTOs)*sizeof(int*));for(i=0;i<(*PCGTOs);++i)(*pqn)[i]=(int*)malloc(3*sizeof(int));
    pcgtos_counter=-1;
    for(i=0;i<GTOs;++i)
    {
        sprintf(buffer,"%s",GTO_data[i]);
        tokenize2(buffer,&elements_out,&depth,&values);
        if(atoi(values[2])==0)
        {
            //s-type
            if(atoi(values[1])==0){
                for(j=1;j<=atoi(values[3]);++j)
                {
                    ++pcgtos_counter;
                    (*pcoeff)[pcgtos_counter]=atof(values[4+1+(j-1)*2-1]);
                    (*palpha)[pcgtos_counter]=atof(values[4+2+(j-1)*2-1]);
                    (*pqn)[pcgtos_counter][0]=0;
                    (*pqn)[pcgtos_counter][1]=0;
                    (*pqn)[pcgtos_counter][2]=0;
                }
            }
            //p-type
            if(atoi(values[1])==1){
                for(j=1;j<=atoi(values[3]);++j)
                {
                    ++pcgtos_counter;
                    (*pcoeff)[pcgtos_counter]=atof(values[4+1+(j-1)*2-1]);
                    (*palpha)[pcgtos_counter]=atof(values[4+2+(j-1)*2-1]);
                    (*pqn)[pcgtos_counter][0]=1;
                    (*pqn)[pcgtos_counter][1]=0;
                    (*pqn)[pcgtos_counter][2]=0;
                }
                for(j=1;j<=atoi(values[3]);++j)
                {
                    ++pcgtos_counter;
                    (*pcoeff)[pcgtos_counter]=atof(values[4+1+(j-1)*2-1]);
                    (*palpha)[pcgtos_counter]=atof(values[4+2+(j-1)*2-1]);
                    (*pqn)[pcgtos_counter][0]=0;
                    (*pqn)[pcgtos_counter][1]=1;
                    (*pqn)[pcgtos_counter][2]=0;
                }
                for(j=1;j<=atoi(values[3]);++j)
                {
                    ++pcgtos_counter;
                    (*pcoeff)[pcgtos_counter]=atof(values[4+1+(j-1)*2-1]);
                    (*palpha)[pcgtos_counter]=atof(values[4+2+(j-1)*2-1]);
                    (*pqn)[pcgtos_counter][0]=0;
                    (*pqn)[pcgtos_counter][1]=0;
                    (*pqn)[pcgtos_counter][2]=1;
                }
            }
            //d-type
            if(atoi(values[1])==2){
                for(j=1;j<=atoi(values[3]);++j)
                {
                    ++pcgtos_counter;
                    (*pcoeff)[pcgtos_counter]=atof(values[4+1+(j-1)*2-1]);
                    (*palpha)[pcgtos_counter]=atof(values[4+2+(j-1)*2-1]);
                    (*pqn)[pcgtos_counter][0]=2;
                    (*pqn)[pcgtos_counter][1]=0;
                    (*pqn)[pcgtos_counter][2]=0;
                }
                for(j=1;j<=atoi(values[3]);++j)
                {
                    ++pcgtos_counter;
                    (*pcoeff)[pcgtos_counter]=atof(values[4+1+(j-1)*2-1]);
                    (*palpha)[pcgtos_counter]=atof(values[4+2+(j-1)*2-1]);
                    (*pqn)[pcgtos_counter][0]=1;
                    (*pqn)[pcgtos_counter][1]=1;
                    (*pqn)[pcgtos_counter][2]=0;
                }
                for(j=1;j<=atoi(values[3]);++j)
                {
                    ++pcgtos_counter;
                    (*pcoeff)[pcgtos_counter]=atof(values[4+1+(j-1)*2-1]);
                    (*palpha)[pcgtos_counter]=atof(values[4+2+(j-1)*2-1]);
                    (*pqn)[pcgtos_counter][0]=1;
                    (*pqn)[pcgtos_counter][1]=0;
                    (*pqn)[pcgtos_counter][2]=1;
                }
                for(j=1;j<=atoi(values[3]);++j)
                {
                    ++pcgtos_counter;
                    (*pcoeff)[pcgtos_counter]=atof(values[4+1+(j-1)*2-1]);
                    (*palpha)[pcgtos_counter]=atof(values[4+2+(j-1)*2-1]);
                    (*pqn)[pcgtos_counter][0]=0;
                    (*pqn)[pcgtos_counter][1]=2;
                    (*pqn)[pcgtos_counter][2]=0;
                }
                for(j=1;j<=atoi(values[3]);++j)
                {
                    ++pcgtos_counter;
                    (*pcoeff)[pcgtos_counter]=atof(values[4+1+(j-1)*2-1]);
                    (*palpha)[pcgtos_counter]=atof(values[4+2+(j-1)*2-1]);
                    (*pqn)[pcgtos_counter][0]=0;
                    (*pqn)[pcgtos_counter][1]=1;
                    (*pqn)[pcgtos_counter][2]=1;
                }
                for(j=1;j<=atoi(values[3]);++j)
                {
                    ++pcgtos_counter;
                    (*pcoeff)[pcgtos_counter]=atof(values[4+1+(j-1)*2-1]);
                    (*palpha)[pcgtos_counter]=atof(values[4+2+(j-1)*2-1]);
                    (*pqn)[pcgtos_counter][0]=0;
                    (*pqn)[pcgtos_counter][1]=0;
                    (*pqn)[pcgtos_counter][2]=2;
                }
            }
        }
        clear_tok(elements_out,&values);
    }
    
    for(i=0;i<GTOs;++i)free(GTO_data[i]);free(GTO_data);
}
