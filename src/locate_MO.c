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
void locate_MO(char *STOproj_MO_file,char *STOproj_MO,int CGTOs,double ***MOcoeffs,int verb);
void locate_MO(char *STOproj_MO_file,char *STOproj_MO,int CGTOs,double ***MOcoeffs,int verb)
{
    FILE *fp;
    char file_path[cmax_length],**MOLog,buffer[cmax_length];
    int lines,i,j,MOcount,*MOpos,*MOcol;
    int elements_out,depth;
    char **values;
    
    if(verb==-1){
        // console output
        printf("Processing %s\n",STOproj_MO_file);
        printf("Target MO: %s\n",STOproj_MO);
    }
    // read file and transfer to memory
    sprintf(file_path,"%s",STOproj_MO_file);
    fp=fopen(file_path,"r");if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
    lines=0;while(fgets(buffer,cmax_length,fp)!=NULL)if(strcmp(buffer,"\n")!=0)++lines;rewind(fp);
    MOLog=(char**)malloc(lines*sizeof(char*));for(i=0;i<lines;++i)MOLog[i]=(char*)malloc(cmax_length*sizeof(char));
    i=-1;while(fgets(buffer,cmax_length,fp)!=NULL)if(strcmp(buffer,"\n")!=0){++i;sprintf(MOLog[i],"%s",buffer);}
    fclose(fp);
    // deduce number of channels
    MOcount=0;
    for(i=0;i<lines;++i)
    {
        sprintf(buffer,"%s",MOLog[i]);
        tokenize2(buffer,&elements_out,&depth,&values);
        if(elements_out<=4)
        {
            for(j=0;j<elements_out;++j)
            {
                if(strcmp(values[j],STOproj_MO)==0)
                {
                    ++MOcount;
                }
            }
        }
        clear_tok(elements_out,&values);
    }
    // locate MO
    MOpos=(int*)malloc(MOcount*sizeof(int));
    MOcol=(int*)malloc(MOcount*sizeof(int));
    MOcount=0;
    for(i=0;i<lines;++i)
    {
        sprintf(buffer,"%s",MOLog[i]);
        tokenize2(buffer,&elements_out,&depth,&values);
        if(elements_out<=4)
        {
            for(j=0;j<elements_out;++j)
            {
                if(strcmp(values[j],STOproj_MO)==0)
                {
                    ++MOcount;
                    MOpos[MOcount-1]=i;
                    MOcol[MOcount-1]=j;
                }
            }
        }
        clear_tok(elements_out,&values);
    }
    // store to memory
    (*MOcoeffs)=(double**)malloc(CGTOs*sizeof(double*));for(i=0;i<CGTOs;++i)(*MOcoeffs)[i]=(double*)malloc(MOcount*sizeof(double));
    for(i=0;i<MOcount;++i)
    {
        for(j=0;j<CGTOs;++j)
        {
            sprintf(buffer,"%s",MOLog[MOpos[i]+3+j]);
            tokenize2(buffer,&elements_out,&depth,&values);
            (*MOcoeffs)[j][i]=atof(values[3+MOcol[i]+1]);
            clear_tok(elements_out,&values);
        }
    }
    if(verb==-1){
        // console output
        printf("MO coefficients:\n");
        for(i=0;i<CGTOs;++i)
        {
            printf("%d\t",i+1);
            for(j=0;j<MOcount;++j)printf("%lf\t",(*MOcoeffs)[i][j]);
            printf("\n");
        }
    }
    if(verb==-1)printf("\n--------------------------------------------------------------------------\n\n");
    // free
    for(i=0;i<lines;++i)free(MOLog[i]);free(MOLog);
    free(MOpos);free(MOcol);
}
