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
void tokenize2(char *buffer,int *elements_out,int *depth,char ***res)
{
    //
    char second_buffer[cmax_length];
    int j,k,l;
    int end;
    int elements;
    
    int i,*positions;
    
    // loop over buffer and tokenize spaces and tabs
    end=-1;
    for(j=0;j<cmax_length;++j)
    {
        end=end+1;
        if(buffer[j]=='\0')break;
        if(buffer[j]==' ' || buffer[j]=='\t' || buffer[j]=='\n')buffer[j]='$';
    }
    // erase multiple dollars
    for(j=0;j<end-1;++j)    // end-1 because we check pairs j,j+1
    {
        if(buffer[j]=='$' && buffer[j+1]=='$')  // if you find two consecutive $s
        {
            l=-1;                               // index for writing to second_buffer
            for(k=0;k<=j;++k)                   // copy the left segment (including the first $) to memory (in second_buffer)
            {
                l=l+1;
                second_buffer[l]=buffer[k];
            }
            for(k=j+2;k<=end;++k)               // copy the remaining segment, without the second $
            {
                l=l+1;
                second_buffer[l]=buffer[k];
            }
            sprintf(buffer,"%s",second_buffer); // overwrite buffer
            j=j-1;                              // retrack to account for multiple $s
            end=end-1;                          // shrink dimension by one position
        }
    }
    // strip leading token
    if(buffer[0]=='$')
    {
        sprintf(second_buffer,"%s",buffer);
        for(j=0;j<cmax_length;++j)
        {
            if(second_buffer[j]=='\0')break;
            buffer[j]=second_buffer[j+1];
        }
    }
    // count the number of final tokens
    elements=0;
    for(j=0;j<cmax_length;++j)
    {
        if(buffer[j]=='\0')break;
        if(buffer[j]=='$')elements=elements+1;
    }
    // replace tokens
    /*
     for(j=0;j<cmax_length;++j)
     {
     if(buffer[j]=='\0')break;
     if(buffer[j]=='$')buffer[j]='\t';
     }
     */
    
    // resolve depth
    positions=(int*)malloc(elements*sizeof(int));
    k=-1;
    for(j=0;j<cmax_length;++j)
    {
        if(buffer[j]=='\0')break;
        if(buffer[j]=='$'){++k;positions[k]=j;}
    }
    
    *depth=positions[0];for(j=1;j<elements;++j)if(positions[j]-positions[j-1]>*depth)*depth=positions[j]-positions[j-1];
    *depth=*depth+1;
    // allocate results array
    *res=(char**)malloc(elements*sizeof(char*));for(j=0;j<elements;++j)(*res)[j]=(char*)malloc(*depth*sizeof(char));
    // populate results array
    i=-1; // column index
    k=0;  // row index
    for(j=0;j<cmax_length;++j)
    {
        if(buffer[j]=='\0')break;
        if(buffer[j]!='$'){++i;(*res)[k][i]=buffer[j];}
        if(buffer[j]=='$'){(*res)[k][i+1]='\0';i=-1;++k;}
    }
    
    //for(i=0;i<elements;++i)printf("%s\n",(*res)[i]);
    
    
    *elements_out=elements;
    
    free(positions);
    
    //
}
void clear_tok(int elements,char ***res);
void clear_tok(int elements,char ***res)
{
    int i;
    for(i=0;i<elements;++i)free((*res)[i]);free(*res);
}
