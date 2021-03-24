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
void read_basis(char *input_cp2k_basis_file,char *input_species,char *input_basis,int verb)
{
    FILE *fp,*fpw;
    char current_folder[cmax_length],file_path[cmax_length],buffer[cmax_length],word[cmax_length],basis[cmax_length];
    char **values;
    int nset,n,lmin,lmax,nexp;
    int i,j,k,elements_out,depth,l,*nshell,m,counter=0,cols,c_counter;
    double **set_matrix;
    
    int max;
    
    getcwd(current_folder,cmax_length);
    
    // log file
    sprintf(file_path,"%s/temp_%s_%s.out",current_folder,input_species,input_basis);
    fpw=fopen(file_path,"w+");
    
    // open cp2k basis file
    sprintf(file_path,"%s",input_cp2k_basis_file);
    fp=fopen(file_path,"r");if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
    // locate species
    while(fgets(buffer,cmax_length,fp)!=NULL)
    {
        sprintf(word,"%s","");
        sscanf(buffer,"%s",word);
        if(strcmp(word,input_species)==0)
        {
            // locate basis type
            sscanf(buffer,"%s\t%s\t%s",word,word,basis);
            if(strcmp(basis,input_basis)==0)
            {
                // read the number of sets
                fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&nset);
                // loop on sets
                for(i=0;i<nset;++i)
                {
                    // read set header
                    fgets(buffer,cmax_length,fp);
                    if(verb==-1){
                        // console out
                        printf("%s",buffer);
                    }
                    tokenize2(buffer,&elements_out,&depth,&values);
                    
                    lmin=atoi(values[1]);
                    lmax=atoi(values[2]);
                    nexp=atoi(values[3]);
                    nshell=(int*)malloc((lmax-lmin+1)*sizeof(int));
                    for(j=4;j<4+(lmax-lmin+1);++j)sscanf(values[j],"%d",&nshell[j-4]);
                    cols=1;
                    for(j=0;j<(lmax-lmin+1);++j)cols=cols+nshell[j];
                    clear_tok(elements_out,&values);
                    //
                    set_matrix=(double**)malloc(nexp*sizeof(double*));for(j=0;j<nexp;++j)set_matrix[j]=(double*)malloc(cols*sizeof(double));
                    // loop on exponents and store matrix in memory
                    for(j=0;j<nexp;++j)
                    {
                        fgets(buffer,cmax_length,fp);if(verb==-1)printf("%s",buffer);
                        tokenize2(buffer,&elements_out,&depth,&values);
                        for(k=0;k<elements_out;++k)sscanf(values[k],"%lf",&set_matrix[j][k]);
                        clear_tok(elements_out,&values);
                    }
                    //
                    //for(j=0;j<nexp;++j){for(k=0;k<(lmax-lmin+1+1);++k)printf("%lf\t",set_matrix[j][k]);printf("\n");}
                    //
                    if(verb==-1)printf("lmin=%d\n",lmin);
                    if(verb==-1)printf("lmax=%d\n",lmax);
                    // loop on l
                    c_counter=0;
                    for(l=lmin;l<=lmax;++l)
                    {
                        //
                        if(verb==-1)printf("l=%d has nshell=%d contractions\n",l,nshell[l-lmin]);
                        
                        for(n=0;n<nshell[l-lmin];++n)
                        {
                            ++c_counter;
                            for(m=-l;m<=l;++m)
                            {
                                ++counter;
                                if(verb==-1)printf("gto_%s_%d = r^%d * ( ",input_species,counter,l);
                                fprintf(fpw,"%s_%d\t%d\t%d\t%d\t",input_species,counter,l,m,nexp);
                                for(k=0;k<nexp;++k)
                                {
                                    if(fabs(set_matrix[k][c_counter])>1.0e-9)
                                    {
                                        if(verb==-1)printf("%lf * exp(-%lf*r^2) + ",set_matrix[k][c_counter],set_matrix[k][0]);
                                    }
                                    fprintf(fpw,"%lf\t%lf\t",set_matrix[k][c_counter],set_matrix[k][0]);
                                }
                                
                                if(verb==-1)printf(" ) * Y(%d,%d)\n",l,m);
                                fprintf(fpw,"\n");
                            }
                        }
                    }
                    
                    
                    //
                    for(j=0;j<nexp;++j)free(set_matrix[j]);free(set_matrix);
                    free(nshell);
                    if(verb==-1)printf("\n");
                }
            }
        }
    }
    fclose(fp);
    
    fclose(fpw);
    
    // restructure log file
    sprintf(file_path,"%s/temp_%s_%s.out",current_folder,input_species,input_basis);
    fp=fopen(file_path,"r");
    sprintf(file_path,"%s/log_%s_%s.dat",current_folder,input_species,input_basis);
    fpw=fopen(file_path,"w+");
    max=1;
    while(fgets(buffer,cmax_length,fp)!=NULL)
    {
        tokenize2(buffer,&elements_out,&depth,&values);
        //for(i=0;i<elements_out;++i)printf("%s\n",values[i]);
        //printf("%s\n",values[3]);
        if(atoi(values[3])>max)max=atoi(values[3]);
        clear_tok(elements_out,&values);
    }
    rewind(fp);
    while(fgets(buffer,cmax_length,fp)!=NULL)
    {
        tokenize2(buffer,&elements_out,&depth,&values);
        for(i=0;i<elements_out;++i)fprintf(fpw,"%s\t",values[i]);
        for(i=0;i<max-atoi(values[3]);++i)fprintf(fpw,"%lf\t%lf\t",0.0,0.0);
        fprintf(fpw,"\n");
        clear_tok(elements_out,&values);
    }
    fclose(fp);fclose(fpw);
    
    sprintf(buffer,"rm temp_%s_%s.out",input_species,input_basis);
    system(buffer);
}
