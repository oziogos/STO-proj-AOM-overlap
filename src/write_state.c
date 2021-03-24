/* -----------------------------------------------------------------------------
 
 Unified STO projection / AOM overlap code
 
 Beta version: 1.0
 1-Mar-2019
 
 Orestis George Ziogos, UCL
 o.ziogos@ucl.ac.uk
 
 For more information, examine the README file in the parent directory.
 
-------------------------------------------------------------------------------- */
#include"general.h"

void write_state(char *current_folder,char *STOproj_name,char *STOproj_MO,
                 int atoms,int STOs,double *px,double *py,double *pz,double **Smatrix,double **STO_matrix,
                 double orb_compl,double *V_array,double *s_array);
void write_state(char *current_folder,char *STOproj_name,char *STOproj_MO,
                 int atoms,int STOs,double *px,double *py,double *pz,double **Smatrix,double **STO_matrix,
                 double orb_compl,double *V_array,double *s_array)
{
    FILE *fp;
    int i,j;
    char file_path[cmax_length];
    
    // write state
    sprintf(file_path,"%s/%s_state.dat",current_folder,STOproj_name);
    fp=fopen(file_path,"w");
    
    fprintf(fp,"{");
                    
    fprintf(fp,"'pvecs': ");
    fprintf(fp,"{");
    fprintf(fp,"'px': ");
    fprintf(fp,"[");
    for(i=0;i<atoms-1;++i)fprintf(fp,"%.16e, ",px[i]);fprintf(fp,"%.16e",px[atoms-1]);
    fprintf(fp,"],");
    fprintf(fp,"'py': ");
    fprintf(fp,"[");
    for(i=0;i<atoms-1;++i)fprintf(fp,"%.16e, ",py[i]);fprintf(fp,"%.16e",py[atoms-1]);
    fprintf(fp,"],");
    fprintf(fp,"'pz': ");
    fprintf(fp,"[");
    for(i=0;i<atoms-1;++i)fprintf(fp,"%.16e, ",pz[i]);fprintf(fp,"%.16e",pz[atoms-1]);
    fprintf(fp,"]");
    fprintf(fp,"}, ");
    
    fprintf(fp,"'S_matrix': ");
    fprintf(fp,"[");
    for(i=0;i<STOs-1;++i){
        fprintf(fp,"[");
        for(j=0;j<STOs-1;++j){fprintf(fp,"%.16e, ",Smatrix[i][j]);}fprintf(fp,"%.16e",Smatrix[i][STOs-1]);
        fprintf(fp,"], ");
    }
    i=STOs-1;
    fprintf(fp,"[");
    for(j=0;j<STOs-1;++j){fprintf(fp,"%.16e, ",Smatrix[i][j]);}fprintf(fp,"%.16e",Smatrix[i][STOs-1]);
    fprintf(fp,"]");
    fprintf(fp,"], ");
    
    fprintf(fp,"'AOM_dict': ");
    fprintf(fp,"{%s: ",STOproj_MO);
    fprintf(fp,"[");
    for(i=0;i<atoms-1;++i)fprintf(fp,"%.16e, ",px[i]*STO_matrix[i][1]+py[i]*STO_matrix[i][2]+pz[i]*STO_matrix[i][3]);
    i=atoms-1;
    fprintf(fp,"%.16e",px[i]*STO_matrix[i][1]+py[i]*STO_matrix[i][2]+pz[i]*STO_matrix[i][3]);
    fprintf(fp,"]");
    fprintf(fp,"}, ");
    
    fprintf(fp,"'V_array': ");
    fprintf(fp,"[");
    for(i=0;i<STOs-1;++i)fprintf(fp,"%.16e, ",V_array[i]);fprintf(fp,"%.16e",V_array[STOs-1]);
    fprintf(fp,"], ");
    
    fprintf(fp,"'compl_dict': ");
    fprintf(fp,"{%s: %.16e",STOproj_MO,orb_compl);
    fprintf(fp,"}, ");
    
    fprintf(fp,"'STO_matrix': ");
    fprintf(fp,"[");
    for(i=0;i<atoms-1;++i)fprintf(fp,"[%.16e, %.16e, %.16e, %.16e], ",STO_matrix[i][0],STO_matrix[i][1],STO_matrix[i][2],STO_matrix[i][3]);
    i=atoms-1;
    fprintf(fp,"[%.16e, %.16e, %.16e, %.16e]",STO_matrix[i][0],STO_matrix[i][1],STO_matrix[i][2],STO_matrix[i][3]);
    fprintf(fp,"], ");
    
    fprintf(fp,"'singular_values': ");
    fprintf(fp,"[");
    for(i=0;i<STOs-1;++i)fprintf(fp,"%.16e, ",s_array[i]);fprintf(fp,"%.16e",s_array[STOs-1]);
    fprintf(fp,"]");
    
    fprintf(fp,"}\n");
    
    fclose(fp);
    
}

