/* -----------------------------------------------------------------------------
 
 Unified STO projection / AOM overlap code
 
 Beta version: 1.0
 1-Mar-2019
 
 Orestis George Ziogos, UCL
 o.ziogos@ucl.ac.uk
 
 For more information, examine the README file in the parent directory.
 
-------------------------------------------------------------------------------- */
#include"general.h"
void resolve_unique_species(int atoms,char **species,int *atom_types,char ***unique_species,int verb);
void resolve_unique_species(int atoms,char **species,int *atom_types,char ***unique_species,int verb)
{
    int i,j,*species_map;
    species_map=(int*)malloc(atoms*sizeof(int));for(i=0;i<atoms;++i)species_map[i]=1;
    for(i=0;i<atoms-1;++i){for(j=i+1;j<atoms;++j){if(strcmp(species[i],species[j])==0)species_map[j]=0;}}
    *atom_types=0;for(i=0;i<atoms;++i)if(species_map[i]==1)++(*atom_types);
    *unique_species=(char**)malloc(*atom_types*sizeof(char*));for(i=0;i<*atom_types;++i)(*unique_species)[i]=(char*)malloc(sub_length*sizeof(char));
    j=-1;for(i=0;i<atoms;++i)if(species_map[i]==1){++j;sprintf((*unique_species)[j],"%s",species[i]);}
    free(species_map);
    if(verb==-1){
        // console output
        printf("Unique species: %d\n",*atom_types);for(i=0;i<*atom_types;++i)printf("%s\t",(*unique_species)[i]);printf("\n");
    }
}
