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
int resolve_atomic_Z(char *species)
{
    int res;
    if(strcmp(species,"H")==0)res=1;
    if(strcmp(species,"C")==0)res=6;
    if(strcmp(species,"N")==0)res=7;
    if(strcmp(species,"O")==0)res=8;
    if(strcmp(species,"F")==0)res=9;
    if(strcmp(species,"S")==0)res=16;
    return res;
}
