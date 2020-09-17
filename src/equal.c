/* -----------------------------------------------------------------------------
 
 Unified STO projection / AOM overlap code
 
 Alpha version: 0.1
 1-Oct-2019
 
 Orestis George Ziogos, UCL
 o.ziogos@ucl.ac.uk
 
 For more information, examine the README file in the parent directory.
 
-------------------------------------------------------------------------------- */
#include"general.h"
int equal(double A,double B);
int equal(double A,double B)
{
    int value;
    if(fabs(A-B)<equal_diff){value=1;}else{value=0;}
    return value;
}
