/***********************************************/
/* Atoms: -llapack -lblas -lm                  */
/*        -lVecMat3 -lVecMat -lScalar -lIO     */
/*                                             */
/* Frequently used Physical Constants, Macros, */
/* and Subroutines for Atomistic Simulation.   */
/*                                             */
/* Dec.12, 1999  Ju Li <liju99@mit.edu>        */
/***********************************************/

/* template for evaluating constants */

#include "NIST.h"
/* #define ULENGTH_IN_M  CM_IN_M                      */
/* #define UMASS_IN_KG   G_IN_KG                      */
/* #define UENERGY_IN_J  G_IN_KG * CM_IN_M * CM_IN_M  */
#include "Atoms.h"

#ifdef CHECK_THIS_CONSTANT
int main (int argc, char *argv[])
{
    printf ("CHECK_THIS_CONSTANT = %lf\n", CHECK_THIS_CONSTANT);
    return (0);
}
#endif
