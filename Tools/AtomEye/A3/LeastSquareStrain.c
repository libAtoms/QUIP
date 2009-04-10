/*******************************************/
/* Threaded Atomistic Configuration Viewer */
/*                                         */
/* Fri Mar 31 2000 Ju Li <liju99@mit.edu>  */
/*******************************************/

#include "A.h"

int ComputeLeastSquareStrain = 0;
IsoAtomicReference ref[1]={0};
char ref_fbasename[MAX_FILENAME_SIZE]={0};


/* Free the auxiliary properties based on least-square strain */
void LeastSquareStrain_Free()
{
    int k;
    if (CONFIG_num_auxiliary >= 11)
    {
        CONFIG_num_auxiliary -= 11;
        for (k=CONFIG_num_auxiliary; k<CONFIG_num_auxiliary+11; k++)
        {
            Free(CONFIG_auxiliary[k]);
            CONFIG_auxiliary_name[k][0] = EOS;
            CONFIG_auxiliary_unit[k][0] = EOS;
        }
    }
    return;
} /* end LeastSquareStrain_Free() */


/* Append least-square strain as auxiliary properties */
void LeastSquareStrain_Append()
{
    M3 *J=NULL, eta;
    int i, j, k;
    char list[3]={'x','y','z'};

    ComputeLeastSquareDeformationGradient (ref, Config_Aapp_to_Alib, N, &J);

    for (k=CONFIG_num_auxiliary; k<CONFIG_num_auxiliary+11; k++)
    {
        REALLOC(LeastSquareStrain_Append, CONFIG_auxiliary[k], np, double);
    }

    strcpy(CONFIG_auxiliary_name[CONFIG_num_auxiliary], "eta_Mises");
    strcpy(CONFIG_auxiliary_unit[CONFIG_num_auxiliary], "reduced unit");
    strcpy(CONFIG_auxiliary_name[CONFIG_num_auxiliary+1], "eta_hydro");
    strcpy(CONFIG_auxiliary_unit[CONFIG_num_auxiliary+1], "reduced unit");
    for (j=0; j<DIMENSION; j++)
        for (k=0; k<DIMENSION; k++)
        {
            sprintf (CONFIG_auxiliary_name[CONFIG_num_auxiliary+2+3*j+k],
                     "J%c%c", list[j], list[k]);
            strcpy (CONFIG_auxiliary_unit[CONFIG_num_auxiliary+2+3*j+k],
                    "reduced unit"); 
        }

    for (i=np; i--;)
    {
        M3MULT (J[i], eta);
        M3SubdiaG (eta, 1);
        M3DividE (eta, 2);
        CONFIG_auxiliary[CONFIG_num_auxiliary][i] =
            SymmetricM3MisesInvariant(eta);
        SET_ZERO_IF_TINY(CONFIG_auxiliary[CONFIG_num_auxiliary][i]);
        CONFIG_auxiliary[CONFIG_num_auxiliary+1][i] =
            SymmetricM3HydroInvariant(eta);
        SET_ZERO_IF_TINY(CONFIG_auxiliary[CONFIG_num_auxiliary+1][i]);
        for (j=0; j<DIMENSION; j++)
            for (k=0; k<DIMENSION; k++)
                CONFIG_auxiliary[CONFIG_num_auxiliary+2+3*j+k][i] = J[i][j][k];
    }
    CONFIG_num_auxiliary += 11;

    Free (J);
    return;
} /* end LeastSquareStrain_Append() */
