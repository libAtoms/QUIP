/*************************************************/
/* Atoms: -llapack -lblas -lm                    */
/*        -lVecMat3 -lVecMat -lScalar -lIO       */
/*                                               */
/* Physical constants, macros, and configuration */
/* operators for atomistic simulations.          */
/*                                               */
/* Dec. 12, 1999  Ju Li <liju99@mit.edu>         */
/*************************************************/

#include "Atoms.h"

/*****************/
/* Atomic Strain */
/*****************/

/* A unique integer that depends on the chemical symbol sequence */
unsigned int ConfigChecksum (Alib_Declare_Config)
{
    unsigned int sum=0;
    unsigned short r=55665, c1=52845, c2=22719;
    char cipher;
    register int i;
    for (i=2*(*np); i--;)
    {
        cipher = (*symbol)[i] ^ (r >> 8);
        r = (cipher + r) * c1 + c2;
        sum += cipher;
    }
    return (sum);
} /* end ConfigChecksum() */


/* Save a reference configuration. Remember to free it some time! */
void IsoAtomicReferenceReImprint (Alib_Declare_Config, IsoAtomicReference *ref)
{
    ref->checksum = ConfigChecksum(Config_Alib_to_Alib);
    M3EQV (H, ref->H);
    VRECLONE (DIMENSION*(*np), *s, ref->s);
    return;
} /* end IsoAtomicReferenceReImprint() */


/* Free the reference configuration */
void IsoAtomicReferenceFree (IsoAtomicReference *ref)
{
    Free (ref->s);
} /* end IsoAtomicReferenceFree() */


/* Compute and save the least-square transformation */
/* matrix  dxij(new) \approx dxij(ref) J.           */
void ComputeLeastSquareDeformationGradient
(IsoAtomicReference *ref, Alib_Declare_Config, Neighborlist *N, M3 **J)
{
    int i, j, k, kmin, kmax;
    M3 *V=NULL, *W=NULL, tmp;
    V3 ds, d0, d;

    if (ConfigChecksum(Config_Alib_to_Alib) != ref->checksum)
        pe ("ComputeLeastSquareDeformationGradient: the configuration\n"
            "is not isoatomic with the imprinted reference.\n");

    REALLOC_ZERO (ComputeLeastSquareDeformationGradient, V,  *np, M3);
    REALLOC_ZERO (ComputeLeastSquareDeformationGradient, W,  *np, M3);
    REALLOC      (ComputeLeastSquareDeformationGradient, *J, *np, M3);
    
    for (i=*np; i--;)
    {
        if ( (N->min_shrinkage == 1) && (N->max_atom_displacement == 0) )
        { /* DISPOSABLE Neighborlist */
            kmin = N->idx[i];
            kmax = N->idx[i+1];
        }
        else
        { /* REUSABLE Neighborlist */
            kmin = N->idx[2*i];
            kmax = N->idx[2*i+1];
        }
        for (k=kmin; k<kmax; k++)
        {
            j = N->list[k];
            V3SUB (&ref->s[DIMENSION*j], &ref->s[DIMENSION*i], ds);
            V3ImagE (ds);
            V3mM3 (ds, ref->H, d0);
            M3ASSIGNV3V3 (d0, d0, tmp);
            M3AdD (tmp, V[i]);
            M3AdD (tmp, V[j]);
            
            V3SUB (&(*s)[DIMENSION*j], &(*s)[DIMENSION*i], ds);
            V3ImagE (ds);
            V3mM3 (ds, H, d);
            M3ASSIGNV3V3 (d0, d, tmp);
            M3AdD(tmp, W[i]);
            M3AdD(tmp, W[j]);
        }
    }
    
    for (i=*np; i--;)
    {
        if (M3VOLUME(V[i]) != 0)
        {
            M3Inv (V[i]);
            M3MUL (V[i], W[i], (*J)[i]);
            M3_SET_ZERO_IF_TINY((*J)[i]);
        }
        else M3IDENTITY((*J)[i]);
    }
    
    Free(W);
    Free(V);
    return;
} /* end ComputeLeastSquareDeformationGradient() */


#ifdef _annotate_atomic_strain
int main (int argc, char *argv[])
{
    Aapp_Define_Config;
    Chemtab ct[1]={{0}};
    Tp *tp=NULL;
    Neighborlist N[1]={{0}};
    char *ref_fname, *cur_fname, *sav_fname, list[3]={'x','y','z'};
    IsoAtomicReference ref[1]={0};
    M3 *J=NULL, eta;
    int i, j, k;
    FILE *fp;

    if (argc == 4)
    {
        ref_fname = argv[1];
        cur_fname = argv[2];
        sav_fname = argv[3];
    }
    else if (argc == 3)
    {
        ref_fname = argv[1];
        cur_fname = argv[2];
        sav_fname = cur_fname;
    }
    else
    {
        printf ("\nPurpose: compute atomic strain and save "
                "as auxiliary properties.\nhttp://alum.mit.edu/www/liju99/"
                "Graphics/A/annotate_atomic_strain/Doc/main.pdf\n\n");
        printf ("Usage: %s ref_cfg cur_cfg (sav_cfg=cur_cfg)\n"
                "       %s ref_cfg cur_cfg sav_cfg\n\n",
                argv[0], argv[0]);
        return (1);
    }

    printf ("Loading reference configuration \"%s\"...\n\n", ref_fname);
    CONFIG_LOAD (ref_fname, Config_Aapp_to_Alib);
    
    IsoAtomicReferenceReImprint (Config_Aapp_to_Alib, ref);

    printf ("Loading current configuration \"%s\"...\n\n", cur_fname);
    CONFIG_LOAD (cur_fname, Config_Aapp_to_Alib);

    rebind_CT (Config_Aapp_to_Alib, "", ct, &tp); cr();
    Neighborlist_Recreate_Form (Config_Aapp_to_Alib, ct, N);
    /* N->max_atom_displacement = 0.01; */
    N->s_overflow_err_handler =
        NEIGHBORLIST_S_OVERFLOW_ERR_HANDLER_FOLD_INTO_PBC;
    Neighborlist_Recreate (Config_Aapp_to_Alib, stdout, ct, &tp, N);

    ComputeLeastSquareDeformationGradient (ref, Config_Aapp_to_Alib, N, &J);
    IsoAtomicReferenceFree (ref);
        
    fp = WOpen(sav_fname);
    fprintf (fp, "Number of particles = %d\n", np);
    fprintf (fp, "H0(1,1) = %.15g A\n", H[0][0]);
    fprintf (fp, "H0(1,2) = %.15g A\n", H[0][1]);
    fprintf (fp, "H0(1,3) = %.15g A\n", H[0][2]);
    fprintf (fp, "H0(2,1) = %.15g A\n", H[1][0]);
    fprintf (fp, "H0(2,2) = %.15g A\n", H[1][1]);
    fprintf (fp, "H0(2,3) = %.15g A\n", H[1][2]);
    fprintf (fp, "H0(3,1) = %.15g A\n", H[2][0]);
    fprintf (fp, "H0(3,2) = %.15g A\n", H[2][1]);
    fprintf (fp, "H0(3,3) = %.15g A\n", H[2][2]);
    fprintf (fp, "entry_count = %d\n", 6+CONFIG_num_auxiliary+11);
    for (k=0; k<CONFIG_num_auxiliary; k++)
        fprintf (fp, "auxiliary[%d] = %s [%s]\n", 
                 k, CONFIG_auxiliary_name[k], CONFIG_auxiliary_unit[k]);

    fprintf (fp, "auxiliary[%d] = eta_Mises [reduced unit]\n",
             CONFIG_num_auxiliary);
    fprintf (fp, "auxiliary[%d] = eta_hydro [reduced unit]\n",
             CONFIG_num_auxiliary+1);
    for (j=0; j<DIMENSION; j++)
        for (k=0; k<DIMENSION; k++)
            fprintf (fp, "auxiliary[%d] = J%c%c [reduced unit]\n", 
                     CONFIG_num_auxiliary+2+3*j+k, list[j], list[k]);
    
    for (i=0; i<np; i++)
    {
        fprintf (fp, "%.5g\n", mass[i] * UMASS_IN_AMU);
        fprintf (fp, "%2s\n",  SYM(i));
        for (j=0; j<DIMENSION; j++) fprintf (fp, " %.10g", s[3*i+j] );
        for (j=0; j<DIMENSION; j++)
            fprintf (fp, " %.10g", s1[3*i+j] / UTIME_IN_NS);
        for (k=0; k<CONFIG_num_auxiliary; k++)
            fprintf (fp, " %.10g", CONFIG_auxiliary[k][i]);
        
        M3MULT (J[i], eta);
        M3SubdiaG (eta, 1);
        M3DividE (eta, 2);
        
        fprintf (fp, " %.5g", SymmetricM3MisesInvariant(eta));
        fprintf (fp, " %.5g", SymmetricM3HydroInvariant(eta));
        for (j=0; j<DIMENSION; j++)
            for (k=0; k<DIMENSION; k++)
                fprintf (fp, " %.5g", J[i][j][k]);

        fcr(fp);
    }
    Zclose (fp, sav_fname);
    printf ("\n\"%s\" / \"%s\" -> \"%s\".\n", cur_fname, ref_fname, sav_fname);

    Free (J);
    return (0);
}
#endif /* _annotate_atomic_strain */
