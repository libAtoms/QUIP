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

/*************************************************/
/* Build Grains by Voronoi Site-Rotation and Cut */
/*************************************************/

/************************************************************/
/* Build n0 x n1 x n2 random rotation Voronoi site lattice. */
/* If st == NULL, then a random position lattice is built.  */
/************************************************************/
VoronoiSpec *VoronoiSpec_rebuild
(int n0, int n1, int n2, Xtal_Abstract *st, VoronoiSpec *V)
{
    register int i;
    Aapp_Define_Config;
    V3 a,b;
    if (st)
    {
        if (st->nkd != 1)
            pe ("VoronoiSpec_rebuild:\n"
                "please, at least use a single-component lattice.\n");
        Config_rebuild_Xtal
            (Config_Aapp_to_Alib, n0, n1, n2, st, "Cu", -1.,
             -1., -1., -1.);
        V->nv = np;
        REALLOC ( VoronoiSpec_rebuild, V->v, V->nv, VoronoiSite );
        for (i=V->nv; i--;)
            V3EQV ( &(s[DIMENSION*i]), V->v[i].s );
        Config_free ( Config_Aapp_to_Alib );
        /* In 1D, 3 points in PBC ensures nearest image rule; otherwise: */
        V->cut_against_all_images = (n0 < 3) || (n1 < 3) || (n2 < 3);
    }
    else  /* st == NULL */
    {
        V->nv = n0 * n1 * n2;
        REALLOC ( VoronoiSpec_rebuild, V->v, V->nv, VoronoiSite );
        for (i=V->nv; i--;)
            V3Frandom ( V->v[i].s );
        V->cut_against_all_images = TRUE;
    }
    for (i=V->nv; i--;)
    {
        V3randomunit (a);
        V3randomunit (b);
        M3geodesic (a, b, V->v[i].R);
    }
    return (V);
} /* end VoronoiSpec_rebuild() */


/* Cut against "ss[]" and all images */
int Config_Voronoi_CUT (V3 s0, int nss, double *ss, Alib_Declare_Config)
{
    int j, k, l;
    double ds0[DIMENSION], ds[DIMENSION], dx[DIMENSION+1], *sm, *SM;
    register int i;
    register double *sn;
    char *taglist;
    MALLOC ( Config_Voronoi_CUT, sm, nss*DIMENSION3*DIMENSION, double );
    SM = sm;
    for (l=nss; l--;)
    {
        V3SUB ( &(ss[DIMENSION*l]), s0, ds0 );
        for (i=-1; i<2; i++)
            for (j=-1; j<2; j++)
                for (k=-1; k<2; k++)
                {
                    ds[0] = ds0[0] + i;
                    ds[1] = ds0[1] + j;
                    ds[2] = ds0[2] + k;
                    if ( V3EQZERO(ds) ) continue;
                    V3mM3( ds, H, dx );
                    dx[DIMENSION] = V3LENGTH2( dx );
                    M3mV3 (H, dx, SM);
                    V3MuL (2./dx[DIMENSION], SM);
                    SM += DIMENSION;
                }
    }
    CALLOC ( Config_Voronoi_CUT, taglist, *np, char );
    for (i=(*np); i--;)
    {
        V3SUB ( &((*s)[DIMENSION*i]), s0, ds );
        for (sn=sm; sn<SM; sn+=DIMENSION)
            if (V3DOT(ds,sn) > 1)
            {
                taglist[i] = 1;
                break;
            }
    }
    Free( sm );
    return(Config_rid_of_TAGLIST(taglist,Config_Alib_to_Alib,i));
} /* end Config_Voronoi_CUT() */


/* Cut only against the nearest image of "ss[]" */
int Config_Voronoi_cut (V3 s0, int nss, double *ss, Alib_Declare_Config)
{
    int l;
    double ds[DIMENSION], dx[DIMENSION+1], *sm, *SM;
    register int i;
    register double *sn;
    char *taglist;
    MALLOC ( Config_Voronoi_cut, sm, nss*DIMENSION, double );
    SM = sm;
    for (l=nss; l--;)
    {
        V3SUB ( &(ss[DIMENSION*l]), s0, ds );
        V3ImagE ( ds );
        if ( V3EQZERO(ds) ) continue;
        V3mM3( ds, H, dx );
        dx[DIMENSION] = V3LENGTH2( dx );
        M3mV3 (H, dx, SM);
        V3MuL (2./dx[DIMENSION], SM);
        SM += DIMENSION;
    }
    CALLOC ( Config_Voronoi_cut, taglist, *np, char );
    for (i=(*np); i--;)
    {
        V3SUB ( &((*s)[DIMENSION*i]), s0, ds );
        for (sn=sm; sn<SM; sn+=DIMENSION)
            if (V3DOT(ds,sn) > 1)
            {
                taglist[i] = 1;
                break;
            }
    }
    Free( sm );
    return(Config_rid_of_TAGLIST(taglist,Config_Alib_to_Alib,i));
} /* end Config_Voronoi_cut() */


/* Generate grains based on specification "V" and the input configuration */
int Config_Voronoirize (VoronoiSpec *V, Alib_Declare_Config, FILE *info)
{
    int i, j, k, m, npold;
    double *ss;
    ConfigStack cs[1], ca[1], cb[1];
    Chemtab ct[1]={{0}};
    Tp *tp=NULL;
    Neighborlist N[1]={{0}};
    char *taglist;
    npold = *np;
    CONFIG_PUSH (cs);
    MALLOC ( Config_Voronoirize, ss, V->nv*DIMENSION, double );
    for (m=0; m<V->nv; m++) V3EQV ( V->v[m].s, &(ss[DIMENSION*m]) );
    Fprintf ( info, "Cut against all images = %s;\n\n",
              BWORD(V->cut_against_all_images) );
    for (m=0; m<V->nv; m++)
    {
        Config_clone (cs, Config_Alib_to_Alib);
        Config_CENTERFOLD (Config_Alib_to_Alib, V->v[m].s);
        Config_rotate_VIA_s (Config_Alib_to_Alib, V->v[m].s, V->v[m].R);
        if ( V->cut_against_all_images )
            Config_Voronoi_CUT  (V->v[m].s, V->nv, ss, Config_Alib_to_Alib);
        else Config_Voronoi_cut (V->v[m].s, V->nv, ss, Config_Alib_to_Alib);
        if (m == 0) Config_push (Config_Alib_to_Alib, cb);
        else
        {
            Config_push (Config_Alib_to_Alib, ca);
            Config_pop (cb, Config_Alib_to_Alib);
            Config_cat (ca, Config_Alib_to_Alib);
            Config_push (Config_Alib_to_Alib, cb);
        }
        Fprintf (info, "Grain %d (0-%d) centered at (%g %g %g) added,\n",
                 m, V->nv-1, V3E(V->v[m].s));
        S3fPR (info, "R = %M.\n ", V->v[m].R);
    }
    Free (ss);
    CONFIG_erase (cs);
    Config_fold_into_PBC (Config_Alib_to_Alib);
    Fprintf (info, "There are total of %d grains, %d atoms (%g%%).\n\n",
             V->nv, *np, 100.*(*np)/npold);
    if (V->min_bond_ratio == 0)
        V->min_bond_ratio = VORONOISPEC_DEF_MIN_BOND_RATIO;
    if (V->min_bond_ratio > 0)
    {
        CALLOC ( Config_Voronoirize, taglist, *np, char );
        rebind_ct (Config_Alib_to_Alib, "", ct, &tp, info); Fcr(info);
        Neighborlist_Recreate_Form (Config_Alib_to_Alib, ct, N);
        /* modify the default cutoffs */
        for (i=0; i<ct->t; i++)
            for (j=0; j<ct->t; j++)
                NEIGHBOR_TABLE(N->rcut,ct,i,j) *=
                    V->min_bond_ratio / NEIGHBORLIST_RCUT_RATIO;
        Neighborlist_Recreate (Config_Alib_to_Alib, info, ct, &tp, N);
        for (i=0; i<*np; i++)
            if (!taglist[i])
                for (k=N->idx[i]; k<N->idx[i+1]; k++)
                {
                    j = N->list[k];
                    if (taglist[j]) continue;
                    if (HALF_DECISION())
                    {
                        taglist[i] = 1;
                        break;
                    }
                    else taglist[j] = 1;
                }
        Neighborlist_Free (N);
        Free (tp);
        Config_rid_of_TAGLIST (taglist, Config_Alib_to_Alib, i);
        Fprintf (info, "Then, %d (%g%%) GB atoms are removed randomly\n"
                 "for bond ratio < %g -> now np = %d (%g%%).\n\n",
                 i, 100.*i/((*np)+i), V->min_bond_ratio,
                 *np, 100.*(*np)/npold);
    }
    return(1);
} /* end Config_Voronoirize() */


#ifdef _voronoirize
int main (int argc, char *argv[])
{
    int nc[DIMENSION];
    Aapp_Define_Config;
    VoronoiSpec V[1] = {{0}};
    if (argc != 7)
    {
        printf ("Purpose: build grains by Voronoi site-rotation and cut.\n");
        printf ("Usage: %s in_file <fcc,bcc,random> n0 n1 n2 out_file\n",
                argv[0]);
        return (1);
    }
    TimeRandomize();
    printf ("Loading \"%s\"...\n\n", argv[1]);
    Config_Load (argv[1], stdout, Config_Aapp_to_Alib);
    nc[0] = atoi(argv[3]);
    nc[1] = atoi(argv[4]);
    nc[2] = atoi(argv[5]);
    if (!strcmp(argv[2], "bcc")) BCC_VoronoiSpec ( nc, V );
    else if (!strcmp(argv[2], "fcc")) FCC_VoronoiSpec ( nc, V );
    else Random_VoronoiSpec ( nc[0]*nc[1]*nc[2], V );
    Config_Voronoirize (V, Config_Aapp_to_Alib, stdout);
    Config_save (Config_Aapp_to_Alib, TRUE, argv[6]);
    printf ("%s %d x %d x %d Voronoi superlattice -> \"%s\".\n",
            argv[2], V3E(nc), argv[6]);
    return (0);
}
#endif /* _voronoirize */
