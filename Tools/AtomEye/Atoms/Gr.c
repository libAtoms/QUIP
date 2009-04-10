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

/************************************************************/
/* Compute g(r)'s for system with multiple chemical species */
/************************************************************/

/* for estimating nearest neighbor distances */
#define ATOM_RADIUS_IN_A(Z)  ATOM_EMPIRICAL_RADIUS_IN_A(Z)

/* Recreate form space and fill in defaults */
void Gr_Recreate_Form (Alib_Declare_Config, Chemtab *ct, Gr *GR)
{
    int i,j;
    REALLOC( Gr_Recreate_Form, GR->rcut, SQUARE(ct->t), double );
    REALLOC( Gr_Recreate_Form, GR->mesh, SQUARE(ct->t), int );
    REALLOC( Gr_Recreate_Form, GR->accuracy, SQUARE(ct->t), double );
    for (i=0; i<DIMENSION; i++) GR->pbc[i] = GR_DEF_PBC;
    for (i=0; i<ct->t; i++)
        for (j=i; j<ct->t; j++)
        { /* against possible first-neighbor distances */
            GR_TABLE(GR->rcut,ct,i,j) = GR_DEF_RCUT_RATIO *
                ( ATOM_RADIUS_IN_A(ct->Z[i]) +
                  ATOM_RADIUS_IN_A(ct->Z[j]) ) / ulength_IN_A;
            GR_TABLE(GR->mesh,ct,i,j) = GR_DEF_MESH;
            GR_TABLE(GR->accuracy,ct,i,j) = GR_DEF_ACCURACY;
        }
    GR->fn_matlab = GR_DEF_FN_MATLAB;
    GR->autosaves = GR_DEF_AUTOSAVES;
    return;
} /* end Gr_Recreate_Form() */


/* Reallocate and clear counters: return the necessary number  */
/* of INDEPENDENT instances to achieve the desired accuracies. */
int Gr_Reset (Alib_Declare_Config, Chemtab *ct, Gr *GR, FILE *info)
{
    double volume, cc;
    int i, j, k, necessary;

    REALLOC( Gr_Reset, GR->rcut2, SQUARE(ct->t), double );

    k = 0;
    for (i=0; i<ct->t; i++)
        for (j=i; j<ct->t; j++)
        { /* user should fill in (i,j), the upper triangular part */
            GR_TABLE(GR->rcut,ct,j,i) = GR_TABLE(GR->rcut,ct,i,j);
            GR_TABLE(GR->rcut2,ct,j,i) = GR_TABLE(GR->rcut2,ct,i,j)
                = SQUARE( GR_TABLE(GR->rcut,ct,i,j) );
            GR_TABLE(GR->mesh,ct,j,i) = GR_TABLE(GR->mesh,ct,i,j);
            k += GR_TABLE(GR->mesh,ct,i,j);
            GR_TABLE(GR->accuracy,ct,j,i) = GR_TABLE(GR->accuracy,ct,i,j);
        }

    REALLOC( Gr_Reset, GR->workspace, k, double );
    VZERO (k, GR->workspace);
    REALLOC( Gr_Reset, GR->count, SQUARE(ct->t), double * );
    REALLOC( Gr_Reset, GR->normalization, SQUARE(ct->t), double );

    k = 0;
    for (i=0; i<ct->t; i++)
        for (j=i; j<ct->t; j++)
        {
            GR_TABLE(GR->count,ct,j,i) = GR_TABLE(GR->count,ct,i,j)
                = GR->workspace + k;
            GR_TABLE(GR->normalization,ct,i,j) = 0;
            k += GR_TABLE(GR->mesh,ct,i,j);
        }
    
    /* Print out the form */
    Fprintf(info, "Prepare to compute g(r)'s with the "
            "following parameters:\n");
    for (i=0; i<DIMENSION; i++)
        Fprintf(info, "pbc[%d] = %s,  ", i, BWORD(GR->pbc[i]));
    Fcr(info);
    if (info) Sfpr(info, "GRcut = %M (reduced),", GR->rcut, ct->t,ct->t);
    if (info) Ifpr(info, "mesh  = %M,", GR->mesh, ct->t,ct->t);
    if (info) Sfpr(info, "accuracy = %M.", GR->accuracy, ct->t,ct->t);

    volume = M3VOLUME(H);
    necessary = 0;
    for (i=0; i<ct->t; i++)
        for (j=i; j<ct->t; j++)
        { /* how many hits we get in a mesh bin at one instance */
            if (i == j)
                cc = ct->count[i] / 2. * (ct->count[i] - 1) / volume *
                    SPHERE_VOLUME(GR_TABLE(GR->rcut,ct,i,j)) /
                    GR_TABLE(GR->mesh,ct,i,j);
            else
                cc = ct->count[i] * ct->count[j] / volume *
                    SPHERE_VOLUME(GR_TABLE(GR->rcut,ct,i,j)) /
                    GR_TABLE(GR->mesh,ct,i,j);
            /* 1/sqrt(N) statistical accuracy: */
            cc = 1. / SQUARE(GR_TABLE(GR->accuracy,ct,i,j)) / cc;
            /* so how many instances are needed to achieve GR->accuracy */
            if (necessary < cc) necessary = cc;
        }
    necessary++;

    GR->save_freq = necessary / GR->autosaves + 1;
    Fprintf(info, "-> need to accumulate %d INDEPENDENT instances.\n",
            necessary);
    Fprintf(info, "Matlab output = \"%s\",  autosaves = %d,\n",
            GR->fn_matlab, GR->autosaves);
    Fprintf(info, "-> automatically save at every %d instances.\n",
            GR->save_freq);
    GR->instances = 0;
    GR->info = info;
    return (necessary);
} /* end Gr_Reset() */


/* Accumulate a configuration instance to counters */
void Gr_Accumulate (Alib_Declare_Config, Chemtab *ct, Tp *tp, Gr *GR)
{
    register int i, m;
    int j, k, di, dj, dk, ii, jj, kk, l, n;
    double cc, thickness[DIMENSION], ds[DIMENSION], dx[DIMENSION+1];

    GR->instances++;
    cc = VMAX (SQUARE(ct->t), GR->rcut);
    M3rowthicknesses (H, thickness);
    GR->nbin[DIMENSION] = 1;
    for (i=0; i<DIMENSION; i++)
    {
        if ( GR->pbc[i] && (thickness[i] < 2*cc) )
            pe ("Gr_Accumulate: losing r<rmax images may happen\n"
                "in direction %d, please increase H[%d][].\n", i, i);
        GR->nbin[i] = thickness[i] / cc;
        if (GR->nbin[i] < 1) GR->nbin[i] = 1;
        GR->nbin[DIMENSION] *= GR->nbin[i];
    }
    /* bin-bin: prepare for maximal bin connectivity */
    GR->binbinlist = Irecreatelist
        (&GR->binbindx, GR->nbin[DIMENSION], DIMENSION3);
    for (i=0; i<GR->nbin[0]; i++)
        for (j=0; j<GR->nbin[1]; j++)
            for (k=0; k<GR->nbin[2]; k++)
            {
                l = GR_BINDEX(GR,i,j,k);
                for (di=-1; di<=1; di++)
                    for (dj=-1; dj<=1; dj++)
                        for (dk=-1; dk<=1; dk++)
                        {
                            ii = i + di;
                            jj = j + dj;
                            kk = k + dk;
                            if (ii >= GR->nbin[0])
                            {
                                if (GR->pbc[0]) ii-=GR->nbin[0];
                                else continue;
                            }
                            if (ii < 0)
                            {
                                if (GR->pbc[0]) ii+=GR->nbin[0];
                                else continue;
                            }
                            if (jj >= GR->nbin[1])
                            {
                                if (GR->pbc[1]) jj-=GR->nbin[1];
                                else continue;
                            }
                            if (jj < 0)
                            {
                                if (GR->pbc[1]) jj+=GR->nbin[1];
                                else continue;
                            }
                            if (kk >= GR->nbin[2])
                            {
                                if (GR->pbc[2]) kk-=GR->nbin[2];
                                else continue;
                            }
                            if (kk < 0)
                            {
                                if (GR->pbc[2]) kk+=GR->nbin[2];
                                else continue;
                            }
                            m = GR_BINDEX(GR,ii,jj,kk);
                            /* make sure it is a new record */
                            for (n=GR->binbindx[2*l];
                                 n<GR->binbindx[2*l+1]; n++)
                                if (m == GR->binbinlist[n]) break;
                            if (n == GR->binbindx[2*l+1])
                                GR->binbinlist[GR->binbindx[2*l+1]++] = m;
                        }
            }

    n = ceil( DOUBLE(*np) / GR->nbin[DIMENSION] * GR_BINATOM_RATIO );
    GR->binlist = Irecreatelist (&GR->bindx, GR->nbin[DIMENSION], n);
    REALLOC( Gr_Accumulate, GR->mybin, *np, int );

    REALLOC( Gr_Accumulate, GR->s, DIMENSION*(*np), double );
    VEQV( DIMENSION*(*np), *s, GR->s );
    /* atoms not under PBC cannot escape: they get stuck on the wall */
    for (j=0; j<DIMENSION; j++)
        if (GR->pbc[j]) for (i=*np; i--;) { Trim(GR->s[DIMENSION*i+j]); }
        else warn_and_correct_if_stuck_on_the_wall (*np, GR->s, j, stderr);

    /* bin-atom: fill atoms into bins */
    for (i=0; i<*np; i++)
        Iappend(GR->bindx, GR->binlist, i, GR->mybin[i]=GR_BINDEXS(GR,GR->s,i),
                0, GR->nbin[DIMENSION]);

    for (i=0; i<*np; i++)
    { /* already compressed bin-bin list */
        ii = GR->binbindx[2*GR->mybin[i]];
        jj = GR->binbindx[2*GR->mybin[i]+1];
        for (j=ii; j<jj; j++)
        { /* neighboring bins */
            l = GR->binbinlist[j];
            for (k=GR->bindx[2*l]; k<GR->bindx[2*l+1]; k++)
            { /* particles in neighboring bins */
                m = GR->binlist[k];
                if (i >= m) continue;
                V3SUB ( &((GR->s)[DIMENSION*m]), &((GR->s)[DIMENSION*i]), ds);
                if (GR->pbc[0]) Image(ds[0]);
                if (GR->pbc[1]) Image(ds[1]);
                if (GR->pbc[2]) Image(ds[2]);
                V3M3LENGTH2 (ds, H, dx);
                /* the meshes are regular in \delta r2 */
                n = dx[DIMENSION] / GR_TABLE(GR->rcut2,ct,tp[i],tp[m]) *
                    GR_TABLE(GR->mesh,ct,tp[i],tp[m]);
                if (n < GR_TABLE(GR->mesh,ct,tp[i],tp[m]))
                    GR_TABLE(GR->count,ct,tp[i],tp[m])[n] ++;
            }
        }
    }

    cc = M3VOLUME(H);
    for (i=0; i<ct->t; i++)
        for (j=i; j<ct->t; j++)
            if (i == j)
                GR_TABLE(GR->normalization,ct,i,j) +=
                    DOUBLE(ct->count[i]) / 2. * (ct->count[i] - 1) / cc;
            else
                GR_TABLE(GR->normalization,ct,i,j) +=
                    DOUBLE(ct->count[i]) * ct->count[j] / cc;

    Free (GR->mybin);
    Ifreelist (&GR->bindx, &GR->binlist);
    Ifreelist (&GR->binbindx, &GR->binbinlist);
    Free (GR->s);

    if (GR->instances % GR->save_freq == 0)
    {
        Gr_Save (Config_Alib_to_Alib, ct, GR, GR->info);
        Fcr (GR->info);
    }
    return;
} /* end Gr_Accumulate() */


/* Save result in counters to GR->fn_matlab */
void Gr_Save (Alib_Declare_Config, Chemtab *ct, Gr *GR, FILE *info)
{
    int i,j,n;
    double r2del,radius,cc,*accuracy;
    char *name;
    FILE *fp_matlab;
    if (GR->instances <= 0) return;
    fp_matlab = wopen (GR->fn_matlab);
    Fprintf(info, "Computing g(r)'s after %d accumulations:\n",
            GR->instances);
    MALLOC ( Gr_Save, accuracy, SQUARE(ct->t), double );
    for (i=0; i<ct->t; i++)
        for (j=i; j<ct->t; j++)
        {
            if (j!=0) fprintf (fp_matlab, "\ninput('press a key... ');\n\n");
            name = str2(blank_advance(SYMBOL(ct->first[i])),
                        blank_advance(SYMBOL(ct->first[j])));
            fprintf (fp_matlab, "%s = [\n", name);
            r2del = GR_TABLE(GR->rcut2,ct,i,j) / GR_TABLE(GR->mesh,ct,i,j);
            for (n=0; n<GR_TABLE(GR->mesh,ct,i,j); n++)
            {
                radius = sqrt((n+0.5)*r2del);
                cc = GR_TABLE(GR->count,ct,i,j)[n] /
                    ( GR_TABLE(GR->normalization,ct,i,j) * (
                        SPHERE_VOLUME(sqrt((n+1)*r2del)) -
                        SPHERE_VOLUME(sqrt(n*r2del)) ) );
                fprintf (fp_matlab, "%g %g %g\n",
                         radius, radius*ulength_IN_A, cc);
            }
            fprintf (fp_matlab, "];\nclf;plot(%s(:,2),%s(:,3));\n",name,name);
            fprintf (fp_matlab, "v = axis; axis([0 %g v(3) v(4)]);\n",
                     GR_TABLE(GR->rcut,ct,i,j)*ulength_IN_A);
            fprintf (fp_matlab, "xlabel('r [A]'); ylabel('g_{%s}(r)');\n",
                     name);
            fprintf (fp_matlab, "title('%s-%s Radial Distribution Function"
                     "');\n", blank_advance(SYMBOL(ct->first[i])),
                     blank_advance(SYMBOL(ct->first[j])));
            GR_TABLE(accuracy,ct,j,i) = GR_TABLE(accuracy,ct,i,j) = pow (
                GR_TABLE(GR->normalization,ct,i,j) *
                SPHERE_VOLUME(GR_TABLE(GR->rcut,ct,i,j)) /
                GR_TABLE(GR->mesh,ct,i,j), -0.5);
        }
    if (info) Sfpr(info, "current accuracy = %M.", accuracy, ct->t,ct->t);
    free(accuracy);
    Fprintf(info, "g(r)'s saved on Matlab script \"%s\".\n", GR->fn_matlab);
    fclose (fp_matlab);
    return;
} /* end Gr_Save() */


/* Free all memory allocations and set NULL */
void Gr_Free (Gr *GR)
{
    Free (GR->mybin);
    Ifreelist (&GR->bindx, &GR->binlist);
    Ifreelist (&GR->binbindx, &GR->binbinlist);
    Free (GR->s);
    Free (GR->rcut2);
    Free (GR->normalization);
    Free (GR->count);
    Free (GR->workspace);
    Free (GR->accuracy);
    Free (GR->mesh);
    Free (GR->rcut);
    return;
} /* end Gr_Free() */


#ifdef _Gr_TEST
#define NC  10
int main (int argc, char *argv[])
{
    register int i,j,instances;
    Aapp_Define_Config;
    Chemtab ct[1];
    Tp *tp=NULL;
    Gr GR[1] = {0};
    Config_build_Xtal (Config_Aapp_to_Alib, NC,NC,NC,
                       Xtal_abstracts+XTAL_ZNS, "Si",-1., "C",-1., -1.);
    rebind_ct (Config_Aapp_to_Alib, NULL, ct, &tp, stdout); cr();
    Gr_Recreate_Form (Config_Aapp_to_Alib, ct, GR);
    GR->autosaves = 1000000;
    instances = Gr_Reset (Config_Aapp_to_Alib, ct, GR, stdout); cr();
    for (i=0; i<instances; i++)
    {
        VFrandom (3*np, s, j);
        Gr_Accumulate (Config_Aapp_to_Alib, ct, tp, GR);
    }
    Gr_Save (Config_Aapp_to_Alib, ct, GR, stdout);
    Gr_Free (GR);
    return (0);
}
#undef NC
#endif /* _Gr_TEST */
