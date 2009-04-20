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

/*************************************************************/
/* bin-bin/bin-atom/atom-atom lists assuming initial [0,1)   */
/* s-bounds under H[][], which does not have to be PBC. If   */
/* min_shrinkage is 1 and max_atom_displacement is 0,        */
/* it is taken to mean a disposable list and all extra       */
/* memory is freed except a final compressed atom-atom list. */
/* Otherwise one registers anchor sa[] and H0[][], frees no  */
/* memory, and provides an uncompressed atom-atom list: at   */
/* each new configuration, the user should call Neighborlist */
/* Maintenance. Under PBC, the user himself should NOT fold  */
/* s[] into [0,1); Maintenance will do that only at the atom */
/* re_anchoring event, as the ANCHORS must be in [0,1). If a */
/* direction is not PBC, and a s[i] exceeds [0,1), it would  */
/* be stuck at the boundary and a warning is issued. The     */
/* atom-atom list is guaranteed to list all pairs <= RCUT    */
/* by making sure that at any atom's re_anchoring event,     */
/* all surrounding anchors <= RLIST are recorded.            */
/*                                                           */
/* Of all libAtoms modules, only this one require [0,1) s-   */
/* bounds. load_from_pdb(), for example, would not assign    */
/* bounding box H[][] if CRYST1 tag exists; yet a lot of PDB */
/* files with CRYST1 tag do exceed [0,1) bounds, and         */
/* load_from_pdb() does not trim them. Same for UUUS format, */
/* although it is recommended that a config which is not     */
/* meant to be non-PBC call Config_fold_into_PBC before it's */
/* saved. Therefore there are three options for this module  */
/* when it is given a non-[0,1) config: QUIT, FOLD_INTO_PBC, */
/* BOUNDING_BOX. There's no guarantee of the correct action, */
/* but a good guess would be that if the config is from PDB, */
/* BOUNDING_BOX should be used; otherwise FOLD_INTO_PBC.     */
/*************************************************************/

#define DEFORMABLE(N)     ( (N)->min_shrinkage  < 1 )
#define MOVABLE(N)        ( (N)->max_atom_displacement > 0 )
#define LIST_REUSABLE(N)  ( MOVABLE(N) || DEFORMABLE(N) )
#define LIST_DISPOSABLE(N) \
  (((N)->min_shrinkage == 1) && ((N)->max_atom_displacement == 0))
#define LIST_COMPRESSED(np,N) ((N)->list-(N)->idx==(np)+1)

/* Recreate form space and fill in defaults */
void Neighborlist_Recreate_Form
(Alib_Declare_Config, Chemtab *ct, Neighborlist *N)
{
    register int i,j,k;
    /* default treatment of s-bounds overflow is folding + warning */
    N->s_overflow_err_handler = NEIGHBORLIST_DEF_S_OVERFLOW_ERR_HANDLER;
    N->small_cell_err_handler = NEIGHBORLIST_DEF_SMALL_CELL_ERR_HANDLER;
    /* use pairwise 1/2 saving? */
    N->pairwise = NEIGHBORLIST_DEF_PAIRWISE;
    /* PBC? */
    for (i=0; i<DIMENSION; i++)
        N->pbc[i] = NEIGHBORLIST_DEF_PBC;
    N->min_shrinkage = NEIGHBORLIST_DEF_MIN_SHRINKAGE;
    N->max_atom_displacement = NEIGHBORLIST_DEF_MAX_ATOM_DISPLACEMENT;
    N->track_dx = NEIGHBORLIST_DEF_TRACK_DX;
    /* create rcut table */
    REALLOC( Neighborlist_Recreate_Form, N->rcut, SQUARE(ct->t), double );
    /* fill in probable first-neighbor distances */
    for (i=0; i<ct->t; i++)
        for (j=0; j<ct->t; j++)
            NEIGHBOR_TABLE(N->rcut,ct,i,j) = NEIGHBORLIST_RCUT_RATIO *
                ( ATOM_RADIUS_IN_A(ct->Z[i]) +
                  ATOM_RADIUS_IN_A(ct->Z[j]) ) / ulength_IN_A;
    /* create neighbor number table */
    REALLOC( Neighborlist_Recreate_Form, N->neighbormax, ct->t, int );
    /* ball filling estimates (before pairwise reduction) */
    for (i=0; i<ct->t; i++)
    {
        N->neighbormax[i] = 0;
        for (j=0; j<ct->t; j++)
        {
            k = CUBE( ATOM_RADIUS_IN_A(ct->Z[i]) /
                      ATOM_RADIUS_IN_A(ct->Z[j]) + 2 ) -
                CUBE( ATOM_RADIUS_IN_A(ct->Z[i]) /
                      ATOM_RADIUS_IN_A(ct->Z[j]) );
            if (k > N->neighbormax[i])
                N->neighbormax[i] = k;
        }
    }
    N->keepcounter = FALSE;
    return;
} /* end Neighborlist_Recreate_Form() */


/* Recreate neighborlist according to submitted form */
void Neighborlist_Recreate
(Alib_Declare_Config, FILE *info, Chemtab *ct, Tp **tp, Neighborlist *N)
{
    register int i, m;
    int j, k, di, dj, dk, ii, jj, kk, l, n, old_n;
    double rmax, thickness[DIMENSION], ds[DIMENSION], dx[DIMENSION];

    /* Verify the only assumption of this program */
    for (i=0; i<*np; i++)
    {
        /* if (IN( (*s)[DIMENSION*i],   -EPS, 0 )) (*s)[DIMENSION*i]   = 0; */
        /* if (IN( (*s)[DIMENSION*i+1], -EPS, 0 )) (*s)[DIMENSION*i+1] = 0; */
        /* if (IN( (*s)[DIMENSION*i+2], -EPS, 0 )) (*s)[DIMENSION*i+2] = 0; */
        if ( V3NEED_TRIM(&((*s)[DIMENSION*i])) )
        {
            Fprintf(info, "** %s atom (1-%d) has\n""** s = (%g %g %g)\n",
                    word_for_order(i+1), *np, (*s)[DIMENSION*i],
                    (*s)[DIMENSION*i+1], (*s)[DIMENSION*i+2]);
            V3mM3( &((*s)[DIMENSION*i]), H, dx );
            V3MuL( ulength_IN_A, dx );
            Fprintf(info, "** or x = (%g %g %g) A.\n\n", dx[0], dx[1], dx[2]);
            Fprintf(info, "This is impermissible for constructing "
                    "neighborlist,\n");
            if (N->s_overflow_err_handler ==
                NEIGHBORLIST_S_OVERFLOW_ERR_HANDLER_QUIT)
                pe ("Neighborlist_Recreate() exiting...\n");
            else if (N->s_overflow_err_handler ==
                     NEIGHBORLIST_S_OVERFLOW_ERR_HANDLER_FOLD_INTO_PBC)
            { /* even though the user may not mean PBC */
                Fprintf(info, "so we have to fold all atoms into [0,1).\n\n");
                Config_fold_into_PBC (Config_Alib_to_Alib);
                break;
            }
            else  /* NEIGHBORLIST_S_OVERFLOW_ERR_HANDLER_BOUNDING_BOX */
            { /* even though the user may mean PBC */
                Fprintf(info, "so we use BOUNDING BOX as H[][] instead:\n\n");
                Config_to_bounding_box_config(Config_Alib_to_Alib,info);
                Fcr(info);
                break;
            }
        }
    }

    /* Also check RCUT symmetry in the case of pairwise interaction */
    if (N->pairwise)
        for (i=0; i<ct->t; i++)
            for (j=i+1; j<ct->t; j++)
                if (NEIGHBOR_TABLE(N->rcut,ct,i,j) !=
                    NEIGHBOR_TABLE(N->rcut,ct,j,i))
                {
                    Fprintf(info, "Warning: rcut[%d][%d] != "
                            "rcut[%d][%d]\n", i,j, j,i);
                    Fprintf(info, "-> %f will be used for both.\n\n",
                            NEIGHBOR_TABLE(N->rcut,ct,i,j)=
                            NEIGHBOR_TABLE(N->rcut,ct,j,i)=
                            MAX(NEIGHBOR_TABLE(N->rcut,ct,i,j),
                                NEIGHBOR_TABLE(N->rcut,ct,j,i)));
                }

    /* Print out the form */
    Fprintf(info, "Create neighborlist with the following parameters:\n");
    Fprintf(info, "pairwise saving = %s,  track_dx = %s\n",
            BWORD(N->pairwise), BWORD(N->track_dx));
    for (i=0; i<DIMENSION; i++)
        Fprintf(info, "pbc[%d] = %s,  ", i, BWORD(N->pbc[i]));
    Fcr(info);
    Fprintf(info, "Strain Session min_shrinkage = %g,\n"
            "All Atom Tether max_atom_displacement = %g (reduced),\n",
            N->min_shrinkage, N->max_atom_displacement);
    if ( (N->min_shrinkage > 1) || (N->min_shrinkage <= 0) )
        pe ("Neighborlist_Recreate: min_shrinkage must be (0,1].\n\n");
    if ( N->max_atom_displacement < 0 )
        pe ("Neighborlist_Recreate: max_atom_displacement must be >= 0.\n\n");
        
    REALLOC( Neighborlist_Recreate_Form, N->rlist, SQUARE(ct->t), double );
    if (LIST_DISPOSABLE(N))
    {
        Fprintf(info, "-> COMPRESSED disposable atom-atom list\n"
                "with Rlist_ij = Rcut_ij;\n\n");
        VEQV( SQUARE(ct->t), N->rcut, N->rlist );
    }
    else
    {
        Fprintf(info, "-> UNCOMPRESSED reusable atom-atom list\n"
                "with %g x Rlist_ij = Rcut_ij + 2 x %g;\n\n",
                N->min_shrinkage, N->max_atom_displacement);
        if (info) Sfpr(info, "Rcut  = %M (reduced)\n ", N->rcut, ct->t,ct->t);
        for (i=0; i<ct->t; i++)
            for (j=0; j<ct->t; j++)
                NEIGHBOR_TABLE(N->rlist,ct,i,j) =
                    (NEIGHBOR_TABLE(N->rcut,ct,i,j) +
                     2 * N->max_atom_displacement) / N->min_shrinkage;
    }

    for (i=0; i<ct->t; i++)
        for (j=0; j<ct->t; j++)
            if (NEIGHBOR_TABLE(N->rcut,ct,i,j) <= 0)
                NEIGHBOR_TABLE(N->rlist,ct,i,j) = 0;
    
    REALLOC( Neighborlist_Recreate_Form, N->rlist2, SQUARE(ct->t), double );
    for (i=0; i<ct->t; i++)
        for (j=0; j<ct->t; j++)
            NEIGHBOR_TABLE(N->rlist2,ct,i,j) =
                SQUARE(NEIGHBOR_TABLE(N->rlist,ct,i,j));

    if (info) Sfpr(info, "rlist = %M (reduced)\n ", N->rlist, ct->t,ct->t);
    rmax = VMAX (SQUARE(ct->t), N->rlist);
    Fprintf(info, "MAX(rlist) = %g  ->  bin thicknesses should >= %g\n",
            rmax, rmax);

    old_n = *np;
  restart:
    M3rowthicknesses (H, thickness);
    if (info) Vfpr(info, "since H thicknesses = %M", thickness, DIMENSION);

    if ( N->small_cell_err_handler ==
         NEIGHBORLIST_SMALL_CELL_ERR_HANDLER_MULTIPLY )
    {
        ii = ceil(2*rmax / thickness[0]);
        jj = ceil(2*rmax / thickness[1]);
        kk = ceil(2*rmax / thickness[2]);
        if ( (ii > 1) || (jj > 1) || (kk > 1) )
        {
            Config_multiply ( ii, jj, kk, Config_Alib_to_Alib );
	    //            rebind_CT (Config_Alib_to_Alib, "", ct, tp);
#ifndef ATOMEYE_LIB
	    rebind_CT (Config_Alib_to_Alib, "", ct, tp);
#else
	    rebind_ct (Config_Alib_to_Alib, "", ct, tp, NULL);
#endif
            goto restart;
        }
    }

    N->nbin[DIMENSION] = 1;
    for (i=0; i<DIMENSION; i++)
    {
        if ( N->pbc[i] && (thickness[i] < 2*rmax) &&
             (N->small_cell_err_handler !=
              NEIGHBORLIST_SMALL_CELL_ERR_HANDLER_NOCHECK) )
            pe ("Neighborlist_Recreate: losing r<rmax images may happen\n"
                "in direction %d, please increase H[%d][] "
                "(rmax=%g,thickness=%g).\n", i, i, rmax,thickness[i]);
        N->nbin[i] = thickness[i] / rmax;
        if (N->nbin[i] < 1) N->nbin[i] = 1;
        if (N->nbin[i] > NEIGHBORLIST_MAXBIN) N->nbin[i]=NEIGHBORLIST_MAXBIN;
        N->nbin[DIMENSION] *= N->nbin[i];
    }
    Fprintf(info, "-> %d x ", N->nbin[0]);
    for (i=1; i<DIMENSION-1; i++)
        Fprintf(info, "%d x ", N->nbin[i]);
    Fprintf(info, "%d = %d bins.\n", N->nbin[DIMENSION-1],
            N->nbin[DIMENSION]);
    
    /* bin-bin: prepare for maximal bin connectivity */
    N->binbinlist = Irecreatelist
        (&N->binbindx, N->nbin[DIMENSION], DIMENSION3);
    for (i=0; i<N->nbin[0]; i++)
        for (j=0; j<N->nbin[1]; j++)
            for (k=0; k<N->nbin[2]; k++)
            {
                l = NEIGHBORLIST_BINDEX(N,i,j,k);
                for (di=-1; di<=1; di++)
                    for (dj=-1; dj<=1; dj++)
                        for (dk=-1; dk<=1; dk++)
                        {
                            ii = i + di;
                            jj = j + dj;
                            kk = k + dk;
                            if (ii >= N->nbin[0])
                            {
                                if (N->pbc[0]) ii-=N->nbin[0]; else continue;
                            }
                            if (ii < 0)
                            {
                                if (N->pbc[0]) ii+=N->nbin[0]; else continue;
                            }
                            if (jj >= N->nbin[1])
                            {
                                if (N->pbc[1]) jj-=N->nbin[1]; else continue;
                            }
                            if (jj < 0)
                            {
                                if (N->pbc[1]) jj+=N->nbin[1]; else continue;
                            }
                            if (kk >= N->nbin[2])
                            {
                                if (N->pbc[2]) kk-=N->nbin[2]; else continue;
                            }
                            if (kk < 0)
                            {
                                if (N->pbc[2]) kk+=N->nbin[2]; else continue;
                            }
                            m = NEIGHBORLIST_BINDEX(N,ii,jj,kk);
                            /* make sure it is a new record */
                            for (n=N->binbindx[2*l]; n<N->binbindx[2*l+1]; n++)
                                if (m == N->binbinlist[n]) break;
                            if (n == N->binbindx[2*l+1])
                                N->binbinlist[N->binbindx[2*l+1]++] = m;
                        }
            }
    Fprintf(info, "bin-bin list created.\n\n");

    /* it is reasonable to assume that bin-connectivity never changes */
    N->binbinlist = Icompresslist (&N->binbindx, N->nbin[DIMENSION]);
    Ilist_statistics ("bin-bin", N->binbindx, N->nbin[DIMENSION], TRUE, info);
    Fcr(info);

    Fprintf(info, "On average, there are %.2f atoms per bin,\n"
            "but we will offer space of %d for each.\n",
            DOUBLE(*np) / N->nbin[DIMENSION], n = ceil
            (DOUBLE(*np) / N->nbin[DIMENSION] * NEIGHBORLIST_BINATOM_RATIO));
    N->binlist = Irecreatelist (&N->bindx, N->nbin[DIMENSION], n);
    REALLOC( Neighborlist_Recreate, N->mybin, *np, int );

    /* bin-atom: fill atoms into bins */
    for (i=0; i<*np; i++)
        Iappend(N->bindx, N->binlist, i, N->mybin[i]=
                NEIGHBORLIST_BINDEXS(N,*s,i), 0, N->nbin[DIMENSION]);
    Fprintf(info, "bin-atom list created.\n\n");

    Ilist_statistics("bin-atom", N->bindx, N->nbin[DIMENSION], FALSE, info);
    Fcr(info);

    /* atom-atom */
    for (n=i=0; i<ct->t; i++)
    {
        if ( N->pairwise && (!N->keepcounter) )
            N->neighbormax[i] = N->neighbormax[i] / 2;
        n += N->neighbormax[i] * ct->count[i];
    }
    if (N->pairwise) Fprintf(info, "After pairwise/own_pair saving,\n");
    if (info)
        ifpr(info, "likely MAX(neighbor) = %M atoms,", N->neighbormax, ct->t);
    REALLOC( Neighborlist_Recreate, N->idx, 2*(*np)+1+n, int );
    N->idx[0] = N->idx[1] = 0;
    for (i=1; i<*np; i++)
        N->idx[2*i] = N->idx[2*i+1] =
            N->idx[2*i-1]+N->neighbormax[(int)(*tp)[i-1]];
    if ( (N->idx[2*(*np)]=n) != N->idx[2*i-1]+N->neighbormax[(int)(*tp)[i-1]] )
        pe ("Neighborlist_Recreate: checksum failed.\n");
    N->list = &N->idx[2*(*np)] + 1;
    for (i=0; i<*np; i++)
    { /* already compressed bin-bin list */
        ii = N->binbindx[N->mybin[i]];
        jj = N->binbindx[N->mybin[i]+1];
        for (j=ii; j<jj; j++)
        { /* neighboring bins */
            l = N->binbinlist[j];
            for (k=N->bindx[2*l]; k<N->bindx[2*l+1]; k++)
            { /* particles in neighboring bins */
                m = N->binlist[k];
                if (N->pairwise && (!own_pair(i,m))) continue;
                if (i == m) continue;
                V3SUB ( &((*s)[DIMENSION*m]), &((*s)[DIMENSION*i]), ds);
                if (N->pbc[0]) Image(ds[0]);
                if (N->pbc[1]) Image(ds[1]);
                if (N->pbc[2]) Image(ds[2]);
                V3mM3 (ds, H, dx);
                if (V3LENGTH2(dx) <=
                    NEIGHBOR_TABLE(N->rlist2,ct,(*tp)[i],(*tp)[m]))
                    Iappend (N->idx, N->list, m, i, 0, *np);
            }
        }
    }
    Fprintf(info, "atom-atom list created.\n\n");
    Ilist_statistics ("atom-atom", N->idx, *np, FALSE, info);
    Fcr(info);

    M3EQV (H, N->H0);
    if (LIST_DISPOSABLE(N))
    {
        N->list = Icompresslist (&N->idx, *np);
        Ilist_statistics ("atom-atom", N->idx, *np, TRUE, info);
        Fcr(info);
        Ifreelist (&N->binbindx, &N->binbinlist);
        Ifreelist (&N->bindx, &N->binlist);
        Free (N->mybin);
        Fprintf(info, "All bin-related allocations freed.\n");
    }
    else  /* LIST_REUSABLE(N)) */
    {
        if (!N->keepcounter) N->remesh_counter = 0;
        REALLOC (Neighborlist_Recreate, N->sa, DIMENSION*(*np), double);
        VEQV (DIMENSION*(*np), (*s), N->sa);
        Fprintf (info, "Particle anchors established.\n");
        REALLOC (Neighborlist_Recreate, N->maxtether2, ct->t, double);
        VMaX (ct->t, N->neighbormax, i, N->stack_max);
        N->stack_max *= NEIGHBORLIST_STACK_MAX_RATIO;
        REALLOC (Neighborlist_Recreate, N->stack, N->stack_max, int);
        if (!N->keepcounter) N->reanchor_counter = 0;
        if (N->track_dx)
        {
            REALLOC (Neighborlist_Recreate, N->dx, 2*DIMENSION*(*np), double);
            N->dxa = N->dx + DIMENSION*(*np);
            if (!N->keepcounter)
            {
                VZERO (2*DIMENSION*(*np), N->dx);
                Fprintf (info, "dx tracking initiated, counters cleared.\n");
            }
        }
        if (!N->keepcounter) N->maintenance_counter = 0;
    }
    N->keepcounter = TRUE;

    // jrk - return to old size
/*     *np = old_n; */
/*     fprintf(info, "Neighbourlist_Recreate: return to size %d atoms\n", *np); */
/*     Config_realloc (Config_Alib_to_Alib); */
/*     fprintf(info, "done\n"); */

    return;

} /* end Neighborlist_Recreate() */


/* Print appropriate allocation and occupation statistics */
void Neighborlist_print_statistics (int np, Neighborlist *N, FILE *info)
{
    if (LIST_REUSABLE(N))
    {
        if (DEFORMABLE(N))
            Fprintf(info, "We are at the %s strain session.\n",
                    word_for_order(N->remesh_counter+1));
        if (MOVABLE(N))
            Fprintf(info, "%d maintenances: total re_anchor count=%d,\n"
                    "or %g per maintenance, %g per atom per maintenance.\n\n",
                    N->maintenance_counter, N->reanchor_counter,
                    safe_avg(N->reanchor_counter,N->maintenance_counter),
                    safe_avg(N->reanchor_counter,N->maintenance_counter)/np);
        Ilist_statistics
            ("bin-bin", N->binbindx, N->nbin[DIMENSION], TRUE, info);
        Fcr(info);
        Ilist_statistics
            ("bin-atom", N->bindx, N->nbin[DIMENSION], FALSE, info);
        Fcr(info);
    }
    Ilist_statistics ("atom-atom", N->idx, np, LIST_DISPOSABLE(N), info);
    return;
} /* end Neighborlist_print_statistics() */


/* Free all memory allocations and set NULL */
void Neighborlist_Free (Neighborlist *N)
{
    Ifreelist (&N->binbindx, &N->binbinlist);
    Ifreelist (&N->bindx, &N->binlist);
    Free (N->mybin);
    Free (N->rlist2);
    Free (N->rlist);
    Free (N->rcut);
    Free (N->neighbormax);
    Ifreelist (&N->idx, &N->list);
    Free (N->sa);
    Free (N->maxtether2);
    Free (N->dx);
    Free (N->stack);
    return;
} /* end Neighborlist_Free() */


/* Transfer particle from one bin to another */
static void transfer_particle (int i, int o, int n, Neighborlist *N)
{
    register int m;
    /* delete this particle from bin o */
    for (m=N->bindx[2*o]; (m<N->bindx[2*o+1]) && (N->binlist[m]!=i); m++);
    /* if (m == N->bindx[2*o+1]) */
    /* pe ("transfer_particle: asked to delete particle %d\n" */
    /* "which does not exist in bin %d.\n", i, o); */
    for (m++; m<N->bindx[2*o+1]; m++) N->binlist[m-1] = N->binlist[m];
    N->bindx[2*o+1]--;
    /* insert this particle at the end of the new N->bin */
    Iappend (N->bindx, N->binlist, i, n, 0, N->nbin[DIMENSION]);
    N->mybin[i] = n;
    return;
} /* end transfer_particle() */


/* Atom refreshes its own and new neighbors' anchor lists */
static void re_anchor
(int np, double H[3][3], Chemtab *ct, Tp **tp, Neighborlist *N, int i)
{  /* better performance */
    register int m, n, stack_top;
    int j, k, l;
    V3 ds, dx;
    /* the bin that this atom (anchor) belongs to */
    m = NEIGHBORLIST_BINDEXS (N, N->sa, i);
    if (m != N->mybin[i]) transfer_particle (i, N->mybin[i], m, N);
    for (stack_top=0,j=N->binbindx[N->mybin[i]];
         j<N->binbindx[N->mybin[i]+1]; j++)
    { /* neighboring bins */
        l = N->binbinlist[j];
        for (k=N->bindx[2*l]; k<N->bindx[2*l+1]; k++)
        { /* particles in neighboring bins */
            m = N->binlist[k];
            if (i==m) continue;
            V3SUB ( &N->sa[DIMENSION*m], &N->sa[DIMENSION*i], ds);
            if (N->pbc[0]) Image(ds[0]);
            if (N->pbc[1]) Image(ds[1]);
            if (N->pbc[2]) Image(ds[2]);
            V3mM3 ( ds, N->H0, dx );
            if (N->pairwise)
            {
                if (V3LENGTH2(dx) <=
                    NEIGHBOR_TABLE(N->rlist2,ct,(*tp)[i],(*tp)[m]))
                { /* this pair interaction "belongs" to either i or m */
                    if (own_pair(i,m))
                        if (stack_top >= N->stack_max)
                            pe ("re_anchor: stack_max = %d exceeded.\n",
                                N->stack_max);
                        else N->stack[stack_top++] = m;
                    else
                    {  /* see if i is already in m's list */
                        for (n=N->idx[2*m];
                             (n<N->idx[2*m+1])&&(N->list[n]!=i); n++);
                        if (n == N->idx[2*m+1]) /* if not */
                            Iappend (N->idx, N->list, i, m, 0, np);
                    }
                }
            }
            else  /* not pairwise */
            {
                ds[0] = V3LENGTH2(dx);
                if ( ds[0] <= NEIGHBOR_TABLE(N->rlist2,ct,(*tp)[i],(*tp)[m]) )
                {
                    if (stack_top >= N->stack_max)
                        pe ("re_anchor: stack_max = %d exceeded.\n",
                            N->stack_max);
                    else N->stack[stack_top++] = m;
                }
                if ( ds[0] <= NEIGHBOR_TABLE(N->rlist2,ct,(*tp)[m],(*tp)[i]) )
                {
                    for (n=N->idx[2*m];
                         (n<N->idx[2*m+1])&&(N->list[n]!=i); n++);
                    if (n == N->idx[2*m+1]) /* if not */
                        Iappend (N->idx, N->list, i, m, 0, np);
                }
            }
        }
    }
    if (N->idx[2*i]+stack_top <= N->idx[2*i+2])
    { /* there is plenty of space */
        memcpy(N->list+N->idx[2*i], N->stack, stack_top*sizeof(int));
        N->idx[2*i+1] = N->idx[2*i] + stack_top;
    }
    else /* collisions do happen */
    { /* fill in some neighbors first */
        memcpy(N->list+N->idx[2*i], N->stack,
               (N->idx[2*i+2]-N->idx[2*i])*sizeof(int));
        N->idx[2*i+1] = N->idx[2*i+2];
        for (n=N->idx[2*i+2]-N->idx[2*i]; n<stack_top; n++)
            Iappend (N->idx, N->list, N->stack[n], i, 0, np);
    }
    return;
} /* end re_anchor() */


/* Warn and correct (to 0 or 1-TINY) if < 0 or >= 1 */
bool warn_and_correct_if_stuck_on_the_wall
(int np, double *s, int dim, FILE *info)
{
    register int i;
    bool stuck_at_0 = FALSE, stuck_at_1 = FALSE;
    for (i=np; i--;)
        if (s[DIMENSION*i+dim] < 0)
        {
            s[DIMENSION*i+dim] = 0;
            stuck_at_0 = TRUE;
        }
        else if (s[DIMENSION*i+dim] >= 1)
        {
            s[DIMENSION*i+dim] = 1-TINY;
            stuck_at_1 = TRUE;
        }
    if (stuck_at_0)
        Fprintf(info, "warn_and_correct_if_stuck_on_the_wall:\n"
                "atom(s) get stuck on s[%d]=0 wall.\n", dim);
    if (stuck_at_1)
        Fprintf(info, "warn_and_correct_if_stuck_on_the_wall:\n"
                "atom(s) get stuck on s[%d]=1 wall.\n", dim);
    return (stuck_at_0 || stuck_at_1);
} /* end warn_and_correct_if_stuck_on_the_wall() */


void Neighborlist_Maintenance
(Alib_Declare_Config, FILE *info, Chemtab *ct, Tp **tp, Neighborlist *N)
{
    register int i;
    int j;
    V3 eigval, ds, dx;
    M3 eta, Q;

    /* Check cell deformation */
    if ( M3NE(H, N->H0) )
    {
        if ( !DEFORMABLE(N) )
        {
            S3fPR (info, "H0 = %M (reduced)", N->H0);
            S3fPR (info, "H  = %M (reduced)", H);
            pe ("Neighborlist_Maintenance: You lied! min_shrinkage\n"
                "was declared to be one in Neighborlist_Recreate() but\n"
                "now H has changed.\n");
        }
        Lagrangian_strain (N->H0, H, eta);
        M3Diag (eta, eigval, Q);
        eigval[0] = sqrt((1+2*eigval[0]));
    }
    else eigval[0] = 1;
    if ( eigval[0] < N->min_shrinkage )
    {
        Fprintf(info, "min_shrinkage reached, re-meshing...\n\n");
        Config_fold_into_PBC (Config_Alib_to_Alib);
        Neighborlist_Recreate (Config_Alib_to_Alib, info, ct, tp, N);
        N->remesh_counter ++;
        N->reanchor_counter += *np;
        goto exit;
    }
    if (!MOVABLE(N)) goto exit;
    for (i=0; i<ct->t; i++)
    {
        N->maxtether2[i] = DOUBLE_PRECISION_INFINITY;
        for (j=0; j<ct->t; j++)
        {
            eigval[1] = (NEIGHBOR_TABLE(N->rlist,ct,i,j) * eigval[0] -
                         NEIGHBOR_TABLE(N->rcut, ct,i,j)) / 2;
            if (eigval[1] < N->maxtether2[i]) N->maxtether2[i] = eigval[1];
            if (!N->pairwise)
            {
                eigval[1] = (NEIGHBOR_TABLE(N->rlist,ct,j,i) * eigval[0] -
                             NEIGHBOR_TABLE(N->rcut, ct,j,i)) / 2;
                if (eigval[1] < N->maxtether2[i]) N->maxtether2[i] = eigval[1];
            }
        }
        N->maxtether2[i] *= N->maxtether2[i];
    }

    /* atoms not under PBC cannot escape: they get stuck on the wall */
    for (i=0; i<DIMENSION; i++)
        if (!N->pbc[i]) warn_and_correct_if_stuck_on_the_wall(*np,*s,i,info);

    for (i=0; i<*np; i++)
    { /* check individual particle drift from anchor */
        V3SUB ( &(*s)[DIMENSION*i], &N->sa[DIMENSION*i], ds);
        /* attention: no PBC is enforced here! */
        V3mM3 ( ds, H, dx );
        if ( N->track_dx )
            V3ADD (dx, &(N->dxa[DIMENSION*i]), &(N->dx[DIMENSION*i]));
        if ( (V3LENGTH2(dx) > N->maxtether2[(int)(*tp)[i]]) ||
             (N->idx[2*i+1] == N->idx[2*i+2]) )
        { /* set anchors */
            V3TriM ( &(*s)[DIMENSION*i] );
            V3EQV ( &(*s)[DIMENSION*i], &N->sa[DIMENSION*i] );
            if ( N->track_dx )
                V3EQV (&(N->dx[DIMENSION*i]), &(N->dxa[DIMENSION*i]));
            re_anchor (*np, H, ct, tp, N, i);
            N->reanchor_counter++;
        }
    }
    
  exit:
    N->maintenance_counter++;
    return;
} /* end Neighborlist_Maintenance() */


/* Self-check the neighborlist using brute force */
void Neighborlist_Check
(Alib_Declare_Config, FILE *info, Chemtab *ct, Tp **tp, Neighborlist *N)
{
    register int i,m,n,count;
    V3 ds, dx;
    for (i=0; i<*np; i++)
        for (m=0; m<*np; m++)
        {
            if (N->pairwise && (!own_pair(i,m))) continue;
            if (i == m) continue;
            V3SUB ( &((*s)[DIMENSION*m]), &((*s)[DIMENSION*i]), ds);
            if (N->pbc[0]) Image(ds[0]);
            if (N->pbc[1]) Image(ds[1]);
            if (N->pbc[2]) Image(ds[2]);
            V3mM3 (ds, H, dx);
            if (V3LENGTH2(dx) <
                SQUARE(NEIGHBOR_TABLE(N->rcut,ct,
                                      (*tp)[i],(*tp)[m])))
            {
                count = 0;
                for (n=N->idx[2*i]; n<N->idx[2*i+1]; n++)
                    if (N->list[n]==m) count++;
                if (count != 1)
                    pe ("Neighborlist_Check: we failed atom %d\n", i);
            }
        }
    return;
} /* end Neighborlist_Check() */


/* Convert a neighborlist with pairwise saving */
/* to full (redundant) list image.             */
void Neighborlist_create_nonpairwise_image
(Alib_Declare_Config, Chemtab *ct, Tp **tp, Neighborlist *N, Neighborlist *M)
{
    register int i,j;
    int n;
    if (!(N->pairwise))
        pe ("Neighborlist_create_nonpairwise_image: "
            "Neighborlist is already non-pairwise.\n");
    *M = *N;
    M->pairwise = FALSE;
    for (n=i=0; i<ct->t; i++)
    { /* N->keepcounter should be true */
        M->neighbormax[i] = N->neighbormax[i] * 2;
        n += M->neighbormax[i] * ct->count[i];
    }
    MALLOC ( Neighborlist_create_nonpairwise_image, M->idx, 2*(*np)+1+n, int );
    M->idx[0] = M->idx[1] = 0;
    for (i=1; i<*np; i++)
        M->idx[2*i] = M->idx[2*i+1] =
            M->idx[2*i-1]+M->neighbormax[(int)(*tp)[i-1]];
    if ( (M->idx[2*(*np)]=n) != M->idx[2*i-1]+M->neighbormax[(int)(*tp)[i-1]] )
        pe ("Neighborlist_create_nonpairwise_image: checksum failed.\n");
    M->list = &M->idx[2*(*np)] + 1;
    if ( LIST_COMPRESSED(*np,N) )
        for (i=0; i<*np; i++)
            for (j=N->idx[i]; j<N->idx[i+1]; j++)
            {
                Iappend (M->idx, M->list, N->list[j], i, 0, *np);
                Iappend (M->idx, M->list, i, N->list[j], 0, *np);
            }
    else
        for (i=0; i<*np; i++)
            for (j=N->idx[2*i]; j<N->idx[2*i+1]; j++)
            {
                Iappend (M->idx, M->list, N->list[j], i, 0, *np);
                Iappend (M->idx, M->list, i, N->list[j], 0, *np);
            }
    return;
} /* end Neighborlist_create_nonpairwise_image() */


/* Convert a neighborlist with pairwise saving */
/* to full (redundant) & compressed list image */
void Neighborlist_create_nonpairwise_compressed_image
(Alib_Declare_Config, Chemtab *ct, Tp **tp, Neighborlist *N, Neighborlist *M)
{
    Neighborlist_create_nonpairwise_image(Config_Alib_to_Alib, ct, tp, N, M);
    M->list = Icompresslist (&M->idx, *np);
    return;
} /* end Neighborlist_create_nonpairwise_compressed_image() */


/* Destroy created list image */
void Neighborlist_free_nonpairwise_image (Neighborlist *M)
{
    Free(M->idx);
} /* end Neighborlist_free_nonpairwise_image() */


/* Maximum number of records in the neighborlist */
int Neighborlist_max_records (Alib_Declare_Config, Neighborlist *N)
{
    register int i;
    int n=0,m;
    if ( LIST_COMPRESSED(*np,N) )
        for (i=0; i<*np; i++)
        {
            m = N->idx[i+1] - N->idx[i];
            if (m > n) n = m;
        }
    else
        for (i=0; i<*np; i++)
        {
            m = N->idx[2*i+1] - N->idx[2*i];
            if (m > n) n = m;
        }
    return(n);
} /* end Neighborlist_max_records() */


/* Maximum number of neighbors in the neighborlist */
int Neighborlist_max_neighbors (Alib_Declare_Config, Neighborlist *N)
{
    register int i,j;
    int n=0;
    short *coordination=NULL;
    CALLOC (Neighborlist_max_neighbors, coordination, *np, short);
    if ( LIST_COMPRESSED(*np,N) )
        for (i=0; i<*np; i++)
            for (j=N->idx[i]; j<N->idx[i+1]; j++)
            {
                coordination[i]++;
                if (N->pairwise) coordination[N->list[j]]++;
            }
    else
        for (i=0; i<*np; i++)
            for (j=N->idx[2*i]; j<N->idx[2*i+1]; j++)
            {
                coordination[i]++;
                if (N->pairwise) coordination[N->list[j]]++;
            }
    for (i=0; i<*np; i++)
        if (coordination[i] > n)
            n = coordination[i];
    Free (coordination);
    return(n);
} /* end Neighborlist_max_neighbors() */


#ifdef _Neighborlist_TEST
#define ULENGTH_IN_A  2.
#define LOAD_FILE  "/Home/Archive/NetApp/Mime/DNA.pdb.bz2"
int main (int argc, char *argv[])
{
    Chemtab ct[1];
    Aapp_Define_Config;
    Tp *tp=NULL;
    Neighborlist N[1] = {0};
    CONFIG_LOAD (LOAD_FILE, Config_Aapp_to_Alib);
    rebind_CT (Config_Aapp_to_Alib, NULL, ct, &tp); cr();
    NEIGHBORLIST_RECREATE (Config_Aapp_to_Alib, ct, tp, N); cr();
    Neighborlist_PRINT_statistics (np, N); cr();
    Neighborlist_Free (N);
    Press_return();
    return (0);
}
#endif /* _Neighborlist_TEST */
