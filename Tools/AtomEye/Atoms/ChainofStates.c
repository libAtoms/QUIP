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

/***************************************************************/
/* Toolkit for chain-of-states methods                         */
/*                                                             */
/* Henkelman, Jóhannesson and Jónsson,                         */
/* Methods for Finding Saddle Points and Minimum Energy Paths, */
/* in Progress on Theoretical Chemistry and Physics,           */
/* Ed. S. D. Schwartz (Kluwer Academic Publishers, 2000)       */
/* pp. 269-300.                                                */
/*                                                             */
/* Henkelman and Jónsson,                                      */
/* J Chem. Phys. 113 (2000) 9978-9985; 9901-9904.              */
/*                                                             */
/* The fundamental assumption is that np, symbol, mass, etc.   */
/* intrinsic properties do not change from chain node to chain */
/* node. Only H[][], s[], s1[] may change.                     */
/***************************************************************/


/* Save current configuration to ChainofStates node i */
void ChainofStates_Push (Alib_Declare_Config, ChainofStates *c, int i)
{
    M3_To_V9 (H, c->H+9*i);
    VEQV (3*(*np), *s, c->s+3*(*np)*i);
    if (c->s1) VEQV(3*(*np), *s1, c->s1+3*(*np)*i);
    return;
} /* end ChainofStates_Push() */


/* Retrieve ChainofStates node i to current configuration */
void ChainofStates_Pop (ChainofStates *c, int i, Alib_Declare_Config)
{
    V9_To_M3 (c->H+9*i, H);
    VEQV (3*(*np),  c->s+3*(*np)*i, *s);
    if (c->s1) VEQV(3*(*np), c->s1+3*(*np)*i, *s1);
    return;
} /* end ChainofStates_Pop() */


/* NULL-safe freeing of memory H, s, s1 */
void ChainofStates_Free (ChainofStates *c)
{
    Free(c->H);
    Free(c->s);
    Free(c->s1);
    c->N = c->checksum = 0;
    return;
} /* end ChainofStates_Free() */


/* Reload chain of states c based on a sorted UNIX glob file pattern */
void ChainofStates_Load (char *GlobPattern, FILE *info, int flags,
                         Alib_Declare_Config, ChainofStates *c)
{
    int i;
    glob_t globbuf;
    if (glob(GlobPattern, GLOB_TILDE | GLOB_ERR, NULL, &globbuf))
        pe("cannot find chain of configuration files \"%s\"\n", GlobPattern);
    for (i=0; i<globbuf.gl_pathc; i++)
    {
        Fprintf (info, "Loading %s%s\n", globbuf.gl_pathv[i],
                 (i==0)? ":" : " ...");
        Config_Load (globbuf.gl_pathv[i], (i==0)?info:NULL,
                     Config_Alib_to_Alib);
        if (i == 0)
        {
            c->checksum = ConfigChecksum(Config_Alib_to_Alib);
            c->N = globbuf.gl_pathc;
            REALLOC( ChainofStates_Load, c->H, 9*c->N, double );
            REALLOC( ChainofStates_Load, c->s, 3*(*np)*c->N, double );
            if (flags | CHAIN_OF_STATES_NO_VELOCITY)
            { /* no need to save velocity information */
                Free(c->s1);
            }
            else
            {
                REALLOC( ChainofStates_Load, c->s1, 3*(*np)*c->N, double );
            }
        }
        else if (ConfigChecksum(Config_Alib_to_Alib) != c->checksum)
            pe ("ChainofStates_Load: the configuration \"%s\"\n"
                "is not isoatomic with \"%s\".\n", globbuf.gl_pathv[i],
                globbuf.gl_pathv[0]);
        ChainofStates_Push (Config_Alib_to_Alib, c, i);
    }
    c->alpha = pow(DOUBLE(*np),1./6);
    /* by default, alpha=N^(1/6) to equalize the stiffness */
    if (i>1)
    {
        ChainofStates_Pop(c,0,Config_Alib_to_Alib);
        Fprintf (info, "\nPopping back to %s configuration.\n\n",
                 globbuf.gl_pathv[0]);
    }
    globfree (&globbuf);
    return;
} /* end ChainofStates_Load() */

#ifdef _ChainofStates_Load_TEST
int main (int argc, char *argv[])
{
    Aapp_Define_Config;
    ChainofStates c[1]={0};
    ChainofStates_LOAD("ChainofStatesDir/FCC*.cfg", Config_Aapp_to_Alib, c);
    ChainofStates_Free (c);
    return (0);
}
#endif /* _ChainofStates_Load_TEST */


#define MAXDIGITS 8
/* Save chain of states c to a directory ("Cfg/") or pattern ("Cfg/#") */
void ChainofStates_Save
(Alib_Declare_Config, ChainofStates *c, FILE *info, char *pattern)
{
    int i, width;
    char *bef, num[MAXDIGITS], *q, *fn;
    double *sold, *s1old, Hold[3][3];

    sold = *s;
    s1old = *s1;
    M3EQV (H, Hold);
    MALLOC(Config_alloc, *s,      DIMENSION*(*np),   double);
    MALLOC(Config_alloc, *s1,     DIMENSION*(*np),   double);

    bef = IOclone(pattern);
    q = strstr(bef, "#");
    if (q)
    {
        *q = EOS;
        q++;
    }
    else q=strend(bef);

    width = sprintf(num, "%d", c->N)+1;
    if (width > MAXDIGITS-1)
        pe ("ChainofStates_Save: increase MAXDIGITS=%d\n", MAXDIGITS);

    for (i=0; i<c->N; i++)
    {
        ChainofStates_Pop (c, i, Config_Alib_to_Alib);
        sprintf(num, "%0"STR(MAXDIGITS)"d", i);
        fn = str3( bef, num+MAXDIGITS-width, q );
        Fprintf (info, "Saving %s ...\n", fn);
        if (c->s1)
            CONFig_SAVE(Config_Alib_to_Alib, fn,
                        "%.15g %.15g %.15g", " %.15g %.15g %.15g");
        else CONFIG_NV_SAVE
                 (Config_Alib_to_Alib, fn, "%.15g %.15g %.15g", 0);
    }

    Free(*s);
    Free(*s1);
    *s = sold;
    *s1 = s1old;
    M3EQV (Hold, H);
    IOfree (bef);
    
    return;
} /* end ChainofStates_Save() */

#ifdef _ChainofStates_Save_TEST
int main (int argc, char *argv[])
{
    Aapp_Define_Config;
    ChainofStates c[1]={0};
    ChainofStates_LOAD("ChainofStatesDir/FCC*.cfg", Config_Aapp_to_Alib, c);
    ChainofStates_SAVE(Config_Aapp_to_Alib, c, "/tmp/a/#.bz");
    ChainofStates_Free (c);
    return (0);
}
#endif /* _ChainofStates_Save_TEST */


/* Compute the difference between node j and i in hyperspace */
/* cd is realloc-ed and should be initialized as {0}.        */
void ChainofStates_Diff (Alib_Declare_Config, ChainofStates *c, int j, int i,
                         ChainofStates_Difference *cd)
{
    register int k;
    double dx[4];
    M3 H0;

    /* difference in [X,Y]=[s*H0,alpha*H] space */
    V9SUB( c->H+9*j, c->H+9*i, cd->dH );
    cd->dR2 = V9LENGTH2( cd->dH ) * SQUARE( c->alpha );

    V9_To_M3 ( c->H, H0 );

    REALLOC (ChainofStates_Diff, cd->ds, DIMENSION*(*np), double);
    for (k=*np; k--;)
    {
        V3SUB (c->s+3*((*np)*j+k), c->s+3*((*np)*i+k), cd->ds+3*k);
        V3ImagE (cd->ds+3*k);
        V3M3LENGTH2 (cd->ds+3*k, H0, dx);
        cd->dR2 += dx[3];
    }
    if (c->s1)
    {
        REALLOC (ChainofStates_Diff, cd->ds1,      DIMENSION*(*np),   double);
        Vsub (3*(*np), c->s1+3*(*np)*j, c->s1+3*(*np)*i, cd->ds1);
    }
    cd->dR = sqrt(cd->dR2);
    return;
} /* end ChainofStates_Diff() */


/* Total path length in [X,Y]=[s*H0,alpha*H] space */
double ChainofStates_TotalPathLength (int N, ChainofStates_Difference *cd)
{
    int n;
    double totalR = 0;
    for (n=0; n<N-1; n++) totalR += cd[n].dR;
    return (totalR);
} /* end ChainofStates_TotalPathLength() */


/* NULL-safe freeing of memory */
void ChainofStates_Difference_Free(ChainofStates_Difference *cd)
{
    V9ZERO(cd->dH);
    Free(cd->ds);
    Free(cd->ds1);
    cd->dR2 = 0;
    cd->dR = 0;
    return;
} /* end ChainofStates_Difference_Free() */


/* Use linear interpolation in hyperspace to add nodes. plan[i] contains */
/* the number of new nodes to be added between original node i and i+1.  */
void ChainofStates_Add_Nodes (int *plan, Alib_Declare_Config, ChainofStates *c)
{
    int i, iend, jend, k, m, sum;
    double ratio;
    ChainofStates_Difference cd[1]={0};

    for (sum=0,i=0; i<c->N-1; i++)
    {
        if (plan[i] < 0)
            pe("Negative node addition (%d:%d) not allowed.\n", i, plan[i]);
        sum += plan[i];
    }

    if (sum > 0)
    {
        jend = c->N + sum;
        REALLOC( ChainofStates_Add_Nodes, c->H, 9*jend, double );
        REALLOC( ChainofStates_Add_Nodes, c->s, 3*(*np)*jend, double );
        if (c->s1)
            REALLOC( ChainofStates_Add_Nodes, c->s1, 3*(*np)*jend, double );
        V9EQV(c->H+9*(c->N-1), c->H+9*(jend-1));
        VEQV(3*(*np), c->s+3*(*np)*(c->N-1), c->s+3*(*np)*(jend-1));
        if (c->s1)
            VEQV(3*(*np), c->s1+3*(*np)*(c->N-1), c->s1+3*(*np)*(jend-1));
    }
    else return;
    
    for (i=c->N-2; i>=0; i--)
    {
        iend = jend-2-plan[i];
        /* printf ("i=%d iend=%d jend=%d\n", i, iend, jend); */
        if (plan[i] > 0)
        {
            ChainofStates_Diff (Config_Alib_to_Alib, c, i+1, i, cd);
            /* printf ("%g\n", cd->dR2); */
            for (k=1; k<=plan[i]; k++)
            {
                ratio = k/(plan[i]+1.);
                V9ADDMUL(c->H+9*i, ratio, cd->dH, c->H+9*(iend+k));
                Vaddmul(3*(*np), c->s+3*(*np)*i, ratio, cd->ds,
                        c->s+3*(*np)*(iend+k));
                if (c->s1) Vaddmul(3*(*np), c->s1+3*(*np)*i, ratio, cd->ds1,
                                   c->s1+3*(*np)*(iend+k));
            }
        }
        if (iend > i)
        {
            V9EQV(c->H+9*i, c->H+9*iend);
            VEQV(3*(*np), c->s+3*(*np)*i, c->s+3*(*np)*iend);
            if (c->s1) VEQV(3*(*np), c->s1+3*(*np)*i, c->s1+3*(*np)*iend);
        }
        jend = iend+1;
    }
    ChainofStates_Difference_Free (cd);
    c->N += sum;
    return;
} /* end ChainofStates_Add_Nodes() */

#ifdef _ChainofStates_Add_Nodes_TEST
int main (int argc, char *argv[])
{
    Aapp_Define_Config;
    ChainofStates c[1]={0};
    /* int plan1[1]={1}, plan2[2]={4,0}; */
    int plan1[3]={10,10,10}, plan2[2]={4,0};
    /* ChainofStates_LOAD("ChainofStatesDir/FCC*.cfg", */
    /* Config_Aapp_to_Alib, c); */
    /* ChainofStates_Add_Nodes (plan1, Config_Aapp_to_Alib, c); */
    /* ChainofStates_Add_Nodes (plan2, Config_Aapp_to_Alib, c); */
    ChainofStates_LOAD("ChainofStatesDir/Travel*.cfg",
                       Config_Aapp_to_Alib, c);
    ChainofStates_Add_Nodes (plan1, Config_Aapp_to_Alib, c);
    ChainofStates_SAVE(Config_Aapp_to_Alib, c, "/tmp/a_saveit#.cfg");
    ChainofStates_Free (c);
    return (0);
}
#endif /* _ChainofStates_Add_Nodes_TEST */


#ifdef _linear_path
int main (int argc, char *argv[])
{
    Aapp_Define_Config;
    ChainofStates c[1]={0};
    char *start_fname, *finish_fname, *saveto;
    int i, plan[1];

    if (argc < 4)
    {
        printf ("\nPurpose: Generate a sequence of atomistic configurations\n"
                "         that linearly interpolate between start and finish."
                "\n\n");
        printf ("Usage: %s start_fname finish_fname total_nodes\n"
                "       %s start_fname finish_fname total_nodes saveto\n\n",
                argv[0], argv[0]);
        return (1);
    }
    start_fname = argv[1];
    finish_fname = argv[2];
    plan[0] = atoi(argv[3])-2;
    if (argc == 4) saveto = str2(start_fname, ".");
    else saveto = argv[4];
    for (i=0; i<2; i++)
    {
        printf ("Loading %s%s\n", argv[i+1], (i==0)? ":" : " ...");
        Config_LOAD (argv[i+1], Config_Aapp_to_Alib);
        if (i == 0)
        {
            c->checksum = ConfigChecksum(Config_Aapp_to_Alib);
            c->N = 2;
            REALLOC( ChainofStates_Load, c->H, 9*c->N, double );
            REALLOC( ChainofStates_Load, c->s, 3*np*c->N, double );
            Free(c->s1);
        }
        else if (ConfigChecksum(Config_Aapp_to_Alib) != c->checksum)
            pe ("ChainofStates_Load: the configuration \"%s\"\n"
                "is not isoatomic with \"%s\".\n", argv[i+1],
                argv[i+0]);
        ChainofStates_Push (Config_Aapp_to_Alib, c, i);
    }
    ChainofStates_Add_Nodes (plan, Config_Aapp_to_Alib, c);
    ChainofStates_SAVE(Config_Aapp_to_Alib, c, saveto);
    ChainofStates_Free (c);
    return (0);
}
#endif /* _linear_path */


/* Get rid of nodes with taglist[i]!=0; return the number of nodes lost */
int ChainofStates_Delete_Nodes
(char taglist[], Alib_Declare_Config, ChainofStates *c)
{
    int i, j;
    for (i=j=0; i<c->N; i++)
        if (taglist[i]==0)
        { /* keep */
            if (j!=i)
            {
                V9EQV(c->H+9*i, c->H+9*j);
                VEQV(3*(*np), c->s+3*(*np)*i, c->s+3*(*np)*j);
                if (c->s1) VEQV(3*(*np), c->s1+3*(*np)*i, c->s1+3*(*np)*j);
            }
            j++;
        }
    j = c->N - j;
    c->N -= j;
    REALLOC( ChainofStates_Delete_Nodes, c->H, 9*c->N, double );
    REALLOC( ChainofStates_Delete_Nodes, c->s, 3*(*np)*c->N, double );
    if (c->s1) REALLOC( ChainofStates_Delete_Nodes, c->s1,
                        3*(*np)*c->N, double );
    return (j);
} /* end ChainofStates_Delete_Nodes() */

#ifdef _ChainofStates_Delete_Nodes_TEST
int main (int argc, char *argv[])
{
    Aapp_Define_Config;
    ChainofStates c[1]={0};
    char taglist1[4]={0,1,1,0};
    int plan1[2]={20};
    ChainofStates_LOAD("ChainofStatesDir/Travel*.cfg",
                       Config_Aapp_to_Alib, c);
    printf ("%d\n", ChainofStates_Delete_Nodes
            (taglist1, Config_Aapp_to_Alib, c));
    ChainofStates_Add_Nodes (plan1, Config_Aapp_to_Alib, c);
    ChainofStates_SAVE(Config_Aapp_to_Alib, c, "/tmp/a");
    ChainofStates_Free (c);
    return (0);
}
#endif /* _ChainofStates_Delete_Nodes_TEST */


/* NULL-safe freeing of NeighborlistPool memory */
void ChainofStates_NeighborlistPool_Free(ChainofStates_NeighborlistPool *Np)
{
    int i;
    for (i=0; i<Np->size; i++) Neighborlist_Free (Np->pool+i);
    Free (Np->pool);
    Np->size = 0;
    return;
} /* end ChainofStates_NeighborlistPool_Free() */


/* Create a neighborlist pool. cramped=0: a neighborlist for */
/* everyone;    cramped=1: everyone shares one neighborlist. */
void ChainofStates_NeighborlistPool_Recreate
(int cramped, ChainofStates *c, ChainofStates_NeighborlistPool *Np)
{
    int i, desired_size;
    desired_size = cramped? 1 : c->N;
    if (desired_size != Np->size)
    {
        ChainofStates_NeighborlistPool_Free (Np);
        Np->size = desired_size;
        REALLOC_ZERO( ChainofStates_NeighborlistPool_Recreate,
                      Np->pool, Np->size, Neighborlist );
    }
    return;
} /* end ChainofStates_NeighborlistPool_Create() */


/* Swap out a member from the pool */
Neighborlist *ChainofStates_NeighborlistPool_Select
(ChainofStates_NeighborlistPool *Np, int i)
{
    if (i<Np->size) return(Np->pool+i);
    else return(Np->pool);
} /* end ChainofStates_NeighborlistPool_Select() */
