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

const Atom_coordination_color ATOM_COORDINATION_COLOR
[ATOM_COORDINATION_MAX+1] =
{
    {"white",             1.,       1.,       1.,     },
    {"turquoise",         64/255.,  224/255., 208/255.},
    {"burlywood",         222/255., 184/255., 135/255.},
    {"ForestGreen",       34/255.,  139/255., 34/255. },
    {"gray",              190/255., 190/255., 190/255.}, /* Si is white? */
    {"red",               1.,       0.,       0.,     },
    {"brown",             165/255., 42/255.,  42/255.,},
    {"purple",            160/255., 32/255.,  240/255.},
    {"LightSteelBlue",    176/255., 196/255., 222/255.}, /* steel is bcc */
    {"magenta",           1.,       0.,       1.,     },
    {"LimeGreen",         50/255.,  205/255., 50/255. },
    {"MediumVioletRed",   199/255., 21/255.,  133/255.},
    {"LightGoldenrod",    238/255., 221/255., 130/255.}, /* gold is fcc */
    {"RoyalBlue",         65/255.,  105/255., 225/255.},
    {"MistyRose",         255/255., 228/255., 225/255.},
    {"cyan",              0.,       1.,       1.,     },
    {"coral",             255/255., 127/255., 80/255. },
    {"green",             0.,       1.,       0.,     },
    {"YellowGreen",       154/255., 205/255., 50/255. },
    {"MediumSpringGreen", 0/255.,   250/255., 154/255.},
    {"yellow",            1.,       1.,       0.,     },
    {"thistle",           216/255., 191/255., 216/255.},
    {"chartreuse",        127/255., 255/255., 0/255.  },
    {"plum",              221/255., 160/255., 221/255.},
    {"DeepPink",          255/255., 20/255.,  147/255.},
};

const struct Mendeleyev MENDELEYEV [MENDELEYEV_MAX+1] = {
#include "Mendeleyev.c"
};

/* From our periodic table, find the atom Z corresponding */
/* to the "symbol" string. If not found, return 0.        */
int search_atom_by_symbol (char *symbol)
{
    register int j;
    for (j=0; j<=MENDELEYEV_MAX; j++)
        if (!strcasecmp(symbol,ATOM_SYMBOL(j)))
            break;
    if (j>MENDELEYEV_MAX) return(0);
    else return(j);
} /* end search_atom_by_symbol() */


/* From our periodic table, find the atom Z corresponding to */
/* the "symbol" string. If not found, print error and exit.  */
int Search_atom_by_symbol (char *symbol)
{
    register int j;
    for (j=0; j<=MENDELEYEV_MAX; j++)
        if (!strcasecmp(symbol,ATOM_SYMBOL(j)))
            break;
    if (j>MENDELEYEV_MAX)
        pe("Search_atom_by_symbol: non-existent atom symbol \"%s\".\n",
           symbol);
    return(j);
} /* end Search_atom_by_symbol() */


/************************************************************************/
/* Given an atomistic configuration, (re)bind sequential chemical index */
/* to each atom, which points to an entry of Z (and count) in "ct".     */
/* When allocating the index, give priority to those that appear in     */
/* "specification" string, for instance "Si C" makes Si=0,C=1, " C O"   */
/* makes C=0,O=1, etc. If specification=NULL, treat it as "". Return    */
/* the total number of chemical species found and print a report.       */
/************************************************************************/
int rebind_ct (Alib_Declare_Config, char *specification,
               Chemtab *ct, Tp **tp, FILE *out)
{
    register int i, j;
    char S[(MENDELEYEV_MAX+1)*SYMBOL_SIZE];
    double totalweight, weight[MENDELEYEV_MAX+1];
    REALLOC( rebind_ct, *tp, *np, Tp );
    if (specification==NULL) specification="";
    i = strlen(specification);
    if (!isfactor(SYMBOL_CHAR,i))
        pe ("rebind_ct: length of \"%s\" is not an integer multiple\n"
            "of std. chemical symbol string length %d.\n", specification,
            SYMBOL_CHAR);
    for (ct->t=0; ct->t<i/SYMBOL_CHAR; ct->t++)
    {
        if (ct->t >= MENDELEYEV_MAX)
            pe ("rebind_ct: you wrote too much.\n");
        for (j=0; j<SYMBOL_CHAR; j++)
            S[ct->t*SYMBOL_SIZE+j] = specification[ct->t*SYMBOL_CHAR+j];
        S[ct->t*SYMBOL_SIZE+j] = EOS;
        /* we may well get a zero - but don't worry */
        ct->Z[ct->t] = search_atom_by_symbol (&S[ct->t*SYMBOL_SIZE]);
        ct->count[ct->t] = 0;
    }
    for (i=0; i<(*np); i++)
    {
        for (j=0; j<ct->t; j++)
            if ( ( (ct->count[j] > 0) &&
                   !strcmp(SYMBOL(i), SYMBOL(ct->first[j])) ) ||
                 ( (ct->count[j] == 0) &&
                   !strcmp(SYMBOL(i), &S[j*SYMBOL_SIZE]) ) )
            {
                (*tp)[i] = j;
                ct->count[j]++;
                if (ct->count[j] == 1) ct->first[j] = i;
                break;
            }
        if (j == ct->t)
        {
            if (j >= MENDELEYEV_MAX)
                pe ("rebind_ct: too many species.\n");
            (*tp)[i] = j;
            ct->Z[j] = search_atom_by_symbol(SYMBOL(i));
            ct->count[j] = 1;
            ct->first[j] = i;
            ct->t++;
        }
    }
    if (out)
    {
        totalweight = 0;
        for (j=0; j<ct->t; j++) weight[j] = 0;
        for (i=0; i<(*np); i++)
        {
            weight[(int)((*tp)[i])] += (*mass)[i] * umass_IN_AMU;
            totalweight += (*mass)[i] * umass_IN_AMU;
        }
        fprintf (out, "------------- Chemical Species Report -----------\n");
        fprintf (out, "Idx Type  Z  Avg.Mass   Count  Abundance  Wt.Pct.\n");
        for (j=0; j<ct->t; j++)
            fprintf (out, "%2d   %2s %3d  %7.3f%9d   %6.2f%%  %6.2f%%\n",
                     j, (ct->count[j] == 0)? &S[j*SYMBOL_SIZE] :
                     SYMBOL(ct->first[j]), ct->Z[j],
                     (ct->count[j]>0)?weight[j]/ct->count[j]:
                     ATOM_MASS_IN_AMU(ct->Z[j]),
                     ct->count[j], ct->count[j]*100./(*np),
                     weight[j]*100./totalweight);
        fprintf (out, "-------------------------------------------------\n");
    }
    return(ct->t);
} /* end rebind_ct() */


/* Assign new symbol[] and mass[] according to new tp[] */
void rematch_ct (Chemtab *ct, Tp **tp, Alib_Declare_Config)
{
    register int i,j,k;
    ConfigStack cs[1];
    CONFIG_PUSH(cs);
    Config_alloc(Config_Alib_to_Alib);
    for (j=0; j<(*np); j++)
    {
        i = ct->first[(int)((*tp)[j])];
        for (k=0; k<SYMBOL_SIZE; k++)
            SYMBOL(j)[k] = Symbol(cs->symBOL,i)[k];
        (*mass)[j] = cs->MASS[i];
    }
    VEQV( DIMENSION*(*np), cs->S,  *s  );
    VEQV( DIMENSION*(*np), cs->S1, *s1 );
    CONFIG_erase (cs);
    return;
} /* end rematch_ct() */


/* Re-index the atoms according to their assigned chemical      */
/* index and mass. This is useful for efficient CONFIG storage. */
void Config_compact_index (Chemtab *ct, Tp **tp, Alib_Declare_Config)
{
    register int i,j,k;
    int *idx,*hook;
    double *tmp;
    Tp *oldtp;
    ConfigStack cs[1];
    MALLOC (Config_compact_index, tmp,   *np, double);
    MALLOC (Config_compact_index, idx,   *np, int);
    CALLOC (Config_compact_index, hook,  *np, int);
    MALLOC (Config_compact_index, oldtp, *np, Tp);
    for (j=0; j<(*np); j++)
    { /* hash values */
        tmp[j] = (double)((*tp)[j]) * 1000. + (*mass)[j] * umass_IN_AMU;
        idx[j] = j;
        oldtp[j] = (*tp)[j];
    }
    for (i=0; i<ct->t; i++) hook[ct->first[i]] = i+1;
    qsort_numerical_recipes (*np, tmp, idx, USE_OLD_IDX);
    free (tmp);
    CONFIG_PUSH(cs);
    Config_alloc(Config_Alib_to_Alib);
    for (j=0; j<(*np);)
    {
        if (hook[idx[j]] > 0) ct->first[hook[idx[j]]-1] = j;
        (*tp)[j] = oldtp[idx[j]];
        Config_RETRIEVE (cs,idx[j],Config_Alib_to_Alib,j,k);
    }
    CONFIG_erase (cs);
    free(idx);
    free(hook);
    free(oldtp);
    return;
} /* end Config_compact_index() */
