/*******************************************/
/* Threaded Atomistic Configuration Viewer */
/*                                         */
/* Fri Mar 31 2000 Ju Li <liju99@mit.edu>  */
/*******************************************/

#include "A.h"

Rcut_patch rcut_patch [RCUT_PATCH_MAX];
int rcut_patching, rcut_patch_top, rcut_patch_item;
char rcut_patch_pairname[TERMSIZE] = {0};

void start_rcut_patch (int iw)
{
    char ini[3], inj[3], *p, *q;
    int i, j, patch_Zi, patch_Zj, k, *count;
    if (rcut_patching)
    {
        printf ("You haven't yet finished patching the last rcut.\n");
        printf ("Press R to submit the last item.\n");
        return;
    }
    if (rcut_patch_pairname[0] == EOS)
    { /* popularity contest */
        CALLOC (start_rcut_patch, count, SQUARE(ct->t), int);
        for (i=0; i<np; i++)
            for (j=N->idx[i]; j<N->idx[i+1]; j++)
                count[ct->t*((int)tp[i])+((int)tp[N->list[j]])]++;
        for (i=0; i<ct->t; i++)
            for (j=i+1; j<ct->t; j++)
                count[ct->t*i+j] = count[ct->t*j+i] =
                    count[ct->t*i+j] + count[ct->t*j+i];
        for (i=j=k=0; i<SQUARE(ct->t); i++)
            if (count[i] > j)
            {
                j = count[i];
                k = i;
            }
        free (count);
        i = k / ct->t;
        j = k % ct->t;
        sprintf (rcut_patch_pairname, "%s %s",
                 COMPACT_SYMBOL(ATOM_SYMBOL(ct->Z[i])),
                 COMPACT_SYMBOL(ATOM_SYMBOL(ct->Z[j])));
    }
    if (ct->t > 1)
    {
        xterm_get_focus(iw); clear_stdin_buffer();
        printf("\nStart patching neighbor distance cutoff between (%s): ",
               rcut_patch_pairname);
        q = eos(rcut_patch_pairname);
        for (p=q; p<rcut_patch_pairname+TERMSIZE-5; p++)
            if ((*p=(char)getc(stdin))=='\n') break;
        if (p==q) q=rcut_patch_pairname;
        *p = EOS;
        printf ("\n"); xterm_release_focus(iw);
        sscanf (q, "%2s %2s", ini, inj);
        SAFE_SYMBOL(ini);
        SAFE_SYMBOL(inj);
        patch_Zi = search_atom_by_symbol(ini);
        patch_Zj = search_atom_by_symbol(inj);
    }
    else
    {
        patch_Zi = patch_Zj = ct->Z[0];
        printf("\nStart patching neighbor distance cutoff between %s %s:\n",
               COMPACT_SYMBOL(ATOM_SYMBOL(patch_Zi)),
               COMPACT_SYMBOL(ATOM_SYMBOL(patch_Zj)));
    }
    sprintf(rcut_patch_pairname, "%s %s",
            COMPACT_SYMBOL(ATOM_SYMBOL(patch_Zi)),
            COMPACT_SYMBOL(ATOM_SYMBOL(patch_Zj)));
    for (k=0; k<rcut_patch_top; k++)
        if ( ( (rcut_patch[k].Zi == patch_Zi) &&
               (rcut_patch[k].Zj == patch_Zj) ) ||
             ( (rcut_patch[k].Zi == patch_Zj) &&
               (rcut_patch[k].Zj == patch_Zi) ) ) break;
    if (k == rcut_patch_top)
    { /* new one: */
        rcut_patch_item = rcut_patch_top++;
        rcut_patch[rcut_patch_item].Zi = MIN(patch_Zi,patch_Zj);
        rcut_patch[rcut_patch_item].Zj = MAX(patch_Zi,patch_Zj);
        rcut_patch[rcut_patch_item].rcut = NEIGHBORLIST_RCUT_RATIO *
            ( ATOM_RADIUS_IN_A(patch_Zi) + ATOM_RADIUS_IN_A(patch_Zj) );
    }
    else
    { /* the pair is already in patch list */
        rcut_patch_item = k;
    }
    if (rcut_patch_top >= RCUT_PATCH_MAX)
    {
        printf ("RCUT_PATCH_MAX = %d reached.\n", RCUT_PATCH_MAX);
        printf ("Cannot add rcut patch no more.\n\n");
        rcut_patch_top--;
        return;
    }
    else rcut_patching = 1;
    printf ("RCUT(%s) = %g.\n", rcut_patch_pairname,
            rcut_patch[rcut_patch_item].rcut);
    return;
} /* end start_rcut_patch() */


bool apply_rcut_patch (int iw)
{
    int i, j, k, rcut_patched = 0;
    for (i=0; i<ct->t; i++)
        for (j=i; j<ct->t; j++)
            for (k=0; k<rcut_patch_top; k++)
                if ( ( ( (rcut_patch[k].Zi == ct->Z[i]) &&
                         (rcut_patch[k].Zj == ct->Z[j]) ) ||
                       ( (rcut_patch[k].Zi == ct->Z[j]) &&
                         (rcut_patch[k].Zj == ct->Z[i]) ) ) &&
                     ( NEIGHBOR_TABLE(N->rcut,ct,i,j) !=
                       rcut_patch[k].rcut ) )
                {
                    rcut_patched = 1;
                    NEIGHBOR_TABLE(N->rcut,ct,i,j) =
                        NEIGHBOR_TABLE(N->rcut,ct,j,i) =
                        rcut_patch[k].rcut;
                }
    if (rcut_patched)
    {
        Neighborlist_Recreate (Config_Aapp_to_Alib, NULL, ct, &tp, N);
        geo_clear_has_evaluated_flags();
        evaluate_geo_measures();
        if ((n[iw].xtal_mode) && (n[iw].color_mode == COLOR_MODE_COORD))
            assign_coordination_color(iw);
        else if ( (n[iw].color_mode == COLOR_MODE_AUXILIARY) &&
                  INW(n[iw].auxiliary_idx-CONFIG_MAX_AUXILIARY,
                      MAX_GEO_MEASURES) )
            color_encode_auxiliary(iw);
        Config_to_3D_Bonds (n[iw].bond_radius);
        if (n[iw].bond_mode)
        {
            bond_xtal_origin_update (iw);
            bond_atom_color_update(iw);
        }
        else
        {
            n[iw].bond_xtal_origin_need_update = TRUE;
            n[iw].bond_atom_color_need_update = TRUE;
        }
    }
    return (rcut_patched);
} /* end apply_rcut_patch() */


void finish_rcut_patch (int iw)
{
    printf ("RCUT(%s) = %g.\n\n",
            rcut_patch_pairname, rcut_patch[rcut_patch_item].rcut);
    Sfpr(stdout, "rlist = %M (reduced)\n ", N->rcut, ct->t,ct->t);
    print_coordination_histogram(); cr();
    rcut_patching = 0;
    return;
} /* end finish_rcut_patch() */
