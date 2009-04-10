/*******************************************/
/* Threaded Atomistic Configuration Viewer */
/*                                         */
/* Fri Mar 31 2000 Ju Li <liju99@mit.edu>  */
/*******************************************/

#include "A.h"

AX_Float *scratch_r=NULL, *scratch_g=NULL, *scratch_b=NULL;
int scratch_np=0, scratch_n[3];


bool scratch_color (int iw)
{
    register int i;
    if (np != scratch_np)
    {
        printf ("scratch_color: current np=%d not equals to that of\n"
                "old scratch=%d, Please rescratch.", np, scratch_np);
        return (FALSE);
    }
    strcpy(AX_title[iw], str4(fbasename, " (", "scratch", ")"));
    AXSetName(iw);
    XStoreName(AX_display[iw],xterm_win,AX_title[iw]);
    XSetIconName(AX_display[iw],xterm_win,AX_title[iw]);
    for (i=0; i<np; i++)
        AX_3D_AssignRGB(B->BALL[i], scratch_r[i], scratch_g[i], scratch_b[i]);
    n[iw].color_mode = COLOR_MODE_SCRATCH;
    if (n[iw].bond_mode)
    {
        bond_xtal_origin_update(iw);
        bond_atom_color_update(iw);
    }
    else
    {
        n[iw].bond_xtal_origin_need_update = TRUE;
        n[iw].bond_atom_color_need_update = TRUE;
    }
    return (TRUE);
} /* end scratch_color() */


/* assign color blocks to serve as scratches */
bool rescratch (int iw)
{
    register int i,j;
    double ds[3],new_s[3];
    char danswer[MAX_FILENAME_SIZE],*buf;
    REALLOC (rescratch, scratch_r, np, AX_Float);
    REALLOC (rescratch, scratch_g, np, AX_Float);
    REALLOC (rescratch, scratch_b, np, AX_Float);
    xterm_get_focus(iw); clear_stdin_buffer();
    sprintf (danswer, "2 2 2");
    buf = readline_gets("\nScratch into n1xn2xn3 color blocks", danswer);
    sscanf (buf, "%d %d %d\n", scratch_n, scratch_n+1, scratch_n+2);
    xterm_release_focus(iw);
    scratch_np = np;
    if (scratch_n[0] < 1) scratch_n[0] = 2;
    if (scratch_n[1] < 1) scratch_n[1] = 2;
    if (scratch_n[2] < 1) scratch_n[2] = 2;
    V3mM3 (n[iw].xtal_origin, HI, ds);
    V3TRIM (ds, ds);
    for (i=0; i<np; i++)
    {
        V3SUB ( &(s[DIMENSION*i]), ds, new_s );
        V3TriM ( new_s );
        j = INT(new_s[0] * scratch_n[0]) +
            INT(new_s[1] * scratch_n[1]) +
            INT(new_s[2] * scratch_n[2]);
        if (ISEVEN(j))
        { /* firebrick4 */
            scratch_r[i] = 159 / 255.;
            scratch_g[i] =  26 / 255.;
            scratch_b[i] =  26 / 255.;
        }
        else
        { /* CadetBlue3 */
            scratch_r[i] = 122 / 255.;
            scratch_g[i] = 197 / 255.;
            scratch_b[i] = 205 / 255.;
        }
    }
    return(scratch_color(iw));
} /* end rescratch() */


/* Release allocated memory */
void scratch_free (int iw)
{
    int i = ((scratch_r!=NULL)+(scratch_g!=NULL)+(scratch_b!=NULL));
    Free (scratch_r);
    Free (scratch_g);
    Free (scratch_b);
    printf ("\nPrevious scratch colors freed (%s bytes).\n",
            strmem((long)i*scratch_np*sizeof(AX_Float)));
    scratch_np = 0;
    return;
} /* end scratch_free() */
