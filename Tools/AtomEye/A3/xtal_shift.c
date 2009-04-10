/*******************************************/
/* Threaded Atomistic Configuration Viewer */
/*                                         */
/* Fri Mar 31 2000 Ju Li <liju99@mit.edu>  */
/*******************************************/

#include "A.h"

void atom_xtal_origin (double xtal_origin[3])
{
    register int i;
    double ds[3], new_s[3];
    V3mM3 (xtal_origin, HI, ds);
    V3TRIM (ds, ds);
    for (i=0; i<np; i++)
    {
        V3SUB ( &(s[DIMENSION*i]), ds, new_s );
        V3TriM ( new_s );
        V3mM3 ( new_s, H, B->BALL[i].x );
    }
    return;
} /* end atom_xtal_origin() */


void bond_xtal_origin_update (int iw)
{
    register int i,j;
    double *new_s, ds[3], tmp[3];
    if (temporary_disable_bond) return;
    MALLOC( bond_xtal_origin_update, new_s, DIMENSION*np, double );
    V3mM3 (n[iw].xtal_origin, HI, ds);
    V3TRIM (ds, ds);
    for (i=0; i<np; i++)
    {
        V3SUB ( &(s[DIMENSION*i]), ds, &(new_s[DIMENSION*i]) );
        V3TriM ( &(new_s[DIMENSION*i]) );
    }
    for (i=0; i<np; i++)
        for (j=N->idx[i]; j<N->idx[i+1]; j++)
        {
            V3SUB (&(new_s[DIMENSION*N->list[j]]), &(new_s[DIMENSION*i]), tmp);
            if ( V3NEED_IMAGE(tmp) )
            { /* g<0 means geometric invisibility */
                C->CYLINDER[j].g = -1;
                C->CYLINDER[j].radius = -1;
            }
            else
            {
                C->CYLINDER[j].g = BOND_G(i, N->list[j]);
                V3EQV (B->BALL[i].x, C->CYLINDER[j].x0);
                if (C->CYLINDER[j].r >= 0)
                    C->CYLINDER[j].radius = n[iw].bond_radius;
            }
        }
    free (new_s);
    n[iw].bond_xtal_origin_need_update = FALSE;
    return;
} /* end bond_xtal_origin_update() */


/* shift the crystal along axis-i by amount d*lengthscale*XTAL_SHIFT_GEAR */
bool xtal_shift (int iw, int i, double d)
{
    if (!n[iw].xtal_mode)
    {
        printf ("Crystal translation is only available under Xtal mode.\n");
        return(FALSE);
    }
    d *= lengthscale * XTAL_SHIFT_GEAR;
    V3ADDmuL (d, AX_3D[iw].V[i], n[iw].xtal_origin);
    atom_xtal_origin (n[iw].xtal_origin);
    if (n[iw].bond_mode) bond_xtal_origin_update (iw);
    else n[iw].bond_xtal_origin_need_update = TRUE;
    return (TRUE);
} /* end xtal_shift() */


/* use pointer device to shift the crystal */
bool pointer_grab_xtal_shift (int iw, int to_x, int to_y)
{
    double z, tmp[3];
    if (!n[iw].xtal_mode)
    {
        printf ("Crystal translation is only available under Xtal mode.\n");
        return(FALSE);
    }
    if ((n[iw].lx == to_x) && (n[iw].ly == to_y)) return (FALSE);
    if (n[iw].anchor >= 0)
        V3SUB (B->BALL[n[iw].anchor].x, AX_3D[iw].x, tmp);
    else V3SUB (n[iw].hook, AX_3D[iw].x, tmp);
    z = V3DOT (AX_3D[iw].V[2], tmp);
    tmp[0] = (n[iw].lx - to_x) / AX_3D[iw].k * z;
    V3ADDmuL (tmp[0], AX_3D[iw].V[0], n[iw].xtal_origin);
    tmp[1] = (n[iw].ly - to_y) / AX_3D[iw].k * z;
    V3ADDmuL (tmp[1], AX_3D[iw].V[1], n[iw].xtal_origin);
    n[iw].lx = to_x;
    n[iw].ly = to_y;
    atom_xtal_origin (n[iw].xtal_origin);
    if (n[iw].bond_mode) bond_xtal_origin_update (iw);
    else n[iw].bond_xtal_origin_need_update = TRUE;
    return (TRUE);
} /* end pointer_grab_xtal_shift() */


/* In viewframe (x0,V), the eyesight is ray (k0,k1,1) * z with z>0. */
/* We would like to know its intersection with box surface s * H,   */
/* where s_i=0 or 1 for certain i (encoded into surface_id 0-5).    */
/* Return z>0 if there is intersection, -1 otherwise.               */
double Eyesight_Intersect_H_Surface
(double H[3][3], int surface_id,
 double x0[3], double V[3][3], double k0, double k1,
 double s[3], double x[3])
{
    int i;
    double b[3], M[3][3], MI[3][3], tmp;
    i = surface_id / 2;
    if (surface_id % 2) V3SUB(x0,H[i],b);
    else V3EQV(x0,b);
    M3EQV(H,M);
    V3NEG(V[2],M[i]);
    V3SUBmuL(M[i],k0,V[0]);
    V3SUBmuL(M[i],k1,V[1]);
    tmp = M3DETERMINANT(M);
    if (ISTINY(tmp)) return(-1);
    M3INV(M,MI,tmp);
    V3mM3(b,MI,s);
    tmp = s[i];
    s[i] = surface_id % 2;
    V3mM3 (s,H,x);
    return (tmp);
} /* end Eyesight_Intersect_H_Surface() */


int Eyesight_Intersect_H_Box (int iw, int to_x, int to_y, double xx[3])
{
    register int i,j;
    int idx[6];
    double k0, k1, s[6][3], x[6][3];
    double z[6];
    k0 = (to_x + 0.5 - AX_3D[iw].wx) / AX_3D[iw].k;
    k1 = (to_y + 0.5 - AX_3D[iw].wy) / AX_3D[iw].k;
    for (j=i=0; i<6; i++)
        if ( ((z[i]=Eyesight_Intersect_H_Surface
               (H, i, AX_3D[iw].x, AX_3D[iw].V, k0, k1, s[i], x[i])) > 0) &&
             V3XIN(s[i],0,1) ) idx[j++] = i;
    if (j > 0)
    {
        qsort_numerical_recipes (j, z, idx, USE_OLD_IDX);
        V3EQV (x[idx[0]], xx);
        return (idx[0]);
    }
    else return (-1);
} /* end Eyesight_Intersect_H_Box() */


/* Use pointer device to shift the crystal */
bool pointer_xtal_shift (int iw, int to_x, int to_y)
{
    int i;
    double s[3], xx_last[3], xx_now[3];
    if (!n[iw].xtal_mode)
    {
        printf ("Crystal translation is only available under Xtal mode.\n");
        return(FALSE);
    }
    if ((n[iw].lx == to_x) && (n[iw].ly == to_y)) return (FALSE);
    V3mM3 (AX_3D[iw].x, HI, s);
    if (V3XIN(s,0,1)) return(pointer_grab_xtal_shift(iw,to_x,to_y));
    i = Eyesight_Intersect_H_Box (iw, n[iw].lx, n[iw].ly, xx_last);
    if (i >= 0) n[iw].last_surface_id = i;
    else if (Eyesight_Intersect_H_Surface
             (H, n[iw].last_surface_id, AX_3D[iw].x, AX_3D[iw].V,
              (n[iw].lx + 0.5 - AX_3D[iw].wx) / AX_3D[iw].k,
              (n[iw].ly + 0.5 - AX_3D[iw].wy) / AX_3D[iw].k, s, xx_last) <= 0)
    {
        n[iw].lx = to_x;
        n[iw].ly = to_y;
        return(FALSE);
    }
    i = Eyesight_Intersect_H_Box (iw, to_x, to_y, xx_now);
    if (i >= 0) n[iw].last_surface_id = i;
    else if (Eyesight_Intersect_H_Surface
             (H, n[iw].last_surface_id, AX_3D[iw].x, AX_3D[iw].V,
              (to_x + 0.5 - AX_3D[iw].wx) / AX_3D[iw].k,
              (to_y + 0.5 - AX_3D[iw].wy) / AX_3D[iw].k, s, xx_now) <= 0)
    {
        n[iw].lx = to_x;
        n[iw].ly = to_y;
        return(FALSE);
    }
    n[iw].lx = to_x;
    n[iw].ly = to_y;
    V3AdD(xx_last,n[iw].xtal_origin);
    V3SuB(n[iw].xtal_origin,xx_now);
    atom_xtal_origin (n[iw].xtal_origin);
    if (n[iw].bond_mode) bond_xtal_origin_update (iw);
    else n[iw].bond_xtal_origin_need_update = TRUE;
    return (TRUE);
} /* end pointer_xtal_shift() */


bool xtal_origin_goto (int iw)
{
    double old_s[3], s[3];
    char danswer[MAX_FILENAME_SIZE],*answer;
    if (!n[iw].xtal_mode)
    {
        printf ("Crystal translation is only available under Xtal mode.\n");
        return(FALSE);
    }
    V3mM3 (AX_3D[iw].x, HI, old_s);
    xterm_get_focus(iw); clear_stdin_buffer();
    sprintf (danswer, "%g %g %g", old_s[0],old_s[1],old_s[2]);
    answer = readline_gets("\nCrystal origin s0,s1,s2",danswer);
    sscanf (answer, "%lf %lf %lf", s, s+1, s+2);
    xterm_release_focus(iw);
    V3TRIM (s,s);
    if (V3EQ(old_s,s)) return(FALSE);
    V3mM3 (s, H, n[iw].xtal_origin);
    atom_xtal_origin (n[iw].xtal_origin);
    if (n[iw].bond_mode) bond_xtal_origin_update (iw);
    else n[iw].bond_xtal_origin_need_update = TRUE;
    return (TRUE);
} /* end xtal_origin_goto() */


bool xtal_origin_zero (int iw)
{
    if (!n[iw].xtal_mode)
    {
        printf ("Crystal translation is only available under Xtal mode.\n");
        return(FALSE);
    }
    if (V3LENGTH(n[iw].xtal_origin) < TINY) return(FALSE);
    V3ZERO(n[iw].xtal_origin);
    atom_xtal_origin (n[iw].xtal_origin);
    if (n[iw].bond_mode) bond_xtal_origin_update (iw);
    else n[iw].bond_xtal_origin_need_update = TRUE;
    return (TRUE);
} /* end xtal_origin_zero() */
