/*******************************************/
/* Threaded Atomistic Configuration Viewer */
/*                                         */
/* Fri Mar 31 2000 Ju Li <liju99@mit.edu>  */
/*******************************************/

#include "A.h"

static char *auxiliary_keycode[16] = {
    "Alt+0",
    "Alt+1",
    "Alt+2",
    "Alt+3",
    "Alt+4",
    "Alt+5",
    "Alt+6",
    "Alt+7",
    "Alt+8",
    "Alt+9",
    "Alt+a",
    "Alt+b",
    "Alt+c",
    "Alt+d",
    "Alt+e",
    "Alt+f"
};

void print_coordination_histogram()
{
    register int i;
    printf ("------------- Coordination Number Statistics --------------\n");
    printf ("Coord.  Count  Percentage    R     G     B    Name\n");
    for (i=0; i<=ATOM_COORDINATION_MAX; i++)
        if (coordination_hist[i] > 0)
            printf("%3d%10d   %6.2f%%    %.3f %.3f %.3f  %s\n",
                   i, coordination_hist[i], 100.*coordination_hist[i]/np,
                   ATOM_COORDINATION_COLOR[i].r, ATOM_COORDINATION_COLOR[i].g,
                   ATOM_COORDINATION_COLOR[i].b,
                   ATOM_COORDINATION_COLOR[i].name );
    printf ("-----------------------------------------------------------\n");
    printf ("average = %g, most populous = %d.\n",
            avg_coordination, coordination_crystal);
    return;
} /* end print_coordination_histogram() */


/* print info about atom "i" */
bool print_atom (int iw, int i)
{
    int j, k;
    M3 J;
    char keycode[3][3]={{'a','b','c'},{'d','e','f'},{'g','h','i'}};

    if (n[iw].suppress_printout) return(FALSE);
    printf ("\n%d atom (0-%d): %d/%s, s=[%g %g %g],\n"
            "x=[%g %g %g] A, coordination number = %d,\n",
            i, np-1, tp[i],
            (ct->Z[(int)tp[i]]>0) ?
            atom_symbol(ct->Z[(int)tp[i]]) :
            "(see table)",
            s[DIMENSION*i], s[DIMENSION*i+1], s[DIMENSION*i+2],
            B->BALL[i].x[0], B->BALL[i].x[1], B->BALL[i].x[2],
            coordination[i]);
    for (j=k=0; k<MAX_GEO_MEASURES; k++)
    {
        if (geolist[k].has_evaluated)
        {
            if (j>0) printf (", ");
            printf ("%s = %g", geolist[k].token, geo[k][i]);
            j++;
        }
    }
    if (j>0) printf (";\n");
    for (k=0; k<CONFIG_num_auxiliary; k++)
    {
        if (strcmp(CONFIG_auxiliary_name[k], "Jxx")==0)
        {
            M3ASSIGN( CONFIG_auxiliary[k  ][i],
                      CONFIG_auxiliary[k+1][i],
                      CONFIG_auxiliary[k+2][i],
                      CONFIG_auxiliary[k+3][i],
                      CONFIG_auxiliary[k+4][i],
                      CONFIG_auxiliary[k+5][i],
                      CONFIG_auxiliary[k+6][i],
                      CONFIG_auxiliary[k+7][i],
                      CONFIG_auxiliary[k+8][i],J);
            for (j=0; j<9; j++)
                keycode[j/3][j%3] = auxiliary_keycode[(k+j)%16][4];
            Mprintf ("%m[%sAlt+%3M]%c %c %c] "
                     " Local J = %3M|| %.8lf %.8lf %.8lf |",
                     (k>=16)?"CapsLockOn+":((k+8>=16)?"(CapsLockOn+)":""),
                     keycode, J);
            k+=8;
        }
        else printf ("[%s%s] %s = %g%s%s.\n",
                     (k>=16)?"CapsLock On + ":"", 
                     auxiliary_keycode[k%16],
                     CONFIG_auxiliary_name[k],
                     CONFIG_auxiliary[k][i],
                     (CONFIG_auxiliary_unit[k][0]!=EOS)?" ":"",
                     CONFIG_auxiliary_unit[k]);
    }
    return (FALSE);
} /* end print_atom() */


bool find_atom (int iw)
{
    static int last_atom = 0;
    int i;
    char question[MAX_FILENAME_SIZE],danswer[MAX_FILENAME_SIZE],*answer;
    xterm_get_focus(iw); clear_stdin_buffer();
    sprintf (question, "\nFind atom [0-%d]", np-1);
    sprintf (danswer, "%d", last_atom);
    answer = readline_gets(question,danswer);
    sscanf (answer, "%d", &i);
    xterm_release_focus(iw);
    if ((i < 0) || (i >= np)) printf("find_atom: illegal index\n");
    else
    {
        n[iw].anchor = i;
        print_atom(iw,i);
        last_atom = i;
    }
    return (FALSE);
} /* end find_atom() */


/* print info about bond "k", return its owner atom */
bool print_bond (int iw, int k)
{
    register int i,j;
    if (n[iw].suppress_printout) return(FALSE);
    for (i=0; i<np; i++) if (k < N->idx[i+1]) break;
    j = N->list[k];
    printf ("=============================================================\n");
    printf ("Bond %d (0-%d) is between\n", k, N->idx[np]-1);
    printf ("-------------------------------------------------------------\n");
    print_atom(iw,i);
    printf ("-------------------------------------------------------------\n");
    print_atom(iw,j);
    printf ("-------------------------------------------------------------\n");
    printf ("direction = [%g %g %g], L = %g A.\n",
            C->CYLINDER[k].axis[0], C->CYLINDER[k].axis[1],
            C->CYLINDER[k].axis[2], C->CYLINDER[k].axis[3]);
    printf ("=============================================================\n");
    return (FALSE);
} /* end print_bond() */


bool print_status (int iw)
{
    int i;
    double x[3], V[3][3];
    SimpleStatistics ss;
    /* xterm_get_focus(iw); */
    for (i=0; i<CONFIG_num_auxiliary; i++)
    {
        CalculateSimpleStatistics
            (np, CHARP(CONFIG_auxiliary[i]), sizeof(double),
             IOVAL_DOUBLE, &ss);
        printf("\nauxiliary[%d]=%s [%s], threshold=[%g, %g]\n", i,
               CONFIG_auxiliary_name[i], CONFIG_auxiliary_unit[i],
               n[iw].auxiliary_threshold[i][0],
               n[iw].auxiliary_threshold[i][1]);
        printf("[%g(%d), %g(%d)], avg=%g, std.dev.=%g\n",
               ss.min, ss.idx_min, ss.max, ss.idx_max,
               ss.average, ss.standard_deviation);
    }
    printf("\n======================= Status of Viewport #%d "
           "=======================\n",iw);
    V3EQV (AX_3D[iw].x, x);
    V3pr ("Viewpoint is at %M A,\n", x);
    M3inv (H, V);
    V3mM3 (AX_3D[iw].x, V, x);
    V3pr("in reduced coordinates it is %M;\n", x);
    M3EQV(AX_3D[iw].V, V); S3PR("viewport axes = %M;\n", V);
    printf ("window width = %d, height = %d pixels,\n",
            AX_size[iw].width, AX_size[iw].height);
    printf ("and conversion factor is %g pixel/radian,\n", AX_3D[iw].k);
    printf ("which converts to %g x %g degrees of field of view.\n",
            RADIAN_TO_DEGREE(2*atan(AX_size[iw].width/2/AX_3D[iw].k)),
            RADIAN_TO_DEGREE(2*atan(AX_size[iw].height/2/AX_3D[iw].k)));
    printf ("The viewport is now anchored to %s",
            (n[iw].anchor>=0)? "atom" : "hook" );
    if (n[iw].anchor >= 0) print_atom(iw,n[iw].anchor);
    else
    {
        M3inv (H, V);
        V3mM3 (n[iw].hook, V, x);
        printf("\nx = [%g %g %g] A, or s = [%g %g %g].\n", n[iw].hook[0],
               n[iw].hook[1], n[iw].hook[2], x[0], x[1], x[2]);
    }
    printf("parallel projection mode is turned %s.\n",
           n[iw].parallel_projection?"ON":"OFF");
    printf("term printout suppression is turned %s.\n",
           n[iw].suppress_printout?"ON":"OFF");
    V3pr ("background color = %M.\n", n[iw].bgcolor);
    printf ("atom r_ratio = %f, bond radius = %f A.\n",
            n[iw].atom_r_ratio, n[iw].bond_radius);
    printf("bond mode is turned %s.\n", n[iw].bond_mode?"ON":"OFF");
    printf("system average IS%s subtracted off from atomistic strains.\n",
           shear_strain_subtract_mean ? "" : "N'T");
    printf("wireframe mode is %s.\n",
           (n[iw].wireframe_mode==WIREFRAME_MODE_CONTRAST)?"CONTRAST":
           (n[iw].wireframe_mode==WIREFRAME_MODE_NONE)?"NONE":
           (n[iw].wireframe_mode==WIREFRAME_MODE_RGBO)?"RGBO":
           (n[iw].wireframe_mode==WIREFRAME_MODE_RGBK)?"RGBK":
           (n[iw].wireframe_mode==WIREFRAME_MODE_RGB)?"RGB":
           "UNKNOWN");
    if (n[iw].xtal_mode)
    {
        printf ("Xtal mode is turned ON:\n");
        V3mM3 (n[iw].xtal_origin, HI, x);
        V3TRIM (x, x);
        V3pr ("xtal_origin = %M.\n", x);
    }
    else printf ("Xtal mode is turned OFF.\n");
    printf ("color mode = %s.\n",
            (n[iw].color_mode==COLOR_MODE_NORMAL)? "NORMAL" :
            (n[iw].color_mode==COLOR_MODE_COORD)? "COORDINATION" :
            (n[iw].color_mode==COLOR_MODE_AUXILIARY)? "Auxiliary Properties" :
            (n[iw].color_mode==COLOR_MODE_SCRATCH)? "SCRATCH" : "UNKNOWN");
    if (n[iw].shell_viewer_mode)
        printf("Shell viewer auto-invoke is turned ON.\n");
    else printf("Shell viewer auto-invoke is turned OFF.\n");
    printf("s[%d]=%d surface is now seen or selected.\n",
           n[iw].last_surface_id/2, n[iw].last_surface_id%2);
    if (rcut_patching)
        printf ("Neighbor distance cutoff between %s = %g.\n",
                rcut_patch_pairname, rcut_patch[rcut_patch_item].rcut);
    printf ("rate of change = %g.\n", n[iw].delta);
    if (n[iw].color_mode==COLOR_MODE_AUXILIARY)
    {
        i = n[iw].auxiliary_idx;
        if (i < CONFIG_num_auxiliary)
            printf("auxiliary[%d] = %s [%s], threshold = [%g, %g],\n", i,
                   CONFIG_auxiliary_name[i], CONFIG_auxiliary_unit[i],
                   n[iw].auxiliary_threshold[i][0],
                   n[iw].auxiliary_threshold[i][1]);
        else printf("auxiliary = %s, threshold = [%g, %g],\n", 
                    geolist[i-CONFIG_MAX_AUXILIARY].token, 
                    n[iw].auxiliary_threshold[i][0],
                    n[iw].auxiliary_threshold[i][1]);
        CalculateSimpleStatistics
            (np, CHARP(INW(n[iw].auxiliary_idx,CONFIG_num_auxiliary) ?
                       CONFIG_auxiliary[i] : geo[i-CONFIG_MAX_AUXILIARY]),
             sizeof(double), IOVAL_DOUBLE, &ss);
        printf("[%g(%d),%g(%d)], avg=%g, std.dev.=%g,\n",
               ss.min, ss.idx_min, ss.max, ss.idx_max,
               ss.average, ss.standard_deviation);
        printf("auxiliaries' colormap = %s \"%s\".\n",
               AX_cmap_funs[n[iw].auxiliary_cmap].name,
               AX_cmap_funs[n[iw].auxiliary_cmap].description);
        printf("invisible outside auxiliary thresholds flag = %s.\n",
               n[iw].auxiliary_thresholds_saturation?"OFF":"ON");
        printf("floating auxiliary thresholds flag = %s.\n",
               n[iw].auxiliary_thresholds_rigid?"OFF":"ON");
    }
    printf ("clicked atoms = [ ");
    for (i=0; i<ATOM_STACK_SIZE; i++) printf ("%d ", n[iw].atom_stack[i]);
    printf ("];\n");
    for (i=0; i<AX_3D_MAX_FILTER_PLANE; i++)
        if (AX_V3NEZERO(AX_3D[iw].fp[i].dx))
            printf("%s fp %d: dx = [%g %g %g], s = [%g %g %g]\n",
                   (n[iw].just_activated_fp==i) ? "*" : " ",
                   i, V3E(AX_3D[iw].fp[i].dx), V3E(n[iw].fp[i].s0));
    printf("=============================================="
           "=======================\n");
    return(FALSE);
} /* end print_status() */


bool print_help (int iw)
{
    /* xterm_get_focus(iw); */
    printf("\n***********************************************************\n");
    printf("(for more goto http://alum.mit.edu/www/liju99/Graphics/A/)\n");
    printf("F1: Help;\n");
    printf("F2: (re)define color tiling blocks (terminal input);\n");
    printf("Ctrl+F2: free allocated tiling colors;\n");
    printf("F3: show previously defined color tiling blocks;\n");
    printf("F4: New a viewport;\n");
    printf("F5: Png screenshot;\n");
    printf("F6: Jpg screenshot;\n");
    printf("F7: Eps screenshot;\n");
    printf("F8: Find an atom (terminal input);\n");
    printf("F9: load in new config file;\n");
    printf("F10: reload original config file;\n");
    printf("F11: load auxiliary properties from file;\n");
    printf("F12: load atom color/radii from file;\n");
    printf("Ctrl+F12: set config list step advance (terminal input);\n");
    printf("TAB: toggle between parallel/perspective projections;\n");
    printf("A: Anchor centering;\n");
    printf("B: toggle show bond;\n");
    printf("C: clone this viewport;\n");
    printf("G: Goto position (terminal input)\n");
    printf("W: Weight (mass) centering;\n");
    printf("Q|Ctrl+C: Quit this viewport;\n");
    printf("ESC: Imprint present configuration as reference state;\n");
    printf("V: toggle auto-invoke shell viewer for screenshot;\n");
    printf("D: background color change (terminal input);\n");
    printf("X: toggle crystal (PBC) mode;\n");
    printf("I: toggle H[][] wireframe: contrast/none/rgbo/rgbk/rgb;\n");
    printf("S: status report;\n");
    printf("Ctrl+S: resize window (terminal input)\n");
    printf("U: upright viewframe;\n");
    printf("Y: run animation script;\n");
    printf("<rate of change gearbox>:\n");
    printf("[1]%.3f [2]%.3f [3]%.3f [4]%.3f [5]%.3f\n", gearbox[1],
           gearbox[2], gearbox[3], gearbox[4], gearbox[5]);
    printf("[6]%.3f [7]%.3f [8]%.3f [9]%.3f [0]%.3f;\n", gearbox[6],
           gearbox[7], gearbox[8], gearbox[9], gearbox[0]);
    printf("PageUp/PageDown: in/decrease radius ratio;\n");
    printf("K: toggle color encode coordination number;\n");
    printf("L: assign normal atomic colors;\n");
    printf("T: xtal_origin (terminal input);\n");
    printf("Z: Zero the xtal_origin;\n");
    printf("R: start/finish changing neighbor distance cutoff (t.i.);\n");
    printf(",: distance between the last two atoms clicked;\n");
    printf(".: bond angle between the last three atoms clicked;\n");
    printf("/: dihedral angle between the last four atoms clicked;\n");
    printf("Shift+Home/End: de/increase view angle;\n");
    printf("Insert/Delete: config list backtrack/advance;\n");
    printf("Ctrl+Insert/Delete: config list first/last;\n");
    printf("SHIFT+R: toggle suppressing term printout;\n");
    printf("Right/Left/Up/Down,Shift+Up/Down: rotate about axes;\n");
    printf("Ctrl+above: translate viewport;\n");
    printf("L-BTN drag to rotate viewport;\n");
    printf("Ctrl+L-BTN select atom%s; drag to translate viewport;\n",
           n[iw].bond_mode?"/bond":"");
    printf("R-BTN selects atom%s; drag (or (Ctrl+)IMWheel) to pull/push;\n",
           n[iw].bond_mode?"/bond":"");
    printf("dble (or CapsLock+) click L-BTN: print atom%s info;\n",
           n[iw].bond_mode?"/bond":"");
    printf("CapsLockOn|Meta+L-click: change an atom's color;\n");
    printf("CapsLockOn|Meta+R-click: make an atom invisible;\n");
    printf("O: assign original atomic colors;\n");
    printf ("input NEGATIVE color value -> make atom invisible.\n");
    printf("CapsLockOn|Meta+[0-9,a-f]: color-encode auxiliaries;\n");
    printf("CapsLockOn|Meta+G: color-encode von Mises strain invariant;\n");
    printf("CapsLockOn|Meta+H: color-encode central symmetry;\n");
    printf("CapsLockOn|Meta+[m-=]: change property colormap;\n");
    printf("Shift+[0-9,a-f]: toggle cutting planes / gain focus;\n");
    printf("Shift+Right/Left: advance/backtrack focused cutting plane;\n");
    printf("Shift+Ctrl+[0-9,a-f]: shift cutting plane to current hook;\n");
    printf("Shift+i: toggle focused cutting plane wireframe mode;\n");
    printf("Shift+p: flip focused cutting plane direction;\n");
    printf("Shift+CapsLockOn|Meta+[0-9,a-f]: delete cutting plane;\n");
    printf("--------------- Context Dependent Commands -----------------\n");
    if (n[iw].xtal_mode)
    {
        printf("CapsLock+above/Shift+(Ctrl+)IMWheel: shift crystal;\n");
        printf("Shift+L-BTN select atom%s; drag to shift crystal;\n",
               n[iw].bond_mode?"/bond":"");
    }
    if (n[iw].bond_mode) printf("Home/End: in/decrease bond radii;\n");
    if (rcut_patching)
        printf("Ctrl+Home/End: in/decrease neighbor distance cutoff;\n");
    if (n[iw].color_mode==COLOR_MODE_NORMAL)
    {
        printf("Ctrl+Shift+L-click: change chemical element coloring;\n");
        printf("Ctrl+Shift+R-click: make chemical element invisible;\n");
    }
    else if (n[iw].color_mode==COLOR_MODE_COORD)
    {
        printf("Ctrl+Shift+L-click: change coordination no. coloring;\n");
        printf("Ctrl+Shift+R-click: make coordination no. invisible;\n");
    }
    else if (n[iw].color_mode == COLOR_MODE_AUXILIARY)
    {
        printf("Ctrl+A: toggle invisibility outside auxiliary thresholds;\n");
        printf("Ctrl+T: toggle floating auxiliary thresholds mode;\n");
        printf("Ctrl+R: reset auxiliary thresholds;\n");
        if (n[iw].auxiliary_idx == CONFIG_MAX_AUXILIARY + GEO_SHEAR_STRAIN) 
            printf("Shift+G: toggle subtract mean from Mises "
                   "strain evaluation;\n");
        if (n[iw].auxiliary_idx == CONFIG_MAX_AUXILIARY + GEO_CENTRAL_SYMM) 
            printf("Shift+H: change central symmetry neighbormax;\n");
    }
    printf("***********************************************************\n");
    return(FALSE);
} /* end print_help() */


/* calculate dxji[]=x_j[]-x_i[] and |dxji|^2 */
void atom_pair (int j, int i, double dxji[4])
{
    double ds[3];
    V3SUB (&s[DIMENSION*j], &s[DIMENSION*i], ds);
    V3ImagE (ds);
    V3M3LENGTH2 (ds, H, dxji);
    return;
} /* end atom_pair() */


void print_atom_pair_info (int iw, int j, int i)
{
    double dxji[4];
    atom_pair (j, i, dxji);
    printf ("x[%d]-x[%d] = [%g %g %g], distance = %g A;\n",
            j, i, V3E(dxji), sqrt(dxji[3]));
    return;
} /* end print_atom_pair_info() */


/* dxkj[], dxij[] and the bond angle in radian [0,pi] */
double atom_triplet (int k, int j, int i, double dxkj[4], double dxij[4])
{
    atom_pair (k, j, dxkj);
    atom_pair (i, j, dxij);
    if ( (dxkj[3]>0) && (dxij[3]>0) )
        return( acos( V3DOT(dxkj, dxij) / sqrt(dxkj[3]) / sqrt(dxij[3]) ) );
    else return (0);
} /* end atom_triplet_angle() */


void print_atom_triplet_info (int iw, int k, int j, int i)
{
    double dxkj[4], dxij[4], angle;
    angle = atom_triplet (k, j, i, dxkj, dxij);
    print_atom_pair_info (iw, i, j);
    printf ("bond angle = %g degrees.\n", RADIAN_TO_DEGREE(angle));
    print_atom_pair_info (iw, k, j);
    return;
} /* end print_atom_triplet_info() */


void print_atom_quartet_info (int iw, int l, int k, int j, int i)
{
    double dxkj[4], dxij[4], angle, normal[4], dxlk[4], dxjk[4], dihedral;
    angle = atom_triplet (k, j, i, dxkj, dxij);
    print_atom_pair_info (iw, i, j);
    printf ("bond angle = %g degrees.\n", RADIAN_TO_DEGREE(angle));
    print_atom_pair_info (iw, k, j);
    V3CROSS (dxkj, dxij, normal);
    normal[3] = V3LENGTH2 (normal);
    angle = atom_triplet (l, k, j, dxlk, dxjk);
    printf ("bond angle = %g degrees.\n", RADIAN_TO_DEGREE(angle));
    print_atom_pair_info (iw, l, k);
    /* right-handed helix gives positive dihedral angle */
    if ( (normal[3]>0) && (dxlk[3]>0) )
        dihedral = acos( V3DOT(normal,dxlk)/sqrt(normal[3])/sqrt(dxlk[3]) ) -
            PI / 2;
    else dihedral = 0;
    printf ("dihedral angle = %g degrees.\n", RADIAN_TO_DEGREE(dihedral));
    return;
} /* end print_atom_quartet_info() */


static int *tag_atoms_in_monoclinic_filter
(V3 origin, M3 HH, double height, double xytolerance,
 int *selected, char *taglist)
{
    register int i;
    M3 HHH, HHHI;
    double dx[4], ds[3];

    M3EQV (HH, HHH);
    V3MuL (height, HHH[2]);
    M3inv (HHH, HHHI);
    selected[0] = 0;
    selected[1] = 0;
    for (i=np; i--;)
    {
        V3SUB (B->BALL[i].x, origin, dx);
        V3mM3 (dx, HHHI, ds);
        if ( XIN(ds[0],-xytolerance,1+xytolerance) &&
             XIN(ds[1],-xytolerance,1+xytolerance) &&
             XIN(ds[2],-0.5,0.5) )
        {
            if (ds[2] < 0)
            {
                selected[0]++;
                taglist[i] = 1;
            }
            else
            {
                selected[1]++;
                taglist[i] = 2;
            }
        }
        else taglist[i] = 0;
    }
    return(selected[0]+selected[1]);
} /* end tag_atoms_in_monoclinic_filter() */


static void save_dipole_indices (int *selected, char *taglist, char *fname)
{
    register int i;
    FILE *out = WOpen(fname);
    fprintf (out, "%d %d\n", selected[0], selected[1]);
    for (i=0; i<np; i++)
        if (taglist[i]==1) fprintf (out, "%d\n", i);
    for (i=0; i<np; i++)
        if (taglist[i]==2) fprintf (out, "%d\n", i);
    Zclose (out, fname);
    return;
} /* end save_dipole_indices() */


/* Save atoms selected in a monoclinic filter to a file */
void save_atoms_in_monoclinic_filter (int iw)
{
    M3 HH;
    double d0, zmargin, xytolerance, origin[3];
    char danswer[MAX_FILENAME_SIZE], *answer, fname[MAX_FILENAME_SIZE];
    char *taglist = NULL;
    int selected[2];
    
    V3SUB( B->BALL[n[iw].atom_stack[0]].x,
           B->BALL[n[iw].atom_stack[1]].x, HH[0] );
    V3SUB( B->BALL[n[iw].atom_stack[2]].x,
           B->BALL[n[iw].atom_stack[1]].x, HH[1] );
    V3CROSS (HH[0], HH[1], HH[2]);
    if (V3ISSMALL(HH[2]))
    {
        printf("The selected parallelogram is ill-conditioned\n"
               "for constructing a monoclinic filter.\n");
        return;
    }
    V3NORMALIZE (HH[2], d0);
    printf ("\"up\" is [%g %g %g]\n"
            "check it agrees with your mirror normal...\n", V3E(HH[2]));
        
    d0 = 2.5;
    zmargin = 0.01; 
    xytolerance = 0.01;
    xterm_get_focus(iw); clear_stdin_buffer();
    REALLOC (save_atoms_in_monoclinic_filter, taglist, np, char);
    while (1)
    {
        sprintf (danswer, "%g %g %g", d0, zmargin, xytolerance);
        answer = readline_gets
            ("Interplanar spacing [A]  z-margin [A]  xy-tolerance", danswer);
        sscanf(answer, "%lf %lf %lf", &d0, &zmargin, &xytolerance);
        V3ADDMUL (B->BALL[n[iw].atom_stack[1]].x, d0/2,HH[2], origin);
        if (tag_atoms_in_monoclinic_filter
            (origin,HH, d0+2*zmargin,xytolerance, selected,taglist)>0) break;
        strcpy (danswer, answer);
    }
    printf("down=%d and up=%d atoms selected in filter.\n",
           selected[0], selected[1]);
    while (1)
    {
        sprintf (danswer, "%s.idx", fbasename);
        answer = readline_gets("Save the selected atoms to", danswer);
        sscanf(answer, "%s", fname);
        if ( tested_to_be_writable(fname) )
        {
            save_dipole_indices (selected, taglist, fname);
            printf ("selected atom indices [0-%d] saved to %s.\n", np-1,fname);
            Free (taglist);
            break;
        }
        else
        {
            printf ("\n** %s: **\n", fname);
            printf ("** This file is unwritable! **\n");
        }
        strcpy(danswer, answer);
    }
    xterm_release_focus(iw);

    return;
} /* end save_atoms_in_monoclinic_filter() */
