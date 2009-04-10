/*******************************************/
/* Threaded Atomistic Configuration Viewer */
/*                                         */
/* Fri Mar 31 2000 Ju Li <liju99@mit.edu>  */
/*******************************************/

#include "A.h"

int temporary_disable_bond = 0;


/* Allocate balls */
void Config_to_3D_Balls (double atom_r_ratio)
{
    register int i;
    AX_3D_Balls_Realloc (B, np);
    for (i=0; i<np; i++)
    {
        V3mM3 ( &(s[DIMENSION*i]), H, B->BALL[i].x );
        B->BALL[i].radius = ATOM_Radius(ct->Z[(int)tp[i]]) * atom_r_ratio;
        AX_3D_AssignRGB (B->BALL[i],
                         ATOM_Color_R(ct->Z[(int)tp[i]]),
                         ATOM_Color_G(ct->Z[(int)tp[i]]),
                         ATOM_Color_B(ct->Z[(int)tp[i]]) );
    }
    return;
} /* end Config_to_3D_Balls() */


/* Allocate bonds */
void Config_to_3D_Bonds (double bond_radius)
{
    register int i,j;
    AX_Float ds[3], DS[3];
    AX_3D_Cylinders_Realloc (C, N->idx[np]);
    for (i=0; i<np; i++)
        for (j=N->idx[i]; j<N->idx[i+1]; j++)
        {
            V3EQV ( B->BALL[i].x, C->CYLINDER[j].x0 );
            V3SUB ( &s[DIMENSION*N->list[j]], &s[DIMENSION*i], ds );
            V3IMAGE ( ds, DS );
            V3mM3 ( DS, H, C->CYLINDER[j].axis );
            AX_V3NORMALIZE ( C->CYLINDER[j].axis, C->CYLINDER[j].axis[3] );
            BONDCOLOR (i, N->list[j], j);
            if (V3NEED_IMAGE(ds))
            {
                C->CYLINDER[j].g = -1;
                C->CYLINDER[j].radius = -1;
            }
            else C->CYLINDER[j].radius = bond_radius;
        }
    return;
} /* end Config_to_3D_Bonds() */


bool change_atom_r_ratio (double atom_r_ratio)
{
    register int i;
    for (i=0; i<B->n_balls; i++)
        B->BALL[i].radius = ATOM_Radius(ct->Z[(int)tp[i]]) * atom_r_ratio;
    return (TRUE);
} /* end change_atom_r_ratio() */


/* Invisibility flag of AX_3D_Cylinders is radius <= 0 */
bool change_bond_radius (double bond_radius)
{
    register int j;
    for (j=0; j<N->idx[np]; j++)
        if (C->CYLINDER[j].radius > 0)
            C->CYLINDER[j].radius = bond_radius;
    return (TRUE);
} /* end change_bond_radius() */


/* hook to the center of the appearing cylinder part */
void hook_to_cylinder (int k, double *hook)
{
    register int i;
    register double tmp;
    for (i=0; i<np; i++)
        if (k < N->idx[i+1]) break;
    tmp = (B->BALL[i].radius + C->CYLINDER[k].axis[3] -
           B->BALL[N->list[k]].radius) / 2.;
    V3ADDMUL (C->CYLINDER[k].x0, tmp, C->CYLINDER[k].axis, hook);
    return;
} /* end hook_to_cylinder() */


/* Reload the configuration but keeping the rendering state */
void reload_config (int iw, bool term_input_filename)
{
    register int i;
    int j, k, old_np;
    char fname[MAX_FILENAME_SIZE], oldfname[MAX_FILENAME_SIZE]; 
    V3 hook_s, tmp, dx;
    char *old_symbol=NULL;
    bool incompatible_config;
    
    strcpy(oldfname, config_fname);
    if (n[iw].anchor >= 0)
    { /* the new configuration may not even have the atom */
        V3EQV (B->BALL[n[iw].anchor].x, n[iw].hook);
        n[iw].anchor = -1;
    }
    /* hook_s[] is what is kept invariant */
    V3mM3 (n[iw].hook, HI, hook_s);
    if (term_input_filename)
    {
        xterm_get_focus(iw); clear_stdin_buffer();
        strcpy(fname,readline_gets("\nLoad configuration",config_fname));
        strcpy(config_fname,fname);
        xterm_release_focus(iw);
    }
    if (!Fexists(config_fname))
    {
        printf ("\n** %s: **\n", config_fname);
        printf ("** There is no such file! **\n");
        strcpy(config_fname, oldfname);
        return;
    }
    if (!Freadable(config_fname))
    {
        printf ("\n** %s: **\n", config_fname);
        printf ("** This file is unreadable! **\n");
        strcpy(config_fname, oldfname);
        return;
    }
    cr();

    old_np = np;
    CLONE(symbol, SYMBOL_SIZE*np, char, old_symbol);
    i = CONFIG_LOAD (config_fname, Config_Aapp_to_Alib);

    for (k=0; k<CONFIG_num_auxiliary; k++)
        if (*blank_advance(CONFIG_auxiliary_name[k])==EOS)
            sprintf(CONFIG_auxiliary_name[k], "auxiliary%d", k);
    rebind_CT (Config_Aapp_to_Alib, "", ct, &tp); cr();
    Neighborlist_Recreate_Form (Config_Aapp_to_Alib, ct, N);
    if (i == CONFIG_CFG_LOADED)
        N->s_overflow_err_handler =
            NEIGHBORLIST_S_OVERFLOW_ERR_HANDLER_FOLD_INTO_PBC;
    else
        N->s_overflow_err_handler =
            NEIGHBORLIST_S_OVERFLOW_ERR_HANDLER_BOUNDING_BOX;
    N->small_cell_err_handler = NEIGHBORLIST_SMALL_CELL_ERR_HANDLER_MULTIPLY;
    for (i=0; i<ct->t; i++)
        for (j=i; j<ct->t; j++)
            for (k=0; k<rcut_patch_top; k++)
                if ( ( ( (rcut_patch[k].Zi == ct->Z[i]) &&
                         (rcut_patch[k].Zj == ct->Z[j]) ) ||
                       ( (rcut_patch[k].Zi == ct->Z[j]) &&
                         (rcut_patch[k].Zj == ct->Z[i]) ) ) )
                    NEIGHBOR_TABLE(N->rcut,ct,i,j) =
                        NEIGHBOR_TABLE(N->rcut,ct,j,i) = rcut_patch[k].rcut;
    Neighborlist_Recreate (Config_Aapp_to_Alib, stdout, ct, &tp, N);
    V3mM3 (hook_s, H, tmp);
    V3SUB (tmp, n[iw].hook, dx);
    V3EQV (tmp, n[iw].hook);
    V3AdD (dx, AX_3D[iw].x);
    M3InV (H, HI, volume);
    lengthscale = cbrt(volume);
    V3ASSIGN (0.5,0.5,0.5,tmp);
    V3mM3 (tmp, H, cm);
    geo_clear_has_evaluated_flags();
    evaluate_geo_measures(); Free(s1); Free(mass);

    if  ( ComputeLeastSquareStrain )
    {
        if (ConfigChecksum(Config_Aapp_to_Alib) != ref->checksum)
            printf ("This configuration is not isoatomic with the imprinted "
                    "reference\n%s. Least-square strain NOT calculated.\n",
                    ref_fbasename);
        else LeastSquareStrain_Append();
    }

    incompatible_config = (np != old_np) ||
        memcmp(symbol, old_symbol, SYMBOL_SIZE*MIN(np,old_np));
    Free(old_symbol);
    
    if (incompatible_config)
        Config_to_3D_Balls (n[iw].atom_r_ratio);
    else for (i=0; i<np; i++) V3mM3 ( &(s[DIMENSION*i]), H, B->BALL[i].x );

    atom_xtal_origin (n[iw].xtal_origin);
    if (!n[iw].auxiliary_thresholds_rigid)
    {
        for (i=0; i<CONFIG_num_auxiliary; i++)
            reset_auxiliary_threshold(iw,i);
        for (i=0; i<MAX_GEO_MEASURES; i++)
            if (geolist[i].has_evaluated)
                reset_auxiliary_threshold(iw,CONFIG_MAX_AUXILIARY+i);
    }
    if (!temporary_disable_bond) Config_to_3D_Bonds (n[iw].bond_radius);
    select_fbasename (config_fname);
    if ((n[iw].xtal_mode) && (n[iw].color_mode == COLOR_MODE_COORD))
        assign_coordination_color(iw);
    else if (n[iw].color_mode == COLOR_MODE_AUXILIARY)
        color_encode_auxiliary(iw);
    else if (n[iw].color_mode == COLOR_MODE_SCRATCH)
        scratch_color (iw);
    else
    {
        strcpy (AX_title[iw],fbasename);
        AXSetName (iw);
        XStoreName(AX_display[iw],xterm_win,AX_title[iw]);
        XSetIconName(AX_display[iw],xterm_win,AX_title[iw]);
        if (!temporary_disable_bond)
        {
            bond_xtal_origin_update (iw);
            bond_atom_color_update (iw);
        }
    }
    return;
} /* end reload_config() */


#define SCRIPT_LINE_SIZE  512
#define SCRIPT_LINE_CHAR (SCRIPT_LINE_SIZE-1)
/* Use script to produce jpeg frames */
void script_animate (int iw)
{
    char buf[SCRIPT_LINE_SIZE],output[SCRIPT_LINE_SIZE],*ptr;
    int i,hascontent,quality,frames;
    FILE *fp;
    glob_t globbuf;
    
    strcpy(buf, "scr_anim");
    if (!(fp=ropen(buf)))
    {
        printf ("\nAnimation script \"%s\" does not exist,\n", buf);
        xterm_get_focus(iw); clear_stdin_buffer();
        if (!strcasecmp("y", readline_gets
                        ("Do you want a default one created (y/n)?","y")))
        {
            numerically_sorted_glob (config_fname, &globbuf);
            fp=wopen(buf);
            fprintf (fp, "%d\n", AX_JPG_DEF_QUALITY);
            for (i=0; i<globbuf.gl_pathc; i++)
                fprintf (fp, "%s Jpg/%05d.jpg\n", globbuf.gl_pathv[i], i);
            globfree (&globbuf);
            fclose(fp);
            fp = ropen(buf);
        }
        else
        {
            xterm_release_focus(iw);
            return;
        }
    }
    if (!(ptr=fgets(buf,SCRIPT_LINE_SIZE,fp)))
    {
        printf ("\nThere is nothing in animation script \"%s\".\n", buf);
        fclose (fp);
        return;
    }
    for (hascontent=i=0; (buf[i]!=EOS) &&
             (ISDIGIT(buf[i]) || ISBLANK(buf[i]) || (buf[i]=='\n')); i++)
        if (ISALNUM(buf[i])) hascontent=1;
    if (!hascontent)
    {
        printf ("\nThere is no content in animation script \"%s\".\n", buf);
        fclose (fp);
        return;
    }
    if (buf[i] == EOS)
    {
        sscanf (buf, "%d", &quality);
        if ((quality<0) || (quality>100))
        {
            printf ("\nquality = %d is out of valid range ([0,100]).\n",
                    quality);
            return;
        }
        else printf ("\nquality = %d\n", quality);
        ptr = fgets(buf,SCRIPT_LINE_SIZE,fp);
    }
    else quality = AX_JPG_DEF_QUALITY;
    frames = 0;
    /* If bonds are not on now, there is no need to refresh */
    temporary_disable_bond = !n[iw].bond_mode;
    /* cylinder data structure during the rendering.        */
    while (ptr)
    {
        buf[SCRIPT_LINE_CHAR] = EOS;
        sscanf (buf, "%s %s", config_fname, output);
        reload_config (iw, FALSE);
        paint_scene(iw);
        AX_dump(iw); AX_show(iw);
        if (str_caseend_with(output,".png"))
            AX_save_pixmap_as_png
                (iw,AX_JPG_QUALITY_TO_PNG_LEVEL(quality),output);
        else if (str_caseend_with(output,".eps"))
            AX_save_pixmap_as_eps(iw,quality,output);
        else AX_save_pixmap_as_jpg(iw,quality,output);
        frames++;
        ptr = fgets(buf,SCRIPT_LINE_SIZE,fp);
    }
    fclose(fp);
    printf ("%d frames saved.\n\n", frames);
    if (temporary_disable_bond)
    {
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
        temporary_disable_bond = 0;
    }
    return;
} /* end script_animate() */
#undef SCRIPT_LINE_CHAR
#undef SCRIPT_LINE_SIZE


bool set_glob_advance (int iw)
{
    char danswer[TERMSIZE];
    xterm_get_focus(iw); clear_stdin_buffer();
    sprintf (danswer, "%d", n[iw].glob_advance);
    sscanf(readline_gets("\nFilelist step advance",danswer),
           "%d", &(n[iw].glob_advance));
    xterm_release_focus(iw);
    n[iw].glob_advance = ABS(n[iw].glob_advance);
    return (FALSE);
} /* set_glob_advance() */


bool config_advance (int iw, int how_much)
{
    char oldfname[MAX_FILENAME_SIZE];
    strcpy (oldfname, config_fname);
    Numerically_sorted_glob_advance (config_fname, how_much);
    if (!Fexists(config_fname))
    {
        printf ("\n** %s: **\n", config_fname);
        printf ("** There is no such file! **\n");
        strcpy (config_fname, oldfname);
        return (FALSE);
    }
    if (!Freadable(config_fname))
    {
        printf ("\n** %s: **\n", config_fname);
        printf ("** This file is unreadable! **\n");
        strcpy (config_fname, oldfname);
        return (FALSE);
    }
    reload_config (iw, FALSE);
    return (TRUE);
} /* end config_advance() */


bool config_advance_first (int iw)
{
    char oldfname[MAX_FILENAME_SIZE];
    strcpy (oldfname, config_fname);
    Numerically_sorted_glob_first (config_fname);
    if (!Fexists(config_fname))
    {
        printf ("\n** %s: **\n", config_fname);
        printf ("** There is no such file! **\n");
        strcpy (config_fname, oldfname);
        return (FALSE);
    }
    if (!Freadable(config_fname))
    {
        printf ("\n** %s: **\n", config_fname);
        printf ("** This file is unreadable! **\n");
        strcpy (config_fname, oldfname);
        return (FALSE);
    }
    reload_config (iw, FALSE);
    return (TRUE);
} /* end config_advance_first() */


bool config_advance_last (int iw)
{
    char oldfname[MAX_FILENAME_SIZE];
    strcpy (oldfname, config_fname);
    Numerically_sorted_glob_last (config_fname);
    if (!Fexists(config_fname))
    {
        printf ("\n** %s: **\n", config_fname);
        printf ("** There is no such file! **\n");
        strcpy (config_fname, oldfname);
        return (FALSE);
    }
    if (!Freadable(config_fname))
    {
        printf ("\n** %s: **\n", config_fname);
        printf ("** This file is unreadable! **\n");
        strcpy (config_fname, oldfname);
        return (FALSE);
    }
    reload_config (iw, FALSE);
    return (TRUE);
} /* end config_advance_last() */


#define COLOR_LINESIZE 1024
/* Load color/radii file for atoms */
bool load_atom_color_from_file (int iw)
{
    char fname[MAX_FILENAME_SIZE], buf[COLOR_LINESIZE];
    FILE *fp;
    int m, j, k, items;
    double c[4];

    xterm_get_focus(iw); clear_stdin_buffer();
    sprintf (buf, "%s.clr", fbasename);
    strcpy(fname,readline_gets("\nLoad color properties from",buf));
    xterm_release_focus(iw);
    if (Freadable(fname))
    {
        fp = ROpen(fname);
        for (m=0; ; m++)
            if (fgets(buf,COLOR_LINESIZE,fp))
            {
                if (m >= np)
                {
                    printf ("\n** %s has more rows than atoms **\n", fname);
                    return (FALSE);
                }
                items = sscanf (buf, "%lf %lf %lf %lf", c, c+1, c+2, c+3);
                if (items == 1)
                {
                    c[3] = c[0];
                    c[0] = B->BALL[m].r;
                    c[1] = B->BALL[m].g;
                    c[2] = B->BALL[m].b;
                }
                else if (items == 3)
                    c[3] = B->BALL[m].radius / n[iw].atom_r_ratio;
                if (c[0]>1) c[0]/=255;
                if (c[1]>1) c[1]/=255;
                if (c[2]>1) c[2]/=255;
                AX_3D_AssignRGB(B->BALL[m], c[0], c[1], c[2]);
                B->BALL[m].radius = c[3] * n[iw].atom_r_ratio;
            }
            else if (m < np) 
            {
                if (m == 0)
                {
                    printf ("\n** %s has no data **\n", fname);
                    return (FALSE);
                }
                else if ( ISFACTOR (m, np) )
                { /* make a bold guess */
                    for (j=m; j<np; j++)
                    {
                        k = j % m;
                        AX_3D_AssignRGB (B->BALL[j], B->BALL[k].r,
                                         B->BALL[k].g, B->BALL[k].b);
                        B->BALL[j].radius = B->BALL[k].radius;
                    }
                }
                else printf ("\n** premature ending of %s **\n", fname);
                bond_atom_color_update (iw);
                return (TRUE);
            }
            else break;
    }
    else
    {
        printf ("\n** %s: **\n", fname);
        printf ("** This file is unreadable! **\n");
        return (FALSE);
    }
    bond_atom_color_update (iw);
    return (TRUE);
} /* end load_atom_color_from_file() */
#undef COLOR_LINESIZE
