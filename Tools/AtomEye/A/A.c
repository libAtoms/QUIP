/*******************************************/
/* Threaded Atomistic Configuration Viewer */
/*                                         */
/* Fri Mar 31 2000 Ju Li <liju99@mit.edu>  */
/*******************************************/

#include "A.h"
#include "Icon/icon.c"

/* from libAtoms: */
Aapp_Define_Config;
bool guess_to_be_PBC;
AX_Float *S=NULL;
double HI[3][3], volume, lengthscale, cm[3];
Chemtab ct[1]={{0}};
Tp *tp=NULL;
Neighborlist N[1]={{0}};
/* from libAX: */
AX_3D_Balls B[1]={{0}}; /* invisibility flag of AX_3D_Balls: red < 0 */
/* invisibility flag of AX_3D_Cylinders: radius<=0 */
AX_3D_Cylinders C[1]={{0}};
/* from this application */
char fbasename[MAX_FILENAME_SIZE];
char config_fname[MAX_FILENAME_SIZE];
const double gearbox[10]={0.15,0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5};
Navigator n[AX_MAXWIN]={{{0}}};
char xterm_identifier[XTERM_IDENTIFIER_SIZE];
Window xterm_win;
Atom_coordination_color CoordColor[ATOM_COORDINATION_MAX+1];
struct Mendeleyev Dmitri[MENDELEYEV_MAX+1];

/* Render scene to shared memory pixmap */
void paint_scene (int iw)
{
    int j;
    AX_Carrier bgcarrier;
    AX_Float HH[3][3],cr,cg,cb,ff,fr,fg,fb;
    AX_3D_Define_PW(L);  /* wireframe */
    bgcarrier = AX_Colorcarrier
        ( n[iw].bgcolor[0], n[iw].bgcolor[1], n[iw].bgcolor[2] );
    AX_bg (iw, bgcarrier);
    AX_clearzmap (iw);
    AX_3D_Balls_Zdet (iw, B);
    /* pr ("%d\n", 0); */
    if (n[iw].bond_mode) AX_3D_Cylinders_Zdet (iw, C);
    /* pr ("%d\n", 1); */
    AX_3D_Balls_Zpaper (iw, B);
    if (n[iw].bond_mode) AX_3D_Cylinders_Zpaper (iw, C);
    AX_3D_Balls_ZBprint (iw, B);
    if (n[iw].bond_mode) AX_3D_Cylinders_ZBprint (iw, C);
    AX_3D_Balls_ZAprint (iw, B);
    if (n[iw].bond_mode) AX_3D_Cylinders_ZAprint (iw, C);
    if (n[iw].wireframe_mode != WIREFRAME_MODE_NONE)
    {
        M3EQV(H,HH);
        cr = 1-n[iw].bgcolor[0];
        cg = 1-n[iw].bgcolor[1];
        cb = 1-n[iw].bgcolor[2];
        for (j=0; j<AX_3D_MAX_FILTER_PLANE; j++)
            if (V3NEZERO(AX_3D[iw].fp[j].dx) &&
                (n[iw].fp[j].wireframe_mode < FILTER_PLANE_WIREFRAME_GRADES))
            {
                ff = ((AX_Float)n[iw].fp[j].wireframe_mode) /
                    FILTER_PLANE_WIREFRAME_GRADES;
                fr = ff * n[iw].bgcolor[0] + (1-ff) * cr;
                fg = ff * n[iw].bgcolor[1] + (1-ff) * cg;
                fb = ff * n[iw].bgcolor[2] + (1-ff) * cb;
                AX_3D_Lines_Zprint
                    (iw, plane_wireframe
                     (AX_3D[iw].fp[j].dx, AX_3D[iw].fp[j].d0, fr,fg,fb));
            }
        AX_3D_PW_Assign (0,0,0, HH, cr, cg, cb, L);
        if (n[iw].wireframe_mode == WIREFRAME_MODE_RGB)
        {
            AX_3D_AssignRGB (L->LINE[0],  1,.35,.35);
            AX_3D_AssignRGB (L->LINE[3],  .35,1,.35);
            AX_3D_AssignRGB (L->LINE[8],  .45,.45,1);
            AX_3D_AssignRGB (L->LINE[1],  .35,1,.35);
            AX_3D_AssignRGB (L->LINE[2],  1,.35,.35);
            AX_3D_AssignRGB (L->LINE[4],  1,.35,.35);
            AX_3D_AssignRGB (L->LINE[5],  .35,1,.35);
            AX_3D_AssignRGB (L->LINE[6],  1,.35,.35);
            AX_3D_AssignRGB (L->LINE[7],  .35,1,.35);
            AX_3D_AssignRGB (L->LINE[9],  .45,.45,1);
            AX_3D_AssignRGB (L->LINE[10], .45,.45,1);
            AX_3D_AssignRGB (L->LINE[11], .45,.45,1);
        }
        else if (n[iw].wireframe_mode == WIREFRAME_MODE_RGBO)
        {
            AX_3D_AssignRGB (L->LINE[0],  1,.35,.35);
            AX_3D_AssignRGB (L->LINE[3],  .35,1,.35);
            AX_3D_AssignRGB (L->LINE[8],  .45,.45,1);
            AX_3D_AssignRGB (L->LINE[1],  -1,1,1);
            AX_3D_AssignRGB (L->LINE[2],  -1,1,1);
            AX_3D_AssignRGB (L->LINE[4],  -1,1,1);
            AX_3D_AssignRGB (L->LINE[5],  -1,1,1);
            AX_3D_AssignRGB (L->LINE[6],  -1,1,1);
            AX_3D_AssignRGB (L->LINE[7],  -1,1,1);
            AX_3D_AssignRGB (L->LINE[9],  -1,1,1);
            AX_3D_AssignRGB (L->LINE[10], -1,1,1);
            AX_3D_AssignRGB (L->LINE[11], -1,1,1);
        }
        else if (n[iw].wireframe_mode == WIREFRAME_MODE_RGBK)
        {
            AX_3D_AssignRGB (L->LINE[0],  1,.35,.35);
            AX_3D_AssignRGB (L->LINE[3],  .35,1,.35);
            AX_3D_AssignRGB (L->LINE[8],  .45,.45,1);
        }
        AX_3D_Lines_Zprint(iw, L);
    }
    return;
} /* end paint_scene() */


/* Manager of key/mouse events */
bool treatevent (int iw)
{
    int i=0;
    double thickness[DIMENSION];
    bool pointer_in_window;
    AXSize newsize;
    pthread_t tid;
    static int WindowUnMapped[AX_MAXWIN] = {0};
    KeySym keysym;
    pointer_in_window = AXQueryPointer(iw);
    switch (AX_event[iw].type)
    {
        case UnmapNotify:
            for (i=0; i<AX_MAXWIN; i++)
                if ((AX_cid[i]!=0) && (i!=iw) && (!WindowUnMapped[i])) break;
            if (i == AX_MAXWIN)
                XIconifyWindow(AX_display[iw], xterm_win, AX_screen[iw]);
            WindowUnMapped[iw] = 1;
            return (FALSE);
        case MapNotify:
            XMapWindow(AX_display[iw], xterm_win);
            AXRaiseWindow(iw);
            WindowUnMapped[iw] = 0;
            return (FALSE);
        case KeyPress:
            if (!pointer_in_window) return (FALSE);
            switch ( (keysym=AXKeysym0(iw)) )
            {
                case XK_F1:
                case XK_h:
                    if ( AXLOCK(iw) || AXMETA(iw) )
                    { /* central symmetry parameter */
                        treat_numeral_event(iw,10+keysym-XK_a);
                        n[iw].auxiliary_threshold
                            [n[iw].auxiliary_idx][0] = 0.00376;
                        n[iw].auxiliary_thresholds_rigid = 1;
                        color_encode_auxiliary(iw);
                        return(TRUE);
                    }
                    else if ( AXSHFT(iw) )
                        return(change_central_symm_neighbormax(iw));
                    return(print_help(iw));
                case XK_F2:
                    if (AXLOCK(iw) || AXMETA(iw) || AXCTRL(iw))
                    {
                        scratch_free (iw);
                        return (FALSE);
                    }
                    else return(rescratch(iw));
                case XK_F3:
                    n[iw].color_mode = COLOR_MODE_SCRATCH;
                    return(scratch_color(iw));
                case XK_F4:
                case XK_n:
                    n[iw].mitosis = -1;
                    pthread_create (&tid, NULL, (void *(*)(void *))
                                    thread_start, (void *)(&(n[iw].mitosis)));
                    return (FALSE);
                case XK_F5:
                case XK_p:
                    if ( AXSHFT(iw) )
                    {
                        i = n[iw].just_activated_fp;
                        if (V3NEZERO(AX_3D[iw].fp[i].dx))
                        {
                            AX_3D[iw].fp[i].d0 = -AX_3D[iw].fp[i].d0;
                            V3NeG ( AX_3D[iw].fp[i].dx );
                            return (TRUE);
                        }
                        return (FALSE);
                    }
                    else return(capture_png(iw));
                case XK_F6:
                case XK_j:
                    return(capture_jpg(iw));
                case XK_F7:
                    return(capture_eps(iw));
                case XK_F8:
                    return(find_atom(iw));
                case XK_minus:
                    if ( ( AXLOCK(iw) || AXMETA(iw) ) &&
                         ( n[iw].color_mode == COLOR_MODE_AUXILIARY ) )
                    {
                        n[iw].auxiliary_cmap =
                            (n[iw].auxiliary_cmap + AX_MAX_CMAP - 1) %
                            AX_MAX_CMAP;
                        return(color_encode_auxiliary(iw));
                    }
                    else return(treat_numeral_event(iw,6));
                    /* French keyboard: see Notes/0001.txt */
                    break;
                case XK_equal:
                    if ( ( AXLOCK(iw) || AXMETA(iw) ) &&
                         ( n[iw].color_mode == COLOR_MODE_AUXILIARY ) )
                    {
                        n[iw].auxiliary_cmap = (n[iw].auxiliary_cmap + 1) %
                            AX_MAX_CMAP;
                        return(color_encode_auxiliary(iw));
                    }
                    break;
                case XK_F9:
                    reload_config(iw,TRUE);
                    return(TRUE);
                case XK_F10:
                    reload_config(iw,FALSE);
                    return(TRUE);
                case XK_F11:
                    return (load_auxiliary_from_file(iw));
                case XK_F12:
                    if (AXCTRL(iw)) return (set_glob_advance(iw));
                    else return (load_atom_color_from_file(iw));
                case XK_Tab:
                    if (n[iw].parallel_projection)
                        return(parallel_to_perspective(iw));
                    else return(perspective_to_parallel(iw));
                case XK_g:
                    if ( AXLOCK(iw) || AXMETA(iw) )
                        return(treat_numeral_event(iw,10+keysym-XK_a));
                    else if ( AXSHFT(iw) )
                        return(change_shear_strain_subtract_mean(iw));
                    return(observer_goto(iw));
                case XK_w:
                    n[iw].anchor = -1;
                    V3EQV(cm, n[iw].hook);
                    return(look_at_the_anchor(iw));
                case XK_Escape:
                    if  ( (!ComputeLeastSquareStrain) ||
                          (strcmp(ref_fbasename, fbasename)!=0) )
                    {
                        ComputeLeastSquareStrain = 1;
                        strcpy (ref_fbasename, fbasename);
                        IsoAtomicReferenceReImprint (Config_Aapp_to_Alib, ref);
                        printf("\n\"%s\" is now reference for least-square "
                               "strain calculation.\n", ref_fbasename);
                    }
                    else
                    {    
                        ComputeLeastSquareStrain = 0;
                        IsoAtomicReferenceFree (ref);
                        LeastSquareStrain_Free();
                        printf("\nLeast-square strain calculation "
                               "turned OFF.\n");
                    }
                    return (FALSE);
                case XK_q:
                    AX_closewindow(iw);
                    pthread_exit ((void *)NULL);
                case XK_v:
                    n[iw].shell_viewer_mode = !n[iw].shell_viewer_mode;
                    if (n[iw].shell_viewer_mode)
                        printf("\nShell viewer auto-invoke is turned ON.\n");
                    else printf("\nShell viewer auto-invoke is turned OFF.\n");
                    return (FALSE);
                case XK_x:
                    n[iw].xtal_mode = !n[iw].xtal_mode;
                    return (TRUE);
                case XK_i:
                    if ( AXSHFT(iw) )
                    {
                        i = n[iw].just_activated_fp;
                        if (V3NEZERO(AX_3D[iw].fp[i].dx))
                        {
                            n[iw].fp[i].wireframe_mode
                                = ( n[iw].fp[i].wireframe_mode + 1 ) %
                                (FILTER_PLANE_WIREFRAME_GRADES + 1 );
                            return (TRUE);
                        }
                        return (FALSE);
                    }
                    else
                    {
                        n[iw].wireframe_mode = (n[iw].wireframe_mode+1) % 5;
                        return (TRUE);
                    }
                case XK_m:
                    if ( ( AXLOCK(iw) || AXMETA(iw) ) &&
                         ( n[iw].color_mode == COLOR_MODE_AUXILIARY ) )
                        return (change_auxiliary_colormap(iw));
                    return(FALSE);
                case XK_s:
                    if (AXCTRL(iw)) return(resize_window(iw));
                    return(print_status(iw));
                case XK_u:
                    M3IDENTITY(AX_3D[iw].V);
                    return(look_at_the_anchor(iw));
                case XK_y:
                    script_animate(iw);
                    return(FALSE);
                case XK_0:
                case XK_1:
                case XK_2:
                case XK_3:
                case XK_4:
                case XK_5:
                case XK_6:
                case XK_7:
                case XK_8:
                case XK_9:
                    return(treat_numeral_event(iw,keysym-XK_0));

                case XK_agrave:     return(treat_numeral_event(iw,0));
                case XK_ampersand:  return(treat_numeral_event(iw,1));
                case XK_eacute:     return(treat_numeral_event(iw,2));
                case XK_quotedbl:   return(treat_numeral_event(iw,3));
                case XK_apostrophe: return(treat_numeral_event(iw,4));
                case XK_parenleft:  return(treat_numeral_event(iw,5));
                    /* case XK_minus: return(treat_numeral_event(iw,6)); */
                case XK_egrave:     return(treat_numeral_event(iw,7));
                case XK_underscore: return(treat_numeral_event(iw,8));
                case XK_ccedilla:   return(treat_numeral_event(iw,9));
                    /* French keyboard: see Notes/0001.txt */

                case XK_a:
                    if ( AXLOCK(iw) || AXMETA(iw) || AXSHFT(iw) )
                        return(treat_numeral_event(iw,10+keysym-XK_a));
                    else if ( AXCTRL(iw) &&
                              (n[iw].color_mode == COLOR_MODE_AUXILIARY) )
                    {
                        n[iw].auxiliary_thresholds_saturation =
                            !n[iw].auxiliary_thresholds_saturation;
                        return (color_encode_auxiliary(iw));
                    }
                    else return(look_at_the_anchor(iw));
                case XK_b:
                    if ( AXLOCK(iw) || AXMETA(iw) || AXSHFT(iw) )
                        return(treat_numeral_event(iw,10+keysym-XK_a));
                    else 
                    {
                        n[iw].bond_mode = !n[iw].bond_mode;
                        if (n[iw].bond_mode)
                        {
                            if (n[iw].bond_xtal_origin_need_update)
                                bond_xtal_origin_update (iw);
                            if (n[iw].bond_atom_color_need_update)
                                bond_atom_color_update (iw);
                        }
                        return (TRUE);
                    }
                case XK_c:
                    if ( AXLOCK(iw) || AXMETA(iw) || AXSHFT(iw) )
                        return(treat_numeral_event(iw,10+keysym-XK_a));
                    else
                    {
                        if (AXCTRL(iw))
                        {
                            AX_closewindow(iw);
                            pthread_exit ((void *)NULL);
                        }
                        n[iw].mitosis = iw;
                        pthread_create (&tid, NULL, (void *(*)(void *))
                                        thread_start,
                                        (void *)(&(n[iw].mitosis)));
                        return (FALSE);
                    }
                case XK_d:
                    if ( AXLOCK(iw) || AXMETA(iw) || AXSHFT(iw) )
                        return(treat_numeral_event(iw,10+keysym-XK_a));
                    else return (bgcolor_change(iw));
                case XK_e:
                    if ( AXLOCK(iw) || AXMETA(iw) || AXSHFT(iw) )
                        return(treat_numeral_event(iw,10+keysym-XK_a));
                    else return(capture_eps(iw));
                case XK_f:
                    if ( AXLOCK(iw) || AXMETA(iw) || AXSHFT(iw) )
                        return(treat_numeral_event(iw,10+keysym-XK_a));
                    else return(find_atom(iw));
                case XK_Page_Up:
                    if ( (n[iw].color_mode == COLOR_MODE_AUXILIARY) &&
                         (AXSHFT(iw) || AXCTRL(iw)) )
                    {
                        n[iw].auxiliary_threshold[n[iw].auxiliary_idx]
                            [AXShft(iw)]
                            += ( n[iw].auxiliary_threshold
                                 [n[iw].auxiliary_idx][1] -
                                 n[iw].auxiliary_threshold
                                 [n[iw].auxiliary_idx][0] ) * n[iw].delta;
                        if (n[iw].auxiliary_idx < CONFIG_num_auxiliary)
                            printf("Thresholds of %s = [%g, %g]%s%s.\n",
                                   CONFIG_auxiliary_name[n[iw].auxiliary_idx],
                                   n[iw].auxiliary_threshold
                                   [n[iw].auxiliary_idx][0],
                                   n[iw].auxiliary_threshold
                                   [n[iw].auxiliary_idx][1],
                                   (*CONFIG_auxiliary_unit
                                    [n[iw].auxiliary_idx]==EOS)?"":" ",
                                   CONFIG_auxiliary_unit[n[iw].auxiliary_idx]);
                        else
                            printf("Thresholds of %s = [%g, %g].\n",
                                   geolist
                                   [n[iw].auxiliary_idx-
                                   CONFIG_MAX_AUXILIARY].token, 
                                   n[iw].auxiliary_threshold
                                   [n[iw].auxiliary_idx][0],
                                   n[iw].auxiliary_threshold
                                   [n[iw].auxiliary_idx][1]);
                        return(color_encode_auxiliary(iw));
                    }
                    change_atom_r_ratio(n[iw].atom_r_ratio*=(1+n[iw].delta));
                    return(TRUE);
                case XK_Page_Down:
                    if ( (n[iw].color_mode == COLOR_MODE_AUXILIARY) &&
                         (AXSHFT(iw) || AXCTRL(iw)) )
                    {
                        n[iw].auxiliary_threshold[n[iw].auxiliary_idx]
                            [AXShft(iw)]
                            -= ( n[iw].auxiliary_threshold
                                 [n[iw].auxiliary_idx][1] -
                                 n[iw].auxiliary_threshold
                                 [n[iw].auxiliary_idx][0] ) * n[iw].delta;
                        if (n[iw].auxiliary_idx < CONFIG_num_auxiliary)
                            printf("Thresholds of %s = [%g, %g]%s%s.\n",
                                   CONFIG_auxiliary_name[n[iw].auxiliary_idx],
                                   n[iw].auxiliary_threshold
                                   [n[iw].auxiliary_idx][0],
                                   n[iw].auxiliary_threshold
                                   [n[iw].auxiliary_idx][1],
                                   (*CONFIG_auxiliary_unit
                                    [n[iw].auxiliary_idx]==EOS)?"":" ",
                                   CONFIG_auxiliary_unit[n[iw].auxiliary_idx]);
                        else
                            printf("Thresholds of %s = [%g, %g].\n",
                                   geolist
                                   [n[iw].auxiliary_idx-
                                   CONFIG_MAX_AUXILIARY].token, 
                                   n[iw].auxiliary_threshold
                                   [n[iw].auxiliary_idx][0],
                                   n[iw].auxiliary_threshold
                                   [n[iw].auxiliary_idx][1]);
                        return (color_encode_auxiliary(iw));
                    }
                    change_atom_r_ratio(n[iw].atom_r_ratio/=(1+n[iw].delta));
                    return (TRUE);
                case XK_k:
                    if (n[iw].color_mode != COLOR_MODE_COORD)
                    {
                        if (n[iw].xtal_mode)
                        {
                            n[iw].last_color_mode = n[iw].color_mode;
                            return(assign_coordination_color(iw));
                        }
                        else return (FALSE);
                    }
                    else if (n[iw].last_color_mode == COLOR_MODE_NORMAL)
                        return(assign_normal_color(iw));
                    else if (n[iw].last_color_mode == COLOR_MODE_AUXILIARY)
                        return(color_encode_auxiliary(iw));
                    else if (n[iw].color_mode == COLOR_MODE_SCRATCH)
                        return(scratch_color(iw));
                    else return (FALSE);
                case XK_l:
                    return(assign_normal_color(iw));
                case XK_o:
                    return(assign_original_normal_color(iw));
                case XK_t:
                    if ( AXCTRL(iw) &&
                         (n[iw].color_mode == COLOR_MODE_AUXILIARY) )
                    {
                        n[iw].auxiliary_thresholds_rigid =
                            !n[iw].auxiliary_thresholds_rigid;
                        printf("floating auxiliary thresholds flag = %s.\n",
                               n[iw].auxiliary_thresholds_rigid?"OFF":"ON");
                        return (FALSE);
                    }
                    else return(xtal_origin_goto(iw));
                case XK_z:
                    return(xtal_origin_zero(iw));
                case XK_Home:
                    if (AXSHFT(iw))
                    {
                        AX_3D[iw].k *= (1+n[iw].delta);
                        return (TRUE);
                    }
                    else if (AXCTRL(iw) && rcut_patching)
                    {
                        M3rowthicknesses (H, thickness);
                        for (i=0; i<DIMENSION; i++)
                            if (thickness[i] < 2 *
                                (rcut_patch[rcut_patch_item].rcut +
                                 n[iw].delta))
                            {
                                printf ("losing r<rcut images may happen "
                                        "in direction %d, request ignored.\n",
                                        i);
                                return (FALSE);
                            }
                        rcut_patch[rcut_patch_item].rcut += n[iw].delta;
                        printf ("rcut(%s) = %g.\n",
                                rcut_patch_pairname,
                                rcut_patch[rcut_patch_item].rcut);
                        return (apply_rcut_patch(iw));
                    }
                    change_bond_radius(n[iw].bond_radius*=(1+n[iw].delta));
                    return (TRUE);
                case XK_End:
                    if (AXSHFT(iw))
                    {
                        AX_3D[iw].k /= (1+n[iw].delta);
                        return (TRUE);
                    }
                    else if (AXCTRL(iw) && rcut_patching)
                    {
                        if ( (rcut_patch[rcut_patch_item].rcut
                              -= n[iw].delta) < 0 )
                            rcut_patch[rcut_patch_item].rcut = 0;
                        printf ("rcut(%s) = %g.\n",
                                rcut_patch_pairname,
                                rcut_patch[rcut_patch_item].rcut);
                        return (apply_rcut_patch(iw));
                    }
                    change_bond_radius(n[iw].bond_radius/=(1+n[iw].delta));
                    return (TRUE);
                case XK_r:
                    if ( AXCTRL(iw) &&
                         (n[iw].color_mode == COLOR_MODE_AUXILIARY) )
                    {
                        i = n[iw].auxiliary_idx;
                        reset_auxiliary_threshold(iw,i);
                        printf("\nAuxiliary[%d] = %s [%s]'s thresholds have "
                               "been reset:\nThresholds now = [%g, %g].\n\n",
                               i, CONFIG_auxiliary_name[i],
                               CONFIG_auxiliary_unit[i],
                               n[iw].auxiliary_threshold[i][0],
                               n[iw].auxiliary_threshold[i][1]);
                        return (color_encode_auxiliary(iw));
                    }
                    else if (AXSHFT(iw))
                    {
                        n[iw].suppress_printout = !n[iw].suppress_printout;
                        return (FALSE);
                    }
                    else
                    {
                        if (!rcut_patching) start_rcut_patch(iw);
                        else finish_rcut_patch(iw);
                        return (FALSE);
                    };
                case XK_comma:
                    if ( (n[iw].atom_stack[0] >= 0) &&
                         (n[iw].atom_stack[1] >= 0) )
                    {
                        printf ("\n");
                        print_atom_pair_info
                            (iw, n[iw].atom_stack[0], n[iw].atom_stack[1]);
                    }
                    return (FALSE);
                case XK_period:
                    if ( (n[iw].atom_stack[0] >= 0) &&
                         (n[iw].atom_stack[1] >= 0) &&
                         (n[iw].atom_stack[2] >= 0) )
                    {
                        printf ("\n");
                        print_atom_triplet_info
                            (iw, n[iw].atom_stack[0], n[iw].atom_stack[1],
                             n[iw].atom_stack[2]);
                    }
                    return (FALSE);
                case XK_colon:
                case XK_semicolon:
                    if ( (n[iw].atom_stack[0] >= 0) &&
                         (n[iw].atom_stack[1] >= 0) &&
                         (n[iw].atom_stack[2] >= 0) )
                    {
                        printf ("\n");
                        save_atoms_in_monoclinic_filter (iw);
                    }
                    return (FALSE);
                case XK_slash:
                    if ( (n[iw].atom_stack[0] >= 0) &&
                         (n[iw].atom_stack[1] >= 0) &&
                         (n[iw].atom_stack[2] >= 0) &&
                         (n[iw].atom_stack[3] >= 0) )
                    {
                        printf ("\n");
                        print_atom_quartet_info
                            (iw, n[iw].atom_stack[0], n[iw].atom_stack[1],
                             n[iw].atom_stack[2], n[iw].atom_stack[3]);
                    }
                    return (FALSE);
                case XK_Return:
                    if (rcut_patching) finish_rcut_patch(iw);
                    return (FALSE);
                case XK_Insert:
                    if (AXCTRL(iw)) return(config_advance_first(iw));
                    return(config_advance(iw,-n[iw].glob_advance));
                case XK_Delete:
                    if (AXCTRL(iw)) return(config_advance_last(iw));
                    return(config_advance(iw,n[iw].glob_advance));
                case XK_Right:
                    if (AXCTRL(iw)) return(translate(iw, 0, n[iw].delta));
                    else if (AXLOCK(iw) && n[iw].xtal_mode)
                        return(xtal_shift(iw, 0, -n[iw].delta));
                    else if (AXSHFT(iw))
                        return(shift_filter_plane(iw,n[iw].delta/3));
                    else return(axis_rotate(iw, 1, -n[iw].delta*PI));
                case XK_Left:
                    if (AXCTRL(iw)) return(translate(iw, 0, -n[iw].delta));
                    else if (AXLOCK(iw) && n[iw].xtal_mode)
                        return(xtal_shift(iw, 0, n[iw].delta));
                    else if (AXSHFT(iw))
                        return(shift_filter_plane(iw,-n[iw].delta/3));
                    else return(axis_rotate(iw, 1, n[iw].delta*PI));
                case XK_Up:
                    if (AXCTRL(iw))
                        if (AXSHFT(iw)) return(translate(iw, 2, n[iw].delta));
                        else return(translate(iw, 1, -n[iw].delta));
                    else if (AXLOCK(iw) && n[iw].xtal_mode)
                        if (AXSHFT(iw))
                            return(xtal_shift (iw, 2, -n[iw].delta));
                        else return(xtal_shift (iw, 1, n[iw].delta));
                    else
                        if (AXSHFT(iw))
                            return(axis_rotate(iw, 2, n[iw].delta*PI));
                        else return(axis_rotate (iw, 0, -n[iw].delta*PI));
                case XK_Down:
                    if (AXCTRL(iw))
                        if (AXSHFT(iw)) return(translate(iw, 2, -n[iw].delta));
                        else return(translate(iw, 1, n[iw].delta));
                    else if (AXLOCK(iw) && n[iw].xtal_mode)
                        if (AXSHFT(iw)) return(xtal_shift(iw, 2, n[iw].delta));
                        else return(xtal_shift(iw, 1, -n[iw].delta));
                    else
                        if (AXSHFT(iw))
                            return(axis_rotate(iw, 2, -n[iw].delta*PI));
                        else return(axis_rotate (iw, 0, n[iw].delta*PI));
                default: return (FALSE);
            }
        case ButtonPress:
            if (AXPRESSBTN(iw) == 4)
            {
                if (AXSHFT(iw) && n[iw].xtal_mode)
                    if (AXCTRL(iw)) return(xtal_shift(iw, 2, -3*n[iw].delta));
                    else return(xtal_shift(iw, 2, -n[iw].delta));
                else
                    if (AXCTRL(iw)) return(advance(iw, 5*n[iw].delta));
                    else return(advance(iw, n[iw].delta));
            }
            else if (AXPRESSBTN(iw) == 5)
            {
                if (AXSHFT(iw) && n[iw].xtal_mode)
                    if (AXCTRL(iw)) return(xtal_shift(iw, 2, 3*n[iw].delta));
                    else return(xtal_shift(iw, 2, n[iw].delta));
                else
                    if (AXCTRL(iw)) return(advance(iw, -5*n[iw].delta));
                    else return(advance(iw, -n[iw].delta));
            }
            else if ( (AXPRESSBTN(iw) == 2) || (AXPRESSBTN(iw) == 3) ||
                      (AXBTNTIME(iw) - n[iw].last_button_press_time <
                       DEFAULT_DBLECLICK_IN_MS) || AXCTRL(iw) ||
                      AXLOCK(iw) || AXMETA(iw) )
            {
                i = AX_3D_Balls_Zgrab (iw, B, AX_win_x[iw], AX_win_y[iw]);
                if (i >= 0)
                {
                    atom_stack_insert (iw, i);
                    /* BwriteLong(((XKeyEvent *)&AX_event[iw])->state);cr(); */
                    if (AXCTRL(iw) && AXSHFT(iw))
                    {
                        if (n[iw].color_mode == COLOR_MODE_NORMAL)
                            return (normal_color_change(iw,tp[i],
                                                        AXPRESSBTN(iw)!=1));
                        else if (n[iw].color_mode == COLOR_MODE_COORD)
                            return(change_coordination_color
                                   (iw,coordination[i],AXPRESSBTN(iw)!=1));
                    }
                    else if ( AXLOCK(iw) || AXMETA(iw) )
                        return(atom_color_change(iw,i,AXPRESSBTN(iw)!=1));
                    else if ( (AXPRESSBTN(iw) != 1) || AXCTRL(iw) )
                    {
                        n[iw].anchor = i;
                        /* atom_stack_insert (iw, i); */
                    }
                    print_atom (iw,i);
                }
                else if (n[iw].bond_mode)
                {
                    i = AX_3D_Cylinders_Zgrab
                        (iw, C, AX_win_x[iw], AX_win_y[iw]);
                    if (i >= 0)
                    {
                        if ( AXLOCK(iw) || AXMETA(iw) )
                            return(bond_color_change(iw,i,AXPRESSBTN(iw)!=1));
                        print_bond(iw,i);
                        if ( (AXPRESSBTN(iw) != 1) || AXCTRL(iw) )
                        {
                            n[iw].anchor = -1;
                            hook_to_cylinder(i, n[iw].hook);
                        }
                    }
                }
                n[iw].lx = AX_win_x[iw];
                n[iw].ly = AX_win_y[iw];
                return (FALSE);
            }
            else
            {
                n[iw].lx = AX_win_x[iw];
                n[iw].ly = AX_win_y[iw];
                n[iw].last_button_press_time = AXBTNTIME(iw);
                return (FALSE);
            }
        case MotionNotify:
            if (AXMOTIONBTN2(iw) || AXMOTIONBTN3(iw))
                return( pointer_advance(iw, AX_win_x[iw], AX_win_y[iw]) );
            else if (AXMOTIONBTN1(iw))
            {
                if (AXSHFT(iw) && n[iw].xtal_mode)
                    return( pointer_xtal_shift
                            (iw, AX_win_x[iw], AX_win_y[iw]) );
                else if (AXCTRL(iw) || AXSHFT(iw))
                    return( pointer_translate
                            (iw, AX_win_x[iw], AX_win_y[iw]) );
                else return( pointer_rotate(iw, AX_win_x[iw], AX_win_y[iw]) );
            }
        case Expose:
            AXGetGeometry (iw, newsize);
            if ( (newsize.width  != AX_size[iw].width) ||
                 (newsize.height != AX_size[iw].height) )
            {
                /* there is a 4-byte alignment requirement for shm pixmap */
                /* in 16-bit graphics, therefore width must be even. */
                AX_size[iw] = newsize;
                if (AX_size[iw].width & 1) AX_size[iw].width++;
                AX_resizewindow(iw,False);
                n[iw].mgs_radius = DEFAULT_MGS_RATIO *
                    sqrt((double)AX_size[iw].width*AX_size[iw].height/PI);
                return (TRUE);
            }
            return (TRUE);
        default: return (FALSE);
    }
} /* end treatevent() */


void thread_start (void *icopy)
{
    int iw,i,j;
    bool nontrivial_event;
    double tmp[3];
    if (*((int *)icopy) >= 0)
    {
        if ((iw=AX_openwindow(getpid(), fbasename,
                              AX_size[*((int *)icopy)].width,
                              AX_size[*((int *)icopy)].height)) < 0)
            pthread_exit ((void *)NULL);
    }
    else
    {
        if ((iw=AX_openwindow(getpid(), fbasename,
                              DEFAULT_WIDTH, DEFAULT_HEIGHT)) < 0)
            pthread_exit ((void *)NULL);
        /* n[iw].bgcolor[0] = n[iw].bgcolor[1] = n[iw].bgcolor[2] = 0; */
        n[iw].bgcolor[0] = n[iw].bgcolor[1] = n[iw].bgcolor[2] = 1;
        n[iw].lx = AX_win_x[iw] / 2;
        n[iw].ly = AX_win_y[iw] / 2;
        n[iw].last_button_press_time = -HUGE_INT;
        n[iw].anchor = -1;
        for (i=0; i<=ATOM_STACK_SIZE; i++) n[iw].atom_stack[i] = -1;
        V3EQV (cm, n[iw].hook);
        n[iw].suppress_printout = FALSE;
        n[iw].delta = gearbox[0];
        n[iw].mgs_radius = DEFAULT_MGS_RATIO *
            sqrt((double)AX_size[iw].width*AX_size[iw].height/PI);
        n[iw].atom_r_ratio = DEFAULT_ATOM_R_RATIO;
        n[iw].color_mode = COLOR_MODE_NORMAL;
        n[iw].wireframe_mode = WIREFRAME_MODE_CONTRAST;
        n[iw].bond_mode = FALSE;
        n[iw].bond_radius = DEFAULT_BOND_RADIUS;
        n[iw].bond_atom_color_need_update = FALSE;
        /* n[iw].shell_viewer_mode = */
        /* commands_exists(NUMBER_RASTER_VIEWERS, raster_viewers) && */
        /* commands_exists(NUMBER_POSTSCRIPT_VIEWERS, postscript_viewers); */
        n[iw].shell_viewer_mode = FALSE;
        n[iw].xtal_mode = guess_to_be_PBC;
        V3ZERO (n[iw].xtal_origin);
        n[iw].bond_xtal_origin_need_update = FALSE;
        n[iw].last_surface_id = 4; /* looking at the z-bottom */
        n[iw].auxiliary_idx = 0;
        n[iw].auxiliary_cmap = AX_CMAP_JET;
        for (i=0; i<CONFIG_num_auxiliary; i++)
            reset_auxiliary_threshold(iw,i);
        for (i=0; i<MAX_GEO_MEASURES; i++)
            if (geolist[i].has_evaluated)
                reset_auxiliary_threshold(iw,CONFIG_MAX_AUXILIARY+i);
        n[iw].parallel_projection = 0;
        n[iw].glob_advance = 1;
    }
    AXSETICON(iw);
    AX_plugin_3D_module(iw);
    AX_plugin_BC_module(iw);
    if (*((int *)icopy) >= 0)
    {
        AX_3D[iw] = AX_3D[*((int *)icopy)];
        n[iw] = n[*((int *)icopy)];
    }
    else
    {
        
        tmp[0] = 0.5 * (MAX(fabs(H[0][0]), fabs(H[1][0])) + fabs(H[2][0]));
        tmp[1] = 0.5 * (MAX(fabs(H[0][1]), fabs(H[1][1])) + fabs(H[2][1]));
        tmp[2] = MAX(tmp[0], tmp[1]) /
            TAN(AX_3D_DEF_VIEW_ANGLE_IN_DEGREE*INITIAL_VIEW_RATIO/2.);
        /* V3CROSS (H[0], H[1], normal); */
        /* V3normalize (normal); */
        /* if (V3DOT(normal, H[2]) > 0) V3NeG(normal); */
        /* V3MUL(tmp[2], normal, AX_3D[iw].x); */
        /* V3ADDmuLmuL(0.5, H[0], 0.5, H[1], AX_3D[iw].x); */
        /* V3ADDmuL(); */
        /* V3ASSIGN(0.5,0.5, */
        /* -0.5*MAX(V3LENGTH(H[0]),V3LENGTH(H[1])) / V3LENGTH(H[2]) / */
        /* TAN(AX_3D_DEF_VIEW_ANGLE_IN_DEGREE*INITIAL_VIEW_RATIO/2.), */
        /* tmp); */
        V3EQV(cm, AX_3D[iw].x);
        AX_3D[iw].x[2] = -tmp[2] +
            ((H[0][2] < 0)?1:0) * H[0][2] +
            ((H[1][2] < 0)?1:0) * H[1][2] +
            ((H[2][2] < 0)?1:0) * H[2][2];
        look_at_the_anchor(iw);
        if (np < 256) perspective_to_parallel(iw);
        xterm_win = AXIdentifyWIN(iw,xterm_identifier);
        /* fprintf (stderr, "xterm_win = 0x%lx\n", (long)xterm_win); */
        XStoreName(AX_display[iw],xterm_win,AX_title[iw]);
        XSetIconName(AX_display[iw],xterm_win,AX_title[iw]);
        /* XReparentWindow(AX_display[iw],xterm_win,AX_win[iw],100,100); */
        /* XMapWindow(AX_display[iw], xterm_win); */
    }
    for (nontrivial_event = TRUE;;)
    { /* blocking is more CPU efficient */
        while (!nontrivial_event)
        {
            AXNextEvent(iw);
            nontrivial_event = treatevent (iw);
            for (j=3; ;j*=2)
            { /* preventing long queue to freeze the control */
                i = AXQLength(iw);
                if (i > 64)
                {
                    for (;i--;) AXNextEvent(iw);
                    break;
                }
                else if (i > j)
                    for (;i--;)
                    {
                        AXNextEvent(iw);
                        nontrivial_event = nontrivial_event || treatevent (iw);
                    }
                else break;
            }
        }
        paint_scene(iw);
        AX_dump(iw); AX_show(iw);
        nontrivial_event = FALSE;
    }
    return;
} /* end thread_start() */


void select_fbasename (char *raw)
{
    char *p;
    if ((p=strrchr(raw,'/'))==NULL) p = raw;
    else p++;
    strcpy (fbasename,p);
    if ((p=strstr(fbasename,".gz")))  *p = EOS;
    if ((p=strstr(fbasename,".bz2"))) *p = EOS;
    if ((p=strstr(fbasename,".z")))   *p = EOS;
    if ((p=strstr(fbasename,".pdb"))) *p = EOS;
    if ((p=strstr(fbasename,".PDB"))) *p = EOS;
    if ((p=strstr(fbasename,".cfg"))) *p = EOS;
    if ((p=strstr(fbasename,".CFG"))) *p = EOS;
    return;
} /* end select_fbasename() */


int main (int argc, char *argv[])
{
    register int i;
    int k,god;
    double tmp[3];
    char command[512], *TTYname;
    memcpy(CoordColor, ATOM_COORDINATION_COLOR,
           (ATOM_COORDINATION_MAX+1)*sizeof(Atom_coordination_color));
    memcpy(Dmitri, MENDELEYEV, (MENDELEYEV_MAX+1)*sizeof(struct Mendeleyev));
    if (argc == 2)
    {
        if (!command_exists("xterm"))
        {
            fprintf(stderr,
                    "\nAtomEye needs X terminal for some input / output,\n"
                    "please make sure \"xterm\" command can be found in\n"
                    "your PATH environment.\n\n");
            return (1);
        }
        TTYname = ttyname(fileno(stderr));
        TimeRandomize();
        RandomBase64String(XTERM_IDENTIFIER_SIZE, xterm_identifier);
        sprintf(command,"xterm -xrm %s -e %s %s %s &",
                xterm_identifier, argv[0], argv[1], xterm_identifier);
        execlp("xterm", "xterm", "-bg", "gray40", "-fg", "white",
               "-sb", "-sl", "3000", "-cr", "yellow", "-fn", "7x13",
               "-xrm", xterm_identifier,
               "-e", argv[0], argv[1], xterm_identifier, TTYname, NULL);
        return (0);
    }
    else if ((argc != 4) || (strlen(argv[2]) != XTERM_IDENTIFIER_SIZE-1))
    {
        fprintf(stderr, "Threaded atomistic configuration viewer V%s "
                "(C) Ju Li %s\nUsage: %s <PDB or CFG file>\n",
                AX_VERSION, JLCPL_YRS, argv[0]);
        return (1);
    }
    strcpy (config_fname, argv[1]);
    strcpy (xterm_identifier, argv[2]);
    redirect_stderr_to (wOpen(argv[3]));
    if (!Fexists(config_fname))
    {
        pr ("\n** %s: **\n", config_fname);
        pr ("** There is no such file! **\n");
        return(1);
    }
    if (!Freadable(config_fname))
    {
        pr ("\n** %s: **\n", config_fname);
        pr ("** This file is unreadable! **\n");
        return(1);
    }
    i = CONFIG_LOAD (config_fname, Config_Aapp_to_Alib);
    for (k=0; k<CONFIG_num_auxiliary; k++)
        if (*blank_advance(CONFIG_auxiliary_name[k])==EOS)
            sprintf(CONFIG_auxiliary_name[k], "auxiliary%d", k);
    /* more functionalities than you may need */
    guess_to_be_PBC = TRUE;  
    M3InV (H, HI, volume);
    lengthscale = cbrt(volume);
    rebind_CT (Config_Aapp_to_Alib, "", ct, &tp); cr();
    Neighborlist_Recreate_Form (Config_Aapp_to_Alib, ct, N);
    /* Everything is assumed to be under PBC except overflow */
    /* error treatment is different for PDB and CFG.         */
    /* CFG would always fold; whereas if the PDB has no      */
    /* CRYST1 tag, a bounding box H[][] would be used that   */
    /* is so large that the PBC is in reality detached. If   */
    if (i == CONFIG_CFG_LOADED)
        N->s_overflow_err_handler =
            NEIGHBORLIST_S_OVERFLOW_ERR_HANDLER_FOLD_INTO_PBC;
    else
        N->s_overflow_err_handler =
            NEIGHBORLIST_S_OVERFLOW_ERR_HANDLER_BOUNDING_BOX;
    if (i == CONFIG_CFG_LOADED)
        N->small_cell_err_handler =
            NEIGHBORLIST_SMALL_CELL_ERR_HANDLER_MULTIPLY;
    else
        N->small_cell_err_handler =
            NEIGHBORLIST_SMALL_CELL_ERR_HANDLER_NOCHECK;
    /* the PDB does have CRYST1 tag but still overflows, the */
    /* above still happens. Only when there is CRYST1 tag    */
    /* and all atoms are rigorously inside the orthogonal    */
    /* box can that really tight PBC be achieved. In effect  */
    /* the s-bounds check for PDB files is a test of whether */
    /* the author means seriously about his CRYST1 tags.     */
    Neighborlist_Recreate (Config_Aapp_to_Alib, stdout, ct, &tp, N);
    /* H[][] may have been modified in Neighborlist_Recreate() */
    M3InV (H, HI, volume);
    lengthscale = cbrt(volume);
    geo_clear_has_evaluated_flags();
    evaluate_geo_measures(); Free(s1); Free(mass);
    print_coordination_histogram(); cr();
    if (guess_to_be_PBC)
    {
        S3PR("avg. M = %M\n ", avg_M);
        printf ("avg. microscopic shear strain = %g\n", avg_shear_strain);
    }
    V3ASSIGN (0.5,0.5,0.5,tmp);
    V3mM3 (tmp,H,cm);
    rcut_patching = rcut_patch_top = 0;
    select_fbasename (config_fname);
    /* pmemusage(); */
    Config_to_3D_Balls (DEFAULT_ATOM_R_RATIO);
    /* pmemusage(); */
    Config_to_3D_Bonds (DEFAULT_BOND_RADIUS);
    pmemusage();
    god = -1;
    thread_start(&god);
    return(0);
} /* end main() */
