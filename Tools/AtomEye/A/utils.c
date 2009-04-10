/*******************************************/
/* Threaded Atomistic Configuration Viewer */
/*                                         */
/* Fri Mar 31 2000 Ju Li <liju99@mit.edu>  */
/*******************************************/

#include "A.h"

static int rgb_url_printed = 0;

/* Given surface_id 0-5 and sn[], d0, calculate possible abscissa  */
/* points; return 1 if there is abscissa line segment and 0 if not */
static int plane_H_surface_abscissa
(V3 sn, double d0, int surface_id, AX_3D_Line *LINE)
{
    int i, a[2], b, c, ia=0;
    double s, ss;
    AX_Float *x[2], tmp[3];
    x[0] = LINE->x0;
    x[1] = LINE->x1;
    s = surface_id / 3;
    surface_id = surface_id % 3;
    if (surface_id != 0) a[0] = 0; else a[0] = 1;
    a[1] = 3 - surface_id - a[0];
    for (i=0; i<4; i++)
    {
        ss = i / 2;
        b = i % 2;
        c = 1 - b;
        if ( (s * sn[surface_id] + ss * sn[a[b]] - d0) *
             (s * sn[surface_id] + ss * sn[a[b]] + sn[a[c]] - d0) < 0 )
        {
            x[ia][surface_id] = s;
            x[ia][a[b]] = ss;
            x[ia][a[c]] =
                - (s * sn[surface_id] + ss * sn[a[b]] - d0) / sn[a[c]];
            V3MM3 (x[ia], H, tmp);
            ia ++;
        }
        if (ia == 2) return(1);
    }
    return(0);
} /* end plane_H_surface_abscissa() */


/* Obtain wireframe representation of plane cutting H[][] box */
AX_3D_Lines *plane_wireframe
(double dx[3], double d0, AX_Float r, AX_Float g, AX_Float b)
{
    int surface_id;
    V3 sn;
    static AX_3D_Lines L[1];
    static AX_3D_Line LINE[6];
    L->n_lines = 0;
    L->LINE = LINE;
    M3mV3 (H, dx, sn);
    for (surface_id=0; surface_id<6; surface_id++)
    {
        AX_3D_AssignRGB (L->LINE[L->n_lines], r, g, b);
        L->n_lines +=
            plane_H_surface_abscissa(sn, d0, surface_id, L->LINE+L->n_lines);
    }
    return (L);
} /* end plane_wireframe() */


/* handles 0-9,a-e numeral key events */
bool treat_numeral_event (int iw, int number)
{
    char danswer[MAX_FILENAME_SIZE], buffer[MAX_FILENAME_SIZE], *answer;
    int j;
    double tmp[3], dxkj[4], dxij[4];
    if ( AXSHFT(iw) )
    { /* filter plane operations */
        if (n[iw].anchor >= 0) V3EQV (B->BALL[n[iw].anchor].x, n[iw].hook);
        if ( V3EQZERO(AX_3D[iw].fp[number].dx) )
        { /* not active but will be activated and focused */
            n[iw].just_activated_fp = number;
            if ( V3EQZERO(n[iw].fp[number].dx_cache) )
            { /* no cached plane normal */
                xterm_get_focus(iw); clear_stdin_buffer();
                if ( (n[iw].atom_stack[0] == n[iw].atom_stack[1]) &&
                     (n[iw].atom_stack[0] >= 0) &&
                     (n[iw].atom_stack[2] != n[iw].atom_stack[0]) &&
                     (n[iw].atom_stack[2] >= 0) )
                {
                    atom_pair (n[iw].atom_stack[0],
                               n[iw].atom_stack[2],
                               n[iw].fp[number].dx_input);
                    V3ADDMULMUL (0.5, &s[DIMENSION*n[iw].atom_stack[0]],
                                 0.5, &s[DIMENSION*n[iw].atom_stack[2]],
                                 n[iw].fp[number].s0);
                }
                else if ( ( IMIN(3, n[iw].atom_stack) >=0 ) &&
                          Iall_different(3, n[iw].atom_stack) )
                {
                    atom_triplet (n[iw].atom_stack[0],
                                  n[iw].atom_stack[1],
                                  n[iw].atom_stack[2],
                                  dxkj, dxij);
                    V3CROSS (dxkj, dxij, n[iw].fp[number].dx_input);
                    V3EQV (&s[DIMENSION*n[iw].atom_stack[0]],
                           n[iw].fp[number].s0);
                }
                else
                {
                    V3mM3 (n[iw].hook, HI, n[iw].fp[number].s0);
                    if ( V3EQZERO(n[iw].fp[number].dx_input) )
                        V3ASSIGN(1,1,1,n[iw].fp[number].dx_input);
                }
                sprintf (danswer, "%g %g %g %g %g %g",
                         V3E(n[iw].fp[number].dx_input),
                         V3E(n[iw].fp[number].s0));
                while (1)
                {
                    sprintf(buffer, "\nCutting plane %d's dx dy dz s0 s1 s2",
                            number);
                    answer = readline_gets(buffer, danswer);
                    j = sscanf(answer, "%lf %lf %lf %lf %lf %lf",
                               V3e(AX_3D[iw].fp[number].dx),
                               V3e(n[iw].fp[number].s0));
                    if ( V3NEZERO(AX_3D[iw].fp[number].dx) ) break;
                    strcpy(danswer, answer);
                    printf("[%g %g %g] is unacceptable as a plane normal!\n",
                           V3E(AX_3D[iw].fp[number].dx));
                }
                V3EQV( AX_3D[iw].fp[number].dx, n[iw].fp[number].dx_input );
                xterm_release_focus(iw);
                V3normalize (AX_3D[iw].fp[number].dx);
                AX_3D[iw].fp[number].d0 =
                    V3ADOT(n[iw].fp[number].s0, H, AX_3D[iw].fp[number].dx,
                           tmp);
                return(TRUE);
            }
            else
            { /* there is cached plane normal */
                if ( AXCTRL(iw) ) V3mM3 (n[iw].hook, HI, n[iw].fp[number].s0);
                else if ( AXLOCK(iw) || AXMETA(iw) )
                {
                    V3ZERO( n[iw].fp[number].dx_cache );
                    for (j=0; j<AX_3D_MAX_FILTER_PLANE; j++)
                        if ( V3NEZERO(AX_3D[iw].fp[j].dx) )
                            n[iw].just_activated_fp = j;
                    return(FALSE);
                }
                V3EQV( n[iw].fp[number].dx_cache, AX_3D[iw].fp[number].dx );
                AX_3D[iw].fp[number].d0 =
                    V3ADOT(n[iw].fp[number].s0, H, AX_3D[iw].fp[number].dx,
                           tmp);
                return(TRUE);
            }
        }
        else if (n[iw].just_activated_fp != number)
        { /* activated but not the focus */
            n[iw].just_activated_fp = number;
            if ( AXCTRL(iw) )
            {
                V3mM3 (n[iw].hook, HI, n[iw].fp[number].s0);
                AX_3D[iw].fp[number].d0 =
                    V3ADOT(n[iw].fp[number].s0, H, AX_3D[iw].fp[number].dx,
                           tmp);
                return(TRUE);
            }
            else if ( AXLOCK(iw) || AXMETA(iw) )
            {
                V3ZERO( n[iw].fp[number].dx_cache );
                V3ZERO( AX_3D[iw].fp[number].dx );
                for (j=0; j<AX_3D_MAX_FILTER_PLANE; j++)
                    if ( V3NEZERO(AX_3D[iw].fp[j].dx) )
                        n[iw].just_activated_fp = j;
                return(TRUE);
            }
            else
            { /* just gain focus */
                V3EQV( AX_3D[iw].fp[number].dx, n[iw].fp[number].dx_cache );
                V3ZERO( AX_3D[iw].fp[number].dx );
                for (j=0; j<AX_3D_MAX_FILTER_PLANE; j++)
                    if ( V3NEZERO(AX_3D[iw].fp[j].dx) )
                        n[iw].just_activated_fp = j;
                return(TRUE);
            }
        }
        else
        { /* activated and is the focus */
            if ( AXCTRL(iw) )
            {
                V3mM3 (n[iw].hook, HI, n[iw].fp[number].s0);
                AX_3D[iw].fp[number].d0 =
                    V3ADOT(n[iw].fp[number].s0, H, AX_3D[iw].fp[number].dx,
                           tmp);
                return(TRUE);
            }
            else if ( AXLOCK(iw) || AXMETA(iw) )
                V3ZERO( n[iw].fp[number].dx_cache );
            else V3EQV( AX_3D[iw].fp[number].dx,
                        n[iw].fp[number].dx_cache );
            V3ZERO( AX_3D[iw].fp[number].dx );
            for (j=0; j<AX_3D_MAX_FILTER_PLANE; j++)
                if ( V3NEZERO(AX_3D[iw].fp[j].dx) )
                    n[iw].just_activated_fp = j;
            return(TRUE);
        }
    }
    else if ( AXLOCK(iw) || AXMETA(iw) )
    {
        if ( (AXLOCK(iw) && AXMETA(iw)) || (number > 15) )
            number += 16;
        n[iw].auxiliary_idx = number;
        if ( OUW(n[iw].auxiliary_idx,CONFIG_num_auxiliary) &&
             OUW(n[iw].auxiliary_idx-CONFIG_MAX_AUXILIARY,MAX_GEO_MEASURES) )
        {
            printf ("Auxiliary %d is not available.\n", n[iw].auxiliary_idx);
            return (FALSE);
        }
        n[iw].color_mode = COLOR_MODE_AUXILIARY;
        return(color_encode_auxiliary(iw));
    }
    else if ( XIN(number,0,9) )
    {
        n[iw].delta = gearbox[number];
        return(FALSE);
    }
    return(FALSE);
} /* end treat_numeral_event() */


/* shift filter plane in the direction of its normal */
bool shift_filter_plane (int iw, double delta)
{
    V3 sn, tmp;
    int i = n[iw].just_activated_fp;
    V3mM3 (AX_3D[iw].fp[i].dx, HI, sn);
    V3normalize (sn);
    if ( AX_V3NEZERO(AX_3D[iw].fp[i].dx) )
    {
        V3ADDmuL (delta, sn, n[iw].fp[i].s0);
        AX_3D[iw].fp[i].d0 =
            V3ADOT(n[iw].fp[i].s0, H, AX_3D[iw].fp[i].dx, tmp);
        return(TRUE);
    }
    else printf ("No filter plane is currently activated **\n");
    return(FALSE);
} /* end shift_filter_plane() */


typedef struct
{
    Pixmap new_pixmap;
    Drawable old_drawable;
    GC old_gc;
    AXSize old_size;
    double scale;
} CaptureSwitch;


static void CaptureResize (int iw, int new_resolution, CaptureSwitch *cs)
{
    cs->old_size = AX_size[iw];
    cs->scale = ((double) new_resolution) /
        MAX(AX_size[iw].width, AX_size[iw].height);
    if (cs->scale == 1) return;
    AX_size[iw].width  = cs->scale * AX_size[iw].width;
    AX_size[iw].height = cs->scale * AX_size[iw].height;
    /* if (AX_size[iw].width & 1) AX_size[iw].width ^= 1; */
    /* since we bypass X-server, there is no 4-byte alignment problem */
    AX_resizewindow(iw,False);
    if (AX_image_videomode)
    { /* I do not want the large snapshot in Window */
        cs->new_pixmap = XCreatePixmap
            ( AX_display[iw],
              AX_root[iw],
              (unsigned int) AX_size[iw].width,
              (unsigned int) AX_size[iw].height,
              AX_depth[iw] );
        cs->old_drawable = AX_drawable[iw];
        cs->old_gc = AX_gc[iw];
        AX_drawable[iw] = cs->new_pixmap;
        AX_gc[iw] = XCreateGC (AX_display[iw], AX_drawable[iw], 0, NULL);
    }
    paint_scene(iw);  /* to client data-memory */
    AX_dump(iw); /* client data-memory -> AX_drawable[iw] */
    /* ... should put Xlib actions here ... */
    /* font_label_scene(iw); */  /* server Xlib -> AX_drawable[iw] */
    /* ... should put Xlib actions here ... */
    return;
} /* end CaptureResize() */


static void CaptureRecover (int iw, CaptureSwitch *cs)
{
    if (cs->scale == 1) return;
    AX_size[iw] = cs->old_size;
    if (AX_image_videomode)
    { /* no, we do not do this often. Recover Window as Drawable */
        XFreeGC (AX_display[iw], AX_gc[iw]);
        XFreePixmap (AX_display[iw], cs->new_pixmap);
        AX_drawable[iw] = cs->old_drawable;
        AX_gc[iw] = cs->old_gc;
    }
    AX_resizewindow(iw,False);
    return;
} /* end CaptureRecover() */


/* Portable Network Graphics screenshot */
bool capture_png (int iw)
{
    int new_resolution;
    char danswer[MAX_FILENAME_SIZE], *answer,
        buffer[MAX_FILENAME_SIZE]={' ', ' ', ' ', ' '}, *fname;
    CaptureSwitch cs[1]={{(Pixmap)NULL,(Pixmap)NULL,NULL}};
    xterm_get_focus(iw); clear_stdin_buffer();
    new_resolution = MAX(AX_size[iw].width, AX_size[iw].height);
    fname = buffer+4;
    while (1)
    {
        sprintf (danswer, "%s.png %d", fbasename, new_resolution);
        answer = readline_gets("\nSave screen on", danswer);
        sscanf(answer, "%s %d", fname, &new_resolution);
        if (new_resolution > 1) break;
        strcpy(danswer, answer);
        printf("new resolution=%d <= 1!\n", new_resolution);
    }
    if (strcasecmp(eos(fname)-4,".png")) strcat(fname, ".png");
    if (n[iw].color_mode == COLOR_MODE_AUXILIARY)
        save_auxiliary_colormap(iw,str2(fname, ".cmap.eps"));
    xterm_release_focus(iw);
    if ( tested_to_be_writable(fname) )
    {
        CaptureResize(iw, new_resolution, cs);
        printf ("image saved on \"%s\" (%d bytes)\n", fname,
                AX_save_pixmap_as_PNG (iw,fname));
        AXQueryPointer (iw);
        sprintf (fname, "-geometry +%d+%d %s",
                 ZERO_RAMP(AX_root_x[iw]-AX_size[iw].width/2),
                 ZERO_RAMP(AX_root_y[iw]-AX_size[iw].height/4),
                 absolute_pathname(fname));
        CaptureRecover (iw, cs);
        if (n[iw].shell_viewer_mode)
            try_to_runbg (NUMBER_RASTER_VIEWERS, raster_viewers, fname);
    }
    else
    {
        printf ("\n** %s: **\n", fname);
        printf ("** This file is unwritable! **\n");
    }
    return (FALSE);
} /* end capture_png() */


/* Joint Photographic Experts Group screenshot */
bool capture_jpg (int iw)
{
    int new_resolution;
    char danswer[MAX_FILENAME_SIZE], *answer,
        buffer[MAX_FILENAME_SIZE]={' ', ' ', ' ', ' '}, *fname;
    CaptureSwitch cs[1]={{(Pixmap)NULL,(Pixmap)NULL,NULL}};
    xterm_get_focus(iw); clear_stdin_buffer();
    new_resolution = MAX(AX_size[iw].width, AX_size[iw].height);
    fname = buffer+4;
    while (1)
    {
        sprintf (danswer, "%s.jpg %d", fbasename, new_resolution);
        answer = readline_gets("\nSave screen on", danswer);
        sscanf(answer, "%s %d", fname, &new_resolution);
        if (new_resolution > 1) break;
        strcpy(danswer, answer);
        printf("new resolution=%d <= 1!\n", new_resolution);
    }
    if (strcasecmp(eos(fname)-4,".jpg")) strcat(fname, ".jpg");
    if (n[iw].color_mode == COLOR_MODE_AUXILIARY)
        save_auxiliary_colormap(iw,str2(fname, ".cmap.eps"));
    xterm_release_focus(iw);
    if ( tested_to_be_writable(fname) )
    {
        CaptureResize(iw, new_resolution, cs);
        printf ("image saved on \"%s\" (%d bytes)\n", fname,
                AX_save_pixmap_as_JPG (iw,fname));
        AXQueryPointer (iw);
        sprintf (fname, "-geometry +%d+%d %s",
                 ZERO_RAMP(AX_root_x[iw]-AX_size[iw].width/2),
                 ZERO_RAMP(AX_root_y[iw]-AX_size[iw].height/4),
                 absolute_pathname(fname));
        CaptureRecover (iw, cs);
        if (n[iw].shell_viewer_mode)
            try_to_runbg (NUMBER_RASTER_VIEWERS, raster_viewers, fname);
    }
    else
    {
        printf ("\n** %s: **\n", fname);
        printf ("** This file is unwritable! **\n");
    }
    return (FALSE);
} /* end capture_jpg() */


/* Encapsulated PostScript screenshot */
bool capture_eps (int iw)
{
    int new_resolution;
    char danswer[MAX_FILENAME_SIZE], *answer,
        buffer[MAX_FILENAME_SIZE]={' ', ' ', ' ', ' '}, *fname;
    CaptureSwitch cs[1]={{(Pixmap)NULL,(Pixmap)NULL,NULL}};
    
    xterm_get_focus(iw); clear_stdin_buffer();
    new_resolution = AX_MAXWIDTH;
    fname = buffer+4;
    while (1)
    {
        sprintf (danswer, "%s.eps %d", fbasename, new_resolution);
        answer = readline_gets("\nSave screen on", danswer);
        sscanf(answer, "%s %d", fname, &new_resolution);
        if (new_resolution > 1) break;
        strcpy(danswer, answer);
        printf("new resolution=%d <= 1!\n", new_resolution);
    }
    if (strcasecmp(eos(fname)-4,".eps")) strcat(fname, ".eps");
    if (n[iw].color_mode == COLOR_MODE_AUXILIARY)
        save_auxiliary_colormap(iw,str2(fname, ".cmap.eps"));
    xterm_release_focus(iw);
    if ( tested_to_be_writable(fname) )
    {
        CaptureResize(iw, new_resolution, cs);
        printf ("image saved on \"%s\" (%d bytes)\n", fname,
                AX_save_pixmap_as_EPS (iw,fname));
        AXQueryPointer (iw);
        sprintf (fname, "-geometry %dx%d+%d+%d %s",
                 2*cs->old_size.width, 2*cs->old_size.height,
                 ZERO_RAMP(AX_root_x[iw]-cs->old_size.width),
                 ZERO_RAMP(AX_root_y[iw]-cs->old_size.height/2),
                 absolute_pathname(fname));
        CaptureRecover (iw, cs);
        if (n[iw].shell_viewer_mode)
            try_to_runbg(NUMBER_POSTSCRIPT_VIEWERS,postscript_viewers,fname);
    }
    else
    {
        printf ("\n** %s: **\n", fname);
        printf ("** This file is unwritable! **\n");
        return (FALSE);
    }
    return (TRUE);
} /* end capture_eps() */


/* Shift viewpoint so the anchor appears at the center */
bool look_at_the_anchor (int iw)
{
    double s[3], tmp[3];
    if (n[iw].anchor >= 0)
        V3SUB(B->BALL[n[iw].anchor].x, AX_3D[iw].x, tmp);
    else V3SUB(n[iw].hook, AX_3D[iw].x, tmp);
    M3mV3 (AX_3D[iw].V, tmp, s);
    if ( (s[0]==0) && (s[1]==0) && (s[2]>=0) ) return(FALSE);
    V3ADDmuL (s[0], AX_3D[iw].V[0], AX_3D[iw].x);
    V3ADDmuL (s[1], AX_3D[iw].V[1], AX_3D[iw].x);
    if (s[2] < 0)
    { /* to maintain determinant > 0 */
        V3NEG(AX_3D[iw].V[1], AX_3D[iw].V[1]);
        V3NEG(AX_3D[iw].V[2], AX_3D[iw].V[2]);
    }
    return(TRUE);
} /* end look_at_the_anchor() */


bool observer_goto (int iw)
{
    double s[3];
    char danswer[MAX_FILENAME_SIZE],*answer;
    V3mM3 (AX_3D[iw].x, HI, s);
    xterm_get_focus(iw); clear_stdin_buffer();
    sprintf (danswer, "%g %g %g", s[0],s[1],s[2]);
    answer = readline_gets("\nObserver goto",danswer);
    xterm_release_focus(iw);
    sscanf(answer, "%lf %lf %lf", s, s+1, s+2);
    V3mM3 (s, H, AX_3D[iw].x);
    return (TRUE);
} /* end observer_goto() */


bool resize_window (int iw)
{
    char danswer[MAX_FILENAME_SIZE],*answer;
    int new_width, new_height;
    xterm_get_focus(iw); clear_stdin_buffer();
    sprintf (danswer, "%d %d", AX_size[iw].width, AX_size[iw].height);
    answer = readline_gets("\nInput new window width height",danswer);
    xterm_release_focus(iw);
    sscanf(answer, "%d %d", &new_width, &new_height);
    if ( (new_width == 0) || (new_width > AX_MAXWIDTH) )
    {
        printf("width = %d, set to AX_MAXWIDTH = %d\n",
               new_width, AX_MAXWIDTH);
        new_width = AX_MAXWIDTH;
    }
    if ( (new_height==0) || (new_height > AX_MAXHEIGHT) )
    {
        printf("height = %d, set to AX_MAXHEIGHT = %d\n",
               new_height, AX_MAXHEIGHT);
        new_height = AX_MAXHEIGHT;
    }
    if ( (new_width < 0) || (new_height < 0) ||
         ( (new_width ==  AX_size[iw].width) &&
           (new_height == AX_size[iw].height) ) ) return (FALSE);
    AX_size[iw].width = new_width;
    AX_size[iw].height = new_height;
    AX_resizewindow(iw,TRUE);
    return (TRUE);
} /* end resize_window() */


/* change the background color */
bool bgcolor_change (int iw)
{
    double c[3];
    char danswer[MAX_FILENAME_SIZE],*answer;
    xterm_get_focus(iw); clear_stdin_buffer();
    if (!rgb_url_printed)
    {
        printf("\nColor choices at:\n%s\n", RGB_URL);
        rgb_url_printed = 1;
    }
    sprintf (danswer, "%.3f %.3f %.3f",
             n[iw].bgcolor[0], n[iw].bgcolor[1], n[iw].bgcolor[2]);
    answer = readline_gets("\nChange background color",danswer);
    sscanf (answer, "%lf %lf %lf", c, c+1, c+2);
    if (c[0]>1) c[0]/=255;
    if (c[1]>1) c[1]/=255;
    if (c[2]>1) c[2]/=255;
    xterm_release_focus(iw);
    if ( V3NE(c, n[iw].bgcolor) )
    {
        V3EQV(c, n[iw].bgcolor);
        return (TRUE);
    }
    else return (FALSE);
} /* end bgcolor_change() */


void bond_atom_color_update (int iw)
{
    register int i,j;
    if (temporary_disable_bond) return;
    for (i=0; i<np; i++)
        for (j=N->idx[i]; j<N->idx[i+1]; j++)
            if ( (B->BALL[i].r >= 0) && (B->BALL[N->list[j]].r >= 0) )
            {
                C->CYLINDER[j].r = BOND_R(i, N->list[j]);
                C->CYLINDER[j].b = BOND_B(i, N->list[j]);
                if (C->CYLINDER[j].g >= 0)
                {
                    C->CYLINDER[j].g = BOND_G(i, N->list[j]);
                    C->CYLINDER[j].radius = n[iw].bond_radius;
                }
            }
            else
            { /* r<0 means atomic invisibility */
                C->CYLINDER[j].r = -1;
                C->CYLINDER[j].radius = -1;
            }
    n[iw].bond_atom_color_need_update = FALSE;
    return;
} /* end bond_atom_color_update() */


/* depends on atom_r[], atom_g[], atom_b[] */
bool assign_normal_color (int iw)
{
    register int i;
    strcpy(AX_title[iw],fbasename);
    AXSetName(iw);
    XStoreName(AX_display[iw],xterm_win,AX_title[iw]);
    XSetIconName(AX_display[iw],xterm_win,AX_title[iw]);
    for (i=0; i<np; i++)
    {
        B->BALL[i].radius = ATOM_Radius(ct->Z[(int)tp[i]])
            * n[iw].atom_r_ratio;
        AX_3D_AssignRGB (B->BALL[i],
                         ATOM_Color_R(ct->Z[(int)tp[i]]),
                         ATOM_Color_G(ct->Z[(int)tp[i]]),
                         ATOM_Color_B(ct->Z[(int)tp[i]]) );
    }
    n[iw].color_mode = COLOR_MODE_NORMAL;
    if (n[iw].bond_mode) bond_atom_color_update(iw);
    else n[iw].bond_atom_color_need_update = TRUE;
    return (TRUE);
} /* end assign_normal_color() */


bool color_encode_auxiliary (int iw)
{
    register int i;
    double x,r,g,b,*measure;
    if ( OUW(n[iw].auxiliary_idx,CONFIG_num_auxiliary) &&
         OUW(n[iw].auxiliary_idx-CONFIG_MAX_AUXILIARY,MAX_GEO_MEASURES) )
    {
        printf ("Auxiliary %d is not available.\n", n[iw].auxiliary_idx);
        return (FALSE);
    }
    n[iw].color_mode = COLOR_MODE_AUXILIARY;
    strcpy(AX_title[iw],
           str4(fbasename," (",
                OUW(n[iw].auxiliary_idx,CONFIG_num_auxiliary) ? 
                geolist[n[iw].auxiliary_idx-CONFIG_MAX_AUXILIARY].token : 
                CONFIG_auxiliary_name[n[iw].auxiliary_idx], ")"));
    AXSetName(iw);
    XStoreName(AX_display[iw],xterm_win,AX_title[iw]);
    XSetIconName(AX_display[iw],xterm_win,AX_title[iw]);
    if ( OUW(n[iw].auxiliary_idx,CONFIG_num_auxiliary) )
    {
        if (!geolist[n[iw].auxiliary_idx-CONFIG_MAX_AUXILIARY].has_evaluated)
        {
            geo_set_should_evaluate_flag
                (n[iw].auxiliary_idx-CONFIG_MAX_AUXILIARY);
            evaluate_geo_measures();
            reset_auxiliary_threshold(iw,n[iw].auxiliary_idx);
        }
        measure = geo[n[iw].auxiliary_idx-CONFIG_MAX_AUXILIARY];
    }
    else measure = CONFIG_auxiliary[n[iw].auxiliary_idx];
    for (i=0; i<np; i++)
    {
        x = measure[i];
        B->BALL[i].radius = ATOM_Radius(ct->Z[(int)tp[i]])
            * n[iw].atom_r_ratio;
        if ( (n[iw].auxiliary_thresholds_saturation == 0) &&
             ( (x > n[iw].auxiliary_threshold[n[iw].auxiliary_idx][1]) ||
               (x < n[iw].auxiliary_threshold[n[iw].auxiliary_idx][0]) ) )
        {
            B->BALL[i].r = -1;
            B->BALL[i].g = 0;
            B->BALL[i].b = 0;
            continue;
        }
        if (x > n[iw].auxiliary_threshold[n[iw].auxiliary_idx][1])
            x = n[iw].auxiliary_threshold[n[iw].auxiliary_idx][1];
        if (x < n[iw].auxiliary_threshold[n[iw].auxiliary_idx][0])
            x = n[iw].auxiliary_threshold[n[iw].auxiliary_idx][0];
        if ( n[iw].auxiliary_threshold[n[iw].auxiliary_idx][1] >
             n[iw].auxiliary_threshold[n[iw].auxiliary_idx][0] )
            x = (x - n[iw].auxiliary_threshold[n[iw].auxiliary_idx][0]) /
                (n[iw].auxiliary_threshold[n[iw].auxiliary_idx][1] -
                 n[iw].auxiliary_threshold[n[iw].auxiliary_idx][0]);
        else x = 0.5;
        AX_cmap (n[iw].auxiliary_cmap, x, &r, &g, &b);
        AX_3D_AssignRGB(B->BALL[i], r, g, b);
    }
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
    return(TRUE);
} /* end color_encode_auxiliary() */


#define AUXILIARY_LINESIZE 1024
/* Load auxiliary file for color encoding purposes */
bool load_auxiliary_from_file (int iw)
{
    char fname[MAX_FILENAME_SIZE],buf[AUXILIARY_LINESIZE];
    FILE *fp;
    int i,j,k,m;
    SimpleStatistics ss;
    xterm_get_focus(iw); clear_stdin_buffer();
    sprintf (buf, "%s.aux", fbasename);
    strcpy(fname,readline_gets("\nLoad auxiliary properties from",buf));
    xterm_release_focus(iw);
    if (Freadable(fname))
    {
        fp = ropen(fname);
        CONFIG_num_auxiliary = CONFIG_MAX_AUXILIARY;
        for (k=0; k<CONFIG_num_auxiliary; k++)
            REALLOC (load_auxiliary_from_file,
                     CONFIG_auxiliary[k], np, double);
        for (m=0; ; m++)
            if (fgets(buf,AUXILIARY_LINESIZE,fp))
            {
                if (m >= np)
                {
                    printf ("\n** %s has more rows than atoms **\n", fname);
                    goto safe_exit;
                }
                for (i=j=0; ;)
                {
                    while ( ISNOTDIGIT(buf[i]) && (buf[i]!='.') &&
                            (buf[i]!=EOS) ) i++;
                    if (buf[i] == EOS)
                    {
                        if (m==0)
                        {
                            CONFIG_num_auxiliary = j;
                            printf ("number of auxiliaries found = %d.\n", j);
                            if (j==0) goto safe_exit;
                            for (k=j; k<CONFIG_MAX_AUXILIARY; k++)
                                Free(CONFIG_auxiliary[k]);
                            for (k=0; k<CONFIG_num_auxiliary; k++)
                            {
                                sprintf(CONFIG_auxiliary_name[k],
                                        "auxiliary%d", k);
                                CONFIG_auxiliary_unit[k][0] = EOS;
                            }
                        }
                        else if (j != CONFIG_num_auxiliary)
                        {
                            printf ("\n** %s corrupted **\n", fname);
                            goto safe_exit;
                        }
                        break;
                    }
                    for (k=i; (buf[k]!='\n') & (buf[k]!=EOS) & (buf[k]!=' ');
                         k++);
                    if (k >= AUXILIARY_LINESIZE-1)
                    {
                        printf ("\nload_auxiliary_from_file: %s"
                                " line too long\n", fname);
                        goto safe_exit;
                    }
                    if (j >= CONFIG_MAX_AUXILIARY)
                    {
                        printf ("\nload_auxiliary_from_file: number of "
                                "entries > CONFIG_MAX_AUXILIARY=%d\n", 
                                CONFIG_MAX_AUXILIARY);
                        goto safe_exit;
                    }
                    if (buf[k]==EOS)
                    {
                        CONFIG_auxiliary[j++][m] = atof(&buf[i]);
                        i = k;
                    }
                    else
                    {
                        buf[k] = EOS;
                        CONFIG_auxiliary[j++][m] = atof(&buf[i]);
                        i = k+1;
                    }
                }
            }
            else if (m < np) 
            {
                printf ("\n** premature ending of %s **\n", fname);
                goto safe_exit;
            }
            else break;
    }
    else
    {
        printf ("\n** %s: **\n", fname);
        printf ("** This file is unreadable! **\n");
        return (FALSE);
    }
    for (i=0; i<CONFIG_num_auxiliary; i++)
    {
        CalculateSimpleStatistics
            (np, CHARP(CONFIG_auxiliary[i]), sizeof(double),
             IOVAL_DOUBLE, &ss);
        n[iw].auxiliary_threshold[i][0] = ss.min;
        n[iw].auxiliary_threshold[i][1] = ss.max;
    }
    n[iw].color_mode = COLOR_MODE_AUXILIARY;
    if (OUW(n[iw].auxiliary_idx,CONFIG_num_auxiliary))
        n[iw].auxiliary_idx=0;
    return (color_encode_auxiliary(iw));
  safe_exit:
    printf ("** all auxiliary properties freed **\n");
    Config_free_auxiliary();
    return (FALSE);
} /* end load_auxiliary_from_file() */


bool change_auxiliary_colormap (int iw)
{
    int i;
    char buf[TERMSIZE];
    xterm_get_focus(iw); clear_stdin_buffer();
    do
    {
        printf("\nPlease choose a colormap index from:\n");
        for (i=0; i<AX_MAX_CMAP; i++)
            printf ("%2d: %s: %s\n",
                    i, AX_cmap_funs[i].name, AX_cmap_funs[i].description);
        sprintf (buf, "%d", n[iw].auxiliary_cmap);
        i = atoi(readline_gets("Colormap index",buf));
    } while (OUW(i,AX_MAX_CMAP));
    xterm_release_focus(iw);
    if ( (i != n[iw].auxiliary_cmap) ||
         (n[iw].color_mode != COLOR_MODE_AUXILIARY) )
    {
        n[iw].auxiliary_cmap = i;
        n[iw].color_mode = COLOR_MODE_AUXILIARY;
        return (color_encode_auxiliary(iw));
    }
    return(FALSE);
} /* end change_auxiliary_colormap() */


#define AUX_COLORMAP_GRADES      320
#define AUX_COLORMAP_BARWIDTH    24
#define AUX_COLORMAP_SEPARATION  4
#define AUX_COLORMAP_FONTWIDTH   72
#define AUX_COLORMAP_FONTRATIO   (1/GOLDEN_RATIO)
#define AUX_COLORMAP_FONTFMT     "%.4g"
#define AUX_COLORMAP_WIDTH       \
  (AUX_COLORMAP_BARWIDTH + AUX_COLORMAP_SEPARATION + AUX_COLORMAP_FONTWIDTH)
#define AUX_COLORMAP_HEIGHT     640

void save_auxiliary_colormap (int iw, char *default_cmap_fname)
{
    int i=0,j;
    char fname[MAX_FILENAME_SIZE],buf[4][TERMSIZE];
    double x,r,g,b,grade,width,height;
    FILE *fp;
    SimpleStatistics ss;
    clear_stdin_buffer();
    if (default_cmap_fname == NULL)
        default_cmap_fname = DEFAULT_CMAP_FNAME;
    strcpy (fname, readline_gets("Save colormap on",default_cmap_fname));
    if (tested_to_be_writable(fname))
    {
        /* pr("here %s %g %g\n", fname, */
        /* n[iw].auxiliary_threshold[n[iw].auxiliary_idx][0], */
        /* n[iw].auxiliary_threshold[n[iw].auxiliary_idx][1]); */
        fp = wopen(fname);
        fprintf (fp, "%%!PS-Adobe-2.0 EPSF-1.2\n");
        fprintf (fp, "%%%%BoundingBox: %d %d %d %d\n", 0, 0,
                 AUX_COLORMAP_WIDTH, AUX_COLORMAP_HEIGHT);
        fprintf (fp, "/b { 0 0 moveto exch dup 0 rlineto exch 0 exch \n"
                 "rlineto neg 0 rlineto closepath } bind def\n");
        sprintf (buf[0], AUX_COLORMAP_FONTFMT,
                 n[iw].auxiliary_threshold[n[iw].auxiliary_idx][0]);
        sprintf (buf[1], AUX_COLORMAP_FONTFMT,
                 n[iw].auxiliary_threshold[n[iw].auxiliary_idx][1]);
        CalculateSimpleStatistics
            (np, CHARP(INW(n[iw].auxiliary_idx,CONFIG_num_auxiliary) ?
                       CONFIG_auxiliary[n[iw].auxiliary_idx] :
                       geo[n[iw].auxiliary_idx-CONFIG_MAX_AUXILIARY]),
             sizeof(double), IOVAL_DOUBLE, &ss);
        if (ss.min > 0.9 * n[iw].auxiliary_threshold[n[iw].auxiliary_idx][0] +
            0.1 * n[iw].auxiliary_threshold[n[iw].auxiliary_idx][1])
            sprintf (buf[2], AUX_COLORMAP_FONTFMT, ss.min);
        else buf[2][0] = EOS;
        if (ss.max < 0.1 * n[iw].auxiliary_threshold[n[iw].auxiliary_idx][0] +
            0.9 * n[iw].auxiliary_threshold[n[iw].auxiliary_idx][1])
            sprintf (buf[3], AUX_COLORMAP_FONTFMT, ss.max);
        else buf[3][0] = EOS;
        j = 0;
        for (i=0; i<4; i++)
            if (strlen(buf[i]) > j) j = strlen(buf[i]);
        grade = DOUBLE(AUX_COLORMAP_HEIGHT) / AUX_COLORMAP_GRADES;
        width = AUX_COLORMAP_FONTWIDTH / j;
        height = width * AUX_COLORMAP_FONTRATIO;
        for (i=0; i<AUX_COLORMAP_GRADES; i++)
        {
            x = (0.5 + i) / AUX_COLORMAP_GRADES;
            AX_cmap (n[iw].auxiliary_cmap, x, &r, &g, &b);
            fprintf (fp, "%.7g %.7g %.7g setrgbcolor\n", r,g,b);
            fprintf (fp, "%.7g %.7g b fill\n", DOUBLE(AUX_COLORMAP_BARWIDTH),
                     grade);
            fprintf (fp, "0 %.7g translate\n", grade);
        }
        fprintf (fp, "0 0 0 setrgbcolor\n");
        fprintf (fp, "/Helvetica findfont [%.7g 0 0 %.7g 0 0] "
                 "makefont setfont\n", height, height);
        fprintf (fp,"%.7g %.7g translate 0 0 moveto (%s) true charpath fill\n",
                 DOUBLE(AUX_COLORMAP_BARWIDTH + AUX_COLORMAP_SEPARATION),
                 2*grade-height, buf[1]);
        fprintf (fp, "0 %.7g translate 0 0 moveto (%s) true charpath fill\n",
                 height-DOUBLE(AUX_COLORMAP_HEIGHT), buf[0]);
        if (buf[2])
            fprintf (fp, "0 %.7g translate 0 0 moveto (%s) true charpath "
                     "fill\n",
                     (ss.min-
                      n[iw].auxiliary_threshold[n[iw].auxiliary_idx][0]) /
                     (n[iw].auxiliary_threshold[n[iw].auxiliary_idx][1] -
                      n[iw].auxiliary_threshold[n[iw].auxiliary_idx][0]) * 
                     AUX_COLORMAP_HEIGHT, buf[2]);
        if (buf[3])
            fprintf (fp, "0 %.7g translate 0 0 moveto (%s) true charpath "
                     "fill\n", (ss.max - ss.min) /
                     (n[iw].auxiliary_threshold[n[iw].auxiliary_idx][1] -
                      n[iw].auxiliary_threshold[n[iw].auxiliary_idx][0]) *
                     AUX_COLORMAP_HEIGHT, buf[3]);
        fprintf (fp, "showpage\n");
        fclose(fp);
        printf ("auxiliary colormap saved on \"%s\" (%d bytes)\n", fname,
                (int)Fsize(fname));
    }
    else
    {
        printf ("\n** %s: **\n", fname);
        printf ("** This file is unwritable! **\n");
    }
    return;
} /* end save_auxiliary_colormap() */


bool assign_coordination_color (int iw)
{
    register int i;
    if (!n[iw].xtal_mode)
    {
        printf ("Coordination color works only in Xtal mode.\n");
        return (FALSE);
    }
    strcpy(AX_title[iw], str4(fbasename, " (", "coord.", ")"));
    AXSetName(iw);
    XStoreName(AX_display[iw],xterm_win,AX_title[iw]);
    XSetIconName(AX_display[iw],xterm_win,AX_title[iw]);
    for (i=0; i<np; i++)
    {
        AX_3D_AssignRGB(B->BALL[i], CoordColor[coordination[i]].r,
                        CoordColor[coordination[i]].g,
                        CoordColor[coordination[i]].b);
        B->BALL[i].radius = ATOM_Radius(ct->Z[(int)tp[i]])
            * n[iw].atom_r_ratio;
    }
    n[iw].color_mode = COLOR_MODE_COORD;
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
} /* end assign_coordination_color() */


/* direct change of bond color and/or radius */
bool bond_color_change (int iw, int i, bool invisible)
{
    int items;
    double c[4];
    char question[MAX_FILENAME_SIZE],danswer[MAX_FILENAME_SIZE],*answer;
    if (invisible)
    {
        C->CYLINDER[i].r = -1;
        C->CYLINDER[i].radius = -1;
    }
    else
    {
        xterm_get_focus(iw); clear_stdin_buffer();
        if (!rgb_url_printed)
        {
            printf("\nColor choices at:\n%s\n", RGB_URL);
            rgb_url_printed = 1;
        }
        sprintf (question, "\nChange color [radius] of bond-%d", i);
        sprintf (danswer, "%.3f %.3f %.3f [%.3f]",
                 C->CYLINDER[i].r, C->CYLINDER[i].g, C->CYLINDER[i].b,
                 C->CYLINDER[i].radius);
        answer = readline_gets(question,danswer);
        items = sscanf (answer, "%lf %lf %lf %lf", c, c+1, c+2, c+3);
        if (items == 1)
        {
            c[3] = c[0];
            c[0] = C->CYLINDER[i].r;
            c[1] = C->CYLINDER[i].g;
            c[2] = C->CYLINDER[i].b;
        }
        else if (items == 3) c[3] = C->CYLINDER[i].radius;
        if (c[0]>1) c[0]/=255;
        if (c[1]>1) c[1]/=255;
        if (c[2]>1) c[2]/=255;
        xterm_release_focus(iw);
    }
    C->CYLINDER[i].r = c[0];
    C->CYLINDER[i].g = c[1];
    C->CYLINDER[i].b = c[2];
    C->CYLINDER[i].radius = c[3];
    return (TRUE);
} /* end bond_color_change() */


/* direct change of atom color and/or radius */
bool atom_color_change (int iw, int i, bool invisible)
{
    int items;
    double c[4];
    char question[MAX_FILENAME_SIZE],danswer[MAX_FILENAME_SIZE],*answer;
    if (invisible)
    {
        c[0] = c[1] = c[2] = -1;
        c[3] = B->BALL[i].radius / n[iw].atom_r_ratio;
    }
    else
    {
        xterm_get_focus(iw); clear_stdin_buffer();
        if (!rgb_url_printed)
        {
            printf("\nColor choices at:\n%s\n", RGB_URL);
            rgb_url_printed = 1;
        }
        sprintf (question, "\nChange color [radius] of atom-%d (\"%s\")",
                 i, COMPACT_SYMBOL(Symbol(symbol,i)));
        sprintf (danswer, "%.3f %.3f %.3f [%.3f]",
                 B->BALL[i].r, B->BALL[i].g, B->BALL[i].b,
                 B->BALL[i].radius / n[iw].atom_r_ratio);
        answer = readline_gets(question,danswer);
        items = sscanf (answer, "%lf %lf %lf %lf", c, c+1, c+2, c+3);
        if (items == 1)
        {
            c[3] = c[0];
            c[0] = B->BALL[i].r;
            c[1] = B->BALL[i].g;
            c[2] = B->BALL[i].b;
        }
        else if (items == 3) c[3] = B->BALL[i].radius / n[iw].atom_r_ratio;
        if (c[0]>1) c[0]/=255;
        if (c[1]>1) c[1]/=255;
        if (c[2]>1) c[2]/=255;
        xterm_release_focus(iw);
    }
    AX_3D_AssignRGB(B->BALL[i], c[0],c[1],c[2]);
    B->BALL[i].radius = c[3] * n[iw].atom_r_ratio;
    bond_atom_color_update (iw);
    return (TRUE);
} /* end atom_color_change() */


bool normal_color_change (int iw, int t, bool invisible)
{
    register int i;
    int items;
    double c[4];
    char question[MAX_FILENAME_SIZE],danswer[MAX_FILENAME_SIZE],*answer;
    /* selection means there must be something there already */
    if (n[iw].color_mode != COLOR_MODE_NORMAL)
    {
        printf("You need to be in NORMAL color mode first.\n");
        return (FALSE);
    }
    if (invisible)
    {
        c[0] = c[1] = c[2] = -1;
        c[3] = ATOM_Radius(ct->Z[t]);
    }
    else
    {
        xterm_get_focus(iw); clear_stdin_buffer();
        if (!rgb_url_printed)
        {
            printf("\nColor choices at:\n%s\n", RGB_URL);
            rgb_url_printed = 1;
        }
        sprintf (question, "\nChange color [radius] of type-%d (\"%s\") atoms",
                 t, COMPACT_SYMBOL(Symbol(symbol,ct->first[t])));
        sprintf (danswer, "%.3f %.3f %.3f [%.3f]",
                 ATOM_Color_R(ct->Z[t]), ATOM_Color_G(ct->Z[t]),
                 ATOM_Color_B(ct->Z[t]), ATOM_Radius(ct->Z[t]));
        answer = readline_gets(question,danswer);
        items = sscanf (answer, "%lf %lf %lf %lf", c, c+1, c+2, c+3);
        if (items == 1)
        {
            c[3] = c[0];
            c[0] = ATOM_Color_R(ct->Z[t]);
            c[1] = ATOM_Color_G(ct->Z[t]);
            c[2] = ATOM_Color_B(ct->Z[t]);
        }
        else if (items == 3) c[3] = ATOM_Radius(ct->Z[t]);
        if (c[0]>1) c[0]/=255;
        if (c[1]>1) c[1]/=255;
        if (c[2]>1) c[2]/=255;
        ATOM_Color_R(ct->Z[t]) = c[0];
        ATOM_Color_G(ct->Z[t]) = c[1];
        ATOM_Color_B(ct->Z[t]) = c[2];
        ATOM_Radius(ct->Z[t])  = c[3];
        xterm_release_focus(iw);
    }
    return (assign_normal_color(iw));
} /* end normal_color_change() */


bool assign_original_normal_color (int iw)
{
    register int i;
    memcpy(Dmitri, MENDELEYEV, (MENDELEYEV_MAX+1)*sizeof(struct Mendeleyev));
    if ((n[iw].xtal_mode) && (n[iw].color_mode == COLOR_MODE_COORD))
    {
        memcpy(CoordColor, ATOM_COORDINATION_COLOR,
               (ATOM_COORDINATION_MAX+1)*sizeof(Atom_coordination_color));
        return(assign_coordination_color(iw));
    }
    return (assign_normal_color(iw));
} /* end assign_original_normal_color() */


bool change_coordination_color (int iw, int coord, bool invisible)
{
    char question[MAX_FILENAME_SIZE],danswer[MAX_FILENAME_SIZE],*answer;
    double c[3];
    if (!n[iw].xtal_mode)
    {
        printf ("Coordination color change works only in Xtal mode.\n");
        return (FALSE);
    }
    /* selection means there must be something there already */
    if (n[iw].color_mode != COLOR_MODE_COORD)
    {
        printf ("You need to be in COORDINATION color mode first.\n");
        return (FALSE);
    }
    if (invisible) c[0] = c[1] = c[2] = -1;
    else
    {
        xterm_get_focus(iw); clear_stdin_buffer();
        if (!rgb_url_printed)
        {
            printf("\nColor choices at:\n%s\n", RGB_URL);
            rgb_url_printed = 1;
        }
        sprintf (question, "\nChange color of %d-coordinated atoms",
                 coord);
        sprintf (danswer, "%.3f %.3f %.3f", CoordColor[coord].r,
                 CoordColor[coord].g, CoordColor[coord].b);
        answer = readline_gets(question,danswer);
        sscanf (answer, "%lf %lf %lf", c, c+1, c+2);
        if (c[0]>1) c[0]/=255;
        if (c[1]>1) c[1]/=255;
        if (c[2]>1) c[2]/=255;
        xterm_release_focus(iw);
    }
    CoordColor[coord].r = c[0];
    CoordColor[coord].g = c[1];
    CoordColor[coord].b = c[2];
    return(assign_coordination_color(iw));
} /* end change_coordination_color() */


/* Change projection method from perspective to parallel */
bool perspective_to_parallel (int iw)
{
    V3 dx,v;
    if (n[iw].parallel_projection)
    {
        printf ("This window is already in parallel projection.\n");
        return (FALSE);
    }
    if (n[iw].anchor >= 0) V3EQV (B->BALL[n[iw].anchor].x, n[iw].hook);
    V3SUB (n[iw].hook, AX_3D[iw].x, dx);
    v[2] = V3DOT(AX_3D[iw].V[2], dx);
    if (v[2] < 0)
    {
        look_at_the_anchor (iw);
        v[2] = -v[2];
    }
    V3ADDmuL ( (1-PARALLEL_AMP)*v[2], AX_3D[iw].V[2], AX_3D[iw].x);
    AX_3D[iw].k *= PARALLEL_AMP;
    n[iw].parallel_projection = 1;
    return (TRUE);
} /* end perspective_to_parallel() */


/* Change projection method from parallel to perspective */
bool parallel_to_perspective (int iw)
{
    V3 dx,v;
    if (!n[iw].parallel_projection)
    {
        printf ("This window is already in perspective projection.\n");
        return (FALSE);
    }
    if (n[iw].anchor >= 0) V3EQV (B->BALL[n[iw].anchor].x, n[iw].hook);
    V3SUB (n[iw].hook, AX_3D[iw].x, dx);
    v[2] = V3DOT(AX_3D[iw].V[2], dx);
    if (v[2] < 0)
    {
        look_at_the_anchor (iw);
        v[2] = -v[2];
    }
    V3ADDmuL ( (1-1./PARALLEL_AMP)*v[2], AX_3D[iw].V[2], AX_3D[iw].x);
    AX_3D[iw].k /= PARALLEL_AMP;
    n[iw].parallel_projection = 0;
    return (TRUE);
} /* end parallel_to_perspective() */


void xterm_get_focus (int iw)
{
    XWindowAttributes WindowAttributes;
    XMapRaised (AX_display[iw], xterm_win);
    AXSYNC(iw);
    do XGetWindowAttributes(AX_display[iw],xterm_win,&WindowAttributes);
    while (WindowAttributes.map_state != IsViewable);
    XSetInputFocus(AX_display[iw], xterm_win,
                   RevertToPointerRoot, CurrentTime);
    AXSYNC(iw);
    return;
} /* end xterm_get_focus() */


void xterm_release_focus (int iw)
{
    AXRaiseWindow(iw);
    XSetInputFocus(AX_display[iw], AX_win[iw], RevertToNone, CurrentTime);
    AXSYNC(iw);
    return;
} /* end xterm_release_focus() */


void reset_auxiliary_threshold (int iw, int i)
{
    SimpleStatistics ss;
    CalculateSimpleStatistics
        (np, CHARP(INW(i,CONFIG_num_auxiliary) ?
                   CONFIG_auxiliary[i] : geo[i-CONFIG_MAX_AUXILIARY]),
         sizeof(double), IOVAL_DOUBLE, &ss);
    n[iw].auxiliary_threshold[i][0] = ss.min;
    n[iw].auxiliary_threshold[i][1] = ss.max;
    return;
} /* end reset_auxiliary_threshold() */


void atom_stack_insert (int iw, int atom)
{
    int i;
    for (i=ATOM_STACK_SIZE-1; i>0; i--)
        n[iw].atom_stack[i] = n[iw].atom_stack[i-1];
    n[iw].atom_stack[0] = atom;
    return;
} /* end atom_stack_insert() */
