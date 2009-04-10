/********************************************/
/* libAX: -lX11 -lXext -lpng -lz -ljpeg -lm */
/*        -lScalar -lIO -lTimer             */
/*                                          */
/* Accelerated low-level graphics library   */
/* with X Window Shared Memory Extension.   */
/*                                          */
/* Jan.13, 2000 Ju Li <liju99@mit.edu>      */
/********************************************/

#include "AX.h"

/*******************************************/
/* Transferable Color and Pixel operations */
/*******************************************/


/** AX_Pixel (unsigned long): interface to XSetForeground() etc. **/

/* color name must exist in "/usr/X11/lib/X11/rgb.txt" */
AX_Pixel AX_color_pixel_via_name (int iw, char *colorname)
{
    XColor screen_def_return, exact_def_return;
    if ( XAllocNamedColor ( AX_display [iw], AX_colormap [iw],
                            colorname, &screen_def_return,
                            &exact_def_return ) == 0 )
        pe ("AX_color_pixel_via_name: color \"%s\" unavailable!\n",
            colorname);
/*     pr ("%s %d %d %d 0x%lx %d %d %d 0x%lx\n", colorname,  */
/* 	screen_def_return.red,  screen_def_return.green, */
/* 	screen_def_return.blue, screen_def_return.pixel, */
/* 	exact_def_return.red,   exact_def_return.green,  */
/* 	exact_def_return.blue,  exact_def_return.pixel); */
    /* most significant byte seems to have value 0xff on Linux x86_64,
       rather than expected value of 0x00, so we & with 0xffffff to clear it. */
    return( (AX_Pixel) screen_def_return.pixel & 0xffffff);
} /* end AX_color_pixel_via_name() */


/************************************************************/
/* MATLAB colors... any self-respecting X server, even the  */
/* 8-bit ones, must have them in "/usr/X11/lib/X11/rgb.txt" */
/* They are considered the birth rights of each window.     */
/************************************************************/
static char *AX_namedpixelname [AX_NAMEDCOLORS] =
{"black", "red", "green", "blue", "cyan", "magenta", "yellow", "white"};

/**************************************************************/
/* Pixels provided by AX gratis. Each call to AX_openwindow() */
/* updates them, in case some applications "steal our colors" */
/**************************************************************/
AX_Pixel AX_namedpixel [AX_NAMEDCOLORS];

/** AX_Carrier: interface to Direct Pixmap Access like AX_set() **/

AX_Carrier AX_namedcarrier [AX_NAMEDCOLORS];

/* caretaker iw updates community named color pixels and carriers */
void AX_update_namedcolors (int iw)
{
    register int i;
    for (i=0; i<AX_NAMEDCOLORS; i++)
    {
        AX_namedpixel[i] = AX_color_pixel_via_name(iw, AX_namedpixelname[i]);
        AX_namedcarrier[i] = AX_PIXEL2CARRIER(AX_namedpixel[i]);
    }
    return;
} /* end AX_update_namedcolors() */


/* Identical to AX_color_pixel_via_name() except we may save some memory */
AX_Carrier AX_color_carrier_via_name (int iw, char *colorname)
{
    XColor screen_def_return, exact_def_return;
    if ( XAllocNamedColor ( AX_display [iw], AX_colormap [iw],
                            colorname, &screen_def_return,
                            &exact_def_return ) == 0 )
        pe ("AX_color_carrier_via_name: color \"%s\" unavailable!\n",
            colorname);
    return( (AX_Carrier) screen_def_return.pixel );
} /* end AX_color_carrier_via_name() */


/* Input: h, s, v in range [0..1]; Outputs: r, g, b in range [0..1] */
static void hsv2rgb
(double h, double s, double v, double *r, double *g, double *b)
{
    int i;
    double aa, bb, cc, f;
    if (s == 0.)
        *r = *g = *b = v;
    else
    {
        if (h == 1.0) h = 0;
        h *= 6.0;
        i = floor (h);
        f = h - i;
        aa = v * (1 - s);
        bb = v * (1 - (s * f));
        cc = v * (1 - (s * (1 - f)));
        switch (i)
        {
            case 0: *r = v;  *g = cc; *b = aa; break;
            case 1: *r = bb; *g = v;  *b = aa; break;
            case 2: *r = aa; *g = v;  *b = cc; break;
            case 3: *r = aa; *g = bb; *b = v;  break;
            case 4: *r = cc; *g = aa; *b = v;  break;
            case 5: *r = v;  *g = aa; *b = bb; break;
        }
    }
    return;
} /* end hsv2rgb() */


/* Input:   r, g, b in range [0..1]; Outputs: h, s, v in range [0..1] */
static void rgb2hsv
(double r, double g, double b, double * h, double * s, double * v)
{
    double max = MAX (r, MAX (g, b)), min = MIN (r, MIN (g, b));
    double delta = max - min;
    *v = max;
    if (max != 0.0)
        *s = delta / max;
    else
        *s = 0.0;
    if (*s == 0.0) *h = -1;
    else {
        if (r == max)
            *h = (g - b) / delta;
        else if (g == max)
            *h = 2 + (b - r) / delta;
        else if (b == max)
            *h = 4 + (r - g) / delta;
        *h *= 60.0;
        if (*h < 0) *h += 360.0;
        *h /= 360.0;
    }
    return;
} /* end rgb2hsv() */


/* NCSA fluid jet image. blue->cyan->yellow->orange->red. */
static void jet (double x, double *r, double *g, double *b)
{
    if (x < 0.125)
    { /* blue */
        *r = 0;
        *g = 0;
        *b = 0.5 + 0.5 * x / 0.125;
    }
    else if (x < 0.125 + 0.25)
    { /* cyan */
        *r = 0;
        *g = (x - 0.125) / 0.25;
        *b = 1;
    }
    else if (x < 0.125 + 0.25 + 0.25)
    { /* yellow */
        *r = (x - 0.375) / 0.25;
        *g = 1;
        *b = 1 - (x - 0.375) / 0.25;
    }
    else if (x < 0.125 + 0.25 + 0.25 + 0.25)
    { /* orange */
        *r = 1;
        *g = 1 - (x - 0.625) / 0.25;
        *b = 0;
    }
    else
    { /* red */
        *r = 1 - 0.5 * (x - 0.875) / 0.125;
        *g = 0;
        *b = 0;
    }
    return;
} /* end jet() */


/* black -> red -> yellow -> white */
static void hot (double x, double *r, double *g, double *b)
{
    if (x < 0.375)
    { /* black -> red */
        *r = x / 0.375;
        *g = 0;
        *b = 0;
    }
    else if (x < 0.75)
    { /* red -> yellow */
        *r = 1;
        *g = (x - 0.375) / 0.375;
        *b = 0;
    }
    else
    { /* yellow -> white */
        *r = 1;
        *g = 1;
        *b = (x - 0.75) / 0.25;
    }
    return;
} /* end hot() */


/* shades of cyan -> magenta */
static void cool (double x, double *r, double *g, double *b)
{
    *r = x;
    *g = 1-x;
    *b = 1;
    return;
} /* end cool() */


/* Linear gray-scale color map. */
static void gray (double x, double *r, double *g, double *b)
{
    *r = x;
    *g = x;
    *b = x;
    return;
} /* end gray() */


/* pastel shades of pink color map */
static void pink (double x, double *r, double *g, double *b)
{
    hot(x, r, g, b);
    *r = sqrt((2*x+*r)/3);
    *g = sqrt((2*x+*g)/3);
    *b = sqrt((2*x+*b)/3);
    return;
} /* end pink() */


/* gray-scale with a tinge of blue color map */
static void bone (double x, double *r, double *g, double *b)
{
    double rh,gh,bh;
    hot(x, &rh, &gh, &bh);
    *r = (7*x + bh)/8;
    *g = (7*x + gh)/8;
    *b = (7*x + rh)/8;
    return;
} /* end bone() */


/* linear copper-tone color map */
static void copper (double x, double *r, double *g, double *b)
{
    *r = MIN(1.25*x,1);
    *g = MIN(0.7812*x,1);
    *b = MIN(0.4975*x,1);
    return;
} /* end copper() */


/* shades of red and yellow color map */
static void autumn (double x, double *r, double *g, double *b)
{
    *r = 1;
    *g = x;
    *b = 0;
    return;
} /* end autumn() */


/* shades of magenta and yellow color map */
static void spring (double x, double *r, double *g, double *b)
{
    *r = 1;
    *g = x;
    *b = 1-x;
    return;
} /* end spring() */


/* shades of blue and green color map */
static void winter (double x, double *r, double *g, double *b)
{
    *r = 0;
    *g = x;
    *b = .5 + .5 * (1-x);
    return;
} /* end winter() */


/* shades of green and yellow color map */
static void summer (double x, double *r, double *g, double *b)
{
    *r = x;
    *g = .5 + .5 * x;
    *b = .4;
    return;
} /* end summer() */


/* hue-saturation-value color map */
static void hsv (double x, double *r, double *g, double *b)
{
    double s=1., v=1.;
    hsv2rgb (x, s, v, r, g, b);
    return;
} /* end hsv() */


const AX_Cmap AX_cmap_funs[AX_MAX_CMAP] =
{
    {"jet", "blue->cyan->yellow->orange->red", &jet},
    {"hot", "black -> red -> yellow -> white", &hot},
    {"cool", "shades of cyan -> magenta", &cool},
    {"gray", "linear gray-scale", &gray},
    {"pink", "pastel shades of pink", &pink},
    {"bone", "gray-scale with a tinge of blue", &bone},
    {"copper", "linear copper-tone", &copper},
    {"autumn", "shades of red -> yellow", &autumn},
    {"spring", "shades of magenta -> yellow", &spring},
    {"winter", "shades of blue -> green", &winter},
    {"summer", "shades of green -> yellow", &summer},
    {"hsv", "linear hue change", &hsv},
};


/* Calculate (r,g,b) for x in [0,1] according to scheme "cmap_idx". */
void AX_cmap (int cmap_idx, double x, double *r, double *g, double *b)
{
    if (OUW(cmap_idx,AX_MAX_CMAP))
    {
        pr ("AX_cmap: cmap_idx=%d >= AX_MAX_CMAP=%d\n",
            cmap_idx, AX_MAX_CMAP);
        cmap_idx = positive_remainder (cmap_idx, AX_MAX_CMAP);
    }
    if (XOU(x,0,1))
    { /* x must be in [0,1] */
        *r = *g = *b = 0;
        return;
    }
    AX_cmap_funs[cmap_idx].fun (x, r, g, b);
    return;
} /* end AX_cmap() */


#ifdef _AX_cmap_TEST
#define AX_CMAP_IDX AX_CMAP_HSV
#define M 256
int main (int argc, char *argv[])
{
    int i;
    double r,g,b;
    printf ("The colormap %d is \"%s\".\n", AX_CMAP_IDX,
            AX_cmap_funs[AX_CMAP_IDX].name);
    printf ("Compare the following printout with Matlab function\n"
            "fprintf(1,'%%5d %%9.4f %%9.4f %%9.4f\\n',"
            "[(1:%d)' %s(%d)]')\n\n",
            M, AX_cmap_funs[AX_CMAP_IDX].name, M);
    for (i=1; i<=M; i++)
    {
        AX_cmap (AX_CMAP_IDX, DOUBLE(i)/M, &r, &g, &b);
        printf ("%5d %9.4f %9.4f %9.4f\n", i, r, g, b);
    }
    return (0);
}
#undef AX_CMAP_IDX
#undef M
#endif /* _AX_cmap_TEST */


#ifdef _cmap
#define CMAP_DEF_WIDTH  24
#define CMAP_DEF_HEIGHT 640
int main (int argc, char *argv[])
{
    int cmap_idx=0, width, height, i, j, iw;
    double r, g, b;
    AX_Carrier carrier;
    char fname[TERMSIZE];
    if (argc <= 1)
    {
        printf ("Purpose: generate jpg colormap for labeling purposes.\n");
        printf ("Usage: %s cmap_idx [height] [width] [jpg_filename]\n",
                argv[0]);
        printf ("cmap_idx:\n");
        for (i=0; i<AX_MAX_CMAP; i++)
            printf ("%2d: %s: %s\n",
                    i, AX_cmap_funs[i].name,
                    AX_cmap_funs[i].description);
        return (1);
    }
    cmap_idx = atoi(argv[1]);
    if (OUW(cmap_idx,AX_MAX_CMAP))
    {
        pr ("cmap_idx=%d >= AX_MAX_CMAP=%d\n",
            cmap_idx, AX_MAX_CMAP);
        cmap_idx = positive_remainder (cmap_idx, AX_MAX_CMAP);
    }
    if (argc == 2)
    {
        height = CMAP_DEF_HEIGHT;
        width = CMAP_DEF_WIDTH;
        sprintf (fname, "%s.jpg", AX_cmap_funs[cmap_idx].name);
    }
    else if (argc == 3)
    {
        height = atoi(argv[2]);
        width = height * DOUBLE(CMAP_DEF_WIDTH) / CMAP_DEF_HEIGHT;
        sprintf (fname, "%s.jpg", AX_cmap_funs[cmap_idx].name);
    }
    else if (argc == 4)
    {
        height = atoi(argv[2]);
        width = atoi(argv[3]);
        sprintf (fname, "%s.jpg", AX_cmap_funs[cmap_idx].name);
    }
    else
    {
        height = atoi(argv[2]);
        width = atoi(argv[3]);
        strcpy (fname, argv[4]);
    }
    AX_WINDOWS_OVERRIDE_REDIRECT = True;
    iw = AX_openwindow (getpid(),NULL,width,height);
    for (j=0; j<height; j++)
    {
        AX_cmap (cmap_idx, 1-DOUBLE(j)/height, &r, &g, &b);
        carrier = AX_Colorcarrier(r,g,b);
        AX_bset(iw, j*width, carrier, width, i);
    }
    printf ("%d x %d %s colormap -> \"%s\" (%d bytes).\n",
            width, height, AX_cmap_funs[cmap_idx].name, fname,
            AX_save_pixmap_as_JPG(0,fname));
    printf ("Compare with Matlab function\n"
            "colormap(%s(64)); surf((1:64)'*ones(1,2));\n",
            AX_cmap_funs[cmap_idx].name);
    AX_dump(iw);
    AX_show(iw);
    Press_return();
    /* try_to_runbg (NUMBER_RASTER_VIEWERS, raster_viewers, fname); */
    AX_closewindow (iw);
    return (iw);
}
#endif /* _cmap */
