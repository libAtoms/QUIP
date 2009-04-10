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

/************************************************************/
/* Scan-convert 2D primitives with floating-point precision */
/************************************************************/

/***************************/
/* Scan Conversion Manager */
/***************************/

/* pixel-based "alias" operation stack */
AX_AP *AX_aop [AX_MAXWIN];
/* pixel-based "block" operation stack */
AX_BP *AX_bop [AX_MAXWIN];
/* "alias" stack memory limit */
AX_Offset AX_maxaop [AX_MAXWIN];
/* "block" stack space limit */
AX_Offset AX_maxbop [AX_MAXWIN];

/* Allocate scan conversion stacks AX_aop[iw] and AX_bop[iw], */
/* sizes of which are based on an estimate using AX_size[iw]. */
int AX_plugin_Scan_module (int iw)
{
    char *p;
    int pixels, aop_bytes, bop_bytes, bytes;
    AX_AssertWinID ("AX_plugin_Scan_module", iw);
    AX_Check_Module_Dependence ("AX_plugin_Scan_module", iw, MODULE_Scan);
    if ( AX_Module_Is_PluggedIn [iw] [MODULE_Scan] )
        pe ("AX_plugin_Scan_module:\n"
            "you are crazy, this module is already plugged in\n");
    if (AX_size[iw].width > AX_MAXWIDTH)
        pe ("AX_plugin_Scan_module: width = %d exceeds AX_MAXWIDTH = %d\n",
            AX_size[iw].width, AX_MAXWIDTH);
    if (AX_size[iw].height > AX_MAXHEIGHT)
        pe ("AX_plugin_Scan_module: height = %d exceeds AX_MAXHEIGHT = %d\n",
            AX_size[iw].height, AX_MAXHEIGHT);
    if (! AX_noneedcolormap)
        pe ("AX_plugin_Scan_module: Direct/TrueColor visuals\n"
            "unsupported by this display. You cannot use any of the\n"
            "anti-aliasing features of the Scan_module.\n");
    if (AX_Colorpixel(1,0,0) != AX_namedpixel[AX_RED])
    {
        pr ("colormap id = 0x%x\n", AX_colormap[iw]);
        pr ("%x %x 0x%lx\n", AX_Colorpixel(1,0,0), AX_namedpixel[AX_RED], AX_rmask);
        pr ("AX_plugin_Scan_module: red mask not working\n");
    }
    if (AX_Colorpixel(0,1,0) != AX_namedpixel[AX_GREEN])
    {
        pr ("colormap id = 0x%x\n", AX_colormap[iw]);
        pr ("%x %x 0x%lx\n", AX_Colorpixel(0,1,0), AX_namedpixel[AX_GREEN], AX_gmask);
        pr ("AX_plugin_Scan_module: green mask not working\n");
    }
    if (AX_Colorpixel(0,0,1) != AX_namedpixel[AX_BLUE])
    {
        pr ("colormap id = 0x%x\n", AX_colormap[iw]);
        pr ("%x %x 0x%lx\n", AX_Colorpixel(0,0,1), AX_namedpixel[AX_BLUE], AX_bmask);
        pr ("AX_plugin_Scan_module: blue mask not working\n");
    }
    pixels = AX_size[iw].width * AX_size[iw].height;
    if (pixels < AX_FREEZE_WIDTH*AX_FREEZE_WIDTH)
        pixels = AX_FREEZE_WIDTH*AX_FREEZE_WIDTH;
    AX_maxaop[iw] = AX_AOP_MAXOBJ + AX_AOP_MAXRATIO * pixels;
    aop_bytes = sizeof(AX_AP) * (AX_maxaop[iw] + AX_AOP_WARNINGZONE);
    AX_maxbop[iw] = AX_BOP_MAXOBJ + AX_BOP_MAXRATIO * pixels;
    bop_bytes = sizeof(AX_BP) * (AX_maxbop[iw] + AX_BOP_WARNINGZONE);
    bytes = aop_bytes + bop_bytes;
    if ( (p=(char *)malloc(bytes)) == NULL )
        pe ("AX_plugin_Scan_module: win %d unable to malloc %d\n"
            "bytes of scan conversion stack memory\n", iw, bytes);
    AX_aop[iw] = (AX_AP *) p;  /* alias stack */
    AX_bop[iw] = (AX_BP *) (p + aop_bytes);  /* block stack */
    AX_Module_Is_PluggedIn [iw] [MODULE_Scan] = 1;
    return (1);
} /* end AX_plugin_Scan_module() */


/* resize scan conversion stacks AX_aop[iw] and AX_bop[iw] */
/* according to current AX_size[iw] using realloc().       */
int AX_resize_Scan_module (int iw)
{
    char *p;
    int pixels, aop_bytes, bop_bytes, bytes;
    AX_AssertWinID ("AX_resize_Scan_module", iw);
    if ( ! AX_Module_Is_PluggedIn [iw] [MODULE_Scan] )
        pe ("AX_resize_Scan_module:\n"
            "you are crazy, this module is not plugged in\n");
    /* if (AX_size[iw].width > AX_MAXWIDTH) */
    /* pe ("AX_resize_Scan_module: width = %d exceeds AX_MAXWIDTH = %d\n", */
    /* AX_size[iw].width, AX_MAXWIDTH); */
    /* if (AX_size[iw].height > AX_MAXHEIGHT) */
    /* pe ("AX_resize_Scan_module: height = %d exceeds AX_MAXHEIGHT = %d\n", */
    /* AX_size[iw].height, AX_MAXHEIGHT); */
    pixels = AX_size[iw].width * AX_size[iw].height;
    if (pixels < AX_FREEZE_WIDTH*AX_FREEZE_WIDTH)
        pixels = AX_FREEZE_WIDTH*AX_FREEZE_WIDTH;
    AX_maxaop[iw] = AX_AOP_MAXOBJ + AX_AOP_MAXRATIO * pixels;
    aop_bytes = sizeof(AX_AP) * (AX_maxaop[iw] + AX_AOP_WARNINGZONE);
    AX_maxbop[iw] = AX_BOP_MAXOBJ + AX_BOP_MAXRATIO * pixels;
    bop_bytes = sizeof(AX_BP) * (AX_maxbop[iw] + AX_BOP_WARNINGZONE);
    bytes = aop_bytes + bop_bytes;
    if ( (p=(char *)realloc(AX_aop[iw],bytes)) == NULL )
        pe ("AX_resize_Scan_module: win %d unable to realloc %d\n"
            "bytes of scan conversion stack memory\n", iw, bytes);
    AX_aop[iw] = (AX_AP *) p;  /* alias stack */
    AX_bop[iw] = (AX_BP *) (p + aop_bytes);  /* block stack */
    return (1);
} /* end AX_resize_Scan_module() */


/* Free scan conversion scratch stack AX_aop[iw] and AX_bop[iw] */
int AX_plugout_Scan_module (int iw)
{
    AX_AssertWinID ("AX_plugout_Scan_module", iw);
    if ( !AX_Module_Is_PluggedIn [iw] [MODULE_Scan] )
        pe ("AX_plugout_Scan_module:\n"
            "you are crazy, this module is not plugged in\n");
    Free(AX_aop[iw]);
    AX_Module_Is_PluggedIn [iw] [MODULE_Scan] = 0;
    return(1);
} /* end AX_plugout_Scan_module() */


/**********************/
/* Clipping Functions */
/**********************/

/* Test if line segment (x1,y1)-(x2,y2) intersects with */
/* (x3,y3)-(x4,y4); if it does, save intersection x, y. */
int AX_SegmentsIntersect ( AX_Float x1, AX_Float y1,
                           AX_Float x2, AX_Float y2,
                           AX_Float x3, AX_Float y3,
                           AX_Float x4, AX_Float y4,
                           AX_Float *x, AX_Float *y )
{
    AX_Float dx1, dy1, dx3, dy3, determinant, dx, dy, a1, a2;
    dx1 = x2 - x1;
    dy1 = y2 - y1;
    dx3 = x4 - x3;
    dy3 = y4 - y3;
    determinant = dx3 * dy1 - dx1 * dy3;
    if (determinant==0) return (AX_LINES_PARALLEL);
    dx = x3 - x1;
    dy = y3 - y1;
    a1 = (dx3 * dy - dx * dy3) / determinant;
    a2 = (dx1 * dy - dx * dy1) / determinant;
    *x = x1 + a1 * dx1;
    *y = y1 + a1 * dy1;
    if ( (a1 < 0) || (a1 > 1) || (a2 < 0) || (a2 > 1) )
        return (AX_LINES_INTERSECT_BUT_NOT_SEGMENTS);
    return (AX_SEGMENTS_INTERSECT);
} /* end AX_SegmentsIntersect() */

#ifdef _SegmentsIntersect_TEST
#define LINE1_x1  0.
#define LINE1_y1  0.
#define LINE1_x2  1.
#define LINE1_y2  1.
#define LINE2_x1  (1./2)
#define LINE2_y1  (1./3)
#define LINE2_x2  1.
#define LINE2_y2  0.
int main (int argc, char *argv[])
{
    AX_Float x, y;
    printf("%d\n", AX_SegmentsIntersect
           (LINE1_x1,LINE1_y1,LINE1_x2,LINE1_y2,
            LINE2_x1,LINE2_y1,LINE2_x2,LINE2_y2, &x, &y));
    printf ("%f %f\n", x, y);
    return (0);
}
#endif /* _SegmentsIntersect_TEST */


/***********************************************************************/
/* Line segment clipped by a rectangular window [0,height] x [0,width] */
/***********************************************************************/

/* The Cohen-Sutherland algorithm */
#define CS_LEFT    1
#define CS_RIGHT   2
#define CS_TOP     4
#define CS_BOTTOM  8
#define CS_XMIN    0
#define CS_YMIN    0
#define CS_XMAX    width
#define CS_YMAX    height
/*
       |      | 
  0101 | 0100 | 0110
 ---------------------
  0001 | 0000 | 0010
 ---------------------
  1001 | 1000 | 1010
       |      |
*/
#define CS_calculate_outcode(x,y,code) code = (x < CS_XMIN)? CS_LEFT : 0; \
  if (x > CS_XMAX) code |= CS_RIGHT; \
  if (y < CS_YMIN) code |= CS_BOTTOM; \
  if (y > CS_YMAX) code |= CS_TOP;

bool AX_CS_LineClipWindow ( AX_Float *x0, AX_Float *y0,
                            AX_Float *x1, AX_Float *y1,
                            AX_Float width, AX_Float height )
{
    int outcode0, outcode1, outcodeOut;
    AX_Float x, y;
    CS_calculate_outcode (*x0, *y0, outcode0);
    CS_calculate_outcode (*x1, *y1, outcode1);
    while (1)
    {
        /* trivially accepted */
        if ( (outcode0 | outcode1) == 0 ) return(TRUE);
        /* trivially rejected */
        if ( outcode0 & outcode1 ) return(FALSE);
        outcodeOut = outcode0 ? outcode0 : outcode1;
        if ( outcodeOut & CS_TOP )
        {
            x = *x0 + ( *x1 - *x0 ) * (CS_YMAX - *y0 ) / ( *y1 - *y0 );
            y = CS_YMAX;
        }
        else if ( outcodeOut & CS_BOTTOM )
        {
            x = *x0 + ( *x1 - *x0 ) * (CS_YMIN - *y0 ) / ( *y1 - *y0 );
            y = CS_YMIN;
        }
        else if ( outcodeOut & CS_RIGHT )
        {
            y = *y0 + ( *y1 - *y0 ) * (CS_XMAX - *x0 ) / ( *x1 - *x0 );
            x = CS_XMAX;
        }
        else
        {
            y = *y0 + ( *y1 - *y0 ) * (CS_XMIN - *x0 ) / ( *x1 - *x0 );
            x = CS_XMIN;
        }
        if (  outcodeOut == outcode0 )
        {
            *x0 = x; *y0 = y;
            CS_calculate_outcode (x, y, outcode0);
        }
        else
        {
            *x1 = x; *y1 = y;
            CS_calculate_outcode (x, y, outcode1);
        }
    }
    return (TRUE);
} /* end AX_CS_LineClipWindow() */
#undef CS_XMIN
#undef CS_YMIN
#undef CS_XMAX
#undef CS_YMAX

#define CS_XMIN   xmin
#define CS_YMIN   ymin
#define CS_XMAX   xmax
#define CS_YMAX   ymax
bool AX_CS_LineClipRectangle ( AX_Float *x0, AX_Float *y0,
                               AX_Float *x1, AX_Float *y1,
                               AX_Float xmin, AX_Float ymin,
                               AX_Float xmax, AX_Float ymax)
{
    int outcode0, outcode1, outcodeOut;
    AX_Float x, y;
    CS_calculate_outcode (*x0, *y0, outcode0);
    CS_calculate_outcode (*x1, *y1, outcode1);
    while (1)
    {
        /* trivially accepted */
        if ( (outcode0 | outcode1) == 0 ) return(TRUE);
        /* trivially rejected */
        if ( outcode0 & outcode1 ) return(FALSE);
        outcodeOut = outcode0 ? outcode0 : outcode1;
        if ( outcodeOut & CS_TOP )
        {
            x = *x0 + ( *x1 - *x0 ) * (CS_YMAX - *y0 ) / ( *y1 - *y0 );
            y = CS_YMAX;
        }
        else if ( outcodeOut & CS_BOTTOM )
        {
            x = *x0 + ( *x1 - *x0 ) * (CS_YMIN - *y0 ) / ( *y1 - *y0 );
            y = CS_YMIN;
        }
        else if ( outcodeOut & CS_RIGHT )
        {
            y = *y0 + ( *y1 - *y0 ) * (CS_XMAX - *x0 ) / ( *x1 - *x0 );
            x = CS_XMAX;
        }
        else
        {
            y = *y0 + ( *y1 - *y0 ) * (CS_XMIN - *x0 ) / ( *x1 - *x0 );
            x = CS_XMIN;
        }
        if (  outcodeOut == outcode0 )
        {
            *x0 = x; *y0 = y;
            CS_calculate_outcode (x, y, outcode0);
        }
        else
        {
            *x1 = x; *y1 = y;
            CS_calculate_outcode (x, y, outcode1);
        }
    }
    return (TRUE);
} /* end AX_CS_LineClipRectangle() */

#ifdef _CS_LineClipRectangle_TEST
#define CS_WIDTH  200
#define CS_HEIGHT 150
#define CS_TRIALS 10
int main (int argc, char *argv[])
{
    int i;
    AX_Float x0,y0,x1,y1,wx,wy;
    AX_openwindow (getpid(),NULL,AX_DEFWIDTH,AX_DEFHEIGHT);
    AXSetNamedForeground(0,AX_GREEN);
    wx = (AX_DEFWIDTH - CS_WIDTH) / 2;
    wy = (AX_DEFHEIGHT - CS_HEIGHT) / 2;
    AXDrawRectangle(0,wx,wy,CS_WIDTH,CS_HEIGHT);
    TimeRandomize();
    for (i=0; i<CS_TRIALS; i++)
    {
        x0 = Frandom() * AX_DEFWIDTH;
        y0 = Frandom() * AX_DEFHEIGHT;
        x1 = Frandom() * AX_DEFWIDTH;
        y1 = Frandom() * AX_DEFHEIGHT;
        AXSetNamedForeground(0,AX_WHITE);
        AXDrawLine(0,x0,y0,x1,y1);
        if ( AX_CS_LineClipRectangle(&x0,&y0,&x1,&y1,wx,wy,
                                     wx+CS_WIDTH,wy+CS_HEIGHT) )
        {
            AXSetNamedForeground(0,AX_RED);
            AXDrawLine(0,x0,y0,x1,y1);
        }
    }
    AX_dump(0); AX_show(0);
    Press_return();
    return (0);
}
#endif /* _CS_LineClipRectangle_TEST */


/* Mark S. Sobkow, Paul Pospisil and Yee-Hong Yang, Computers & */
/* Graphics Vol. 11, No. 4, pp. 459-467, 1987, adopted by John  */
/* Schultz in Notes/clip2d.readme, clip2d.mod.                  */

#define ClipPBottom() \
  (*x0 = (*x1 - *x0) * (SPY_YMAX - *y0) / (*y1 - *y0) + *x0, *y0 = SPY_YMAX)
#define ClipPTop() \
  (*x0 = (*x1 - *x0) * (SPY_YMIN - *y0) / (*y1 - *y0) + *x0, *y0 = SPY_YMIN)
#define ClipPRight() \
  (*y0 = (*y1 - *y0) * (SPY_XMAX - *x0) / (*x1 - *x0) + *y0, *x0 = SPY_XMAX)
#define ClipPLeft() \
  (*y0 = (*y1 - *y0) * (SPY_XMIN - *x0) / (*x1 - *x0) + *y0, *x0 = SPY_XMIN)
#define ClipQBottom() \
  (*x1 = (*x0 - *x1) * (SPY_YMAX - *y1) / (*y0 - *y1) + *x1, *y1 = SPY_YMAX)
#define ClipQTop() \
  (*x1 = (*x0 - *x1) * (SPY_YMIN - *y1) / (*y0 - *y1) + *x1, *y1 = SPY_YMIN)
#define ClipQRight() \
  (*y1 = (*y0 - *y1) * (SPY_XMAX - *x1) / (*x0 - *x1) + *y1, *x1 = SPY_XMAX)
#define ClipQLeft() \
  (*y1 = (*y0 - *y1) * (SPY_XMIN - *x1) / (*x0 - *x1) + *y1, *x1 = SPY_XMIN)

#define SPY_XMIN    0
#define SPY_YMIN    0
#define SPY_XMAX    width
#define SPY_YMAX    height
bool AX_SPY_LineClipWindow ( AX_Float *x0, AX_Float *y0,
                             AX_Float *x1, AX_Float *y1,
                             AX_Float width, AX_Float height )
{
    int code = 0;

    if (*y1 > SPY_YMAX)  code |= 8;
    else if (*y1 < SPY_YMIN)  code |= 4;

    if (*x1 > SPY_XMAX)  code |= 2;
    else if (*x1 < SPY_XMIN)  code |= 1;

    if (*y0 > SPY_YMAX)  code |= 128;
    else if (*y0 < SPY_YMIN)  code |= 64;

    if (*x0 > SPY_XMAX)  code |= 32;
    else if (*x0 < SPY_XMIN)  code |= 16;

    switch (code)
    {
        case 0x00: return TRUE;
        case 0x01: ClipQLeft(); return TRUE;
        case 0x02: ClipQRight(); return TRUE;
        case 0x04: ClipQTop(); return TRUE;
        case 0x05: ClipQLeft(); if (*y1 < SPY_YMIN) ClipQTop(); return TRUE;
        case 0x06: ClipQRight(); if (*y1 < SPY_YMIN) ClipQTop(); return TRUE;
        case 0x08: ClipQBottom(); return TRUE;
        case 0x09: ClipQLeft(); if (*y1 > SPY_YMAX) ClipQBottom(); return TRUE;
        case 0x0A: ClipQRight(); if (*y1 > SPY_YMAX) ClipQBottom();
            return TRUE;
        case 0x10: ClipPLeft(); return TRUE;
        case 0x11: return FALSE;
        case 0x12: ClipPLeft(); ClipQRight(); return TRUE;
        case 0x14: ClipPLeft(); if (*y0 < SPY_YMIN) return FALSE; else
        { ClipQTop(); return TRUE; }
        case 0x15: return FALSE;
        case 0x16: ClipPLeft(); if (*y0 < SPY_YMIN) return FALSE; else
        { ClipQTop(); if (*x1 > SPY_XMAX) ClipQRight(); return TRUE; }
        case 0x18: ClipPLeft(); if (*y0 > SPY_YMAX) return FALSE; else
        { ClipQBottom(); return TRUE; }
        case 0x19: return FALSE;
        case 0x1A: ClipPLeft(); if (*y0 > SPY_YMAX) return FALSE; else
        { ClipQBottom(); if (*x1 > SPY_XMAX) ClipQRight(); return TRUE; }
        case 0x20: ClipPRight(); return TRUE;
        case 0x21: ClipPRight(); ClipQLeft(); return TRUE;
        case 0x22: return FALSE;
        case 0x24: ClipPRight(); if (*y0 < SPY_YMIN) return FALSE; else
        { ClipQTop(); return TRUE; }
        case 0x25: ClipPRight(); if (*y0 < SPY_YMIN) return FALSE; else
        { ClipQTop(); if (*x1 < SPY_XMIN) ClipQLeft(); return TRUE; }
        case 0x26: return FALSE;
        case 0x28: ClipPRight(); if (*y0 > SPY_YMAX) return FALSE; else
        { ClipQBottom(); return TRUE; }
        case 0x29: ClipPRight(); if (*y0 > SPY_YMAX) return FALSE; else
        { ClipQBottom(); if (*x1 < SPY_XMIN) ClipQLeft(); return TRUE; }
        case 0x2A: return FALSE;
        case 0x40: ClipPTop(); return TRUE;
        case 0x41: ClipPTop(); if (*x0 < SPY_XMIN) return FALSE; else
        { ClipQLeft(); if (*y1 < SPY_YMIN) ClipQTop(); return TRUE; }
        case 0x42: ClipPTop(); if (*x0 > SPY_XMAX) return FALSE; else
        { ClipQRight(); return TRUE; }
        case 0x44: return FALSE;
        case 0x45: return FALSE;
        case 0x46: return FALSE;
        case 0x48: ClipPTop(); ClipQBottom(); return TRUE;
        case 0x49: ClipPTop(); if (*x0 < SPY_XMIN) return FALSE; else
          { ClipQLeft(); if (*y1 > SPY_YMAX) ClipQBottom(); return TRUE; }
        case 0x4A: ClipPTop(); if (*x0 > SPY_XMAX) return FALSE; else
        { ClipQRight(); if (*y1 > SPY_YMAX) ClipQBottom(); return TRUE; }
        case 0x50: ClipPLeft(); if (*y0 < SPY_YMIN) ClipPTop(); return TRUE;
        case 0x51: return FALSE;
        case 0x52: ClipQRight(); if (*y1 < SPY_YMIN) return FALSE; else
        { ClipPTop(); if (*x0 < SPY_XMIN) ClipPLeft(); return TRUE; }
        case 0x54: return FALSE;
        case 0x55: return FALSE;
        case 0x56: return FALSE;
        case 0x58: ClipQBottom(); if (*x1 < SPY_XMIN) return FALSE; else
        { ClipPTop(); if (*x0 < SPY_XMIN) ClipPLeft(); return TRUE; }
        case 0x59: return FALSE;
        case 0x5A: ClipPLeft(); if (*y0 > SPY_YMAX) return FALSE; else
        { ClipQRight(); if (*y1 < SPY_YMIN) return FALSE; else
        { if (*y0 < SPY_YMIN) ClipPTop(); if (*y1 > SPY_YMAX) ClipQBottom();
        return TRUE; } }
        case 0x60: ClipPRight(); if (*y0 < SPY_YMIN) ClipPTop(); return TRUE;
        case 0x61: ClipQLeft(); if (*y1 < SPY_YMIN) return FALSE; else
        { ClipPTop(); if (*x0 > SPY_XMAX) ClipPRight(); return TRUE; }
        case 0x62: return FALSE;
        case 0x64: return FALSE;
        case 0x65: return FALSE;
        case 0x66: return FALSE;
        case 0x68: ClipQBottom(); if (*x1 > SPY_XMAX) return FALSE; else
        { ClipPRight(); if (*y0 < SPY_YMIN) ClipPTop(); return TRUE; }
        case 0x69: ClipQLeft(); if (*y1 < SPY_YMIN) return FALSE; else
        { ClipPRight(); if (*y0 > SPY_YMAX) return FALSE; else
        { if (*y1 > SPY_YMAX) ClipQBottom(); if (*y0 < SPY_YMIN) ClipPTop();
        return TRUE; } }
        case 0x6A: return FALSE;
        case 0x80: ClipPBottom(); return TRUE;
        case 0x81: ClipPBottom(); if (*x0 < SPY_XMIN) return FALSE; else
        { ClipQLeft(); return TRUE; }
        case 0x82: ClipPBottom(); if (*x0 > SPY_XMAX) return FALSE; else
        { ClipQRight(); return TRUE; }
        case 0x84: ClipPBottom(); ClipQTop(); return TRUE;
        case 0x85: ClipPBottom(); if (*x0 < SPY_XMIN) return FALSE; else
        { ClipQLeft(); if (*y1 < SPY_YMIN) ClipQTop(); return TRUE; }
        case 0x86: ClipPBottom(); if (*x0 > SPY_XMAX) return FALSE; else
        { ClipQRight(); if (*y1 < SPY_YMIN) ClipQTop(); return TRUE; }
        case 0x88: return FALSE;
        case 0x89: return FALSE;
        case 0x8A: return FALSE;
        case 0x90: ClipPLeft(); if (*y0 > SPY_YMAX) ClipPBottom();
            return TRUE;
        case 0x91: return FALSE;
        case 0x92: ClipQRight(); if (*y1 > SPY_YMAX) return FALSE; else
        { ClipPBottom(); if (*x0 < SPY_XMIN) ClipPLeft(); return TRUE; }
        case 0x94: ClipQTop(); if (*x1 < SPY_XMIN) return FALSE; else
        { ClipPLeft(); if (*y0 > SPY_YMAX) ClipPBottom(); return TRUE; }
        case 0x95: return FALSE;
        case 0x96: ClipPLeft(); if (*y0 < SPY_YMIN) return FALSE; else
        { ClipQRight(); if (*y1 > SPY_YMAX) return FALSE; else
        { if (*y0 > SPY_YMAX) ClipPBottom(); if (*y1 < SPY_YMIN) ClipQTop();
        return TRUE; } }
        case 0x98: return FALSE;
        case 0x99: return FALSE;
        case 0x9A: return FALSE;
        case 0xA0: ClipPRight(); if (*y0 > SPY_YMAX) ClipPBottom();
            return TRUE;
        case 0xA1: ClipQLeft(); if (*y1 > SPY_YMAX) return FALSE; else
        { ClipPBottom(); if (*x0 > SPY_XMAX) ClipPRight(); return TRUE; }
        case 0xA2: return FALSE;
        case 0xA4: ClipQTop(); if (*x1 > SPY_XMAX) return FALSE; else
        { ClipPRight(); if (*y0 > SPY_YMAX) ClipPBottom(); return TRUE; }
        case 0xA5: ClipQLeft(); if (*y1 > SPY_YMAX) return FALSE; else
        { ClipPRight(); if (*y0 < SPY_YMIN) return FALSE; else
        { if (*y1 < SPY_YMIN) ClipQTop(); if ( *y0 > SPY_YMAX) ClipPBottom();
        return TRUE; } }
        case 0xA6: return FALSE;
        case 0xA8: return FALSE;
        case 0xA9: return FALSE;
        case 0xAA: return FALSE;
        default: break;
    }
    return FALSE;
}  /* end AX_SPY_LineClipWindow() */
#undef SPY_XMIN
#undef SPY_YMIN
#undef SPY_XMAX
#undef SPY_YMAX

#define SPY_XMIN  xmin
#define SPY_YMIN  ymin
#define SPY_XMAX  xmax
#define SPY_YMAX  ymax
bool AX_SPY_LineClipRectangle ( AX_Float *x0, AX_Float *y0,
                                AX_Float *x1, AX_Float *y1,
                                AX_Float xmin, AX_Float ymin,
                                AX_Float xmax, AX_Float ymax )
{
    int code = 0;

    if (*y1 > SPY_YMAX)  code |= 8;
    else if (*y1 < SPY_YMIN)  code |= 4;

    if (*x1 > SPY_XMAX)  code |= 2;
    else if (*x1 < SPY_XMIN)  code |= 1;

    if (*y0 > SPY_YMAX)  code |= 128;
    else if (*y0 < SPY_YMIN)  code |= 64;

    if (*x0 > SPY_XMAX)  code |= 32;
    else if (*x0 < SPY_XMIN)  code |= 16;

    switch (code)
    {
        case 0x00: return TRUE;
        case 0x01: ClipQLeft(); return TRUE;
        case 0x02: ClipQRight(); return TRUE;
        case 0x04: ClipQTop(); return TRUE;
        case 0x05: ClipQLeft(); if (*y1 < SPY_YMIN) ClipQTop(); return TRUE;
        case 0x06: ClipQRight(); if (*y1 < SPY_YMIN) ClipQTop(); return TRUE;
        case 0x08: ClipQBottom(); return TRUE;
        case 0x09: ClipQLeft(); if (*y1 > SPY_YMAX) ClipQBottom(); return TRUE;
        case 0x0A: ClipQRight(); if (*y1 > SPY_YMAX) ClipQBottom();
            return TRUE;
        case 0x10: ClipPLeft(); return TRUE;
        case 0x11: return FALSE;
        case 0x12: ClipPLeft(); ClipQRight(); return TRUE;
        case 0x14: ClipPLeft(); if (*y0 < SPY_YMIN) return FALSE; else
        { ClipQTop(); return TRUE; }
        case 0x15: return FALSE;
        case 0x16: ClipPLeft(); if (*y0 < SPY_YMIN) return FALSE; else
        { ClipQTop(); if (*x1 > SPY_XMAX) ClipQRight(); return TRUE; }
        case 0x18: ClipPLeft(); if (*y0 > SPY_YMAX) return FALSE; else
        { ClipQBottom(); return TRUE; }
        case 0x19: return FALSE;
        case 0x1A: ClipPLeft(); if (*y0 > SPY_YMAX) return FALSE; else
        { ClipQBottom(); if (*x1 > SPY_XMAX) ClipQRight(); return TRUE; }
        case 0x20: ClipPRight(); return TRUE;
        case 0x21: ClipPRight(); ClipQLeft(); return TRUE;
        case 0x22: return FALSE;
        case 0x24: ClipPRight(); if (*y0 < SPY_YMIN) return FALSE; else
        { ClipQTop(); return TRUE; }
        case 0x25: ClipPRight(); if (*y0 < SPY_YMIN) return FALSE; else
        { ClipQTop(); if (*x1 < SPY_XMIN) ClipQLeft(); return TRUE; }
        case 0x26: return FALSE;
        case 0x28: ClipPRight(); if (*y0 > SPY_YMAX) return FALSE; else
        { ClipQBottom(); return TRUE; }
        case 0x29: ClipPRight(); if (*y0 > SPY_YMAX) return FALSE; else
        { ClipQBottom(); if (*x1 < SPY_XMIN) ClipQLeft(); return TRUE; }
        case 0x2A: return FALSE;
        case 0x40: ClipPTop(); return TRUE;
        case 0x41: ClipPTop(); if (*x0 < SPY_XMIN) return FALSE; else
        { ClipQLeft(); if (*y1 < SPY_YMIN) ClipQTop(); return TRUE; }
        case 0x42: ClipPTop(); if (*x0 > SPY_XMAX) return FALSE; else
        { ClipQRight(); return TRUE; }
        case 0x44: return FALSE;
        case 0x45: return FALSE;
        case 0x46: return FALSE;
        case 0x48: ClipPTop(); ClipQBottom(); return TRUE;
        case 0x49: ClipPTop(); if (*x0 < SPY_XMIN) return FALSE; else
          { ClipQLeft(); if (*y1 > SPY_YMAX) ClipQBottom(); return TRUE; }
        case 0x4A: ClipPTop(); if (*x0 > SPY_XMAX) return FALSE; else
        { ClipQRight(); if (*y1 > SPY_YMAX) ClipQBottom(); return TRUE; }
        case 0x50: ClipPLeft(); if (*y0 < SPY_YMIN) ClipPTop(); return TRUE;
        case 0x51: return FALSE;
        case 0x52: ClipQRight(); if (*y1 < SPY_YMIN) return FALSE; else
        { ClipPTop(); if (*x0 < SPY_XMIN) ClipPLeft(); return TRUE; }
        case 0x54: return FALSE;
        case 0x55: return FALSE;
        case 0x56: return FALSE;
        case 0x58: ClipQBottom(); if (*x1 < SPY_XMIN) return FALSE; else
        { ClipPTop(); if (*x0 < SPY_XMIN) ClipPLeft(); return TRUE; }
        case 0x59: return FALSE;
        case 0x5A: ClipPLeft(); if (*y0 > SPY_YMAX) return FALSE; else
        { ClipQRight(); if (*y1 < SPY_YMIN) return FALSE; else
        { if (*y0 < SPY_YMIN) ClipPTop(); if (*y1 > SPY_YMAX) ClipQBottom();
        return TRUE; } }
        case 0x60: ClipPRight(); if (*y0 < SPY_YMIN) ClipPTop(); return TRUE;
        case 0x61: ClipQLeft(); if (*y1 < SPY_YMIN) return FALSE; else
        { ClipPTop(); if (*x0 > SPY_XMAX) ClipPRight(); return TRUE; }
        case 0x62: return FALSE;
        case 0x64: return FALSE;
        case 0x65: return FALSE;
        case 0x66: return FALSE;
        case 0x68: ClipQBottom(); if (*x1 > SPY_XMAX) return FALSE; else
        { ClipPRight(); if (*y0 < SPY_YMIN) ClipPTop(); return TRUE; }
        case 0x69: ClipQLeft(); if (*y1 < SPY_YMIN) return FALSE; else
        { ClipPRight(); if (*y0 > SPY_YMAX) return FALSE; else
        { if (*y1 > SPY_YMAX) ClipQBottom(); if (*y0 < SPY_YMIN) ClipPTop();
        return TRUE; } }
        case 0x6A: return FALSE;
        case 0x80: ClipPBottom(); return TRUE;
        case 0x81: ClipPBottom(); if (*x0 < SPY_XMIN) return FALSE; else
        { ClipQLeft(); return TRUE; }
        case 0x82: ClipPBottom(); if (*x0 > SPY_XMAX) return FALSE; else
        { ClipQRight(); return TRUE; }
        case 0x84: ClipPBottom(); ClipQTop(); return TRUE;
        case 0x85: ClipPBottom(); if (*x0 < SPY_XMIN) return FALSE; else
        { ClipQLeft(); if (*y1 < SPY_YMIN) ClipQTop(); return TRUE; }
        case 0x86: ClipPBottom(); if (*x0 > SPY_XMAX) return FALSE; else
        { ClipQRight(); if (*y1 < SPY_YMIN) ClipQTop(); return TRUE; }
        case 0x88: return FALSE;
        case 0x89: return FALSE;
        case 0x8A: return FALSE;
        case 0x90: ClipPLeft(); if (*y0 > SPY_YMAX) ClipPBottom();
            return TRUE;
        case 0x91: return FALSE;
        case 0x92: ClipQRight(); if (*y1 > SPY_YMAX) return FALSE; else
        { ClipPBottom(); if (*x0 < SPY_XMIN) ClipPLeft(); return TRUE; }
        case 0x94: ClipQTop(); if (*x1 < SPY_XMIN) return FALSE; else
        { ClipPLeft(); if (*y0 > SPY_YMAX) ClipPBottom(); return TRUE; }
        case 0x95: return FALSE;
        case 0x96: ClipPLeft(); if (*y0 < SPY_YMIN) return FALSE; else
        { ClipQRight(); if (*y1 > SPY_YMAX) return FALSE; else
        { if (*y0 > SPY_YMAX) ClipPBottom(); if (*y1 < SPY_YMIN) ClipQTop();
        return TRUE; } }
        case 0x98: return FALSE;
        case 0x99: return FALSE;
        case 0x9A: return FALSE;
        case 0xA0: ClipPRight(); if (*y0 > SPY_YMAX) ClipPBottom();
            return TRUE;
        case 0xA1: ClipQLeft(); if (*y1 > SPY_YMAX) return FALSE; else
        { ClipPBottom(); if (*x0 > SPY_XMAX) ClipPRight(); return TRUE; }
        case 0xA2: return FALSE;
        case 0xA4: ClipQTop(); if (*x1 > SPY_XMAX) return FALSE; else
        { ClipPRight(); if (*y0 > SPY_YMAX) ClipPBottom(); return TRUE; }
        case 0xA5: ClipQLeft(); if (*y1 > SPY_YMAX) return FALSE; else
        { ClipPRight(); if (*y0 < SPY_YMIN) return FALSE; else
        { if (*y1 < SPY_YMIN) ClipQTop(); if ( *y0 > SPY_YMAX) ClipPBottom();
        return TRUE; } }
        case 0xA6: return FALSE;
        case 0xA8: return FALSE;
        case 0xA9: return FALSE;
        case 0xAA: return FALSE;
        default: break;
    }
    return FALSE;
}  /* end AX_SPY_LineClipRectangle() */

#ifdef _SPY_LineClipRectangle_TEST
#define SPY_WIDTH  200
#define SPY_HEIGHT 150
#define SPY_TRIALS 10
int main (int argc, char *argv[])
{
    int i;
    AX_Float x0,y0,x1,y1,wx,wy;
    AX_openwindow (getpid(),NULL,AX_DEFWIDTH,AX_DEFHEIGHT);
    AXSetNamedForeground(0,AX_GREEN);
    wx = (AX_DEFWIDTH - SPY_WIDTH) / 2;
    wy = (AX_DEFHEIGHT - SPY_HEIGHT) / 2;
    AXDrawRectangle(0,wx,wy,SPY_WIDTH,SPY_HEIGHT);
    TimeRandomize();
    for (i=0; i<SPY_TRIALS; i++)
    {
        x0 = Frandom() * AX_DEFWIDTH;
        y0 = Frandom() * AX_DEFHEIGHT;
        x1 = Frandom() * AX_DEFWIDTH;
        y1 = Frandom() * AX_DEFHEIGHT;
        AXSetNamedForeground(0,AX_WHITE);
        AXDrawLine(0,x0,y0,x1,y1);
        if ( AX_SPY_LineClipRectangle(&x0,&y0,&x1,&y1,wx,wy,
                                      wx+SPY_WIDTH,wy+SPY_HEIGHT) )
        {
            AXSetNamedForeground(0,AX_RED);
            AXDrawLine(0,x0,y0,x1,y1);
        }
    }
    AX_dump(0); AX_show(0);
    Press_return();
    return (0);
}
#endif /* _SPY_LineClipRectangle_TEST */


#ifdef _LineClipContest_TEST
#define LineClipContest_WIDTH   (AX_DEFWIDTH / 2)
#define LineClipContest_HEIGHT  (AX_DEFHEIGHT / 2)
#define LineClipContest_TRIALS  3000000
int main (int argc, char *argv[])
{
    int i, j;
    AX_Float X0,Y0,X1,Y1,wx,wy,x0[LineClipContest_TRIALS],
        y0[LineClipContest_TRIALS],x1[LineClipContest_TRIALS],
        y1[LineClipContest_TRIALS];
    struct List
    {
	char *name;
	bool (*fun) (AX_Float *x0, AX_Float *y0, AX_Float *x1, AX_Float *y1,
                     AX_Float width, AX_Float height );
    } funs[] = { 
        {"Cohen-Sutherland", &AX_CS_LineClipWindow},
        {"Sobkow-Pospisil-Yang", &AX_SPY_LineClipWindow},
    };
    wx = (AX_DEFWIDTH  - LineClipContest_WIDTH)  / 2;
    wy = (AX_DEFHEIGHT - LineClipContest_HEIGHT) / 2;
    start_chronometer();
    for (j=0; j<LineClipContest_TRIALS; j++)
    {
        x0[j] = Frandom() * AX_DEFWIDTH  - wx;
        y0[j] = Frandom() * AX_DEFHEIGHT - wy;
        x1[j] = Frandom() * AX_DEFWIDTH  - wx;
        y1[j] = Frandom() * AX_DEFHEIGHT - wy;
    }
    stop_chronometer();
    printf ("preparation of %d line segments takes %lf seconds\n",
            LineClipContest_TRIALS, stopped_usertime() );
    for (i=0; i<sizeof(funs)/sizeof(struct List); i++)
    {
        printf ("Testing %s..", funs[i].name);
        start_chronometer();
        for (j=0; j<LineClipContest_TRIALS; j++)
        {
            X0 = x0[j]; Y0 = y0[j]; X1 = x1[j]; Y1 = y1[j];
            funs[i].fun ( &X0, &Y0, &X1, &Y1, 
                          LineClipContest_WIDTH, LineClipContest_HEIGHT);
        }
        stop_chronometer();
        printf (" %lf seconds\n", stopped_usertime());
    }
    return (0);
}
#endif /* _LineClipContest_TEST */


/* Convert directed edge v0->v1 to an infinite clipping boundary */
void AX_EdgeToClipBoundary (AX_Vertex *v0, AX_Vertex *v1, AX_ClipBoundary *c)
{
    if ( AX_SameVertex( *v0, *v1 ) )
        pe ("AX_EdgeToClipBoundary: edge is of zero length\n"
            "v0.x = v1.x = %f, v0.y = v1.y = %f\n", v0->x, v0->y);
    c->outNormal_x = v1->y - v0->y;
    c->outNormal_y = v0->x - v1->x;
    c->outNormal_self = AX_VertexOutwardDistance (*v0, *c);
    return;
} /* end AX_EdgeToClipBoundary() */


/* Given that "v0" is inside ClipBoundary "c" with "v1" */
/* outside, or vice versa, get the intersection point.  */
void AX_IntersectClipBoundary
(AX_Vertex *v0, AX_Vertex *v1, AX_ClipBoundary *c, AX_Vertex *intersection)
{
    AX_Float dx,dy,d0,d1,retraction;
    dx = v1->x - v0->x;
    dy = v1->y - v0->y;
    d0 = AX_VertexOutwardDistance (*v0, *c);
    d1 = AX_VertexOutwardDistance (*v1, *c);
    retraction = (d0-c->outNormal_self) / (d0-d1);
    intersection->x = v0->x + retraction * dx;
    intersection->y = v0->y + retraction * dy;
    return;
} /* end AX_IntersectClipBoundary() */


/* Calculate the intersection of two clip boundaries */
void AX_ClipIntersectClip
(AX_ClipBoundary *c0, AX_ClipBoundary *c1, AX_Vertex *intersection)
{
    AX_Float determinant;
    determinant = c0->outNormal_x * c1->outNormal_y -
        c1->outNormal_x * c0->outNormal_y;
    intersection->x = (c1->outNormal_y * c0->outNormal_self -
                       c0->outNormal_y * c1->outNormal_self) / determinant;
    intersection->y = (c0->outNormal_x * c1->outNormal_self -
                       c1->outNormal_x * c0->outNormal_self) / determinant;
    return;
} /* end AX_ClipIntersectClip() */


/* Allocate memory for a polygon (should be freed */
/* later!) and initialize its vertex coordinates. */
AX_Polygon *AX_NewPolygon (int nVertex, ...)
{
    AX_Polygon *p;
    int i;
    va_list ap;
    /* if (nVertex < 3) */
    /* pe ("AX_NewPolygon: nVertex = %d < 3\n", nVertex); */
    if (nVertex > AX_PolygonMaxVertex)
        pe ("AX_NewPolygon: nVertex = %d exceeds "
            "AX_PolygonMaxVertex = %d\n", nVertex,
            AX_PolygonMaxVertex);
    p = (AX_Polygon *) malloc(sizeof(AX_Polygon));
    p->nVertex = nVertex;
    va_start(ap, nVertex);
    for (i=0; i<nVertex; i++)
    {
        p->Vertex[i].x = va_arg(ap, double);
        p->Vertex[i].y = va_arg(ap, double);
    }
    va_end (ap);
    return (p);
} /* end AX_NewPolygon() */


/* (Re)assign polygon vertices */
AX_Polygon *AX_PolygonAssign (AX_Polygon *p, int nVertex, ...)
{
    int i;
    va_list ap;
    /* if (nVertex < 3) */
    /* pe ("AX_PolygonAssign: nVertex = %d < 3\n", nVertex); */
    if (nVertex > AX_PolygonMaxVertex)
        pe ("AX_PolygonAssign: nVertex = %d exceeds "
            "AX_PolygonMaxVertex = %d\n", nVertex,
            AX_PolygonMaxVertex);
    p->nVertex = nVertex;
    va_start(ap, nVertex);
    for (i=0; i<nVertex; i++)
    {
        p->Vertex[i].x = va_arg(ap, double);
        p->Vertex[i].y = va_arg(ap, double);
    }
    va_end (ap);
    return (p);
} /* end AX_PolygonAssign() */


AX_Polygon *AX_PolygonASSIGN (AX_Polygon *p, int nVertex, XPoint *V)
{
    int i;
    if (nVertex > AX_PolygonMaxVertex)
        pe ("AX_PolygonASSIGN: nVertex = %d exceeds "
            "AX_PolygonMaxVertex = %d\n", nVertex,
            AX_PolygonMaxVertex);
    p->nVertex = nVertex;
    for (i=0; i<nVertex; i++)
    {
        p->Vertex[i].x = V[i].x;
        p->Vertex[i].y = V[i].y;
    }
    return (p);
} /* end AX_PolygonASSIGN() */


/* from polygon vertices, get the set of clip boundaries */
void AX_PolygonToClipper (AX_Polygon *p, AX_Clipper *c)
{
    int i;
    c->nClipBoundary = p->nVertex;
    for (i=0; i<p->nVertex-1; i++)
        AX_EdgeToClipBoundary( p->Vertex + i,
                               p->Vertex + i + 1,
                               c->ClipBoundary + i);
    AX_EdgeToClipBoundary( p->Vertex + i,
                           p->Vertex,
                           c->ClipBoundary + i);
    return;
}  /* end AX_PolygonToClipper() */


/* from a set of clip boundaries, get polygon vertices */
void AX_ClipperToPolygon (AX_Clipper *c, AX_Polygon *p)
{
    int i;
    p->nVertex = c->nClipBoundary;
    for (i=0; i<c->nClipBoundary-1; i++)
        AX_ClipIntersectClip ( c->ClipBoundary + i,
                               c->ClipBoundary + i + 1,
                               p->Vertex + i + 1 );
    AX_ClipIntersectClip ( c->ClipBoundary + i,
                           c->ClipBoundary,
                           p->Vertex );
    return;
}  /* end AX_ClipperToPolygon() */


/* Allocate memory for a Clipper (should be freed */
/* later!) and initialize by vertex coordinates.  */
AX_Clipper *AX_NewClipper (int nVertex, ...)
{
    AX_Polygon p[1];
    AX_Clipper *c;
    int i;
    va_list ap;
    if (nVertex > AX_PolygonMaxVertex)
        pe ("AX_NewClipper: nVertex = %d exceeds "
            "AX_PolygonMaxVertex = %d\n", nVertex,
            AX_PolygonMaxVertex);
    p->nVertex = nVertex;
    va_start(ap, nVertex);
    for (i=0; i<nVertex; i++)
    {
        p->Vertex[i].x = va_arg(ap, double);
        p->Vertex[i].y = va_arg(ap, double);
    }
    va_end (ap);
    c = (AX_Clipper *) malloc(sizeof(AX_Clipper));
    AX_PolygonToClipper(p,c);
    return (c);
} /* end AX_NewClipper() */


/* (Re)assign Clipper by vertices */
AX_Clipper *AX_ClipperAssign (AX_Clipper *c, int nVertex, ...)
{
    AX_Polygon p[1];
    int i;
    va_list ap;
    if (nVertex > AX_PolygonMaxVertex)
        pe ("AX_ClipperAssign: nVertex = %d exceeds "
            "AX_PolygonMaxVertex = %d\n", nVertex,
            AX_PolygonMaxVertex);
    p->nVertex = nVertex;
    va_start(ap, nVertex);
    for (i=0; i<nVertex; i++)
    {
        p->Vertex[i].x = va_arg(ap, double);
        p->Vertex[i].y = va_arg(ap, double);
    }
    va_end (ap);
    AX_PolygonToClipper(p,c);
    return (c);
} /* end AX_ClipperAssign() */


AX_Clipper *AX_ClipperASSIGN (AX_Clipper *c, int nVertex, XPoint *V)
{
    AX_Polygon p[1];
    int i;
    if (nVertex > AX_PolygonMaxVertex)
        pe ("AX_ClipperASSIGN: nVertex = %d exceeds "
            "AX_PolygonMaxVertex = %d\n", nVertex,
            AX_PolygonMaxVertex);
    p->nVertex = nVertex;
    for (i=0; i<nVertex; i++)
    {
        p->Vertex[i].x = V[i].x;
        p->Vertex[i].y = V[i].y;
    }
    AX_PolygonToClipper(p,c);
    return (c);
} /* end AX_ClipperASSIGN() */


/* Test whether a vertex is inside Clipper */
bool AX_InsideClipper (AX_Vertex *v, AX_Clipper *c)
{
    int i;
    for (i=0; i<c->nClipBoundary; i++)
        if (!AX_VertexInsideClipBoundary(*v, c->ClipBoundary[i]))
            return FALSE;
    return TRUE;
}  /* end AX_InsideClipper() */

bool AX_INSIDEClipper (XPoint *V, AX_Clipper *c)
{
    int i;
    AX_Vertex v;
    v.x = V->x;
    v.y = V->y;
    for (i=0; i<c->nClipBoundary; i++)
        if (!AX_VertexInsideClipBoundary(v, c->ClipBoundary[i]))
            return FALSE;
    return TRUE;
}  /* end AX_INSIDEClipper() */


#ifdef _InsideClipper_TEST
#define InsideClipper_WIDTH  200
#define InsideClipper_HEIGHT 150
#define InsideClipper_TRIALS 100
#define InsideClipper_CC     4
int main (int argc, char *argv[])
{
    int i;
    AX_Float wx, wy;
    XPoint cc[InsideClipper_CC]; AX_Clipper c;
    XPoint pp[InsideClipper_TRIALS];

    AX_openwindow (getpid(),NULL,AX_DEFWIDTH,AX_DEFHEIGHT);

    wx = (AX_DEFWIDTH  - InsideClipper_WIDTH)  / 2;
    wy = (AX_DEFHEIGHT - InsideClipper_HEIGHT) / 2;
    cc[0].x = wx;            cc[0].y = wy;
    cc[1].x = wx+InsideClipper_WIDTH; cc[1].y = wy;
    cc[2].x = wx+InsideClipper_WIDTH; cc[2].y = wy+InsideClipper_HEIGHT;
    cc[3].x = wx;            cc[3].y = wy+InsideClipper_HEIGHT;
    AX_ClipperASSIGN (&c, InsideClipper_CC, cc);
    AXSetNamedForeground (0,AX_GREEN);
    AXDrawPolygon (0,cc,InsideClipper_CC);
    
    for (i=0; i<InsideClipper_TRIALS; i++)
    {
        pp[i].x = Frandom() * AX_DEFWIDTH;
        pp[i].y = Frandom() * AX_DEFHEIGHT;
        if (AX_INSIDEClipper (pp+i, &c))
        {
            AXSetNamedForeground (0,AX_GREEN);
            AXDrawPoint(0,pp[i].x,pp[i].y);
        }
        else
        {
            AXSetNamedForeground (0,AX_RED);
            AXDrawPoint(0,pp[i].x,pp[i].y);
        }
    }
    AX_dump(0); AX_show(0);
    Press_return();
    return (0);
}
#endif /* _InsideClipper_TEST */


/* Check if a polygon Clipper region is "feasible" */
bool AX_feasibleClipper (AX_Clipper *c)
{
    int i,j;
    AX_Polygon p[1];
    AX_ClipperToPolygon (c, p);
    for (i=0; i<c->nClipBoundary; i++)
        for (j=0; j<c->nClipBoundary; j++)
            if ( ! AX_VertexIncidentToEdge(j, i, c->nClipBoundary) )
                if ( ! AX_VertexInsideClipBoundary
                     (p->Vertex[j], c->ClipBoundary[i]) )
                    return FALSE;
    return TRUE;
}  /* end AX_InsideClipper() */


/* Sutherland-Hodgman polygon-polygon clipping algorithm */
void AX_SH_PolygonPolygonClip
(AX_Polygon *p, AX_Clipper *c, AX_Polygon *result)
{
    int i, j;
    AX_Polygon pp;
    AX_Polygon *from, *to, *tmp;
    AX_Vertex *v0, *v1, intersection[1];
    from = &pp;
    AX_PolygonCopy (*p, *from, i);
    to = result;
    for (i=0; i<c->nClipBoundary; i++)
    { /* pipelined operation: hardware implementation */
        to->nVertex = 0;
        /* start with the last vertex of the polygon */
        v0 = from->Vertex + from->nVertex - 1;
        for (j=0; j<from->nVertex; j++)
        {
            v1 = from->Vertex + j;
            if  ( AX_VertexInsideClipBoundary (*v1, c->ClipBoundary[i]) )
            { /* case 1 and 4 */
                if ( AX_VertexInsideClipBoundary (*v0, c->ClipBoundary[i]) )
                    AX_PolygonAddVertex (*v1, *to);  /* case 1 */ 
                else
                {  /* case 4 */
                    AX_IntersectClipBoundary
                        (v0, v1, c->ClipBoundary+i, intersection);
                    AX_PolygonAddVertex (*intersection, *to);
                    AX_PolygonAddVertex (*v1, *to);
                }
            }
            else  /* case 2 and 3 */
                if ( AX_VertexInsideClipBoundary (*v0, c->ClipBoundary[i]) ) 
                {  /* case 2 */
                    AX_IntersectClipBoundary
                        (v0, v1, c->ClipBoundary+i, intersection);
                    AX_PolygonAddVertex (*intersection, *to);
                }  /* no action for case 3 */
            v0 = v1;
        } /* for j */
        SWAP(from, to, tmp);
    } /* for i */
    if (from == &pp) AX_PolygonCopy(pp, *result, i);
    return;
} /* end AX_SH_PolygonPolygonClip() */

#ifdef _SH_PolygonPolygonClip_TEST
#define SH_CC_WIDTH  200
#define SH_CC_HEIGHT 150
#define SH_CC 3
#define SH_PP 3
int main (int argc, char *argv[])
{
    int i;
    AX_Float wx, wy;
    XPoint cc[AX_PolygonMaxVertex]; AX_Clipper c;
    XPoint pp[AX_PolygonMaxVertex]; AX_Polygon p;
    XPoint rr[AX_PolygonMaxVertex]; AX_Polygon result;

    AX_openwindow (getpid(),NULL,AX_DEFWIDTH,AX_DEFHEIGHT);
    TimeRandomize();

    wx = (AX_DEFWIDTH  - SH_CC_WIDTH)  / 2;
    wy = (AX_DEFHEIGHT - SH_CC_HEIGHT) / 2;
    do
    {
        for (i=0; i<SH_CC; i++)
        {
            AX_RandomVertex(cc[i], AX_DEFWIDTH, AX_DEFHEIGHT);
            printf ("%d %d %d\n", i, cc[i].x, cc[i].y);
        }
        AX_ClipperASSIGN (&c, SH_CC, cc);        
    } while ( ! AX_feasibleClipper(&c) );
    AXSetNamedForeground (0, AX_GREEN);
    AXDrawPolygon (0, cc, SH_CC);

    for (i=0; i<SH_PP; i++) AX_RandomVertex(pp[i], AX_DEFWIDTH, AX_DEFHEIGHT);
    AX_PolygonASSIGN (&p, SH_PP, pp);
    AXSetNamedForeground (0, AX_WHITE);
    AXDrawPolygon (0, pp, SH_PP);

    AX_SH_PolygonPolygonClip (&p, &c, &result);
    AX_PolygonToXPoints (result, rr, i);
    AXSetNamedForeground (0,AX_RED);
    AXDrawPolygon (0, rr, result.nVertex);

    AX_dump(0); AX_show(0);
    Press_return();
    return (0);
}
#endif /* _SH_PolygonPolygonClip_TEST */


#define SH_XMIN    0
#define SH_XMAX    width
#define SH_YMIN    0
#define SH_YMAX    height
#define SH_InsideXMIN(v) ( (v).x >= SH_XMIN )
#define SH_ClipXMIN(v0,v1,intersection) \
  ( (intersection).y = ((v1).y - (v0).y) * (SH_XMIN - (v0).x) / \
  ((v1).x - (v0).x) + (v0).y, (intersection).x = SH_XMIN)
#define SH_InsideXMAX(v) ( (v).x <= SH_XMAX )
#define SH_ClipXMAX(v0,v1,intersection) \
  ( (intersection).y = ((v1).y - (v0).y) * (SH_XMAX - (v0).x) / \
  ((v1).x - (v0).x) + (v0).y, (intersection).x = SH_XMAX)
#define SH_InsideYMIN(v) ( (v).y >= SH_YMIN )
#define SH_ClipYMIN(v0,v1,intersection) \
  ( (intersection).x = ((v1).x - (v0).x) * (SH_YMIN - (v0).y) / \
  ((v1).y - (v0).y) + (v0).x, (intersection).y = SH_YMIN)
#define SH_InsideYMAX(v) ( (v).y <= SH_YMAX )
#define SH_ClipYMAX(v0,v1,intersection) \
  ( (intersection).x = ((v1).x - (v0).x) * (SH_YMAX - (v0).y) / \
  ((v1).y - (v0).y) + (v0).x, (intersection).y = SH_YMAX)

/* Sutherland-Hodgman polygon-window clipping */
void AX_SH_PolygonWindowClip
(AX_Polygon *p, AX_Float width, AX_Float height, AX_Polygon *result)
{
    int j;
    AX_Polygon pp[1];
    AX_Vertex *v0, *v1, intersection[1];
    v0 = p->Vertex + p->nVertex - 1;
    pp->nVertex = 0;
    for (j=0; j<p->nVertex; j++)
    {
        v1 = p->Vertex + j;
        if  ( SH_InsideYMIN (*v1) )
        { /* case 1 and 4 */
            if ( SH_InsideYMIN (*v0) )
                AX_PolygonAddVertex (*v1, *pp);  /* case 1 */ 
            else
            {  /* case 4 */
                SH_ClipYMIN (*v0, *v1, *intersection);
                AX_PolygonAddVertex (*intersection, *pp);
                AX_PolygonAddVertex (*v1, *pp);
            }
        }
        else  /* case 2 and 3 */
            if ( SH_InsideYMIN (*v0) ) 
            {  /* case 2 */
                SH_ClipYMIN (*v0, *v1, *intersection);
                AX_PolygonAddVertex (*intersection, *pp);
            }  /* no action for case 3 */
        v0 = v1;
    }  /* bottom */

    v0 = pp->Vertex + pp->nVertex - 1;
    result->nVertex = 0;
    for (j=0; j<pp->nVertex; j++)
    {
        v1 = pp->Vertex + j;
        if  ( SH_InsideXMAX (*v1) )
        { /* case 1 and 4 */
            if ( SH_InsideXMAX (*v0) )
                AX_PolygonAddVertex (*v1, *result);  /* case 1 */ 
            else
            {  /* case 4 */
                SH_ClipXMAX (*v0, *v1, *intersection);
                AX_PolygonAddVertex (*intersection, *result);
                AX_PolygonAddVertex (*v1, *result);
            }
        }
        else  /* case 2 and 3 */
            if ( SH_InsideXMAX (*v0) ) 
            {  /* case 2 */
                SH_ClipXMAX (*v0, *v1, *intersection);
                AX_PolygonAddVertex (*intersection, *result);
            }  /* no action for case 3 */
        v0 = v1;
    }  /* right */

    v0 = result->Vertex + result->nVertex - 1;
    pp->nVertex = 0;
    for (j=0; j<result->nVertex; j++)
    {
        v1 = result->Vertex + j;
        if  ( SH_InsideYMAX (*v1) )
        { /* case 1 and 4 */
            if ( SH_InsideYMAX (*v0) )
                AX_PolygonAddVertex (*v1, *pp);  /* case 1 */ 
            else
            {  /* case 4 */
                SH_ClipYMAX (*v0, *v1, *intersection);
                AX_PolygonAddVertex (*intersection, *pp);
                AX_PolygonAddVertex (*v1, *pp);
            }
        }
        else  /* case 2 and 3 */
            if ( SH_InsideYMAX (*v0) ) 
            {  /* case 2 */
                SH_ClipYMAX (*v0, *v1, *intersection);
                AX_PolygonAddVertex (*intersection, *pp);
            }  /* no action for case 3 */
        v0 = v1;
    }  /* top */

    v0 = pp->Vertex + pp->nVertex - 1;
    result->nVertex = 0;
    for (j=0; j<pp->nVertex; j++)
    {
        v1 = pp->Vertex + j;
        if  ( SH_InsideXMIN (*v1) )
        { /* case 1 and 4 */
            if ( SH_InsideXMIN (*v0) )
                AX_PolygonAddVertex (*v1, *result);  /* case 1 */ 
            else
            {  /* case 4 */
                SH_ClipXMIN (*v0, *v1, *intersection);
                AX_PolygonAddVertex (*intersection, *result);
                AX_PolygonAddVertex (*v1, *result);
            }
        }
        else  /* case 2 and 3 */
            if ( SH_InsideXMIN (*v0) ) 
            {  /* case 2 */
                SH_ClipXMIN (*v0, *v1, *intersection);
                AX_PolygonAddVertex (*intersection, *result);
            }  /* no action for case 3 */
        v0 = v1;
    }  /* left */
    return;
} /* end AX_SH_PolygonWindowClip() */

#ifdef _SH_PolygonWindowClip_TEST
#define SHW_CC_WIDTH  200
#define SHW_CC_HEIGHT 150
#define SHW_CC 4
#define SHW_PP 4
int main (int argc, char *argv[])
{
    int i;
    AX_Float wx, wy;
    XPoint cc[AX_PolygonMaxVertex];
    XPoint pp[AX_PolygonMaxVertex]; AX_Polygon p;
    XPoint rr[AX_PolygonMaxVertex]; AX_Polygon result;
    AX_openwindow (getpid(),NULL,AX_DEFWIDTH,AX_DEFHEIGHT);
    wx = (AX_DEFWIDTH  - SHW_CC_WIDTH)  / 2;
    wy = (AX_DEFHEIGHT - SHW_CC_HEIGHT) / 2;
    AX_VertexAssign ( cc[0], wx, wy );
    AX_VertexAssign ( cc[1], wx+SHW_CC_WIDTH, wy );
    AX_VertexAssign ( cc[2], wx+SHW_CC_WIDTH, wy+SHW_CC_HEIGHT );
    AX_VertexAssign ( cc[3], wx, wy+SHW_CC_HEIGHT );
    AXSetNamedForeground (0, AX_GREEN);
    AXDrawPolygon (0, cc, SHW_CC);
    TimeRandomize();
    for (i=0; i<SHW_PP; i++) AX_RandomVertex(pp[i], AX_DEFWIDTH, AX_DEFHEIGHT);
    AX_PolygonASSIGN (&p, SHW_PP, pp);
    AXSetNamedForeground (0, AX_WHITE);
    AXDrawPolygon (0, pp, SHW_PP);
    for (i=0; i<SHW_PP; i++) AX_VertexShift(p.Vertex[i], -wx, -wy);
    AX_SH_PolygonWindowClip (&p, SHW_CC_WIDTH, SHW_CC_HEIGHT, &result);
    AX_PolygonToXPoints (result, rr, i);
    for (i=0; i<result.nVertex; i++) AX_VertexShift(rr[i], wx, wy);
    AXSetNamedForeground (0,AX_RED);
    AXDrawPolygon (0, rr, result.nVertex);
    AX_dump(0); AX_show(0);
    Press_return();
    return (0);
}
#endif /* _SH_PolygonWindowClip_TEST */


/* Check if segment (xmin,y)-(xmax,y) has finite */
/* length intersection with circle (x0,y0,r).    */
bool AX_HorizontalSegmentIntersectCircle
(AX_Float y, AX_Float xmin, AX_Float xmax,
 AX_Float x0, AX_Float y0, AX_Float radius)
{
    register AX_Float dy, Xmin, Xmax;
    dy = y - y0;
    if ( (dy <= -radius) || (dy >= radius) ) return (FALSE);
    Xmax = sqrt(radius*radius-dy*dy);
    Xmin = x0 - Xmax;
    Xmax += x0;
    if (Xmax > xmax) Xmax = xmax;
    if (Xmin < xmin) Xmin = xmin;
    return (Xmax > Xmin);
} /* end AX_HorizontalSegmentIntersectCircle() */


/* Check if segment (x,ymin)-(x,ymax) has finite */
/* length intersection with circle (x0,y0,r).    */
bool AX_VerticalSegmentIntersectCircle
(AX_Float x, AX_Float ymin, AX_Float ymax,
 AX_Float x0, AX_Float y0, AX_Float radius)
{
    register AX_Float dx, Ymin, Ymax;
    dx = x - x0;
    if ( (dx <= -radius) || (dx >= radius) ) return (FALSE);
    Ymax = sqrt(radius*radius-dx*dx);
    Ymin = y0 - Ymax;
    Ymax += y0;
    if (Ymax > ymax) Ymax = ymax;
    if (Ymin < ymin) Ymin = ymin;
    return (Ymax > Ymin);
} /* end AX_VerticalSegmentIntersectCircle() */


#ifdef _CircleIntersectWindow_TEST
#define x0     200
#define y0     200
#define radius 142
#define width  100
#define height 100
int main (int argc, char *argv[])
{
    printf ("%d\n",AX_CircleIntersectWindow(x0,y0,radius,width,height));
    return (0);
}
#endif /* _CircleIntersectWindow_TEST */


/* Calculate length of line segment (xmin,y)-(xmax,y) in circle (x0,y0,r) */
AX_Float AX_HorizontalSegmentInCircle
(AX_Float y, AX_Float xmin, AX_Float xmax,
 AX_Float x0, AX_Float y0, AX_Float radius2)
{
    register AX_Float dy2, Xmin, Xmax;
    dy2 = y - y0;
    dy2 *= dy2;
    if ( dy2 >= radius2 ) return(0.);
    Xmax = sqrt(radius2-dy2);
    Xmin = x0 - Xmax;
    Xmax += x0;
    if (Xmax > xmax) Xmax = xmax;
    if (Xmin < xmin) Xmin = xmin;
    if (Xmax > Xmin) return (Xmax-Xmin);
    else return(0.);
} /* end AX_HorizontalSegmentInCircle() */


/* Calculate length of line segment (x,ymin)-(x,ymax) in circle (x0,y0,r) */
AX_Float AX_VerticalSegmentInCircle
(AX_Float x, AX_Float ymin, AX_Float ymax,
 AX_Float x0, AX_Float y0, AX_Float radius2)
{
    register AX_Float dx2, Ymin, Ymax;
    dx2 = x - x0;
    dx2 *= dx2;
    if ( dx2 >= radius2 ) return(0.);
    Ymax = sqrt(radius2-dx2);
    Ymin = y0 - Ymax;
    Ymax += y0;
    if (Ymax > ymax) Ymax = ymax;
    if (Ymin < ymin) Ymin = ymin;
    if (Ymax > Ymin) return (Ymax-Ymin);
    else return(0.);
} /* end AX_VerticalSegmentInCircle() */


/* If segment (*xmin,y)-(*xmax,y) intersects circle (x0,y0,r), */
/* revise *xmin,*xmax and return TRUE; otherwise return FALSE. */
bool AX_HorizontalSegmentClipCircle
(AX_Float y, AX_Float *xmin, AX_Float *xmax,
 AX_Float x0, AX_Float y0, AX_Float radius2)
{
    register AX_Float dy2, Xmin, Xmax;
    dy2 = y - y0;
    dy2 *= dy2;
    if ( dy2 > radius2 ) return (FALSE);
    Xmax = sqrt(radius2-dy2);
    Xmin = x0 - Xmax;
    Xmax += x0;
    if (*xmax > Xmax) *xmax = Xmax;
    if (*xmin < Xmin) *xmin = Xmin; 
    return (*xmax >= *xmin);
} /* end AX_HorizontalSegmentClipCircle() */


/* If segment (x,*ymin)-(x,*ymax) intersects circle (x0,y0,r), */
/* revise *ymin,*ymax and return TRUE; otherwise return FALSE. */
bool AX_VerticalSegmentClipCircle
(AX_Float x, AX_Float *ymin, AX_Float *ymax,
 AX_Float x0, AX_Float y0, AX_Float radius2)
{
    register AX_Float dx2, Ymin, Ymax;
    dx2 = x - x0;
    dx2 *= dx2;
    if ( dx2 > radius2 ) return (FALSE);
    Ymax = sqrt(radius2-dx2);
    Ymin = y0 - Ymax;
    Ymax += y0;
    if (*ymax > Ymax) *ymax = Ymax;
    if (*ymin < Ymin) *ymin = Ymin; 
    return (*ymax >= *ymin);
} /* end AX_VerticalSegmentClipCircle() */


/* same as AX_HorizontalSegmentClipCircle() except ensured to cut */
static void AX_HorizontalSegmentCutCircle
(AX_Float y, AX_Float *xmin, AX_Float *xmax,
 AX_Float x0, AX_Float y0, AX_Float radius2)
{
    register AX_Float dy2, Xmin, Xmax;
    dy2 = y - y0;
    dy2 *= dy2;
    Xmax = sqrt(radius2-dy2);
    Xmin = x0 - Xmax;
    Xmax += x0;
    if (*xmax > Xmax) *xmax = Xmax;
    if (*xmin < Xmin) *xmin = Xmin; 
    return;
} /* end AX_HorizontalSegmentCutCircle() */


/* same as AX_VerticalSegmentClipCircle() except ensured to cut */
static void AX_VerticalSegmentCutCircle
(AX_Float x, AX_Float *ymin, AX_Float *ymax,
 AX_Float x0, AX_Float y0, AX_Float radius2)
{
    register AX_Float dx2, Ymin, Ymax;
    dx2 = x - x0;
    dx2 *= dx2;
    Ymax = sqrt(radius2-dx2);
    Ymin = y0 - Ymax;
    Ymax += y0;
    if (*ymax > Ymax) *ymax = Ymax;
    if (*ymin < Ymin) *ymin = Ymin; 
    return;
} /* end AX_VerticalSegmentCutCircle() */


/* Ellipse */

#define AX_EllipseEquation(e) ( (e).D0  = (e).a2 * (e).b2, \
  (e).Dxx = (e).b2 * (e).ax * (e).ax + (e).a2 * (e).ay * (e).ay, \
  (e).Dxy = ( (e).b2 - (e).a2 ) * (e).ax * (e).ay, \
  (e).Dyy = (e).b2 * (e).ay * (e).ay + (e).a2 * (e).ax * (e).ax, \
  (e).ymax = sqrt( (e).Dxx * (e).D0 / ((e).Dxx*(e).Dyy - (e).Dxy*(e).Dxy) ), \
  (e).ymin_x = (e).Dxy * (e).ymax / (e).Dxx, \
  (e).ymax_x = (e).x0 - (e).ymin_x, (e).ymin_x += (e).x0, \
  (e).ymin = (e).y0 - (e).ymax, (e).ymax += (e).y0 )

/* Specify an ellipse by major axis displacements and minor axis length */
AX_Ellipse *AX_EllipseAssign (AX_Float x0, AX_Float y0, AX_Float ax,
                              AX_Float ay, AX_Float b, AX_Ellipse *e)
{
    e->x0 = x0;
    e->y0 = y0;
    e->a2 = ax*ax + ay*ay;
    e->a  = sqrt(e->a2);
    e->ax = ax / e->a;
    e->ay = ay / e->a;
    e->b  = b;
    e->b2 = b * b;
    AX_EllipseEquation (*e);
    return (e);
} /* end AX_EllipseAssign() */


/* Specify an ellipse by major axis angle (in radians) */
/* with x-axis and major, minor axis lengths.          */
AX_Ellipse *AX_EllipseASSIGN (AX_Float x0, AX_Float y0, AX_Float angle,
                              AX_Float a, AX_Float b, AX_Ellipse *e)
{
    e->x0 = x0;
    e->y0 = y0;
    e->ax = cos (angle);
    e->ay = sin (angle);
    e->a  = a;
    e->a2 = a * a;
    e->b  = b;
    e->b2 = b * b;
    AX_EllipseEquation (*e);
    return (e);
} /* end AX_EllipseASSIGN() */


/* Specify an ellipse by A*(x-x0)^2+2*B*(x-x0)*(y-y0)+C*(y-y0)^2=1 */
AX_Ellipse *AX_Ellipseassign (AX_Float x0, AX_Float y0, AX_Float A,
                              AX_Float B, AX_Float C, AX_Ellipse *e)
{
    register AX_Float D;
    e->x0 = x0;
    e->y0 = y0;
    if (B == 0)
    {
        e->ax = 1;
        e->ay = 0;
        e->a2 = 1 / A;
        e->a  = sqrt(e->a2);
        e->b2 = 1 / C;
        e->b  = sqrt(e->b2);
    }
    else
    {
        D = A-C;
        D = sqrt(D*D+4*B*B);
        e->b2 = 2/(A+C-D);
        e->b = sqrt(e->b2);
        D = (A+C+D)/2;
        e->a2 = 1/D;
        e->a  = sqrt(e->a2);
        e->ay = (D - A) / B;
        e->ax = 1 / sqrt(1 + e->ay * e->ay);
        e->ay *= e->ax;
    }
    AX_EllipseEquation (*e);
    return (e);
} /* end AX_Ellipseassign() */


#define A  (e->Dxx)
/* Check whether segment (xmin,y)-(xmax,y) has */
/* finite length intersection with ellipse.    */
bool AX_HorizontalSegmentIntersectEllipse
(AX_Float y, AX_Float xmin, AX_Float xmax, AX_Ellipse *e)
{
    register AX_Float dy, B, C, D, Xmin, Xmax;
    dy = y - e->y0;
    B = e->Dxy * dy;
    C = e->Dyy * dy * dy - e->D0;
    D = B * B - A * C;
    if (D <= 0) return (FALSE);
    D = sqrt(D) / A;
    C = e->x0 - B / A;
    Xmax = C + D;
    Xmin = C - D;
    if (Xmax > xmax) Xmax = xmax;
    if (Xmin < xmin) Xmin = xmin;
    return (Xmax > Xmin);
} /* end AX_HorizontalSegmentIntersectEllipse() */


/* Calculate length of line segment (xmin,y)-(xmax,y) in ellipse */
AX_Float AX_HorizontalSegmentInEllipse
(AX_Float y, AX_Float xmin, AX_Float xmax, AX_Ellipse *e)
{
    register AX_Float dy, B, C, D, Xmin, Xmax;
    dy = y - e->y0;
    B = e->Dxy * dy;
    C = e->Dyy * dy * dy - e->D0;
    D = B * B - A * C;
    if (D <= 0) return(0.);
    D = sqrt(D) / A;
    C = e->x0 - B / A;
    Xmax = C + D;
    Xmin = C - D;
    if (Xmax > xmax) Xmax = xmax;
    if (Xmin < xmin) Xmin = xmin;
    if (Xmax > Xmin) return (Xmax-Xmin);
    else return(0.);
} /* end AX_HorizontalSegmentInEllipse() */


/* If segment (*xmin,y)-(*xmax,y) intersects ellipse, revise */
/* *xmin, *xmax and return TRUE; otherwise return FALSE.     */
bool AX_HorizontalSegmentClipEllipse
(AX_Float y, AX_Float *xmin, AX_Float *xmax, AX_Ellipse *e)
{
    register AX_Float dy, B, C, D, Xmin, Xmax;
    dy = y - e->y0;
    B = e->Dxy * dy;
    C = e->Dyy * dy * dy - e->D0;
    D = B * B - A * C;
    if (D < 0) return (FALSE);
    D = sqrt(D) / A;
    C = e->x0 - B / A;
    Xmax = C + D;
    Xmin = C - D;
    if (*xmax > Xmax) *xmax = Xmax;
    if (*xmin < Xmin) *xmin = Xmin;
    return (*xmax >= *xmin);
} /* end AX_HorizontalSegmentClipEllipse() */


/* same as AX_HorizontalSegmentClipEllipse() except ensured to cut */
static void AX_HorizontalSegmentCutEllipse
(AX_Float y, AX_Float *xmin, AX_Float *xmax, AX_Ellipse *e)
{
    register AX_Float dy, B, C, D, Xmin, Xmax;
    dy = y - e->y0;
    B = e->Dxy * dy;
    C = e->Dyy * dy * dy - e->D0;
    D = B * B - A * C;
    if (D > 0) D = sqrt(D) / A; else D = 0;
    C = e->x0 - B / A;
    Xmax = C + D;
    Xmin = C - D;
    if (*xmax > Xmax) *xmax = Xmax;
    if (*xmin < Xmin) *xmin = Xmin;
    return;
} /* end AX_HorizontalSegmentCutEllipse() */
#undef A


#ifdef _HorizontalSegmentClipEllipse_TEST
#define x0     0
#define y0     0
#define angle  90
#define a      2
#define b      1
#define y      0
#define XMIN  -2
#define XMAX   0
int main (int argc, char *argv[])
{
    AX_Float xmin=XMIN, xmax=XMAX;
    AX_Ellipse e[1];
    AX_ELLIPSEASSIGN(x0,y0,angle,a,b,e);
    printf ("%f %f %f %f\n", e->ymin, e->ymax, e->ymin_x, e->ymax_x);
    if ( AX_HorizontalSegmentIntersectEllipse (y, XMIN, XMAX, e) )
    {
        AX_HorizontalSegmentClipEllipse (y, &xmin, &xmax, e);
        printf ("%f %f\n", xmin, xmax);
        printf ("%f\n", AX_HorizontalSegmentInEllipse (y, XMIN, XMAX, e));
    }
    return (0);
}
#undef x0
#undef y0
#undef angle
#undef a
#undef b
#undef y
#undef XMIN
#undef XMAX
#endif /* _HorizontalSegmentClipEllipse_TEST */


#define A  (e->Dyy)
/* Check whether segment (x,ymin)-(x,ymax) has */
/* finite length intersection with ellipse.    */
bool AX_VerticalSegmentIntersectEllipse
(AX_Float x, AX_Float ymin, AX_Float ymax, AX_Ellipse *e)
{
    register AX_Float dx, B, C, D, Ymin, Ymax;
    dx = x - e->x0;
    B = e->Dxy * dx;
    C = e->Dxx * dx * dx - e->D0;
    D = B * B - A * C;
    if (D <= 0) return (FALSE);
    D = sqrt(D) / A;
    C = e->y0 - B / A;
    Ymax = C + D;
    Ymin = C - D;
    if (Ymax > ymax) Ymax = ymax;
    if (Ymin < ymin) Ymin = ymin;
    return (Ymax > Ymin);
} /* end AX_VerticalSegmentIntersectEllipse() */


/* Calculate length of line segment (x,ymin)-(x,ymax) in ellipse */
AX_Float AX_VerticalSegmentInEllipse
(AX_Float x, AX_Float ymin, AX_Float ymax, AX_Ellipse *e)
{
    register AX_Float dx, B, C, D, Ymin, Ymax;
    dx = x - e->x0;
    B = e->Dxy * dx;
    C = e->Dxx * dx * dx - e->D0;
    D = B * B - A * C;
    if (D <= 0) return(0.);
    D = sqrt(D) / A;
    C = e->y0 - B / A;
    Ymax = C + D;
    Ymin = C - D;
    if (Ymax > ymax) Ymax = ymax;
    if (Ymin < ymin) Ymin = ymin;
    if (Ymax > Ymin) return (Ymax-Ymin);
    else return(0.);
} /* end AX_VerticalSegmentInEllipse() */


/* If segment (x,*ymin)-(x,*ymax) intersects ellipse, revise */
/* *ymin, *ymax and return TRUE; otherwise return FALSE.     */
bool AX_VerticalSegmentClipEllipse
(AX_Float x, AX_Float *ymin, AX_Float *ymax, AX_Ellipse *e)
{
    register AX_Float dx, B, C, D, Ymin, Ymax;
    dx = x - e->x0;
    B = e->Dxy * dx;
    C = e->Dxx * dx * dx - e->D0;
    D = B * B - A * C;
    if (D < 0) return (FALSE);
    D = sqrt(D) / A;
    C = e->y0 - B / A;
    Ymax = C + D;
    Ymin = C - D;
    if (*ymax > Ymax) *ymax = Ymax;
    if (*ymin < Ymin) *ymin = Ymin;
    return (*ymax >= *ymin);
} /* end AX_VerticalSegmentClipEllipse() */


/* same as AX_VerticalSegmentClipEllipse() except ensured to cut */
static void AX_VerticalSegmentCutEllipse
(AX_Float x, AX_Float *ymin, AX_Float *ymax, AX_Ellipse *e)
{
    register AX_Float dx, B, C, D, Ymin, Ymax;
    dx = x - e->x0;
    B = e->Dxy * dx;
    C = e->Dxx * dx * dx - e->D0;
    D = B * B - A * C;
    if (D > 0) D = sqrt(D) / A; else D = 0;
    C = e->y0 - B / A;
    Ymax = C + D;
    Ymin = C - D;
    if (*ymax > Ymax) *ymax = Ymax;
    if (*ymin < Ymin) *ymin = Ymin;
    return;
} /* end AX_VerticalSegmentCutEllipse() */
#undef A


/***************/
/* Color Field */
/***************/

/* Assert [0,1] limit on color field brightnesses */
AX_Carrier AX_safe_gradientcolor (AX_GradientColor *c, AX_Float x, AX_Float y)
{
    register AX_Float r,g,b;
    r = AX_extrapolate (c->r0, c->rx, c->ry, x, y);
    if (r < 0) r = 0; if (r > 1) r = 1;
    g = AX_extrapolate (c->g0, c->gx, c->gy, x, y);
    if (g < 0) g = 0; if (g > 1) g = 1;
    b = AX_extrapolate (c->b0, c->bx, c->by, x, y);
    if (b < 0) b = 0; if (b > 1) b = 1;
    return (AX_Colorcarrier(r,g,b));
} /* end AX_safe_gradientcolor() */


/* Deduce color field given two points: it assumes */
/* there is no variation in the normal direction.  */
void AX_LinearGradientColor
(AX_Float x0, AX_Float y0, AX_Float r0, AX_Float g0, AX_Float b0,
 AX_Float x1, AX_Float y1, AX_Float r1, AX_Float g1, AX_Float b1,
 AX_GradientColor *c)
{
    register AX_Float dx1, dy1, d, dc;
    if ( (x0 == x1) && (y0 == y1) )
    {
        AX_ZeroGradientColor(*c,(r0+r1)/2,(g0+g1)/2,(b0+b1)/2);
        return;
    }
    dx1 = x1 - x0;
    dy1 = y1 - y0;
    d = dx1*dx1 + dy1*dy1;
    dc = r1 - r0;
    c->rx = dc * dx1 / d;
    c->ry = dc * dy1 / d;
    c->r0 = r0 - c->rx * x0 - c->ry * y0;
    dc = g1 - g0;
    c->gx = dc * dx1 / d;
    c->gy = dc * dy1 / d;
    c->g0 = g0 - c->gx * x0 - c->gy * y0;
    dc = b1 - b0;
    c->bx = dc * dx1 / d;
    c->by = dc * dy1 / d;
    c->b0 = b0 - c->bx * x0 - c->by * y0;
    return;
} /* end AX_LinearGradientColor() */


/* deduce color field given three points */
void AX_TriangulateGradientColor
(AX_Float x0, AX_Float y0, AX_Float r0, AX_Float g0, AX_Float b0,
 AX_Float x1, AX_Float y1, AX_Float r1, AX_Float g1, AX_Float b1,
 AX_Float x2, AX_Float y2, AX_Float r2, AX_Float g2, AX_Float b2,
 AX_GradientColor *c)
{
    register AX_Float dx1, dy1, dx2, dy2, determinant, dc1, dc2;
    dx1 = x1 - x0;
    dy1 = y1 - y0;
    dx2 = x2 - x0;
    dy2 = y2 - y0;
    determinant = dx1 * dy2 - dx2 * dy1;
    if (determinant == 0)
        AX_LinearGradientColor (x0,y0,r0,g0,b0, (x1+x2)/2, (y1+y2)/2,
                                (r1+r2)/2, (g1+g2)/2, (b1+b2)/2, c);
    dc1 = r1 - r0;
    dc2 = r2 - r0;
    c->rx = (dy2 * dc1 - dy1 * dc2) / determinant;
    c->ry = (dx1 * dc2 - dx2 * dc1) / determinant;
    c->r0 = r0 - c->rx * x0 - c->ry * y0;
    dc1 = g1 - g0;
    dc2 = g2 - g0;
    c->gx = (dy2 * dc1 - dy1 * dc2) / determinant;
    c->gy = (dx1 * dc2 - dx2 * dc1) / determinant;
    c->g0 = g0 - c->gx * x0 - c->gy * y0;
    dc1 = b1 - b0;
    dc2 = b2 - b0;
    c->bx = (dy2 * dc1 - dy1 * dc2) / determinant;
    c->by = (dx1 * dc2 - dx2 * dc1) / determinant;
    c->b0 = b0 - c->bx * x0 - c->by * y0;
    return;
} /* end AX_TriangulateGradientColor() */


/**********************/
/* Various Primitives */
/**********************/

/* read ScanModuleAPI.txt */

/*****************/
/* straight line */
/*****************/

/* Simple line demo of the ScanModule API. The valid domain */
/* of x,y for AX_SimpleLineScan() is [0,height) x [0,width) */
void AX_SimpleLineScan
(AX_Float x0, AX_Float y0, AX_Float x1, AX_Float y1, AX_BP *bop)
{
    register AX_Float  x, y, increment;
    register AX_IJ  LP, RP;
    AX_clearbop();
    x = x1 - x0;
    y = y1 - y0;
    if ( ABS(x) >= ABS(y) )  /* scan x */
    {
        increment = y / x;
        if ( x0 <= x1 )
        {
            LP = x0;
            if (x0 != LP)
            {
                AX_addbop (LP, y0);
                LP++;
            }
            y = y0 + (LP + HALF - x0) * increment;
            RP = x1;
            while (LP < RP)
            {
                AX_addbop (LP, y);
                LP++;
                y += increment;
            }
            if (x1 > RP) AX_addbop (RP, y1);
        }
        else
        {
            LP = x1;
            if (x1 != LP)
            {
                AX_addbop (LP, y1);
                LP++;
            }
            y = y1 + (LP + HALF - x1) * increment;
            RP = x0;
            while (LP < RP)
            {
                AX_addbop (LP, y);
                LP++;
                y += increment;
            }
            if (x0 > RP) AX_addbop (RP, y0);
        }
    }
    else  /* scan y */
    {
        increment = x / y;
        if ( y0 <= y1 )
        {
            LP = y0;
            if (y0 != LP)
            {
                AX_addbop (x0, LP);
                LP++;
            }
            x = x0 + (LP + HALF - y0) * increment;
            RP = y1;
            while (LP < RP)
            {
                AX_addbop (x, LP);
                LP++;
                x += increment;
            }
            if (y1 > RP) AX_addbop (x1, RP);
        }
        else
        {
            LP = y1;
            if (y1 != LP)
            {
                AX_addbop (x1, LP);
                LP++;
            }
            x = x1 + (LP + HALF - y1) * increment;
            RP = y0;
            while (LP < RP)
            {
                AX_addbop (x, LP);
                LP++;
                x += increment;
            }
            if (y0 > RP) AX_addbop (x0, RP);
        }
    }
    return;
} /* end AX_SimpleLineScan() */


/* Because line is a 1D instead of 2D object, its scan conversion is  */
/* not perfectly defined. The drawing board is [0,width) x [0,height) */
/* for AX_SimpleLine(), which is also a superset of XDrawLine().      */
void AX_SimpleLine (int iw, AX_Float x0, AX_Float y0,
                    AX_Float x1, AX_Float y1, AX_Carrier c)
{
    register int m;
    /* filter out sick primitives */
    if ( (x0 == x1) && (y0 == y1) ) return;
    /* clip against window -> valid floating-point coordination */
    if (!AX_LineClipWindow(&x0, &y0, &x1, &y1, AX_DEFWIDTH*(1-AX_TINY),
                           AX_DEFHEIGHT*(1-AX_TINY))) return;
    /* scan convert this line to BOP command stack */
    AX_SimpleLineScan (x0, y0, x1, y1, AX_BOP_top(iw));
    /* draw this line to pixmap */
    if (AX_32b)
        for (m=1; m<AX_BOP_TOP(iw,0).offset; m++)
            AX_mem[iw].i4
                [AX_OFF(iw,AX_BOP_TOP(iw,m).p.i,AX_BOP_TOP(iw,m).p.j)] = c;
    else  /* AX_16b */
        for (m=1; m<AX_BOP_TOP(iw,0).offset; m++)
            AX_mem[iw].i2
                [AX_OFF(iw,AX_BOP_TOP(iw,m).p.i,AX_BOP_TOP(iw,m).p.j)] = c;
    return;
}  /* end AX_SimpleLine() */

#ifdef _SimpleLine_TEST
int main (int argc, char *argv[])
{
    int i;
    AX_Float ang;
    AX_openwindow (getpid(),NULL,AX_DEFWIDTH,AX_DEFHEIGHT);
    for (ang=0,i=3; ang<=PI/2; ang+=PI/(4*20),i+=5)
    {
        AX_SimpleLine (0, 3, i, 170*cos(ang)+3, 170*sin(ang)+i,
                       AX_namedcarrier[AX_WHITE]);
        AXDrawLine (0, 203, i, 170*cos(ang)+203, 170*sin(ang)+i);
    }
    AX_dump(0); AX_show(0);
    Press_return();
    return (0);
}
#endif /* _SimpleLine_TEST */


void AX_SimpleLineGradientColor
(int iw, AX_Float x0, AX_Float y0, AX_Carrier c0,
 AX_Float x1, AX_Float y1, AX_Carrier c1)
{
    register int m;
    register AX_Float Dx, Dy, D2, dx, dy, f;
    /* filter out sick primitives */
    if ( (x0 == x1) && (y0 == y1) ) return;
    /* clip against window -> valid floating-point coordination */
    if (!AX_LineClipWindow(&x0, &y0, &x1, &y1, AX_DEFWIDTH*(1-AX_TINY),
                           AX_DEFHEIGHT*(1-AX_TINY))) return;
    /* scan convert this line to BOP command stack */
    AX_SimpleLineScan (x0, y0, x1, y1, AX_BOP_top(iw));
    /* draw this line to pixmap */
    Dx = x1 - x0;
    Dy = y1 - y0;
    D2 = Dx * Dx + Dy * Dy;
    Dx /= D2;
    Dy /= D2;
    if (AX_32b)
        for (m=1; m<AX_BOP_TOP(iw,0).offset; m++)
        {
            dx = AX_BOP_TOP(iw,m).p.i + HALF - x0;
            dy = AX_BOP_TOP(iw,m).p.j + HALF - y0;
            f  = dx * Dx + dy * Dy;
            if (f>=1)
                AX_mem[iw].i4
                    [AX_OFF(iw,AX_BOP_TOP(iw,m).p.i,AX_BOP_TOP(iw,m).p.j)]=c1;
            else if (f<=0)
                AX_mem[iw].i4
                    [AX_OFF(iw,AX_BOP_TOP(iw,m).p.i,AX_BOP_TOP(iw,m).p.j)]=c0;
            else
                AX_mem[iw].i4
                    [AX_OFF(iw,AX_BOP_TOP(iw,m).p.i,AX_BOP_TOP(iw,m).p.j)]=
                    AX_mix(c0, c1, f);
        }
    else  /* AX_16b */
        for (m=1; m<AX_BOP_TOP(iw,0).offset; m++)
        {
            dx = AX_BOP_TOP(iw,m).p.i + HALF - x0;
            dy = AX_BOP_TOP(iw,m).p.j + HALF - y0;
            f  = dx * Dx + dy * Dy;
            if (f>=1)
                AX_mem[iw].i2
                    [AX_OFF(iw,AX_BOP_TOP(iw,m).p.i,AX_BOP_TOP(iw,m).p.j)]=c1;
            else if (f<=0)
                AX_mem[iw].i2
                    [AX_OFF(iw,AX_BOP_TOP(iw,m).p.i,AX_BOP_TOP(iw,m).p.j)]=c0;
            else
                AX_mem[iw].i2
                    [AX_OFF(iw,AX_BOP_TOP(iw,m).p.i,AX_BOP_TOP(iw,m).p.j)]=
                    AX_mix(c0, c1, f);
        }
    return;
}  /* end AX_SimpleLineGradientColor() */

#ifdef _SimpleLineGradientColor_TEST
int main (int argc, char *argv[])
{
    AX_openwindow (getpid(),NULL,AX_DEFWIDTH,AX_DEFHEIGHT);
    AX_SimpleLineGradientColor
        (0,0,0,AX_namedcarrier[AX_BLACK],
         AX_DEFWIDTH,AX_DEFHEIGHT,AX_namedcarrier[AX_WHITE]);
    AX_SimpleLineGradientColor
        (0,0,AX_DEFHEIGHT,AX_namedcarrier[AX_RED],
         AX_DEFWIDTH,0,AX_namedcarrier[AX_GREEN]);
    AX_dump(0); AX_show(0);
    Press_return();
    return (0);
}
#endif /* _SimpleLineGradientColor_TEST */


/* Gupta-Sproull anti-aliased scan conversion of straight line */

#define add_ycenter(ix,iy) ( \
  d = dy * (ix + (HALF - x0)) - dx * (iy + (HALF - y0)), \
  AX_addaop (ix, iy-1, AX_GS_intensity(d+dx)), \
  AX_addaop (ix, iy,   AX_GS_intensity(d)), \
  AX_addaop (ix, iy+1, AX_GS_intensity(d-dx)) )
#define add_ycenter_s(ix,iy,s) ( \
  d = dy * (ix + (HALF - x0)) - dx * (iy + (HALF - y0)), \
  AX_addaop (ix, iy-1, AX_GS_intensity(d+dx)*(s)), \
  AX_addaop (ix, iy,   AX_GS_intensity(d)*(s)), \
  AX_addaop (ix, iy+1, AX_GS_intensity(d-dx)*(s)) )

#define add_xcenter(iy,ix) ( \
  d = dy * (ix + (HALF - x0)) - dx * (iy + (HALF - y0)), \
  AX_addaop (ix-1, iy, AX_GS_intensity(d-dy)), \
  AX_addaop (ix,   iy, AX_GS_intensity(d)), \
  AX_addaop (ix+1, iy, AX_GS_intensity(d+dy)) )
#define add_xcenter_s(iy,ix,s) ( \
  d = dy * (ix + (HALF - x0)) - dx * (iy + (HALF - y0)), \
  AX_addaop (ix-1, iy, AX_GS_intensity(d-dy)*(s)), \
  AX_addaop (ix,   iy, AX_GS_intensity(d)*(s)), \
  AX_addaop (ix+1, iy, AX_GS_intensity(d+dy)*(s)) )

/* Because it involves three pixels, the valid domain of */
/* x,y for AX_GSLineScan() is [1,height-1) x [1,width-1) */
void AX_GSLineScan(AX_Float x0,AX_Float y0,AX_Float x1,AX_Float y1,AX_AP *aop)
{
    register AX_Float  dx, dy, slope, d;
    register AX_IJ  LP, RP, ix, iy;
    AX_clearaop();
    dx = x1 - x0;
    dy = y1 - y0;
    d = sqrt(dx*dx + dy*dy);
    dx /= d;
    dy /= d;
    if ( ABS(dx) >= ABS(dy) )  /* scan x */
    {
        slope = dy / dx;
        if ( x0 <= x1 )
        {
            LP = x0;
            RP = x1;
            if (RP > LP)
            {
                if (x0 != LP)
                {
                    iy = y0;
                    add_ycenter_s (LP, iy, LP+1-x0);
                    LP++;
                }
                while (LP < RP) 
                {
                    iy = y0 + (LP + (HALF - x0)) * slope;
                    add_ycenter (LP, iy);
                    LP++;
                }
                if (x1 > RP)
                {
                    iy = y1;
                    add_ycenter_s (RP, iy, x1-RP);
                }
            }
            else
            {
                iy = (y1 + y0) / 2;
                add_ycenter_s (RP, iy, x1-x0);
            }
        }
        else
        {
            LP = x1;
            RP = x0;
            if (RP > LP)
            {
                if (x1 != LP)
                {
                    iy = y1;
                    add_ycenter_s (LP, iy, LP+1-x1);
                    LP++;
                }
                while (LP < RP) 
                {
                    iy = y0 + (LP + (HALF - x0)) * slope;
                    add_ycenter (LP, iy);
                    LP++;
                }
                if (x0 > RP)
                {
                    iy = y0;
                    add_ycenter_s (RP, iy, x0-RP);
                }
            }
            else
            {
                iy = (y0 + y1) / 2;
                add_ycenter_s (RP, iy, x0-x1);
            }
        }
    }
    else  /* scan y */
    {
        slope = dx / dy;
        if ( y0 <= y1 )
        {
            LP = y0;
            RP = y1;
            if (RP > LP)
            {
                if (y0 != LP)
                {
                    ix = x0;
                    add_xcenter_s (LP, ix, LP+1-y0);
                    LP++;
                }
                while (LP < RP) 
                {
                    ix = x0 + (LP + (HALF - y0)) * slope;
                    add_xcenter (LP, ix);
                    LP++;
                }
                if (y1 > RP)
                {
                    ix = x1;
                    add_xcenter_s (RP, ix, y1-RP);
                }
            }
            else
            {
                ix = (x1 + x0) / 2;
                add_xcenter_s (RP, ix, y1-y0);
            }
        }
        else
        {
            LP = y1;
            RP = y0;
            if (RP > LP)
            {
                if (y1 != LP)
                {
                    ix = x1;
                    add_xcenter_s (LP, ix, LP+1-y1);
                    LP++;
                }
                while (LP < RP) 
                {
                    ix = x0 + (LP + (HALF - y0)) * slope;
                    add_xcenter (LP, ix);
                    LP++;
                }
                if (y0 > RP)
                {
                    ix = x0;
                    add_xcenter_s (RP, ix, y0-RP);
                }
            }
            else
            {
                ix = (x0 + x1) / 2;
                add_xcenter_s (RP, ix, y0-y1);
            }
        }
    }
    return;
} /* end AX_GSLineScan() */


/* Because line is a 1D instead of 2D object, its scan conversion is  */
/* not perfectly defined. The drawing board is [1,width) x [1,height) */
void AX_GSLine (int iw, AX_Float x0, AX_Float y0,
                AX_Float x1, AX_Float y1, AX_Carrier c)
{
    register int m, offset;
    if ( (x0 == x1) && (y0 == y1) ) return;
    if ( ! AX_LineClipRectangle
         ( &x0, &y0, &x1, &y1, 1, 1, (AX_size[iw].width-1)*(1-AX_TINY),
           (AX_size[iw].height-1)*(1-AX_TINY) ) )  return;
    AX_GSLineScan (x0, y0, x1, y1, AX_AOP_top(iw));
    if (AX_32b)
        for (m=1; m<AX_AOP_TOP(iw,0).b.offset; m++)
        {
            offset = AX_OFF(iw,AX_AOP_TOP(iw,m).b.p.i,AX_AOP_TOP(iw,m).b.p.j);
            AX_mem[iw].i4[offset] = AX_MIX
                ( AX_mem[iw].i4[offset], c, AX_AOP_TOP(iw,m).c.area );
        }
    else  /* AX_16b */
        for (m=1; m<AX_AOP_TOP(iw,0).b.offset; m++)
        {
            offset = AX_OFF(iw,AX_AOP_TOP(iw,m).b.p.i,AX_AOP_TOP(iw,m).b.p.j);
            AX_mem[iw].i2[offset] = AX_MIX
                ( AX_mem[iw].i2[offset], c, AX_AOP_TOP(iw,m).c.area );
        }
    return;
}  /* end AX_GSLine() */

#ifdef _GSLine_TEST
int main (int argc, char *argv[])
{
    int i;
    AX_Float ang;
    AX_openwindow (getpid(),NULL,AX_DEFWIDTH,AX_DEFHEIGHT);
    /* AX_namedmop(0,AX_BLUE); */
    /* AXSetLineAttributes(0,2,LineSolid,CapButt,JoinMiter); */
    AX_GSLine (0,1000,1000,-1000,-1000,AX_namedcarrier[AX_RED]);
    AX_GSLine (0,0,200,200,0,AX_namedcarrier[AX_RED]);
    AX_GSLine (0,0,200,1000,200,AX_namedcarrier[AX_RED]);
    for(ang=0,i=3; ang<=PI/2; ang+=PI/(4*20),i+=5)
    {
        AX_GSLine (0, 3, i, 170*cos(ang)+3, 170*sin(ang)+i,
                   AX_namedcarrier[AX_WHITE]);
        AXDrawLine (0, 203, i, 170*cos(ang)+203, 170*sin(ang)+i);
    }
    AX_dump(0); AX_show(0);
    Press_return();
    return(0);
}
#endif /* _GSLine_TEST */


void AX_GSLineGradientColor
(int iw, AX_Float x0, AX_Float y0, AX_Carrier c0,
 AX_Float x1, AX_Float y1, AX_Carrier c1)
{
    register int m, offset;
    register AX_Carrier carrier;
    register AX_Float Dx, Dy, D2, dx, dy, f;
    if ( (x0 == x1) && (y0 == y1) ) return;
    if ( ! AX_LineClipRectangle
         ( &x0, &y0, &x1, &y1, 1, 1, (AX_size[iw].width-1)*(1-AX_TINY),
           (AX_size[iw].height-1)*(1-AX_TINY) ) )  return;
    AX_GSLineScan (x0, y0, x1, y1, AX_AOP_top(iw));
    Dx = x1 - x0;
    Dy = y1 - y0;
    D2 = Dx * Dx + Dy * Dy;
    Dx /= D2;
    Dy /= D2;
    if (AX_32b)
        for (m=1; m<AX_AOP_TOP(iw,0).b.offset; m++)
        {
            dx = AX_AOP_TOP(iw,m).b.p.i + 0.5 - x0;
            dy = AX_AOP_TOP(iw,m).b.p.j + 0.5 - y0;
            f  = dx * Dx + dy * Dy;
            if (f>1) f=1; else if (f<0) f=0;
            carrier = AX_mix (c0, c1, f);
            offset = AX_OFF(iw,AX_AOP_TOP(iw,m).b.p.i,AX_AOP_TOP(iw,m).b.p.j);
            AX_mem[iw].i4[offset] = AX_MIX
                ( AX_mem[iw].i4[offset], carrier, AX_AOP_TOP(iw,m).c.area );
        }
    else  /* AX_16b */
        for (m=1; m<AX_AOP_TOP(iw,0).b.offset; m++)
        {
            dx = AX_AOP_TOP(iw,m).b.p.i + 0.5 - x0;
            dy = AX_AOP_TOP(iw,m).b.p.j + 0.5 - y0;
            f  = dx * Dx + dy * Dy;
            if (f>1) f=1; else if (f<0) f=0;
            carrier = AX_mix (c0, c1, f);
            offset = AX_OFF(iw,AX_AOP_TOP(iw,m).b.p.i,AX_AOP_TOP(iw,m).b.p.j);
            AX_mem[iw].i2[offset] = AX_MIX
                ( AX_mem[iw].i2[offset], carrier, AX_AOP_TOP(iw,m).c.area );
        }
    return;
}  /* end AX_GSLineGradientColor() */

#ifdef _GSLineGradientColor_TEST
int main (int argc, char *argv[])
{
    AX_openwindow (getpid(),NULL,AX_DEFWIDTH,AX_DEFHEIGHT);
    AX_namedmop(0,AX_WHITE);
    AX_GSLineGradientColor
        (0,0,0,AX_namedcarrier[AX_BLACK],
         AX_DEFWIDTH,AX_DEFHEIGHT,AX_namedcarrier[AX_WHITE]);
    AX_GSLineGradientColor
        (0,0,AX_DEFHEIGHT,AX_namedcarrier[AX_RED],
         AX_DEFWIDTH,0,AX_namedcarrier[AX_BLUE]);
    AX_dump(0); AX_show(0);
    Press_return();
    return (0);
}
#endif /* _GSLineGradientColor_TEST */


/************************************************************/
/* Scan convert polygons with EXACT area anti-aliasing:     */
/*                                                          */
/* A singly connected linear shape can often be described   */
/* by a tandem of horizontal cuts 0<=y0<..<y_j<..<=HEIGHT,  */
/* each y_j associated with two feet 0<=l_j<=r_j<=WIDTH. If */
/* y_j is integral, it is called a basecut; otherwise it is */
/* called a nonbasecut. We stipulate that our cut tandem    */
/* cannot own consecutive nonbasecuts that span integer y.  */
/* Once a shaped is clipped and "parsed" as cuts, we can    */
/* compute its area in each pixel exactly: pixels that are  */
/* covered 100% go to "block"-stack "bop", otherwise to an  */
/* "alias"-stack "aop" with exact area.  During cutting     */
/* there is a temporary horizontal line buffer which tags   */
/* pixels that are cut open but cannot yet be confirmed to  */
/* be "block" or not. This buffer gets cleared whenever we  */
/* come to a basecut, or to the ending cut.                 */
/************************************************************/

typedef struct
{
    int LE;
    int IB;
    int LP;
    int RP;
    int IE;
    int RE;
    AX_Float lw;
    AX_Float rw;
} AX_State;

static void update_state (AX_Float y, AX_Float l, AX_Float r, AX_State *s)
{
    AX_LP (l, s->LP);
    AX_RP (r, s->RP);
    /* if ((l<0) || (r>AX_MAXWIDTH)) */
    /* printf ("hopho %f %f %d %d\n", l, r, s->LP, s->RP); */
    s->LE = s->LP - 1;  /* loose bound */
    s->RE = s->RP + 1;
    if ( s->LE < s->RE )
    { /* two different residue pixels */
        s->lw = s->LP - l;
        s->rw = r - s->RE;
    }
    else
    { /* same residue pixel */
        s->lw = 0;
        s->rw = r - l;
    }
    AX_IB (l, s->IB);   /* tight bound */
    AX_IE (r, s->IE);
    return;
} /* end update_state() */

/********************************/
/* dy=0->Dy integration events: */
/*                              */
/* Left foot crosses integer x  */
/* Right foot crosses integer x */
/* Line crosses pending cut y   */
/* Line crosses integer y       */
/********************************/

#define next_intx_event_dy(x,slope,dy) \
  if ((slope) > 0) (dy) = (INT(x)+1-(x)) / (slope); \
  else if ((slope) < 0) (dy) = (((x)==INT(x))?-1.:(INT(x)-(x))) / (slope); \
  else (dy) = AX_INFINITY;
/* event with next pending dy >= AX_INFINITY never means to happen */
#define adv_loser(Ldy,Wdy) if ((Ldy) < AX_INFINITY) (Ldy) -= (Wdy);
#define adv_winner(Wdy,slope) if (slope < 0) Wdy=-1./slope; else Wdy=1./slope;

/* increase area[] by those covered by trapezoidal slice 0 < Dy < 1 */
void AX_sliceConversion ( AX_Float Dy, AX_Float l, AX_Float lslope,
                          AX_Float r, AX_Float rslope, AX_Float *area )
{
    register int i;
    /* state machine maintains pixel beneficiary list */
    AX_State state[2], *s0 = state, *s1 = state+1;
    /* event model as we integrate y from 0 to Dy: between two consecutive */
    /* events there is no change in the list of pixel beneficiaries. */
    AX_Float y=0, dy;
    AX_Float ldy, rdy;
    /* initialize state */
    update_state (y, l, r, s0);
    /* what is the next event: if slope=0 dy would be AX_INFINITY */
    next_intx_event_dy (l, lslope, ldy);
    next_intx_event_dy (r, rslope, rdy);
    do
    {
        FIND3 ( ldy, rdy, Dy-y, >, dy );
        y += dy;
        l += lslope * dy;
        r += rslope * dy;
        update_state (y, l, r, s1);
        if ( s0->IB <= s0->IE )  /* add up the area: trapezoidal rule */
        {
            for (i=s0->LP; i<=s0->RP; i++)
            {
                area [i] += dy / 2;
                /* if ((i<0) || (i>=AX_MAXWIDTH)) */
                /* printf ("fdsaf %d\n", i); */
            }
            area [s0->IB] += s0->lw * dy / 2;
            area [s0->IE] += s0->rw * dy / 2;
            /* if ((s0->IB<0) || (s0->IB>=AX_MAXWIDTH)) */
            /* printf ("s0->IB vfdv %d\n", s0->IB); */
            /* if ((s0->IE<0) || (s0->IE>=AX_MAXWIDTH)) */
            /* printf ("s0->IE vfdv %d\n", s0->IB); */
        }
        if (s1->IB <= s1->IE)
        {
            for (i=s1->LP; i<=s1->RP; i++)
            {
                area [i] += dy / 2;
                /* if ((i<0) || (i>=AX_MAXWIDTH)) */
                /* printf ("vfdv %d\n", i); */
            }
            area [s1->IB] += s1->lw * dy / 2;
            area [s1->IE] += s1->rw * dy / 2;
            /* if ((s1->IB<0) || (s1->IB>=AX_MAXWIDTH)) */
            /* printf ("s1->IB vfdv %d\n", s1->IB); */
            /* if ((s1->IE<0) || (s1->IE>=AX_MAXWIDTH)) */
            /* printf ("s1->IE vfdv %d\n", s1->IB); */
        }
        if ( dy == ldy )
        {
            adv_loser (rdy, dy);
            adv_winner (ldy, lslope);
        }
        else if ( dy == rdy )
        {
            adv_loser (ldy, dy);
            adv_winner (rdy, rslope);
        }
        else break;
        SWITCHSTATE (s0, s1, state);
    } while(1);
    return;
}  /* end AX_sliceConversion() */

#ifdef _sliceConversion_TEST
#define SLICEWIDTH  10
#define SLICEDY     0.5
#define SLICEL      0
#define SLICER      0
#define SLICELL     0
#define SLICERR     10
int main (int argc, char *argv[])
{
    int i;
    AX_Float area[SLICEWIDTH] = {0};
    AX_Float lslope = (SLICELL - SLICEL) / SLICEDY;
    AX_Float rslope = (SLICERR - SLICER) / SLICEDY;
    AX_sliceConversion (SLICEDY, SLICEL, lslope, SLICER, rslope, area);
    for (i=0; i<SLICEWIDTH; i++) printf("%d %f\n", i, area[i]);
    return (0);
}
#endif /* _sliceConversion_TEST */


typedef struct
{
    AX_Float l;
    AX_Float r;
    int IB;
    int LP;
    int RP;
    int IE;
} AX_Section;

/* process a thick trapezoidal slice that lies between two basecuts */
void AX_basecutsConversion
(int j0, int Dj, AX_Float l0, AX_Float lslope, AX_Float r0, AX_Float rslope,
 AX_AP *aop, AX_BP *bop)
{
    register int i, j, IB, LP, RP, IE;
    AX_Float area [AX_MAXWIDTH+2];
    AX_Section section [AX_MAXHEIGHT+2];
    /* carve up the trapezoid */
    for (j=0; j<=Dj; j++)  /* start from 0 to enhance cache hit */
    {
        section[j].l = l0 + lslope * j;
        section[j].r = r0 + rslope * j;
        AX_IB (section[j].l, section[j].IB);
        AX_LP (section[j].l, section[j].LP);
        AX_RP (section[j].r, section[j].RP);
        AX_IE (section[j].r, section[j].IE);
    }
    for (j=0; j<Dj; j++)
    {
        IB = MIN ( section[j].IB, section[j+1].IB );
        LP = MAX ( section[j].LP, section[j+1].LP );
        RP = MIN ( section[j].RP, section[j+1].RP );
        IE = MAX ( section[j].IE, section[j+1].IE );
        if ( LP <= RP )  /* peace zone exists */
        {
            for (i=LP; i<=RP; i++)  AX_addbop(i, j0+j);
            if (IB < LP)
            {
                for (i=IB; i<LP; i++) area[i] = 0;
                AX_sliceConversion (1.0, section[j].l, lslope, LP, 0., area);
                for (i=IB; i<LP; i++) AX_addaop(i, j0+j, area[i]);
            }
            if (RP+1 <= IE)
            {
                for (i=RP+1; i<=IE; i++) area[i] = 0;
                AX_sliceConversion (1.0, RP+1, 0., section[j].r, rslope, area);
                for (i=RP+1; i<=IE; i++) AX_addaop(i, j0+j, area[i]);
            }
        }
        else  /* there is no peace zone */
        {
            for (i=IB; i<=IE; i++) area[i] = 0;
            AX_sliceConversion (1.0, section[j].l, lslope,
                                section[j].r, rslope, area);
            for (i=IB; i<=IE; i++) AX_addaop(i, j0+j, area[i]);
        }
    }
    return;
}  /* end AX_basecutsConversion() */

#ifdef _basecutsConversion_TEST
#define BASEDY     1
#define BASEL      0
#define BASER      1.095445
#define BASELL     0
#define BASERR     0.632456
int main (int argc, char *argv[])
{
    register int m;
    AX_AP aop[1024];
    AX_BP bop[1024];
    AX_Float lslope = (BASELL - BASEL) / BASEDY;
    AX_Float rslope = (BASERR - BASER) / BASEDY;
    /* Conversion class intrinsics does not */
    AX_clearaop();
    AX_clearbop();
    /* initialize header, only accumulates. */
    AX_basecutsConversion(0,BASEDY,BASEL,lslope,BASER,rslope,aop,bop);
    AX_printaop(m);
    AX_printbop(m);
    return (0);
}
#endif /* _basecutsConversion_TEST */


#define clear_area() \
  bzero((void *)(area+t->IB), (t->IE-t->IB+1)*sizeof(AX_Float))
#define collect_area() for (i=t->IB; i<=t->IE; i++) if (area[i] > 0) \
  { if (area[i] > 1-AX_TINY) AX_addbop(i,t->cut[k].j); \
  else AX_addaop(i,t->cut[k].j,area[i]); }
#define add_basecut(k,J,i) ( i++, t->cut[i].y=t->cut[i].j=(J), \
  t->cut[i].l=t0->cut[k].l+(t->cut[i].y-t0->cut[k].y)*t0->cut[k].lslope, \
  t->cut[i].r=t0->cut[k].r+(t->cut[i].y-t0->cut[k].y)*t0->cut[k].rslope, \
  t->cut[i].is_basecut=TRUE, t->cut[i].no_basecut=FALSE, \
  t->cut[i].lslope=t0->cut[k].lslope, t->cut[i].rslope=t0->cut[k].rslope )

/* Scan conversion of shapes parsed into tandem horizontal cuts */
void AX_TandemConversion (AX_Tandem *t0, AX_AP *aop, AX_BP *bop)
{
    register int i,j,k;
    AX_Tandem t[1];
    AX_Float area[AX_MAXWIDTH+2];
    AX_IB (t0->cut[0].l, t->IB);
    AX_IE (t0->cut[0].r, t->IE);
    for (k=0; k<t0->ncut-1; k++)
    {
        AX_IB (t0->cut[k+1].l, i);
        if (i < t->IB) t->IB = i;
        AX_IE (t0->cut[k+1].r, j);
        if (j > t->IE) t->IE = j;
        t0->cut[k].lslope =
            (t0->cut[k+1].l - t0->cut[k].l) / (t0->cut[k+1].y - t0->cut[k].y);
        t0->cut[k].rslope =
            (t0->cut[k+1].r - t0->cut[k].r) / (t0->cut[k+1].y - t0->cut[k].y);
        t0->cut[k].j = INT(t0->cut[k].y);
        t0->cut[k].is_basecut = (t0->cut[k].y == t0->cut[k].j);
        t0->cut[k].no_basecut = !t0->cut[k].is_basecut;
    }
    t0->cut[k].j = INT(t0->cut[k].y);
    t0->cut[k].is_basecut = (t0->cut[k].y == t0->cut[k].j);
    t0->cut[k].no_basecut = !t0->cut[k].is_basecut;
    for (i=k=0; k<t0->ncut-1; i++,k++)
    {
        t->cut[i] = t0->cut[k];
        if ( t0->cut[k].no_basecut )
        {
            if ( t0->cut[k+1].no_basecut )
            {
                if ( t0->cut[k].j+1 < t0->cut[k+1].j )
                {
                    add_basecut(k, t0->cut[k].j+1, i);
                    add_basecut(k, t0->cut[k+1].j, i);
                }
                else if ( t0->cut[k].j+1 == t0->cut[k+1].j )
                    add_basecut(k, t0->cut[k+1].j, i);
            }
            else if ( t0->cut[k].j+1 < t0->cut[k+1].j )
                add_basecut(k, t0->cut[k].j+1, i);
        }
        else if ( t0->cut[k+1].no_basecut )
            if ( t0->cut[k].j < t0->cut[k+1].j )
                add_basecut(k, t0->cut[k+1].j, i);
    }
    t->cut[i] = t0->cut[k];
    t->ncut = i+1;
    if (t->cut[0].no_basecut) clear_area();
    /* for (k=0; k<t->ncut-1; k++) */
    /* if (t->cut[k+1].y-t->cut[k].y < 0) */
    /* printf ("ha %f %f\n", t->cut[k+1].y, t->cut[k].y); */
    for (k=0; k<t->ncut-1; k++)
    {
        if ( t->cut[k].no_basecut )  /* k is not basecut */
        {
            /* if (t->cut[k+1].y-t->cut[k].y < 0) */
            /* printf ("mama %f %f\n", t->cut[k+1].y, t->cut[k].y); */
            AX_sliceConversion ( t->cut[k+1].y-t->cut[k].y,
                                 t->cut[k].l, t->cut[k].lslope,
                                 t->cut[k].r, t->cut[k].rslope, area );
            /* if ((t->IB < -1) || (t->IB > 10000) || */
            /* (t->IE < -1) || (t->IE > 10000)) */
            /* printf ("%d %d\n", t->IB, t->IE); */
            if ( t->cut[k+1].is_basecut ) collect_area();
        }
        else  /* k is basecut */
        {
            if ( t->cut[k+1].is_basecut )
                AX_basecutsConversion ( t->cut[k].j, t->cut[k+1].j-t->cut[k].j,
                                        t->cut[k].l, t->cut[k].lslope,
                                        t->cut[k].r, t->cut[k].rslope,
                                        aop, bop );
            else  /* k+1 is not basecut */
            {
                clear_area();
                AX_sliceConversion ( t->cut[k+1].y - t->cut[k].y, t->cut[k].l,
                        t->cut[k].lslope, t->cut[k].r, t->cut[k].rslope, area);
            }
        }
    }
    if ( t->cut[k].no_basecut ) collect_area();
    return;
}  /* end AX_TandemConversion() */

#ifdef _TandemConversion_TEST
int main (int argc, char *argv[])
{
    int m;
    AX_AP aop[1024];
    AX_BP bop[1024];
    AX_Tandem t[1];
    t->ncut = 3;
    t->cut[0].y = 0;
    t->cut[0].l = 0;
    t->cut[0].r = 1.095445;
    t->cut[1].y = 1;
    t->cut[1].l = 0;
    t->cut[1].r = 0.632456;
    t->cut[2].y = 1.2;
    t->cut[2].l = 0;
    t->cut[2].r = 0;
    /* Conversion class intrinsics does not */
    AX_clearaop();
    AX_clearbop();
    /* initialize header, only accumulates. */
    AX_TandemConversion (t,aop,bop);
    AX_printaop(m);
    AX_printbop(m);
    return (0);
}
#endif /* _TandemConversion_TEST */


/* calculate the leftmost and rightmost intersections of a */
/* horizontal cut with the polygon, if they do intersect.  */
bool AX_PolygonHorizontalCut (AX_Polygon *p, AX_Float y, AX_Cut *c)
{
    int i, j, nABOVE, nBELOW, nCross;
    int ABOVE[AX_PolygonMaxVertex], Touch[AX_PolygonMaxVertex],
        BELOW[AX_PolygonMaxVertex], idx[AX_PolygonMaxVertex];
    AX_Float Cross[AX_PolygonMaxVertex];
    
    nABOVE = 0;
    nBELOW = 0;
    for (i=0; i<p->nVertex; i++)
    {
        ABOVE[i] = (p->Vertex[i].y  > y);
        Touch[i] = (p->Vertex[i].y == y);
        BELOW[i] = (p->Vertex[i].y  < y);
        if (ABOVE[i]) nABOVE++;
        else if (BELOW[i]) nBELOW++;
    }
    if ((nABOVE==p->nVertex)||(nBELOW==p->nVertex)) return(FALSE);
    nCross = 0;
    for (i=0; i<p->nVertex; i++)
        if (Touch[i]) Cross[nCross++] = p->Vertex[i].x;
        else
        {
            j = i + 1;
            if (j == p->nVertex) j = 0;
            if ( (ABOVE[i] && BELOW[j]) || (BELOW[i] && ABOVE[j]) )
                Cross[nCross++] = p->Vertex[i].x + (y - p->Vertex[i].y) *
                    (p->Vertex[j].x - p->Vertex[i].x) /
                    (p->Vertex[j].y - p->Vertex[i].y);
        }
    AX_qsort_glibc (nCross, Cross, idx, AX_USE_NEW_IDX);
    c->y = y;
    c->l = Cross[idx[0]];
    c->r = Cross[idx[nCross-1]];
    return(TRUE);
} /* end AX_PolygonHorizontalCut() */


/* Scan convert polygons with exact anti-aliasing */
void AX_PolygonScan (AX_Polygon *p, AX_AP *aop, AX_BP *bop)
{
    int i, idx [AX_PolygonMaxVertex];
    AX_Float cuty [AX_PolygonMaxVertex], prevy;
    AX_Tandem t[1];
    AX_clearaop();
    AX_clearbop();
    for (i=0; i<p->nVertex; i++) cuty[i] = p->Vertex[i].y;
    AX_qsort_glibc (p->nVertex, cuty, idx, AX_USE_NEW_IDX);
    prevy = cuty[idx[0]] - 1;
    t->ncut = 0;
    for (i=0; i<p->nVertex; i++)
        if (cuty[idx[i]] > prevy)
        {
            AX_PolygonHorizontalCut (p, cuty[idx[i]], t->cut+t->ncut);
            prevy = cuty[idx[i]];
            t->ncut++;
        }
    if (t->ncut > 1)
    {
        /* for (k=0; k<t->ncut-1; k++) */
        AX_TandemConversion (t, aop, bop);
    }
    return;
} /* end AX_PolygonScan() */

#ifdef _PolygonScan_TEST
int main (int argc, char *argv[])
{
    register int m;
    AX_AP aop[1024];
    AX_BP bop[1024];
    AX_Polygon p[1];
    AX_PolygonAssign(p, 5, 0.,0., 1.,1., 1.,2., 2.,2., 2.,0.);
    AX_PolygonScan (p, aop, bop);
    AX_printaop(m);
    AX_printbop(m);
    return (0);
}
#endif /* _PolygonScan_TEST */


/* Draw filled polygon in single color "c" with exact anti-aliasing */
void AX_POLYGON (int iw, AX_Polygon *p, AX_Carrier c)
{
    register int m, offset;
    AX_Polygon result[1];
    AX_SH_PolygonWindowClip (p, (AX_Float)AX_size[iw].width,
                             (AX_Float)AX_size[iw].height, result);
    if (result->nVertex < 3) return;
    AX_PolygonScan (result, AX_AOP_top(iw), AX_BOP_top(iw));
    if (AX_32b)
    {
        for (m=1; m<AX_AOP_TOP(iw,0).b.offset; m++)
        {
            offset = AX_OFF(iw,AX_AOP_TOP(iw,m).b.p.i,AX_AOP_TOP(iw,m).b.p.j);
            AX_mem[iw].i4[offset] = AX_MIX
                ( AX_mem[iw].i4[offset], c, AX_AOP_TOP(iw,m).c.area );
        }
        for (m=1; m<AX_BOP_TOP(iw,0).offset; m++)
            AX_mem[iw].i4
                [AX_OFF(iw,AX_BOP_TOP(iw,m).p.i,AX_BOP_TOP(iw,m).p.j)] = c;
    }
    else  /* AX_16b */
    {
        for (m=1; m<AX_AOP_TOP(iw,0).b.offset; m++)
        {
            offset = AX_OFF(iw,AX_AOP_TOP(iw,m).b.p.i,AX_AOP_TOP(iw,m).b.p.j);
            AX_mem[iw].i2[offset] = AX_MIX
                ( AX_mem[iw].i2[offset], c, AX_AOP_TOP(iw,m).c.area );
        }
        for (m=1; m<AX_BOP_TOP(iw,0).offset; m++)
            AX_mem[iw].i2
                [AX_OFF(iw,AX_BOP_TOP(iw,m).p.i,AX_BOP_TOP(iw,m).p.j)] = c;
    }
    return;
} /* end AX_POLYGON() */

#ifdef _POLYGON_TEST
int main (int argc, char *argv[])
{
    AX_Polygon p[1];
    AX_openwindow (getpid(), NULL, AX_DEFWIDTH, AX_DEFHEIGHT);
    AX_POLYGON (0, AX_PolygonAssign
                (p, 5, 0.,0., 100.,100., 100.,200., 200.,200., 200.,0.),
                AX_namedcarrier[AX_RED]);
    AX_dump(0); AX_show(0);
    Press_return();
    AX_POLYGON (0, AX_PolygonAssign
                (p, 3, 0.,-100., 200.,-100., 100.,200.),
                AX_namedcarrier[AX_GREEN]);
    AX_dump(0); AX_show(0);
    Press_return();
    AX_POLYGON (0, AX_PolygonAssign
                (p, 3, 0.,-100., 200.,-100., 600.,200.),
                AX_namedcarrier[AX_BLUE]);
    AX_dump(0); AX_show(0);
    Press_return();
    AX_POLYGON (0, AX_PolygonAssign
                (p, 3, -10.,300., 300.,100., 150.,400.),
                AX_namedcarrier[AX_WHITE]);
    AX_dump(0); AX_show(0);
    Press_return();
    return (0);
}
#endif /* _POLYGON_TEST */


/* Draw filled triangle in single color "c" with exact anti-aliasing */
void AX_Triangle (int iw, AX_Float x0, AX_Float y0, AX_Float x1,
                  AX_Float y1, AX_Float x2, AX_Float y2, AX_Carrier c)
{
    AX_Polygon p[1];
    AX_PolygonAssign (p, 3, x0,y0, x1,y1, x2,y2);
    AX_POLYGON (iw, p, c);
    return;
} /* end AX_Triangle() */

#ifdef _Triangle_TEST
int main (int argc, char *argv[])
{
    int m;
    AX_openwindow (getpid(), NULL, AX_DEFWIDTH, AX_DEFHEIGHT);
    AX_Triangle(0,0,0,AX_DEFWIDTH,0,0,AX_DEFHEIGHT,
                AX_namedcarrier[AX_RED]);
    AX_dump(0); AX_show(0);
    Press_return();
    AX_Triangle(0,AX_DEFWIDTH,0,AX_DEFWIDTH,AX_DEFHEIGHT,0,AX_DEFHEIGHT,
                AX_namedcarrier[AX_RED]);
    AX_dump(0); AX_show(0);
    m = AX_GET(0, 0, AX_DEFHEIGHT-2); Bwriteint (&m); cr();
    m = AX_GET(0, 0, AX_DEFHEIGHT-1); Bwriteint (&m); cr();
    Press_return();
    return (0);
}
#endif /* _Triangle_TEST */


/* Draw rectangular bar with centerline (x0,y0)-(x1,y1) and thickness */
/* +-"barwidth"/2 in a single color "c" with exact anti-aliasing.     */
void AX_Bar (int iw, AX_Float x0, AX_Float y0, AX_Float x1,
             AX_Float y1, AX_Float barwidth, AX_Carrier c)
{
    register AX_Float dx, dy, d;
    AX_Polygon p [1];
    if ( (x0==x1) && (y0==y1) ) return;
    dx = x1 - x0;
    dy = y1 - y0;
    d = barwidth / 2 / sqrt(dx*dx+dy*dy);
    AX_PolygonAssign (p, 4, x0-d*dy, y0+d*dx, x1-d*dy, y1+d*dx,
                      x1+d*dy, y1-d*dx, x0+d*dy, y0-d*dx);
    AX_POLYGON (iw, p, c);
} /* end AX_Bar() */

#ifdef _Bar_TEST
#define BarWidth 1
int main (int argc, char *argv[])
{
    int i;
    AX_Float ang;
    AX_openwindow (getpid(),NULL,AX_DEFWIDTH,AX_DEFHEIGHT);
    /* AX_namedmop(0,AX_BLUE); */
    AX_GSLine (0,1000,1000,-1000,-1000,AX_namedcarrier[AX_RED]);
    AX_GSLine (0,0,200,200,0,AX_namedcarrier[AX_RED]);
    AX_GSLine (0,0,200,1000,200,AX_namedcarrier[AX_RED]);
    for(ang=0,i=3; ang<=PI/2; ang+=PI/(4*20),i+=5)
    {
        AX_GSLine (0, 3, i, 170*cos(ang)+3, 170*sin(ang)+i,
                   AX_namedcarrier[AX_WHITE]);
        AX_Bar (0, 203, i, 170*cos(ang)+203, 170*sin(ang)+i,
                BarWidth, AX_namedcarrier[AX_WHITE]);
    }
    AX_dump(0); AX_show(0);
    Press_return();
    AX_Bar (0,0,AX_DEFHEIGHT,AX_DEFWIDTH,0,BarWidth,AX_namedcarrier[AX_RED]);
    AX_dump(0); AX_show(0);
    Press_return();
    return (0);
}
#endif /* _Bar_TEST */


/* Draw filled polygon in color field "c" with exact anti-aliasing */
void AX_PolygonGradientColor(int iw, AX_Polygon *p, AX_GradientColor *c)
{
    register int m, offset;
    register AX_Carrier carrier;
    AX_Polygon result[1];
    AX_SH_PolygonWindowClip (p, (AX_Float)AX_size[iw].width,
                             (AX_Float)AX_size[iw].height, result);
    if (result->nVertex < 3) return;
    AX_PolygonScan (result, AX_AOP_top(iw), AX_BOP_top(iw));
    if (AX_32b)
    {
        for (m=1; m<AX_AOP_TOP(iw,0).b.offset; m++)
        {
            carrier = AX_safe_gradientcolor
                (c, AX_AOP_TOP(iw,m).b.p.i+0.5, AX_AOP_TOP(iw,m).b.p.j+0.5);
            offset = AX_OFF(iw,AX_AOP_TOP(iw,m).b.p.i,AX_AOP_TOP(iw,m).b.p.j);
            AX_mem[iw].i4[offset] = AX_MIX
                ( AX_mem[iw].i4[offset], carrier, AX_AOP_TOP(iw,m).c.area );
        }
        for (m=1; m<AX_BOP_TOP(iw,0).offset; m++)
            AX_mem[iw].i4
                [AX_OFF(iw,AX_BOP_TOP(iw,m).p.i,AX_BOP_TOP(iw,m).p.j)] =
                AX_gradientcolor
                (*c, AX_BOP_TOP(iw,m).p.i+0.5, AX_BOP_TOP(iw,m).p.j+0.5);
    }
    else  /* AX_16b */
    {
        for (m=1; m<AX_AOP_TOP(iw,0).b.offset; m++)
        {
            carrier = AX_safe_gradientcolor
                (c, AX_AOP_TOP(iw,m).b.p.i+0.5, AX_AOP_TOP(iw,m).b.p.j+0.5);
            offset = AX_OFF(iw,AX_AOP_TOP(iw,m).b.p.i,AX_AOP_TOP(iw,m).b.p.j);
            AX_mem[iw].i2[offset] = AX_MIX
                ( AX_mem[iw].i2[offset], carrier, AX_AOP_TOP(iw,m).c.area );
        }
        for (m=1; m<AX_BOP_TOP(iw,0).offset; m++)
            AX_mem[iw].i2
                [AX_OFF(iw,AX_BOP_TOP(iw,m).p.i,AX_BOP_TOP(iw,m).p.j)] =
                AX_gradientcolor
                (*c, AX_BOP_TOP(iw,m).p.i+0.5, AX_BOP_TOP(iw,m).p.j+0.5);
    }
    return;
} /* end AX_PolygonGradientColor() */


/* Draw filled color-bar with exact anti-aliasing */
void AX_BarGradientColor
(int iw, AX_Float x0, AX_Float y0, AX_Float r0, AX_Float g0, AX_Float b0,
 AX_Float x1, AX_Float y1, AX_Float r1, AX_Float g1, AX_Float b1,
 AX_Float barwidth)
{
    register AX_Float dx, dy, d;
    AX_Polygon p [1];
    AX_GradientColor c [1];
    if ((x0==x1) && (y0==y1)) return;
    dx = x1 - x0;
    dy = y1 - y0;
    d = barwidth / 2. / sqrt(dx*dx+dy*dy);
    AX_PolygonAssign (p, 4, x0-d*dy, y0+d*dx, x1-d*dy, y1+d*dx,
                      x1+d*dy, y1-d*dx, x0+d*dy, y0-d*dx);
    AX_LinearGradientColor (x0, y0, r0, g0, b0, x1, y1, r1, g1, b1, c);
    AX_PolygonGradientColor (iw, p, c);
    return;
} /* end AX_BarGradientColor() */

#ifdef _BarGradientColor_TEST
#define BarWidth  100
#define BarHeight 400
int main (int argc, char *argv[])
{
    AX_Float wx, wy;
    AX_openwindow (getpid(), NULL, AX_DEFWIDTH, AX_DEFHEIGHT);
    wx = AX_DEFWIDTH  / 2;
    wy = AX_DEFHEIGHT / 2;
    AX_BarGradientColor (0, wx, wy-BarHeight/2, 1, 0.4, 0.3,
                         wx, wy+BarHeight/2, 0.45, 1, 0.9, BarWidth);
    AX_dump(0); AX_show(0);
    Press_return();
    return (0);
}
#endif /* _BarGradientColor_TEST */


/* Draw filled color-triangle with exact anti-aliasing */
void AX_TriangleGradientColor
(int iw, AX_Float x0, AX_Float y0, AX_Float r0, AX_Float g0, AX_Float b0,
 AX_Float x1, AX_Float y1, AX_Float r1, AX_Float g1, AX_Float b1,
 AX_Float x2, AX_Float y2, AX_Float r2, AX_Float g2, AX_Float b2)
{
    AX_Polygon p [1];
    AX_GradientColor c [1];
    AX_PolygonAssign (p, 3, x0,y0, x1,y1, x2,y2);
    AX_TriangulateGradientColor
        (x0, y0, r0, g0, b0, x1, y1, r1, g1, b1, x2, y2, r2, g2, b2, c);
    AX_PolygonGradientColor (iw, p, c);
    return;
} /* end AX_TriangleGradientColor() */

#ifdef _TriangleGradientColor_TEST
int main (int argc, char *argv[])
{
    AX_Float side, high, wx, wy;
    AX_openwindow (getpid(), NULL, AX_DEFWIDTH, AX_DEFHEIGHT);
    side = MIN(AX_DEFWIDTH, AX_DEFHEIGHT) * 0.9;
    high = side * sqrt(3.) / 2;
    wx = AX_DEFWIDTH  / 2;
    wy = AX_DEFHEIGHT / 2;
    while(1)
    {
        AX_bzero(0);
        AX_TriangleGradientColor (0, wx-side/2, wy-high/2, 1, 0, 0,
                                  wx+side/2, wy-high/2, 0, 1, 0,
                                  wx,        wy+high/2, 0, 0, 1);
        AX_dump(0); AX_show(0);
        Press_return();
        AX_bzero(0);
        AX_TriangleGradientColor (0, wx-side/2, wy-high/2, 1, 0, 0,
                                  wx+side/2, wy-high/2, 1, 1, 0,
                                  wx,        wy+high/2, 0, 0, 1);
        AX_dump(0); AX_show(0);
        Press_return();
        AX_bzero(0);
        AX_TriangleGradientColor (0, wx-side/2, wy-high/2, 1, 1, 0,
                                  wx+side/2, wy-high/2, 0, 1, 1,
                                  wx,        wy+high/2, 1, 0, 1);
        AX_dump(0); AX_show(0);
        Press_return();
        AX_bzero(0);
        AX_TriangleGradientColor (0, wx-side/2, wy-high/2, 1, 0, 0,
                                  wx+side/2, wy-high/2, 0, 0, 0,
                                  wx,        wy+high/2, 1, 1, 1);
        AX_dump(0); AX_show(0);
        Press_return();
    }
    return (0);
}
#endif /* _TriangleGradientColor_TEST */


/* Perform brute-force integration: using VerticalSegment is better than */
/* using HorizontalSegment near the circle cap where dy / dx is small.   */
#define MESH  16
#define AX_PixelAreaInCircle(i,j,x0,y0,radius2) (( \
  AX_VerticalSegmentInCircle((i)+0.5/MESH,  (j),(j)+1, x0,y0,radius2) + \
  AX_VerticalSegmentInCircle((i)+1.5/MESH,  (j),(j)+1, x0,y0,radius2) + \
  AX_VerticalSegmentInCircle((i)+2.5/MESH,  (j),(j)+1, x0,y0,radius2) + \
  AX_VerticalSegmentInCircle((i)+3.5/MESH,  (j),(j)+1, x0,y0,radius2) + \
  AX_VerticalSegmentInCircle((i)+4.5/MESH,  (j),(j)+1, x0,y0,radius2) + \
  AX_VerticalSegmentInCircle((i)+5.5/MESH,  (j),(j)+1, x0,y0,radius2) + \
  AX_VerticalSegmentInCircle((i)+6.5/MESH,  (j),(j)+1, x0,y0,radius2) + \
  AX_VerticalSegmentInCircle((i)+7.5/MESH,  (j),(j)+1, x0,y0,radius2) + \
  AX_VerticalSegmentInCircle((i)+8.5/MESH,  (j),(j)+1, x0,y0,radius2) + \
  AX_VerticalSegmentInCircle((i)+9.5/MESH,  (j),(j)+1, x0,y0,radius2) + \
  AX_VerticalSegmentInCircle((i)+10.5/MESH, (j),(j)+1, x0,y0,radius2) + \
  AX_VerticalSegmentInCircle((i)+11.5/MESH, (j),(j)+1, x0,y0,radius2) + \
  AX_VerticalSegmentInCircle((i)+12.5/MESH, (j),(j)+1, x0,y0,radius2) + \
  AX_VerticalSegmentInCircle((i)+13.5/MESH, (j),(j)+1, x0,y0,radius2) + \
  AX_VerticalSegmentInCircle((i)+14.5/MESH, (j),(j)+1, x0,y0,radius2) + \
  AX_VerticalSegmentInCircle((i)+15.5/MESH, (j),(j)+1, x0,y0,radius2) ) / MESH)

#define XMIN  0
#define XMAX  width
#define YMIN  0
#define YMAX  height

#define add_cut(Y) { t->cut[i].y=Y; t->cut[i].r=sqrt(radius-SQUARE((Y)-y0)); \
  t->cut[i].l = x0 - t->cut[i].r; t->cut[i].r += x0; \
  if (t->cut[i].l < XMIN) t->cut[i].l = XMIN; \
  if (t->cut[i].r > XMAX) t->cut[i].r = XMAX; i++; }

/***********************************************************/
/* Scan convert circle clipped by window with almost exact */
/* anti-aliasing. The circle must have finite area in the  */
/* window (use AX_CircleIntersectWindow() beforehand).     */
/***********************************************************/
void AX_CircleWindowScan
(AX_Float x0, AX_Float y0, AX_Float radius, AX_Float width, AX_Float height,
 AX_AP *aop, AX_BP *bop)
{
    register int i, j, LP, RP;
    register AX_Float ybtm, ytop, ymid;
    AX_Float min, max;
    AX_Tandem t[1];
    ytop = y0 + radius;
    ybtm = y0 - radius;
    AX_clearaop();
    AX_clearbop();
    radius *= radius;
    if (ytop > YMAX)
    {
        min = XMIN; max = XMAX;
        if (! AX_HorizontalSegmentClipCircle (YMAX, &min,&max, x0,y0,radius) )
            if ( max == XMAX )
            {
                min = YMIN; max = YMAX;
                AX_VerticalSegmentCutCircle (XMAX, &min,&max, x0,y0,radius);
                ytop = max;
            }
            else
            {
                min = YMIN; max = YMAX;
                AX_VerticalSegmentCutCircle (XMIN, &min,&max, x0,y0,radius);
                ytop = max;
            }
        else ytop = YMAX;
    }
    else if (x0 > XMAX)
    {
        min = YMIN; max = YMAX;
        AX_VerticalSegmentCutCircle (XMAX, &min,&max, x0,y0,radius);
        ytop = max;
    }
    else if (x0 < XMIN)
    {
        min = YMIN; max = YMAX;
        AX_VerticalSegmentCutCircle (XMIN, &min,&max, x0,y0,radius);
        ytop = max;
    }
    if (ybtm < YMIN)
    {
        min = XMIN; max = XMAX;
        if (! AX_HorizontalSegmentClipCircle (YMIN, &min,&max, x0,y0,radius) )
            if ( max == XMAX )
            {
                min = YMIN; max = YMAX;
                AX_VerticalSegmentCutCircle (XMAX, &min,&max, x0,y0,radius);
                ybtm = min;
            }
            else
            {
                min = YMIN; max = YMAX;
                AX_VerticalSegmentCutCircle (XMIN, &min,&max, x0,y0,radius);
                ybtm = min;
            }
        else ybtm = YMIN;
    }
    else if (x0 > XMAX)
    {
        min = YMIN; max = YMAX;
        AX_VerticalSegmentCutCircle (XMAX, &min,&max, x0,y0,radius);
        ybtm = min;
    }
    else if (x0 < XMIN)
    {
        min = YMIN; max = YMAX;
        AX_VerticalSegmentCutCircle (XMIN, &min,&max, x0,y0,radius);
        ybtm = min;
    }
    /* prevent NaN exception due to floating-point error */
    radius *= (1 + 2 * AX_TINY);
    /* determine block pixels and proscribe alias pixels */
    LP = ybtm;
    i = 0;
    if (ybtm != LP)
    {
        add_cut (ybtm);
        LP++;
    }
    if ( radius < 1 )
    { /* add cut at middle to guarantee at least one cut at waist */
        ymid = (ybtm + ytop) / 2;
        RP = ymid;
        for (j=LP; j<=RP; j++) add_cut(j);
        if (ymid != RP) add_cut(ymid);
        LP = RP + 1;
    }
    RP = ytop;
    for (j=LP; j<=RP; j++) add_cut(j);
    if (ytop != RP) add_cut(ytop);
    t->ncut = i;
    AX_TandemConversion (t, aop, bop);
    /* refine alias pixel areas close to the bottom and top */
    for (i=1; i<aop[0].b.offset; i++)
        if (aop[i].b.p.j <= ybtm + (AX_CAPHEIGHT-0.5))
            aop[i].c.area = AX_PixelAreaInCircle
                (aop[i].b.p.i, aop[i].b.p.j, x0,y0,radius);
        else break;
    for (j=aop[0].b.offset-1; j>=i; j--)
        if (aop[j].b.p.j >= ytop - (AX_CAPHEIGHT+0.5))
            aop[j].c.area = AX_PixelAreaInCircle
                (aop[j].b.p.i, aop[j].b.p.j, x0,y0,radius);
        else break;
    return;
} /* end AX_CircleWindowScan() */
#undef XMIN
#undef XMAX
#undef YMIN
#undef YMAX
#undef add_cut

#ifdef _CircleWindowScan_TEST
#define X0     0
#define Y0     0.1
#define RADIUS 1.1
#define WIDTH  10
#define HEIGHT 10
int main (int argc, char *argv[])
{
    register int m;
    AX_AP aop[WIDTH*HEIGHT+1];
    AX_BP bop[WIDTH*HEIGHT+1];
    AX_CircleWindowScan (X0, Y0, RADIUS, WIDTH, HEIGHT, aop, bop);
    AX_printaop(m);
    AX_printbop(m);
    return (0);
}
#endif /* _CircleWindowScan_TEST */


/* Draw filled circle in single color "c" */
void AX_Circle
(int iw, AX_Float x0, AX_Float y0, AX_Float radius, AX_Carrier c)
{
    register int m, offset;
    if (radius <= 0) return;
    if (!AX_CircleIntersectWindow
        (x0,y0,radius,AX_size[iw].width,AX_size[iw].height)) return;
    AX_CircleWindowScan (x0, y0, radius,
                         AX_size[iw].width, AX_size[iw].height,
                         AX_AOP_top(iw), AX_BOP_top(iw));
    if (AX_32b)
    {
        for (m=1; m<AX_AOP_TOP(iw,0).b.offset; m++)
        {
            offset = AX_OFF(iw,AX_AOP_TOP(iw,m).b.p.i,AX_AOP_TOP(iw,m).b.p.j);
            AX_mem[iw].i4[offset] = AX_MIX
                ( AX_mem[iw].i4[offset], c, AX_AOP_TOP(iw,m).c.area );
        }
        for (m=1; m<AX_BOP_TOP(iw,0).offset; m++)
            AX_mem[iw].i4
                [AX_OFF(iw,AX_BOP_TOP(iw,m).p.i,AX_BOP_TOP(iw,m).p.j)] = c;
    }
    else  /* AX_16b */
    {
        for (m=1; m<AX_AOP_TOP(iw,0).b.offset; m++)
        {
            offset = AX_OFF(iw,AX_AOP_TOP(iw,m).b.p.i,AX_AOP_TOP(iw,m).b.p.j);
            AX_mem[iw].i2[offset] = AX_MIX
                ( AX_mem[iw].i2[offset], c, AX_AOP_TOP(iw,m).c.area );
        }
        for (m=1; m<AX_BOP_TOP(iw,0).offset; m++)
            AX_mem[iw].i2
                [AX_OFF(iw,AX_BOP_TOP(iw,m).p.i,AX_BOP_TOP(iw,m).p.j)] = c;
    }
    return;
} /* end AX_Circle() */

#ifdef _Circle_TEST
#define r1  100
#define r2  80
#define r3  2.5
int main (int argc, char *argv[])
{
    AX_Float wx,wy;
    wx = AX_DEFWIDTH/4;
    wy = AX_DEFHEIGHT/4;
    AX_openwindow (getpid(), NULL, AX_DEFWIDTH, AX_DEFHEIGHT);
    /* AX_mop(0,AX_namedcarrier[AX_WHITE]); */
    AXSetNamedForeground(0,AX_RED);
    /* AX_Circle(0,0,wy*5,wy*4,AX_namedcarrier[AX_RED]); */
    AXFillArc(0,wx-r1,wy-r1,2*r1,2*r1,0,360*64);
    AX_Circle(0,wx*3,wy*3,r1,AX_namedcarrier[AX_RED]);
    AXFillArc(0,wx-r2,wy*3-r2,2*r2,2*r2,0,360*64);
    AX_Circle(0,wx*3,wy,r2,AX_namedcarrier[AX_RED]);
    AXFillArc(0,2*wx-10-r3,2*wy-r3,2*r3,2*r3,0,360*64);
    AX_Circle(0,2*wx+10,2*wy,r3,AX_namedcarrier[AX_RED]);
    AX_Circle(0,3*wx,3*wy,20,AX_namedcarrier[AX_BLUE]);
    AX_dump(0); AX_show(0);
    Press_return();
    return (0);
}
#endif /* _Circle_TEST */


#define AX_PixelAreaInEllipse(i,j,e) ( ( \
  AX_VerticalSegmentInEllipse((i)+0.5/MESH,  (j),(j)+1, e) + \
  AX_VerticalSegmentInEllipse((i)+1.5/MESH,  (j),(j)+1, e) + \
  AX_VerticalSegmentInEllipse((i)+2.5/MESH,  (j),(j)+1, e) + \
  AX_VerticalSegmentInEllipse((i)+3.5/MESH,  (j),(j)+1, e) + \
  AX_VerticalSegmentInEllipse((i)+4.5/MESH,  (j),(j)+1, e) + \
  AX_VerticalSegmentInEllipse((i)+5.5/MESH,  (j),(j)+1, e) + \
  AX_VerticalSegmentInEllipse((i)+6.5/MESH,  (j),(j)+1, e) + \
  AX_VerticalSegmentInEllipse((i)+7.5/MESH,  (j),(j)+1, e) + \
  AX_VerticalSegmentInEllipse((i)+8.5/MESH,  (j),(j)+1, e) + \
  AX_VerticalSegmentInEllipse((i)+9.5/MESH,  (j),(j)+1, e) + \
  AX_VerticalSegmentInEllipse((i)+10.5/MESH, (j),(j)+1, e) + \
  AX_VerticalSegmentInEllipse((i)+11.5/MESH, (j),(j)+1, e) + \
  AX_VerticalSegmentInEllipse((i)+12.5/MESH, (j),(j)+1, e) + \
  AX_VerticalSegmentInEllipse((i)+13.5/MESH, (j),(j)+1, e) + \
  AX_VerticalSegmentInEllipse((i)+14.5/MESH, (j),(j)+1, e) + \
  AX_VerticalSegmentInEllipse((i)+15.5/MESH, (j),(j)+1, e) ) / MESH )

#define XMIN  0
#define XMAX  width
#define YMIN  0
#define YMAX  height
#define add_cut(Y) { t->cut[i].y=(Y); t->cut[i].l=XMIN; t->cut[i].r=XMAX; \
  AX_HorizontalSegmentCutEllipse(Y, &t->cut[i].l, &t->cut[i].r, e); i++; }

/************************************************************/
/* Scan convert ellipse clipped by window with almost exact */
/* anti-aliasing. The ellipse must have finite area in the  */
/* window (use AX_EllipseIntersectWindow() beforehand).     */
/************************************************************/
void AX_EllipseWindowScan
(AX_Ellipse *e, AX_Float width, AX_Float height, AX_AP *aop, AX_BP *bop)
{
    register int i, j, LP, RP;
    register AX_Float ybtm, ytop, ymid;
    AX_Float min, max;
    AX_Tandem t[1];
    ytop = e->ymax;
    ybtm = e->ymin;
    AX_clearaop();
    AX_clearbop();
    if (ytop > YMAX)
    {
        min = XMIN; max = XMAX;
        if (! AX_HorizontalSegmentClipEllipse (YMAX, &min,&max, e))
            if ( max == XMAX )
            {
                min = YMIN; max = YMAX;
                AX_VerticalSegmentCutEllipse (XMAX, &min,&max, e);
                ytop = max;
            }
            else
            {
                min = YMIN; max = YMAX;
                AX_VerticalSegmentCutEllipse (XMIN, &min,&max, e);
                ytop = max;
            }
        else ytop = YMAX;
    }
    else if (e->ymax_x > XMAX)
    {
        min = YMIN; max = YMAX;
        AX_VerticalSegmentCutEllipse (XMAX, &min,&max, e);
        ytop = max;
    }
    else if (e->ymax_x < XMIN)
    {
        min = YMIN; max = YMAX;
        AX_VerticalSegmentCutEllipse (XMIN, &min,&max, e);
        ytop = max;
    }
    if (ybtm < YMIN)
    {
        min = XMIN; max = XMAX;
        if (! AX_HorizontalSegmentClipEllipse (YMIN, &min,&max, e) )
            if ( max == XMAX )
            {
                min = YMIN; max = YMAX;
                AX_VerticalSegmentCutEllipse (XMAX, &min,&max, e);
                ybtm = min;
            }
            else
            {
                min = YMIN; max = YMAX;
                AX_VerticalSegmentCutEllipse (XMIN, &min,&max, e);
                ybtm = min;
            }
        else ybtm = YMIN;
    }
    else if (e->ymin_x > XMAX)
    {
        min = YMIN; max = YMAX;
        AX_VerticalSegmentCutEllipse (XMAX, &min,&max, e);
        ybtm = min;
    }
    else if (e->ymin_x < XMIN)
    {
        min = YMIN; max = YMAX;
        AX_VerticalSegmentCutEllipse (XMIN, &min,&max, e);
        ybtm = min;
    }
    /* determine block pixels and proscribe alias pixels */
    LP = ybtm;
    i = 0;
    if (ybtm != LP)
    {
        add_cut (ybtm);
        LP++;
    }
    if ( e->ymax - e->ymin < 2 )
    { /* add cut at middle to guarantee at least one cut at waist */
        ymid = (ybtm + ytop) / 2;
        RP = ymid;
        for (j=LP; j<=RP; j++) add_cut(j);
        if (ymid != RP) add_cut(ymid);
        LP = RP + 1;
    }
    RP = ytop;
    for (j=LP; j<=RP; j++) add_cut(j);
    if (ytop != RP) add_cut(ytop);
    t->ncut = i;
    AX_TandemConversion (t, aop, bop);
    /* refine alias pixels close to the bottom and top */
    for (i=1; i<aop[0].b.offset; i++)
        if (aop[i].b.p.j <= ybtm + (AX_CAPHEIGHT-0.5))
            aop[i].c.area = AX_PixelAreaInEllipse(aop[i].b.p.i,aop[i].b.p.j,e);
        else break;
    for (j=aop[0].b.offset-1; j>=i; j--)
        if (aop[j].b.p.j >= ytop - (AX_CAPHEIGHT+0.5))
            aop[j].c.area = AX_PixelAreaInEllipse(aop[j].b.p.i,aop[j].b.p.j,e);
        else break;
    return;
} /* end AX_EllipseWindowScan() */
#undef XMIN
#undef XMAX
#undef YMIN
#undef YMAX
#undef add_cut

#ifdef _EllipseWindowScan_TEST
#define WIDTH  10
#define HEIGHT 10
int main (int argc, char *argv[])
{
    register int m;
    AX_Ellipse e[1];
    AX_Float x0=0, y0=0, angle=90, a=2, b=1;
    AX_AP aop[WIDTH*HEIGHT+1];
    AX_BP bop[WIDTH*HEIGHT+1];
    AX_EllipseWindowScan
        (AX_ELLIPSEASSIGN(x0,y0,angle,a,b,e),WIDTH,HEIGHT,aop,bop);
    AX_printaop(m);
    AX_printbop(m);
    return (0);
}
#endif /* _EllipseWindowScan_TEST */


/* Draw filled ellipse in single color "c" */
void AX_ELLIPSE (int iw, AX_Ellipse *e, AX_Carrier c)
{
    register int m, offset;
    if (!AX_EllipseIntersectWindow(e,AX_size[iw].width,AX_size[iw].height))
        return;
    AX_EllipseWindowScan (e, AX_size[iw].width, AX_size[iw].height,
                          AX_AOP_top(iw), AX_BOP_top(iw));
    if (AX_32b)
    {
        for (m=1; m<AX_AOP_TOP(iw,0).b.offset; m++)
        {
            offset = AX_OFF(iw,AX_AOP_TOP(iw,m).b.p.i,AX_AOP_TOP(iw,m).b.p.j);
            AX_mem[iw].i4[offset] = AX_MIX
                ( AX_mem[iw].i4[offset], c, AX_AOP_TOP(iw,m).c.area );
        }
        for (m=1; m<AX_BOP_TOP(iw,0).offset; m++)
            AX_mem[iw].i4
                [AX_OFF(iw,AX_BOP_TOP(iw,m).p.i,AX_BOP_TOP(iw,m).p.j)] = c;
    }
    else  /* AX_16b */
    {
        for (m=1; m<AX_AOP_TOP(iw,0).b.offset; m++)
        {
            offset = AX_OFF(iw,AX_AOP_TOP(iw,m).b.p.i,AX_AOP_TOP(iw,m).b.p.j);
            AX_mem[iw].i2[offset] = AX_MIX
                ( AX_mem[iw].i2[offset], c, AX_AOP_TOP(iw,m).c.area );
        }
        for (m=1; m<AX_BOP_TOP(iw,0).offset; m++)
            AX_mem[iw].i2
                [AX_OFF(iw,AX_BOP_TOP(iw,m).p.i,AX_BOP_TOP(iw,m).p.j)] = c;
    }
    return;
} /* end AX_ELLIPSE() */

#ifdef _ELLIPSE_TEST
int main (int argc, char *argv[])
{
    AX_Float wx, wy, i;
    AX_Ellipse e[1];
    AX_Float x0, y0, angle, a, b;
    wx = AX_DEFWIDTH/4;
    wy = AX_DEFHEIGHT/4;
    AX_openwindow (getpid(), NULL, AX_DEFWIDTH, AX_DEFHEIGHT);
    for (angle=0,x0=3; x0<AX_DEFWIDTH; angle+=5,x0+=15)
        AX_ELLIPSE (0, AX_ELLIPSEASSIGN (x0, 20, angle, 4, 2, e),
                    AX_namedcarrier[AX_RED]);
    for (angle=90,x0=0; x0<AX_DEFWIDTH; angle+=10,x0+=30)
        AX_ELLIPSE (0, AX_ELLIPSEASSIGN (x0, 40, angle, 10, 6, e),
                    AX_namedcarrier[AX_GREEN]);
    AX_ELLIPSE (0, AX_ELLIPSEASSIGN (2*wx, 2*wy, 30, 100, 70, e),
                AX_namedcarrier[AX_RED]);
    AX_ELLIPSE (0, AX_ELLIPSEASSIGN (wx-10, 3*wy-5, 90, 100, 70, e),
                AX_namedcarrier[AX_RED]);
    AX_ELLIPSE (0, AX_ELLIPSEASSIGN (2*wx, 5*wy, 120, 300, 150, e),
                AX_namedcarrier[AX_RED]);
    AX_ELLIPSE (0, AX_ELLIPSEASSIGN (3.3*wx, 5*wy, 120, 300, 150, e),
                AX_namedcarrier[AX_BLUE]);
    AX_dump(0); AX_show(0);
    Press_return();
    return (0);
}
#endif /* _ELLIPSE_TEST */
