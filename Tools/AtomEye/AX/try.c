/*********/
/* try.c */
/*********/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <sys/types.h>


typedef struct
{
    Pixmap new_pixmap;
    Drawable old_drawable;
    GC old_gc;
    AXSize old_size;
} CaptureSwitch;


static void CaptureResize (int iw, int new_resolution, CaptureSwitch *cs)
{
    double scale_x, scale_y;
    cs->old_size = AX_size[iw];
    scale_x = ((double) AX_MAXWIDTH)  / AX_size[iw].width;
    scale_y = ((double) AX_MAXHEIGHT) / AX_size[iw].height;
    if ( scale_x > scale_y )
    {
        AX_size[iw].width = scale_y * AX_size[iw].width;
        AX_size[iw].height = AX_MAXHEIGHT;
    }
    else /* scale_x <= scale_y */
    {
        AX_size[iw].width = AX_MAXWIDTH;
        AX_size[iw].height = scale_x * AX_size[iw].height;
    }
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
