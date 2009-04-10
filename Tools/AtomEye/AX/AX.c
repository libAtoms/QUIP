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


/********************/
/* Global variables */
/********************/

/** Shared Public Resources **/

int AX_videomode;

/* "TrueColor" or "DirectColor" visuals; read Notes/color.txt */
bool AX_noneedcolormap;
/* defined only when AX_noneedcolormap */
AX_Carrier AX_rmask, AX_gmask, AX_bmask;
/* how many bytes of the pixmap memory belongs to a pixel */
int AX_bytes, AX_8b, AX_16b, AX_32b;
/* should all fits into one cacheline */


/** Individual Public Resources **/
/** accessible via a type int key often called "iw" **/

/* display is a connection, so can only be used by a single thread at a time */
Display *AX_display [AX_MAXWIN];

AXSize AX_size [AX_MAXWIN];

Window AX_win [AX_MAXWIN];

Pixmap AX_pixmap [AX_MAXWIN];

Drawable AX_drawable [AX_MAXWIN];

/* if the server is 8-bit, you will probably need this a lot. */
Colormap AX_colormap [AX_MAXWIN];/* colormaps from the same server / screen should be identical, */
/* unfortunately the initial process may die early and free it. */

/* shared memory: direct pixmap access like in mmap() */
AX_Shm AX_mem [AX_MAXWIN];

XImage *AX_img [AX_MAXWIN];

/* graphics context built on shared memory pixmap */
GC AX_gc [AX_MAXWIN];

/* Default interfaces via which AX routines communicate */
/* XEvent interface */
XEvent AX_event [AX_MAXWIN];
/* XGetGeometry interface */
unsigned int AX_depth [AX_MAXWIN];
int AX_borderwidth [AX_MAXWIN];
XSizeHints AX_hint [AX_MAXWIN];
/* XQueryPointer interface */
Window AX_root [AX_MAXWIN], AX_child [AX_MAXWIN];
/* mouse coordinates with respect to root */
int AX_root_x [AX_MAXWIN], AX_root_y [AX_MAXWIN];
/* mouse coordinates with respect to this window */
int AX_win_x [AX_MAXWIN], AX_win_y [AX_MAXWIN];
/* modifier keys and pointer buttons */
unsigned int AX_pointermask [AX_MAXWIN];
/* XpmCreatePixmapFromData and XSetWMHints interface */
Pixmap AX_icon_pixmap[AX_MAXWIN];
XWMHints* AX_win_hints[AX_MAXWIN];
char AX_title [AX_MAXWIN] [AX_MAXSTRSIZE];

/** Individual Private Resources **/

/******************************************************************/
/* Each AX caller, be it a forked process or thread, must supply  */
/* a unique nonzero number to identify itself. If the same caller */
/* calls twice (asking for two windows to draw non-concurrently), */
/* the X-server may save some resources like display connection.  */
/******************************************************************/
int AX_cid [AX_MAXWIN] = {0};
#define invalid(iw) (((iw)<0)||((iw)>=AX_MAXWIN))
#define vacant(iw)  (AX_cid[iw]==0)
#define nonvacant(iw)  (AX_cid[iw]!=0)
int AX_screen [AX_MAXWIN];
static Visual *AX_visual [AX_MAXWIN];
static XShmSegmentInfo AX_shminfo [AX_MAXWIN];
/* pop-up window toggle */
bool AX_WINDOWS_OVERRIDE_REDIRECT = False;

/*********************/
/* Window Management */
/*********************/

/* Check if this WinID is in use; it is supposed to */
void AX_AssertWinID (char *who, int iw)
{
    if ( invalid(iw) )
        pe ("%s: you are crazy, invalid AX_WinID %d\n", who, iw);
    if ( vacant(iw) )
        pe ("%s: you are crazy, vacant AX_WinID %d\n", who, iw);
    return;
} /* end AX_AssertWinID() */

char *AX_Module_Name [MODULE_MAX] =
{"Shared Memory GC", "Scan Conversion", "3D", "Ball Cache"};

int AX_Module_Is_PluggedIn [AX_MAXWIN] [MODULE_MAX] = {{0}};

static struct
{
    int (*plugin)  (int iw);
    int (*plugout) (int iw);
    int (*resize)  (int iw);
} AX_modules [MODULE_MAX] =
{{&AX_plugin_ShmGC_module, &AX_plugout_ShmGC_module, &AX_resize_ShmGC_module},
 {&AX_plugin_Scan_module, &AX_plugout_Scan_module, &AX_resize_Scan_module},
 {&AX_plugin_3D_module, &AX_plugout_3D_module, &AX_resize_3D_module},
 {&AX_plugin_BC_module, &AX_plugout_BC_module, &AX_resize_BC_module}};


/* Check if module dependence is satisfied before plugging in a new one */
int AX_Check_Module_Dependence (char *who, int iw, int module_idx)
{
    register int i;
    for (i=0; i<module_idx; i++)
        if (! AX_Module_Is_PluggedIn[iw][i])
            pe ("%s: in order to plug in %s Module\n"
                "you have to plug in %s Module first.\n",
                who, AX_Module_Name[module_idx], AX_Module_Name[i]);
    return(1);
} /* end AX_Check_Module_Dependence() */


/* Create shared memory graphics context that is accessible */
/* both via AX_mem[] and standard Xlib GC pixmap draws. The */
/* size of AX_mem[] created is dependent on AX_size[iw].    */
int AX_plugin_ShmGC_module (int iw)
{
    AX_AssertWinID ("AX_plugin_ShmGC_module", iw);
    AX_Check_Module_Dependence ("AX_plugin_ShmGC_module", iw, MODULE_ShmGC);
    if ( AX_Module_Is_PluggedIn [iw] [MODULE_ShmGC] )
        pe ("AX_plugin_ShmGC_module:\n"
            "you are crazy, this module is already plugged in\n");
    if (AX_videomode)
    {
#ifdef _AX_USE_SHM
        AX_img[iw] = XShmCreateImage
            ( AX_display[iw],
              AX_visual[iw],
              AX_depth[iw],
              ZPixmap,
              NULL,
              &AX_shminfo[iw],
              AX_size[iw].width,
              AX_size[iw].height );
#endif
    }
    else /* Shared Memory Video not supported */
    {
        switch (AX_depth[iw])
        {
            case 8:  AX_bytes = 1; break;
            case 16: AX_bytes = 2; break;
            case 24: case 32: AX_bytes = 4; break;
            default: AX_bytes = 4;
        }
        /* non-shared client data-memory */
        MALLOC ( AX_plugin_ShmGC_module, AX_mem[iw].uc,
                 AX_size[iw].width * AX_size[iw].height * AX_bytes,
                 unsigned char );
        AX_img[iw] = XCreateImage
            ( AX_display[iw],
              AX_visual[iw],
              AX_depth[iw],
              ZPixmap,
              0,
              (char *) AX_mem[iw].uc,
              AX_size[iw].width,
              AX_size[iw].height,
              ( (AX_size[iw].width * AX_bytes) % 4 == 0 ) ? 32 :
              ( ( (AX_size[iw].width * AX_bytes) % 2 == 0 ) ? 16 : 8 ),
              0 );
#if ( \
  defined(_IRIX)   || \
  defined(_IRIX64) || \
  defined(_SunOS)  || \
  defined(_HPUX))
  //  defined(_Darwin) )
        AX_img[iw]->byte_order = MSBFirst;
        AX_img[iw]->bitmap_bit_order = MSBFirst;
#else
        AX_img[iw]->byte_order = LSBFirst;
        AX_img[iw]->bitmap_bit_order = LSBFirst;
#endif
    }
    if ( AX_img[iw] == NULL )
        pe ("AX_plugin_ShmGC_module: cannot create image.\n");
    /* how many bytes of an image shared memory belong to a pixel */
    AX_bytes = AX_img[iw] -> bytes_per_line / AX_size[iw].width;
/*     pr ("bytes per pixel = %d\n", AX_bytes); */
    if (AX_bytes > sizeof(AX_Carrier))
        pe ("AX_plugin_ShmGC_module: AX_Carrier unable to carry\n"
            "one pixmap unit of %d bytes\n", AX_bytes);
    AX_8b  = AX_bytes & 1;
    AX_16b = AX_bytes & 2;
    /* calling it 24-bit is inconsistent because no one calls above 15-bit */
    AX_32b = AX_bytes & 4;
    if ( ! (AX_8b || AX_16b || AX_32b) )
        pe ("AX_plugin_ShmGC_module: AX_bytes = %d, neither 8,\n16, or 32-bit "
            "graphics. Small revisions required.\n", AX_bytes);
    if (AX_videomode)
    {
        AX_shminfo[iw].shmid =
            shmget ( IPC_PRIVATE, AX_img[iw]->bytes_per_line *
                     AX_size[iw].height,
                     IPC_CREAT | 0777 );
        if ( AX_shminfo[iw].shmid < 0)
            pe ("AX_plugin_ShmGC_module: cannot allocate shared memory\n");
        AX_shminfo[iw].shmaddr = (char *) shmat(AX_shminfo[iw].shmid, 0, 0);
        if ( AX_shminfo[iw].shmaddr != (char *) -1 )
        {
            if (AX_videomode == AX_VIDEOMODE_SHM_IMAGE)
                AX_img[iw]->data = AX_shminfo[iw].shmaddr; /* arm the XImage */
            AX_mem[iw].uc = (unsigned char *) AX_shminfo[iw].shmaddr;
        }
        else pe ("AX_plugin_ShmGC_module: failed to attach shared memory\n");
        AX_shminfo[iw].readOnly = False;

#ifdef _AX_USE_SHM
        /* ask X server to attach shm segment for read/write */
        if ( XShmAttach(AX_display[iw], &AX_shminfo[iw]) == 0 )
            pe ("AX_plugin_ShmGC_module:\n"
                "X server fails to attach shm for read/write\n");
#endif
        
#if ( ! ( \
  defined(_IRIX64) || \
  defined(_IRIX)   || \
  defined(_SunOS)  || \
  defined(_HPUX)   || \
  defined(_OSF1)   || \
  defined(_Darwin) ) )
        /* Unlink shared memory id such that block is automatically freed */
        /* when this process terminates (since after X server disconnects */
        /* the link count of this allocation becomes 0). */
        shmctl (AX_shminfo[iw].shmid, IPC_RMID, 0);
#endif
    }
    if (AX_videomode == AX_VIDEOMODE_SHM_PIXMAP)
    {
#ifdef _AX_USE_SHM
        AX_pixmap[iw] = XShmCreatePixmap
            ( AX_display[iw],
              AX_root[iw],
              (char *) AX_mem[iw].uc,
              &AX_shminfo[iw],
              (unsigned int) AX_size[iw].width,
              (unsigned int) AX_size[iw].height,
              AX_depth[iw] );
#endif
        XSetWindowBackgroundPixmap(AX_display[iw], AX_win[iw], AX_pixmap[iw]);
        /* server acts on Pixmap instead of Window for backbuffer drawing */
        AX_drawable[iw] = AX_pixmap[iw];
    }
    else /* AX_VIDEOMODE_NO_SHM || AX_VIDEOMODE_SHM_IMAGE */
    { /* by our AX.h discussion, server's default Drawable is Window */
        AX_drawable[iw] = AX_win[iw];
    }
    AX_gc[iw] = XCreateGC (AX_display[iw], AX_drawable[iw], 0, NULL);
    AX_namedbg(iw,AX_DEF_BACKGROUND);
    if (AX_gc[iw] == 0)
        pe ("AX_plugin_ShmGC_module:\n"
            "cannot create GC based on shm pixmap.\n");
    AX_Module_Is_PluggedIn [iw] [MODULE_ShmGC] = 1;
    return (1);
} /* end AX_plugin_ShmGC_module() */


/* Free GC and shared memory: you are left with an empty shell of a window */
int AX_plugout_ShmGC_module (int iw)
{
    AX_AssertWinID ("AX_plugout_ShmGC_module", iw);
    if ( ! AX_Module_Is_PluggedIn [iw] [MODULE_ShmGC] )
        pe ("AX_plugout_ShmGC_module:\n"
            "you are crazy, this module is not plugged in\n");
    XFreeGC (AX_display[iw], AX_gc[iw]); /* free GC structure */
    if (AX_videomode == AX_VIDEOMODE_SHM_PIXMAP)  /* free Pixmap structure */
        XFreePixmap (AX_display[iw], AX_pixmap[iw]);
#ifdef _AX_USE_SHM
    if (AX_videomode)    /* ask X server to detach itself from shm segment */
        XShmDetach (AX_display[iw], &AX_shminfo[iw]);
#endif
    if (AX_image_videomode) /* free XImage struct */
        XDestroyImage (AX_img[iw]); /* also frees client data-memory */
#if ( \
  defined(_IRIX64) || \
  defined(_IRIX)   || \
  defined(_SunOS)  || \
  defined(_HPUX)   || \
  defined(_OSF1)   || \
  defined(_Darwin) )
    shmctl (AX_shminfo[iw].shmid, IPC_RMID, 0);
#endif
    AXSYNC(iw);
    AX_Module_Is_PluggedIn [iw] [MODULE_ShmGC] = 0;
    return (1);
} /* end AX_plugout_ShmGC_module() */


/* There is no way to realloc shared memory, so plugout and plugin */
int AX_resize_ShmGC_module (int iw)
{
    AX_plugout_ShmGC_module (iw);
    AX_plugin_ShmGC_module (iw);
    return (1);
} /* end AX_resize_ShmGC_module() */


/* Determine whether the X server is running on local machine */
Bool AX_is_local (int iw)
{
    char * ptr;
    ptr = DisplayString (AX_display[iw]);
    /* printf ("%s\n", ptr); */
    return ( !ptr ||
             (*ptr==':') ||
             !strncmp(ptr,"localhost:0",11) ||
             !strncmp(ptr,"unix:",5) ||
             !strncmp(ptr,"local:",6) );
} /* AX_is_local */


/***********************************************************************/
/* Open on local X-server a window of width x height (= 0 means deft.  */
/* value) with title. If title = NULL, AX will assign a unique name.   */
/* Returns AX window ID granting the sole write permission to this     */
/* window. The caller must supply a unique nonzero ident. number "cid" */
/* (such as getpid() of forked processes) and from now on becomes the  */
/* care provider of this window, ensuring sufficient frame-rate etc.   */
/* One cid can own multiple window ID's, on condition of never writing */
/* to them at the same time, as they share the same display.           */
/***********************************************************************/

/* The policy is that we never depend on data from other, possibly  */
/* alive threads: the synchronization job is too hard. Therefore    */
/* although a lot of the info like the color depth are the same, we */
/* nevertheless query and save it again and again. The only savings */
/* we have in mind are those with same cid, with no synchronization */
/* problem. Most significantly, we'll save server resource DISPLAY. */

int AX_openwindow
(int cid, char *title, unsigned int width, unsigned int height)
{
    int iw, jw;
    int major, minor;
    unsigned int depth;
    XSetWindowAttributes WindowAttributes;
    Bool support_pixmap;

    if (cid == 0) pe ("AX_openwindow: invalid cid\n");
    
    /* find a vacant window slot */
    for ( iw = 0; (iw < AX_MAXWIN) && nonvacant(iw); iw++ );

    if ( iw == AX_MAXWIN )
    {
        fprintf(stderr, "AX_openwindow: AX_MAXWIN = %d reached\n",
                AX_MAXWIN);
        return(-1);
    }

    /* I vow no one could seg-fault my libraries unintentionally */
    if (title == NULL)
        snprintf ( AX_title[iw], AX_MAXSTRSIZE, "AXwin %d", iw);
    else 
        strncpy ( AX_title[iw], title, AX_MAXSTRSIZE );

    AX_hint[iw].width = (width == 0) ? AX_DEFWIDTH : width;
    AX_hint[iw].height = (height == 0) ? AX_DEFHEIGHT : height;
    AX_hint[iw].flags = PSize;

#ifdef _Linux
    for ( jw = 0; jw < AX_MAXWIN; jw++ )
        if ( ( jw != iw ) && (AX_cid[jw] == cid) )
        { /* there is an officemate */
            AX_display[iw]  =  AX_display[jw];
            AX_screen[iw]   =  AX_screen[jw];
            AX_visual[iw]   =  AX_visual[jw];
            AX_root[iw]     =  RootWindow (AX_display[iw], AX_screen[iw]);
            AX_colormap[iw] =  AX_colormap[jw];
            /* printf ("officemate %d found for %d.\n", jw, iw); */
            break;
        }
#else
    jw = AX_MAXWIN;
#endif

    if (jw == AX_MAXWIN)    /* there is no officemate */
    {
        AX_display[iw] = XOpenDisplay(NULL);
        if (AX_display[iw] == NULL)
            pe ("AX_openwindow: cannot open X display.\n");

        /* default screen */
        AX_screen[iw] = DefaultScreen (AX_display[iw]);

        /* root window */
        AX_root[iw] = RootWindow (AX_display[iw], AX_screen[iw]);

        /* visual type and color depth */
        if (!(AX_visual[iw] = AXBestVisual(iw)))
            pe ("AX_openwindow: no visuals available.\n");

        /* XCreateColormap() dynamically allocates ... */
        AX_colormap[iw] = XCreateColormap
            (AX_display[iw], AX_root[iw], AX_visual[iw], AllocNone);
        XInstallColormap(AX_display[iw], AX_colormap[iw]);
        AXSYNC(iw);

        /* Get those named colors: their pixel and carrier values */
        AX_update_namedcolors (iw);

        AX_noneedcolormap =
            ( AX_visual[iw]->class == DirectColor ) ||
            ( AX_visual[iw]->class == TrueColor );
        if ( AX_noneedcolormap )
        {
            AX_rmask = AX_visual[iw] -> red_mask;
            AX_gmask = AX_visual[iw] -> green_mask;
            AX_bmask = AX_visual[iw] -> blue_mask;
        }

#ifdef _AX_USE_SHM
        if ( ( !AX_is_local(iw) ) ||
             ( !XQueryExtension(AX_display[iw], "MIT-SHM",
                                &minor, &minor, &minor) ) ||
             ( !XShmQueryVersion(AX_display[iw], &major,
                                 &minor, &support_pixmap) ) )
        {
            AX_videomode = AX_VIDEOMODE_NO_SHM;
            /* pe ("AX_openwindow: your X server does not support\n" */
            /* "MIT shared memory extension. Sorry.\n"); */
        }
        else if ( (!support_pixmap) ||
                  (XShmPixmapFormat(AX_display[iw]) != ZPixmap) )
        {
#if ( defined(_Darwin) )
            AX_videomode = AX_VIDEOMODE_NO_SHM;
#else            
            pr ("AX_openwindow: your X server does not support\n"
                "shared memory pixmap in Z format (%d %d).\n"
                "Shared memory XShmPutImage is used instead.\n",
                support_pixmap, XShmPixmapFormat(AX_display[iw]));
            AX_videomode = AX_VIDEOMODE_SHM_IMAGE;
#endif
        }
        else AX_videomode = AX_VIDEOMODE_SHM_PIXMAP;
#else
        AX_videomode = AX_VIDEOMODE_NO_SHM;
#endif

    }  /* no officemate */

    /* window border width */
    AX_borderwidth[iw] = 0;

    WindowAttributes.background_pixel = AX_namedpixel [AX_BLACK];
    WindowAttributes.border_pixel = AX_namedpixel [AX_RED];
    WindowAttributes.colormap = AX_colormap[iw];
    WindowAttributes.override_redirect = AX_WINDOWS_OVERRIDE_REDIRECT;
    /* create window */
    AX_win [iw] = XCreateWindow
        ( AX_display[iw],
          AX_root[iw],
          (AXDisplayWidth(iw)-AX_hint[iw].width)/2,
          (AXDisplayHeight(iw)-AX_hint[iw].height)/2,
          AX_hint[iw].width, AX_hint[iw].height,
          AX_borderwidth [iw],
          AX_depth[iw],
          InputOutput,
          AX_visual[iw],
          CWBackPixel | CWBorderPixel | CWColormap | CWOverrideRedirect,
          &WindowAttributes);

    if ( AX_win[iw] == (Window)NULL )
        pe ("AX_openwindow: cannot open window titled \"%s\"\n",
            AX_title[iw]);

    /* event filter */
    XSelectInput ( AX_display[iw], AX_win[iw],
		   ButtonMotionMask |
                   /* Button1MotionMask | */
                   /* Button2MotionMask | */
                   /* Button3MotionMask | */
                   /* Button4MotionMask | */
                   /* Button5MotionMask | */
                   ButtonPressMask |
                   ButtonReleaseMask |
                   /* ColormapChangeMask | */
                   /* EnterWindowMask | */
                   /* LeaveWindowMask | */
                   ExposureMask |
                   /* FocusChangeMask | */
                   /* KeymapStateMask | */
                   KeyPressMask |
                   /* KeyReleaseMask | */
                   /* OwnerGrabButtonMask | */
                   /* PointerMotionMask | */
                   /* PointerMotionHintMask | */
                   /* PropertyChangeMask | */
                   /* ResizeRedirectMask | */
                   StructureNotifyMask
                   /* SubstructureNotifyMask | */
                   /* SubstructureRedirectMask | */
                   /* VisibilityChangeMask */
                   );

    /* title, icon, etc. */
    XSetStandardProperties
	( AX_display [iw], AX_win [iw], AX_title [iw],
	  AX_title [iw], None, NULL, 0, &AX_hint[iw] );
    
    /* pop this window up on the screen */
    XMapRaised (AX_display[iw], AX_win[iw]);

    /* confirm that we obtained the geometry */
    XGetGeometry
        ( AX_display [iw],
          AX_win [iw],
          &AX_root[iw],
          &AX_hint[iw].x,
          &AX_hint[iw].y,
          (unsigned int *) &AX_size[iw].width,
          (unsigned int *) &AX_size[iw].height,
          (unsigned int *) &AX_borderwidth [iw],
          &depth );
    
    if ( AX_size[iw].width != width )
        pe ("AX_openwindow: cannot obtain width %d\n"
            "in window titled \"%s\"\n", width, AX_title[iw]);
    
    if ( AX_size[iw].height != height )
        pe ("AX_openwindow: cannot obtain height %d\n"
            "in window titled \"%s\"\n", height, AX_title[iw]);

    if ( AX_depth[iw] != depth )
        pe ("AX_openwindow: cannot obtain depth %d\n"
            "in window titled \"%s\"\n", AX_depth[iw], AX_title[iw]);

    /* finally register this care provider */
    AX_cid[iw] = cid;

    /* plug in modules */
    for (minor=0; minor<=AX_DEF_Module_Support_Level; minor++)
        AX_modules[minor].plugin(iw);

    /* default pen color for X-server (non-direct) draws */
    AXSetNamedForeground (iw, AX_DEF_FOREGROUND);

    return (iw);
} /* end AX_openwindow() */


/* Resize all modules that are already plugged in, to AX_size[iw] */
void AX_resizewindow (int iw, Bool do_window_too)
{
    int module;
    AX_AssertWinID ("AX_resizewindow", iw);
    for (module=0; module<MODULE_MAX; module++)
        if (AX_Module_Is_PluggedIn[iw][module])
            AX_modules[module].resize(iw);
    if (do_window_too)
    {
        AXResizeWindow(iw);
        AXSYNC(iw);
    }
    return;
} /* end AX_resizewindow() */


/* close the window and free all resources */
void AX_closewindow (int iw)
{
    int jw;
    AX_AssertWinID ("AX_closewindow", iw);
    for (jw=0; jw<MODULE_MAX; jw++)
        if (AX_Module_Is_PluggedIn [iw] [jw])
            AX_modules[jw].plugout(iw);
    /* close the window */
    AXDestroyWindow(iw);
    AXSYNC(iw);
    for (jw = 0; jw < AX_MAXWIN; jw++)
        if ( (jw != iw) && (AX_cid[jw] == AX_cid[iw]) )
            break;
    /* if no officemate, also pull the phone line */
    if (jw == AX_MAXWIN)
    {
        XFreeColormap (AX_display[iw], AX_colormap[iw]);
        AXCloseDisplay(iw);
    }
    /* unregister this care provider */
    AX_cid [iw] = 0;
    return;
} /* end AX_closewindow() */


/* single care provider test */
#ifdef _single_TEST
int main (int argc, char *argv[])
{
    register int i;
    int iw, jw;
    iw = AX_openwindow (getpid(),NULL,AX_DEFWIDTH,AX_DEFHEIGHT);
    AX_namedbg(iw,AX_BLUE);
    AX_dump(iw);
    AXSetNamedForeground (iw, AX_RED);
    AXFillRectangle (iw,AX_DEFWIDTH/4,AX_DEFHEIGHT/4,
                     AX_DEFWIDTH/2, AX_DEFHEIGHT/2);
    AXDrawLine(iw,0,0,AX_DEFWIDTH-1,AX_DEFWIDTH-1);
    AX_show(iw);
    dumph (16, AX_mem[iw].uc);
    Press_return();
    AX_Bset(iw,0,AX_Colorcarrier(1.,1.,0.),AX_DEFWIDTH*10);
    dumph (16, AX_mem[iw].uc);
    Bwriteint (AX_mem[iw].uc); cr(); cr();
    AX_dump(iw);
    AX_show(iw);
    Press_return();
    jw = AX_openwindow (getpid(),NULL,AX_DEFWIDTH,AX_DEFHEIGHT);
    AX_closewindow(iw);
    AX_bset(jw,0,AX_Colorcarrier(0.,1.,1.),AX_DEFWIDTH*10,i);
    AX_dump(jw);
    AX_show(jw);
    Press_return();
    AX_closewindow(jw);
    return (0);
}
#endif /* _single_TEST */


/******************************************************/
/* Save what's in the pixmap to a .png file with zlib */
/* compression level: 0 (no compression) - 9 (maximal */
/* compression): Z_NO_COMPRESSION=0, Z_BEST_SPEED=1,  */
/* Z_DEFAULT_COMPRESSION=6, Z_BEST_COMPRESSION=9.     */
/* Consult http://www.cdrom.com/pub/png/libpng.html   */
/* and /usr/doc/libpng-1.0.5/example.c for details.   */
/******************************************************/
#define AX_PNG_TEXT_CHUNK 0
int AX_save_pixmap_as_png (int iw, int compression_level,
                           char *file_name)
{
    FILE *fp;
    png_structp png_ptr;
    png_infop info_ptr;
    png_colorp palette=NULL;
    png_bytep *row_pointers;
    XColor color;
    char *buffer=NULL;
    int i, rshift, gshift, bshift, rmag, gmag, bmag;
    AX_capture(iw);
    fp = wOpen(file_name);
    if (fp == NULL)
        pe ("AX_save_pixmap_as_png: cannot create write handle "
            "for \"%s\".\n", file_name);
    png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,
                                      NULL, NULL, NULL);
    if (png_ptr == NULL)
    {
        fclose(fp);
        pe ("AX_save_pixmap_as_png: png_create_write_struct failed.\n");
    }

    info_ptr = png_create_info_struct(png_ptr);
    if (info_ptr == NULL)
    {
        fclose(fp);
        png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
        pe ("AX_save_pixmap_as_png: png_create_info_struct failed.\n");
    }

    if (setjmp(png_ptr->jmpbuf))
    {
        fclose(fp);
        png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
        pe ("AX_save_pixmap_as_png: error using libpng/zlib "
            "for \"%s\".\n", file_name);
    }

    png_init_io(png_ptr, fp);

    AX_C( {
        png_set_IHDR (png_ptr, info_ptr,
                      AX_size[iw].width, AX_size[iw].height, 8,
                      PNG_COLOR_TYPE_PALETTE, AX_PNG_INTERLACE_TYPE,
                      PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
        /* allocate palette */
        palette = (png_colorp) png_malloc(png_ptr, 256*sizeof(png_color));
        for (i=0; i<256; i++)
        {
            color.pixel = i;
            XQueryColor(AX_display[iw], AX_colormap[iw], &color);
            palette[i].red = SAFE_RSHIFT
                ( INT(color.red),
                  (sizeof(unsigned short)-sizeof(png_byte)) * 8 );
            palette[i].green = SAFE_RSHIFT
                ( INT(color.green),
                  (sizeof(unsigned short)-sizeof(png_byte)) * 8 );
            palette[i].blue = SAFE_RSHIFT
                ( INT(color.blue),
                  (sizeof(unsigned short)-sizeof(png_byte)) * 8 ); 
        }
        png_set_PLTE(png_ptr, info_ptr, palette, 256);
    },
          png_set_IHDR (png_ptr, info_ptr,
                        AX_size[iw].width, AX_size[iw].height, 8, 
                        PNG_COLOR_TYPE_RGB, AX_PNG_INTERLACE_TYPE,
                        PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE),
          png_set_IHDR (png_ptr, info_ptr,
                        AX_size[iw].width, AX_size[iw].height, 8, 
                        PNG_COLOR_TYPE_RGB, AX_PNG_INTERLACE_TYPE,
                        PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE)
          );
    png_write_info (png_ptr, info_ptr);

    /* allocate row_pointers */
    row_pointers = (png_bytep *)malloc(AX_size[iw].height*sizeof(png_bytep));

    AX_C(
    {
        for (i=0; i<AX_size[iw].height; i++)
            row_pointers[i] = (png_bytep)
                (AX_mem[iw].uc + i * AX_size[iw].width);
    },
    { /* allocate buffer */
        for (rshift=0; !((AX_rmask>>rshift)&1); rshift++);
        rmag = 256 / ( (AX_rmask>>rshift) + 1 );
        for (gshift=0; !((AX_gmask>>gshift)&1); gshift++);
        gmag = 256 / ( (AX_gmask>>gshift) + 1 );
        for (bshift=0; !((AX_bmask>>bshift)&1); bshift++);
        bmag = 256 / ( (AX_bmask>>bshift) + 1 );
        buffer = (char *)malloc(3*AX_size[iw].width*AX_size[iw].height);
        for (i=0; i<AX_size[iw].width*AX_size[iw].height; i++)
        {
            buffer[3*i]   = ((AX_mem[iw].i2[i] & AX_rmask) >> rshift) * rmag;
            buffer[3*i+1] = ((AX_mem[iw].i2[i] & AX_gmask) >> gshift) * gmag;
            buffer[3*i+2] = ((AX_mem[iw].i2[i] & AX_bmask) >> bshift) * bmag;
        }
        for (i=0; i<AX_size[iw].height; i++)
            row_pointers[i] = (png_bytep)
                (buffer + i * 3 * AX_size[iw].width);
    },
    { /* allocate buffer */
        for (rshift=0; !((AX_rmask>>rshift)&1); rshift++);
        for (gshift=0; !((AX_gmask>>gshift)&1); gshift++);
        for (bshift=0; !((AX_bmask>>bshift)&1); bshift++);
        buffer = (char *)malloc(3*AX_size[iw].width*AX_size[iw].height);
        for (i=0; i<AX_size[iw].width*AX_size[iw].height; i++)
        {
            buffer[3*i]   = (AX_mem[iw].i4[i] & AX_rmask) >> rshift;
            buffer[3*i+1] = (AX_mem[iw].i4[i] & AX_gmask) >> gshift;
            buffer[3*i+2] = (AX_mem[iw].i4[i] & AX_bmask) >> bshift;
        }
        for (i=0; i<AX_size[iw].height; i++)
            row_pointers[i] = (png_bytep)
                (buffer + i * 3 * AX_size[iw].width);
    } );
    png_write_image(png_ptr, row_pointers);
    png_write_end(png_ptr, info_ptr);
    AX_C( free(palette), free(buffer), free(buffer) );
    free(row_pointers);
    png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
    i = ftell(fp);
    fclose(fp);
    return(i);
} /* end AX_save_pixmap_as_png() */


#ifdef _save_pixmap_as_png_TEST
#define FILENAME "/tmp/try.png"
int main (int argc, char *argv[])
{
    AX_openwindow (getpid(),NULL,AX_DEFWIDTH,AX_DEFHEIGHT);
    AX_namedbg(0, AX_BLUE);
    AX_dump(0);
    AXSetNamedForeground (0, AX_RED);
    AXFillRectangle (0,AX_DEFWIDTH/4,AX_DEFHEIGHT/4,
                     AX_DEFWIDTH/2, AX_DEFHEIGHT/2);
    AXDrawLine(0,0,0,AX_DEFWIDTH-1,AX_DEFWIDTH-1);
    AX_show(0);
    Press_return();
    printf ("Image saved on \"%s\" (%d bytes).\n", FILENAME,
            AX_save_pixmap_as_PNG(0, FILENAME));
    return (0);
}
#endif /* _save_pixmap_as_png_TEST */


/*****************************************************************/
/* Save what's in the pixmap to a .jpg file with quality 0 (bad) */
/* to 100 (no loss). Consult http://www.ijg.org/ for details.    */
/*****************************************************************/
int AX_save_pixmap_as_jpg (int iw, int quality, char *file_name)
{
    struct jpeg_compress_struct cinfo;
    struct jpeg_error_mgr jerr;
    FILE *outfile;
    JSAMPROW row_pointer[1];
    int i, rshift, gshift, bshift, rmag, gmag, bmag;
    register AX_Carrier c;
    char *palette;
    XColor color;
    AX_capture(iw);
    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_compress(&cinfo);
    if ((outfile = wOpen(file_name)) == NULL)
        pe("AX_save_pixmap_as_jpg: cannot open \"%s\"\n", file_name);
    jpeg_stdio_dest(&cinfo, outfile);
    cinfo.image_width = AX_size[iw].width;
    cinfo.image_height = AX_size[iw].height;
    cinfo.input_components = 3;
    cinfo.in_color_space = JCS_RGB;
    jpeg_set_defaults (&cinfo);
    /* enable maximal entropy optimization */
    cinfo.optimize_coding = TRUE;
    jpeg_set_quality(&cinfo, quality, TRUE);
    jpeg_start_compress(&cinfo, TRUE);
    row_pointer[0] = (JSAMPROW) malloc(cinfo.image_width*3);
    AX_C (
    {
        palette = (char *) malloc(256*3);
        for (i=0; i<256; i++)
        {
            color.pixel = i;
            XQueryColor(AX_display[iw], AX_colormap[iw], &color);
            palette[3*i]   = color.red   >> ((sizeof(unsigned short)-1)*8);
            palette[3*i+1] = color.green >> ((sizeof(unsigned short)-1)*8);
            palette[3*i+2] = color.blue  >> ((sizeof(unsigned short)-1)*8);
        }
        while (cinfo.next_scanline < cinfo.image_height)
        {
            for (i=0; i<cinfo.image_width; i++)
            {
                c = AX_mem[iw].uc[cinfo.next_scanline*cinfo.image_width+i];
                row_pointer[0][3*i]   = palette[3*c];
                row_pointer[0][3*i+1] = palette[3*c+1];
                row_pointer[0][3*i+2] = palette[3*c+2];
            }
            jpeg_write_scanlines (&cinfo, row_pointer, 1);
        }
        free(palette);
    },
    {
        for (rshift=0; !((AX_rmask>>rshift)&1); rshift++);
        rmag = 256 / ( (AX_rmask>>rshift) + 1 );
        for (gshift=0; !((AX_gmask>>gshift)&1); gshift++);
        gmag = 256 / ( (AX_gmask>>gshift) + 1 );
        for (bshift=0; !((AX_bmask>>bshift)&1); bshift++);
        bmag = 256 / ( (AX_bmask>>bshift) + 1 );
        while (cinfo.next_scanline < cinfo.image_height)
        {
            for (i=0; i<cinfo.image_width; i++)
            {
                c = AX_mem[iw].i2[cinfo.next_scanline*cinfo.image_width+i];
                row_pointer[0][3*i]   = ((c & AX_rmask) >> rshift) * rmag;
                row_pointer[0][3*i+1] = ((c & AX_gmask) >> gshift) * gmag;
                row_pointer[0][3*i+2] = ((c & AX_bmask) >> bshift) * bmag;
            }
            jpeg_write_scanlines (&cinfo, row_pointer, 1);
        }
    },
    {
        for (rshift=0; !((AX_rmask>>rshift)&1); rshift++);
        for (gshift=0; !((AX_gmask>>gshift)&1); gshift++);
        for (bshift=0; !((AX_bmask>>bshift)&1); bshift++);
        while (cinfo.next_scanline < cinfo.image_height)
        {
            for (i=0; i<cinfo.image_width; i++)
            {
                c = AX_mem[iw].i4[cinfo.next_scanline*cinfo.image_width+i];
                row_pointer[0][3*i]   = (c & AX_rmask) >> rshift;
                row_pointer[0][3*i+1] = (c & AX_gmask) >> gshift;
                row_pointer[0][3*i+2] = (c & AX_bmask) >> bshift;
            }
            jpeg_write_scanlines (&cinfo, row_pointer, 1);
        }
    } );
    free(row_pointer[0]);
    jpeg_finish_compress(&cinfo);
    i = ftell(outfile);
    fclose(outfile);
    jpeg_destroy_compress(&cinfo);
    return(i);
} /* end AX_save_pixmap_as_jpg() */


#ifdef _save_pixmap_as_jpg_TEST
#define FILENAME "/tmp/try.jpg"
int main (int argc, char *argv[])
{
    AX_openwindow (getpid(),NULL,AX_DEFWIDTH,AX_DEFHEIGHT);
    AX_namedbg(0, AX_BLUE);
    AX_dump(0);
    AXSetNamedForeground (0, AX_RED);
    AXFillRectangle (0,AX_DEFWIDTH/4,AX_DEFHEIGHT/4,
                     AX_DEFWIDTH/2, AX_DEFHEIGHT/2);
    AXDrawLine(0,0,0,AX_DEFWIDTH-1,AX_DEFWIDTH-1);
    AX_show(0);
    Press_return();
    printf ("Image saved on \"%s\" (%d bytes).\n", FILENAME,
            AX_save_pixmap_as_JPG(0, FILENAME));
    return (0);
}
#endif /* _save_pixmap_as_jpg_TEST */


/*****************************************************************/
/* Save what's in the pixmap to a .eps file with quality 0 (bad) */
/* to 100 (no loss). Consult http://www.ijg.org/ for details.    */
/*****************************************************************/
int AX_save_pixmap_as_eps (int iw, int quality, char *file_name)
{
    char *filename_jpg = TMPNAM(NULL);
    int i;
    AX_save_pixmap_as_jpg (iw, quality, filename_jpg);
    i = JPG2EPS (filename_jpg, file_name);
    remove (filename_jpg);
    return (i);
} /* end AX_save_pixmap_as_eps() */


#ifdef _save_pixmap_as_eps_TEST
#define FILENAME "/tmp/try.eps"
int main (int argc, char *argv[])
{
    AX_openwindow (getpid(),NULL,AX_DEFWIDTH,AX_DEFHEIGHT);
    AX_namedbg(0, AX_BLUE);
    AX_dump(0);
    AXSetNamedForeground (0, AX_GREEN);
    AXFillRectangle (0,AX_DEFWIDTH/4,AX_DEFHEIGHT/4,
                     AX_DEFWIDTH/2, AX_DEFHEIGHT/2);
    AXDrawLine(0,0,0,AX_DEFWIDTH-1,AX_DEFWIDTH-1);
    AX_show(0);
    Press_return();
    printf ("Image saved on \"%s\" (%d bytes).\n", FILENAME,
            AX_save_pixmap_as_EPS(0,FILENAME));
    return (0);
}
#endif /* _save_pixmap_as_eps_TEST */


/* Set the largest available icon from an n-line stream of .xpm data */
int AXSetIcon (int iw, int n, char *icon[])
{
    int i,j,width,height,palette,size;
    double ratio;
    XIconSize *size_list_return;
    XpmAttributes attrib={0};

    attrib.valuemask |= XpmCloseness;
    attrib.closeness = 40000;
    attrib.valuemask |= XpmRGBCloseness;
    attrib.red_closeness   = 40000;
    attrib.green_closeness = 40000;
    attrib.blue_closeness  = 40000;
    
    if (!XGetIconSizes(AX_display[iw],AX_root[iw],&size_list_return,&j))
    {
        /* printf ("Warning: no WM or which doesn't allow icon.\n\n"); */
        return(AX_WM_GIVES_NO_ICON_SIZES);
    }
    else
    {
        /* pr ("Icon info: %d %d %d %d %d.\n\n", j, */
        /* size_list_return[0].min_width, */
        /* size_list_return[0].min_height, */
        /* size_list_return[0].max_width, */
        /* size_list_return[0].max_height); */
    }
    ratio = 1000.; size = 1;
    for (i=0; i<n; i+=palette+height+1)
    {
        sscanf(icon[i], "%d %d %d", &width, &height, &palette);
        if ( (width  >= size_list_return[0].min_width) &&
             (height >= size_list_return[0].min_height) &&
             (width  <= size_list_return[0].max_width) &&
             (height <= size_list_return[0].max_height) &&
             ( (   width*height > size ) ||
               ( ( width*height == size ) &&
                 ( (double)width/height < ratio ) &&
                 ( (double)width/height > 1/ratio ) ) ) )
        {
            j = i;
            size = width*height;
            ratio = MAX( (double)width/height,
                         (double)height/width );
        }
    }
    if (size == 1)
    {
        printf ("Warning: minimum icon size of this WM is %d,\n",
                size_list_return[0].min_width);
        printf ("maximum icon size of this WM is %d.\n\n",
                size_list_return[0].max_width);
        return (AX_NO_SUITABLE_ICON_SIZES);
    }
    if ((i=XpmCreatePixmapFromData
         (AX_display[iw],AX_win[iw],&icon[j],&AX_icon_pixmap[iw],NULL,&attrib))
        == XpmSuccess )
    {
        AX_win_hints[iw] = XAllocWMHints();
        AX_win_hints[iw]->flags = IconPixmapHint;
        AX_win_hints[iw]->icon_pixmap = AX_icon_pixmap[iw];
        XSetWMHints(AX_display[iw], AX_win[iw], AX_win_hints[iw]);
        XFree(AX_win_hints[iw]);
    }
    else
    {
        pr("AXpmCreateIcon: icon allocated but insertion failed (%d).\n", i); 
    }
    XFree(size_list_return);
    return(0);
} /* end AXSetIcon() */


/* Check "win" and all its children for a window with "identifier" string in */
/* WM_COMMAND (set after '-xrm' and use xprop to see). Return NULL if fails. */
Window AXIdentifyWin (int iw, Window win, char *identifier)
{
    Window root_return, parent_return, *children_return;
    unsigned int i, nchildren_return;
    char **argv_return;
    int argc_return;
    if (XGetCommand (AX_display[iw], win, &argv_return, &argc_return))
    {
        for (i=0; i<argc_return; i++)
        {
            if (!strcmp(argv_return[i],identifier))
            {
                XFreeStringList(argv_return);
                return(win);
            }
        }
        if (argc_return>0) XFreeStringList(argv_return);
    }
    if (XQueryTree (AX_display[iw], win, &root_return,
                    &parent_return, &children_return, &nchildren_return))
    {
        for (i=0; i<nchildren_return; i++)
            if ((win=AXIdentifyWin(iw,children_return[i],identifier)))
            {
                XFree(children_return);
                return (win);
            }
        if (nchildren_return>0) XFree(children_return);
    }
    return((Window)NULL);
} /* end AXIdentifyWin() */


/* Obtain the visual with the largest color depth on the current display */
Visual *AXBestVisual (int iw)
{
    int i, nitems_return;
    XVisualInfo *VisualInfoList;
    Visual *BestVisual = NULL;
    AX_depth[iw] = 0;
    if ((VisualInfoList =
         XGetVisualInfo(AX_display[iw], VisualNoMask, NULL, &nitems_return)))
    {
        for (i=0; i<nitems_return; i++)
            if ( (VisualInfoList[i].depth > AX_depth[iw]) ||
                 ( (VisualInfoList[i].depth == AX_depth[iw]) &&
                   (VisualInfoList[i].class == TrueColor) ) )
            {
                AX_depth[iw] = VisualInfoList[i].depth;
                BestVisual = VisualInfoList[i].visual;
            }
        if (nitems_return > 0) XFree(VisualInfoList);
        /* pr ("visual id = 0x%lx, depth = %d\n", */
        /* BestVisual->visualid, AX_depth[iw]); */
        return(BestVisual);
    }
    return(BestVisual);
} /* end AXBestVisual() */
