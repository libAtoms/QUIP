/********************************************/
/* libAX: -lX11 -lXext -lpng -lz -ljpeg -lm */
/*        -lScalar -lIO -lTimer             */
/*                                          */
/* Accelerated low-level graphics library   */
/* with X Window Shared Memory Extension.   */
/*                                          */
/* Jan.13, 2000 Ju Li <liju99@mit.edu>      */
/********************************************/

/** AX, pixel, jpg2eps.c: **/
/* - multiple thread-safe resizable windows */
/* - 8, 16 and 32-bit colors */
/* - direct video access */
/* - .png .jpg .eps saves */

/** Scan, Sorter.c: **/
/* - floating-point drawing board */
/* - line/polygon/conics clipping & scan conversion */
/* - exact area-weight anti-aliasing */

/** 3D.c: **/
/* - ball caching */
/* - line/polygons */
/* - ellipsoid/cylinders */
/* - z-buffer & recombinant anti-aliasing */

#ifndef _AX_h
#define _AX_h

#define AX_VERSION "2.2"
/* use single-precision "float" instead of "double" */
#define AX_MEMORY_EFFICIENT
/* review alignment by "grep aligned *.c *.h" */

#ifndef _CYGWIN
#define _AX_USE_SHM
#endif

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/keysym.h>
#include <X11/extensions/XShm.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/types.h>
#ifdef _SunOS
#include <SunOS/xpm.h>
#else
#include <X11/xpm.h>
#endif
#include <png.h>
#include <jpeglib.h>
/* Ju Li's C performance libraries */
#include <Scalar.h>
#include <IO.h>
#include <Timer.h>


/* AX Window ID: an application wide system resource like UNIX fileno */
/* typedef int AX_WinID; */
#define AX_MAXWIN       8
#define AX_MAXSTRSIZE   80

/* "DEF" stands for default */
#define AX_DEFWIDTH  (64*10)
/* width is compatible with PII cacheline */
#define AX_DEFHEIGHT  AX_DEFWIDTH

/* named colors: provided by AX gratis */
#define AX_NAMEDCOLORS  8
#define AX_BLACK        0
#define AX_RED          1
#define AX_GREEN        2
#define AX_BLUE         3
#define AX_CYAN         4
#define AX_MAGENTA      5
#define AX_YELLOW       6
#define AX_WHITE        7
#define AX_DEF_BACKGROUND  AX_BLACK
#define AX_DEF_FOREGROUND  AX_WHITE

/* AX_Pixel (unsigned long): interface to XSetForeground() etc. */
typedef unsigned long AX_Pixel;
#define AX_PIXEL(c) ((AX_Pixel)(c))

/* AX_Carrier (unsigned int): interface to Direct Pixmap Access */
/* like AX_set(). sizeof(AX_Carrier) must >= possible AX_bytes. */
typedef unsigned int AX_Carrier;
#define AX_CARRIER(p)  ((AX_Carrier)(p))
#define AX_CARRIER2PIXEL(c)  AX_PIXEL(c)
#define AX_PIXEL2CARRIER(p)  AX_CARRIER(p)

/* AX_Offset (int): interface to raster pixel numbering like in a */
/* C array of [0..height-1][0..width-1] with upper left as origin */
typedef int AX_Offset;


/** AX.c: **/

/** Shared Public Resources **/

/************************************************************************/
/* By convention, Pixmap is allocated and acted on by the server, Image */
/* is allocated & acted on by the client (XCreateImage do not malloc    */
/* client data-memory but XDestroyImage free client datamemory). Server */
/* can draw on either Pixmap or Window called Drawable, and how to draw */
/* on a Drawable is in server-side GC. Pixmap is used for backbuffering */
/* because human perceives animation as small differences between       */
/* finished snapshots. Client/server communication of large block of    */
/* data is via IPC or TCP/IP XPutImage() and XGetSubImage(). This is    */
/* often the faster strategy for image (finished or almost finished in  */
/* client data/memory, then send to server) (vs. font or simple graph)  */
/* if a pixel may be redrawn several times, or if the client needs to   */
/* manipulate the image a lot; and also because sending a block of data */
/* to the server has smaller communication and execution overhead.      */
/************************************************************************/

extern int AX_videomode;
#define AX_VIDEOMODE_NO_SHM      0  /* conventional XPutImage */
#define AX_VIDEOMODE_SHM_IMAGE   1  /* faster XShmPutImage */
#define AX_VIDEOMODE_SHM_PIXMAP  2  /* act directly on Drawable: the best */

#define AX_image_videomode  ( (AX_videomode == AX_VIDEOMODE_NO_SHM) || \
  (AX_videomode == AX_VIDEOMODE_SHM_IMAGE) )

/*********************************************************************/
/* With shared memory Pixmap, the client can directly act on server  */
/* Pixmap while the server still can. This is the best way since     */
/* client and server's commands can now mix with appropriate sync.   */
/* (if the client actions are instantaneous, we just need to XSync() */
/* after every chunk of server commands to commit to Pixmap). The    */
/* finished Pixmap is an appropriate backbuffer (XCopyArea to show). */
/*                                                                   */
/* With only shared memory Image however, we are NOT dramatically    */
/* different from AX_VIDEOMODE_NO_SHM except XShmPutImage is faster  */
/* than XPutImage as the server can directly access client memory    */
/* rather than from IPC. But one is still limited by rectangle-wise  */
/* (often whole-scene) communication than pixel-wise teamwork. On    */
/* the issue of backbuffering, we PutImage to where? If to Window,   */
/* then the server commands must also be applied to the Window which */
/* are then not back-buffered. If to certain Pixmap, the server must */
/* allocate Pixmap memory and data traffic's also increased (Image-> */
/* Pixmap->Window rather than Image->Window). Because AX is very     */
/* powerful, Xlib actions are lightweight, e.g. font labeling. We    */
/* thus choose the first method: Image -> Window.                    */
/*                                                                   */
/* For configuration-independent application, you should AX-draw on  */
/* client data-memory AX_mem first, then AX_dump(), then use light   */
/* weight Xlib routines, and finally AX_show(). If Xlib actions are  */
/* not lightweight or if they must intermix with AX, you can either  */
/* bet on AX_VIDEOMODE_SHM_PIXMAP is supported, or not use AX_dump() */
/* AX_show() sequence. You are free and encouraged to do anything    */
/* including creating Pixmap,Image on your own. AX just does not     */
/* provide ready-to-use solution for such particular cases.          */
/*********************************************************************/

/* "TrueColor" or "DirectColor" visuals; read Notes/color.txt */
extern bool AX_noneedcolormap;
/* defined only when AX_noneedcolormap */
extern AX_Carrier AX_rmask, AX_gmask, AX_bmask;
/* how many bytes of the pixmap memory belongs to a pixel */
extern int AX_bytes, AX_8b, AX_16b, AX_32b;

/* prioritize jobs with respect to likely graphics types */
#define AX_c(j8,j16,j32) (AX_32b?(j32):(AX_16b?(j16):(j8)))
#define AX_C(j8,j16,j32) if (AX_32b) {j32;} else if (AX_16b) {j16;} else {j8;};
/* this is only for debugging jobs in the ground layer */

/** Individual Public Resources **/
/** accessible via a type-AX_WinID key often called "iw" **/

/* display is a connection, so can only be used by a single thread at a time */
extern Display *AX_display [AX_MAXWIN];

typedef struct
{
    int width;
    int height;
} AXSize;
extern AXSize AX_size [AX_MAXWIN];

extern Window AX_win [AX_MAXWIN];

extern Pixmap AX_pixmap [AX_MAXWIN];

extern Drawable AX_drawable [AX_MAXWIN];

/* If your server is 8-bit, you will probably need this a lot. */
extern Colormap AX_colormap [AX_MAXWIN];
/* colormaps from the same server / screen should be identical, */
/* unfortunately the initial process may die early and free it. */

/* shared memory: direct pixmap access like in mmap() */
typedef union
{
    unsigned char *uc;
    int2 *i2;
    int4 *i4;
    void *vd;
} AX_Shm;
extern AX_Shm AX_mem [AX_MAXWIN];

extern XImage *AX_img [AX_MAXWIN];

/* graphics context built on shared memory pixmap */
extern GC AX_gc [AX_MAXWIN];

/* Default interfaces via which AX routines communicate */
/* XEvent interface */
extern XEvent AX_event [AX_MAXWIN];
/* XGetGeometry interface */
extern unsigned int AX_depth [AX_MAXWIN];
extern int AX_borderwidth [AX_MAXWIN];
extern XSizeHints AX_hint [AX_MAXWIN];
/* XQueryPointer interface */
extern Window AX_root [AX_MAXWIN], AX_child [AX_MAXWIN];
/* mouse coordinates with respect to root */
extern int AX_root_x [AX_MAXWIN], AX_root_y [AX_MAXWIN];
/* mouse coordinates with respect to this window */
extern int AX_win_x [AX_MAXWIN], AX_win_y [AX_MAXWIN];
/* modifier keys and pointer buttons */
extern unsigned int AX_pointermask [AX_MAXWIN];
/* color icon creation (XpmCreatePixmapFromData,XSetWMHints) interface */
extern Pixmap AX_icon_pixmap[AX_MAXWIN];
extern XWMHints* AX_win_hints[AX_MAXWIN];
extern char AX_title [AX_MAXWIN] [AX_MAXSTRSIZE];
extern int AX_cid [AX_MAXWIN];
extern int AX_screen [AX_MAXWIN];
/* pop-up window toggle */
extern bool AX_WINDOWS_OVERRIDE_REDIRECT;

/*********************/
/* Window Management */
/*********************/

/* Check if this WinID is in use; it is supposed to */
void AX_AssertWinID (char *who, int iw);

/* Module List */
#define MODULE_ShmGC 0  /* shared memory GC instead of window GC */
#define MODULE_Scan  1  /* Scan conversion of 2D primitives */
#define MODULE_3D    2  /* 3D shapes rendering */
#define MODULE_BC    3  /* ball caching */
#define MODULE_MAX   4
/* module level that AX_openwindow() supports by default */
#define AX_DEF_Module_Support_Level MODULE_Scan

extern char *AX_Module_Name [MODULE_MAX];
extern int AX_Module_Is_PluggedIn [AX_MAXWIN] [MODULE_MAX];

/* Check whether module dependence is satisfied */
int AX_Check_Module_Dependence (char *who, int iw, int module_idx);

/* Create shared memory graphics context that is accessible */
/* both via AX_mem[] and standard Xlib GC pixmap draws. The */
/* size of AX_mem[] created is dependent on AX_size[iw].    */
int AX_plugin_ShmGC_module (int iw);

/* Free GC and shared memory: you are left with an empty shell of a window */
int AX_plugout_ShmGC_module (int iw);

/* There is no way to realloc shared memory, so plugout and plugin */
int AX_resize_ShmGC_module (int iw);

/* Determine whether the X server is running on local machine */
Bool AX_is_local (int iw);

/***********************************************************************/
/* Open on local X-server a window of width x height (<= 0 means def.  */
/* value) with title. If title = NULL, AX will assign a unique name.   */
/* Returns AX window ID granting the sole write permission to this     */
/* window. The caller must supply a unique nonzero ident. number "cid" */
/* (such as getpid() of forked processes) and from now on becomes the  */
/* care provider of this window, ensuring sufficient frame-rate etc.   */
/* One cid can own multiple window ID's, on condition of never writing */
/* to them at the same time, as they share the same display.           */
/***********************************************************************/
int AX_openwindow
(int cid, char *title, unsigned int width, unsigned int height);

/* Resize all modules that are already plugged in, to AX_size[iw] */
void AX_resizewindow (int iw, Bool do_window_too);

/* close the window and free all resources */
void AX_closewindow (int iw);

/******************************************************/
/* Save what's in the pixmap to a .png file with zlib */
/* compression level: 0 (no compression) - 9 (maximal */
/* compression): Z_NO_COMPRESSION=0, Z_BEST_SPEED=1,  */
/* Z_DEFAULT_COMPRESSION=6, Z_BEST_COMPRESSION=9.     */
/* Consult http://www.cdrom.com/pub/png/libpng.html   */
/* and /usr/doc/libpng-1.0.5/example.c for details.   */
/******************************************************/
#define AX_PNG_INTERLACE_TYPE  PNG_INTERLACE_NONE
#define AX_PNG_DEF_COMPRESSION_LEVEL  Z_DEFAULT_COMPRESSION
int AX_save_pixmap_as_png
(int iw, int compression_level, char *file_name);
#define AX_save_pixmap_as_PNG(iw,file_name) \
  AX_save_pixmap_as_png(iw,AX_PNG_DEF_COMPRESSION_LEVEL,file_name)

/*****************************************************************/
/* Save what's in the pixmap to a .jpg file with quality 0 (bad) */
/* to 100 (no loss). Consult http://www.ijg.org/ for details.    */
/*****************************************************************/
#define AX_JPG_DEF_QUALITY 95
int AX_save_pixmap_as_jpg (int iw, int quality, char *file_name);
#define AX_save_pixmap_as_JPG(iw,file_name) \
  AX_save_pixmap_as_jpg(iw,AX_JPG_DEF_QUALITY,file_name)

/* conversion from jpeg quality (0..100) to png compression level (0..9) */
#define AX_JPG_QUALITY_TO_PNG_LEVEL(jpg) \
  INT( (1.-jpg/100.) * Z_BEST_COMPRESSION )

/* conversion from png compression level (0..9) to jpeg quality (0..100) */
#define AX_PNG_LEVEL_TO_JPG_QUALITY(png) \
  INT( (1.-DOUBLE(png)/Z_BEST_COMPRESSION) * 100. )

/*****************************************************************/
/* Save what's in the pixmap to a .eps file with quality 0 (bad) */
/* to 100 (no loss). Consult http://www.ijg.org/ for details.    */
/*****************************************************************/
#define AX_EPS_DEF_QUALITY AX_JPG_DEF_QUALITY  /* same scale as JPG */
int AX_save_pixmap_as_eps (int iw, int quality, char *file_name);
#define AX_save_pixmap_as_EPS(iw,file_name) \
  AX_save_pixmap_as_eps(iw,AX_EPS_DEF_QUALITY,file_name)

#define AX_WM_GIVES_NO_ICON_SIZES 1
#define AX_NO_SUITABLE_ICON_SIZES 2
/* Set the largest available icon from an n-line stream of .xpm data */
int AXSetIcon (int iw, int n, char *icon[]);
#define AXSetICON(iw,icon) AXSetIcon(iw,sizeof(icon)/sizeof(char *),icon)
#define AXSETICON(iw) AXSetICON(iw,icon)

/* Check "win" and all its children for a window with "identifier" string in */
/* WM_COMMAND (set after '-xrm' and use xprop to see). Return NULL if fails. */
Window AXIdentifyWin (int iw, Window win, char *identifier);
#define AXIdentifyWIN(iw,identifier) AXIdentifyWin(iw,AX_root[iw],identifier)

/* Obtain the visual with the largest color depth on the current display */
Visual *AXBestVisual (int iw);


/* jpg2eps.c: */

/**********************************************/
/* Bundled jpeg2ps (C) 1994-1999 Thomas Merz  */
/* http://www.pdflib.com/jpeg2ps/. See also   */
/* /usr/share/ghostscript/6.0/lib/viewjpeg.ps */
/**********************************************/
int jpg2eps (char *filename_jpg, char *filename_eps, char *options);
#define AX_EPS_DEF_DPI 300
#define JPG2EPS(filename_jpg,filename_eps) \
  jpg2eps(filename_jpg,filename_eps,"-r "STR(AX_EPS_DEF_DPI))

/* pixel.c: */

/*******************************************/
/* Transferable Color and Pixel operations */
/*******************************************/

/** AX_Pixel (unsigned long): interface to XSetForeground() etc. **/

/* color name must exist in "/usr/X11/lib/X11/rgb.txt" */
AX_Pixel AX_color_pixel_via_name (int iw, char *colorname);

/**************************************************************/
/* Pixels provided by AX gratis. Each call to AX_openwindow() */
/* updates them, in case some applications "steal our colors" */
/**************************************************************/
extern AX_Pixel AX_namedpixel [AX_NAMEDCOLORS];

/* recombining floats */
#define AX_pixel(r,g,b) ( \
  (AX_PIXEL(r)&AX_rmask) | \
  (AX_PIXEL(g)&AX_gmask) | \
  (AX_PIXEL(b)&AX_bmask) )

/* real number interface:  0.0 <= r,g,b <= 1.0 */
#define AX_Colorpixel(r,g,b) AX_pixel((r)*AX_rmask, (g)*AX_gmask, (b)*AX_bmask)

/* integer interface: r,g,b is in 0..MAX */
#define AX_ColorPixel(r,g,b,MAX) \
  AX_Colorpixel(iw, (double)(r)/(MAX), (double)(g)/(MAX), (double)(b)/(MAX))

/* 0xFFFFFF interface */
#define AX_COLORPIXEL(c) \
  AX_ColorPixel(iw, (c)>>16, ((c)>>8)&0xFF, (c)&0xFF, 0xFF)

/** AX_Carrier: interface to Direct Pixmap Access like AX_set() **/
extern AX_Carrier AX_namedcarrier [AX_NAMEDCOLORS];

/* caretaker iw updates community named color pixels and carriers */
void AX_update_namedcolors (int iw);

/* Identical to AX_color_pixel_via_name() except we may save some memory */
AX_Carrier AX_color_carrier_via_name (int iw, char *colorname);

#define AX_CMAP_JET     0
#define AX_CMAP_HOT     1
#define AX_CMAP_COOL    2
#define AX_CMAP_GRAY    3
#define AX_CMAP_PINK    4
#define AX_CMAP_BONE    5
#define AX_CMAP_COPPER  6
#define AX_CMAP_AUTUMN  7
#define AX_CMAP_SPRING  8
#define AX_CMAP_WINTER  9
#define AX_CMAP_SUMMER  10
#define AX_CMAP_HSV     11
#define AX_MAX_CMAP     12

typedef struct
{
    char *name;
    char *description;
    void (*fun) (double x, double *r, double *g, double *b);
}  AX_Cmap;
extern const AX_Cmap AX_cmap_funs[AX_MAX_CMAP];

/* Calculate (r,g,b) for x in [0,1] according to scheme "cmap_idx". */
void AX_cmap (int cmap_idx, double x, double *r, double *g, double *b);


/* recombining floats */
#define AX_carrier(r,g,b) ( \
  (AX_CARRIER(r)&AX_rmask) | \
  (AX_CARRIER(g)&AX_gmask) | \
  (AX_CARRIER(b)&AX_bmask) )

/* real number interface: 0.0 <= r,g,b <= 1.0 */
#define AX_Colorcarrier(r,g,b) \
  AX_carrier((r)*AX_rmask, (g)*AX_gmask, (b)*AX_bmask)

/* integer interface: r,g,b is in 0..MAX */
#define AX_ColorCarrier(iw,r,g,b,MAX) \
  AX_Colorcarrier(iw, (double)(r)/(MAX), (double)(g)/(MAX), (double)(b)/(MAX))

/* 0xFFFFFF interface */
#define AX_COLORCARRIER(iw,c) \
  AX_ColorCarrier(iw, (c)>>16, ((c)>>8)&0xFF, (c)&0xFF, 0xFF)

/* Mixing two color carriers */
#define AX_mix(c0,c1,a1) AX_carrier( \
  (1-(a1))*((c0)&AX_rmask)+(a1)*((c1)&AX_rmask), \
  (1-(a1))*((c0)&AX_gmask)+(a1)*((c1)&AX_gmask), \
  (1-(a1))*((c0)&AX_bmask)+(a1)*((c1)&AX_bmask) )

/* #define AX_YSC(iw,j)    (AX_size[iw].height-1-j) */
#define AX_YSC(iw,j)  (j)
#define AX_OFF(iw,i,j)  (AX_size[iw].width*INT(AX_YSC(iw,j))+INT(i))

/************************/
/* Direct Pixmap Access */
/************************/

#define AX_uc(iw,offset) (AX_mem[iw].uc+(offset))
#define AX_i2(iw,offset) (AX_mem[iw].i2+(offset))
#define AX_i4(iw,offset) (AX_mem[iw].i4+(offset))
#define AX_UC(iw,offset)  AX_mem[iw].uc[offset]
#define AX_I2(iw,offset)  AX_mem[iw].i2[offset]
#define AX_I4(iw,offset)  AX_mem[iw].i4[offset]

/* get the color of a pixel (type AX_Pixel) */
#define AX_getpixel(iw,offset) \
  AX_PIXEL(AX_c(AX_UC(iw,offset),AX_I2(iw,offset),AX_I4(iw,offset)))

/* get the color of a pixel (type AX_Carrier) */
#define AX_get(iw,offset) \
  AX_CARRIER(AX_c(AX_UC(iw,offset),AX_I2(iw,offset),AX_I4(iw,offset)))
#define AX_GET(iw,i,j) AX_get(iw,AX_OFF(iw,i,j))

/* set the color of a pixel */
#define AX_set(iw,offset,carrier) AX_C(AX_UC(iw,offset)=carrier, \
  AX_I2(iw,offset)=carrier, AX_I4(iw,offset)=carrier)
#define AX_SET(iw,i,j) AX_set(iw,AX_OFF(iw,i,j))

/* setting a single color to a small chunk of shared memory */
#define AX_bset(iw,offset,carrier,count,i) AX_C( \
  for (i=(offset); i<(offset)+(count); i++) AX_UC(iw,i)=carrier, \
  for (i=(offset); i<(offset)+(count); i++) AX_I2(iw,i)=carrier, \
  for (i=(offset); i<(offset)+(count); i++) AX_I4(iw,i)=carrier );

/* setting a single color to a large chunk of shared memory */
#define AX_Bset(iw,offset,carrier,count) AX_C( \
  ucharset(AX_uc(iw,offset), carrier, count), \
  int2set(AX_i2(iw,offset), carrier, count), \
  int4set(AX_i4(iw,offset), carrier, count))

/* set the entire pixmap to a certain color */
#define AX_mop(iw,carrier) \
  AX_Bset(iw,0,carrier,AX_size[iw].width*AX_size[iw].height)
#define AX_namedmop(iw,carrier_idx) AX_mop(iw,AX_namedcarrier[carrier_idx])

/* set the entire pixmap to 0 */
#define AX_bzero(iw) \
  bzero(AX_mem[iw].vd,AX_size[iw].width*AX_size[iw].height*AX_bytes)

#define AX_bg(iw,carrier) \
  { if ((carrier)==0) AX_bzero(iw); else AX_mop(iw,carrier); }
#define AX_namedbg(iw,carrier_idx) AX_bg(iw,AX_namedcarrier[carrier_idx])

/* set the entire pixmap to certain byte value */
#define AX_byteset(iw,c) \
  memset(AX_mem[iw].vd,(int)(c),AX_size[iw].width*AX_size[iw].height*AX_bytes)


/******************************************************/
/* Simplified X protocols using default AX interfaces */
/******************************************************/

#ifdef _AX_USE_SHM
/* For AX_VIDEOMODE_NO_SHM and AX_VIDEOMODE_SHM_IMAGE video modes, */
#define AX_dump(iw) { if (AX_videomode == AX_VIDEOMODE_NO_SHM) \
  XPutImage(AX_display[iw], AX_drawable[iw], AX_gc[iw], AX_img[iw], \
  0, 0, 0, 0, AX_size[iw].width, AX_size[iw].height); \
  else if (AX_videomode == AX_VIDEOMODE_SHM_IMAGE) \
  XShmPutImage(AX_display[iw], AX_drawable[iw], AX_gc[iw], AX_img[iw], \
  0, 0, 0, 0, AX_size[iw].width, AX_size[iw].height, False); }
#else
#define AX_dump(iw) { \
  XPutImage(AX_display[iw], AX_drawable[iw], AX_gc[iw], AX_img[iw], \
  0, 0, 0, 0, AX_size[iw].width, AX_size[iw].height); }
#endif

/* dump what is in client datamemory to server Drawable (DEF=Window) */

/* For AX_VIDEOMODE_SHM_PIXMAP video mode, copy Pixmap to Window; */
#define AX_show(iw) { if (AX_videomode == AX_VIDEOMODE_SHM_PIXMAP) \
  XCopyArea(AX_display[iw], AX_pixmap[iw], AX_win[iw], AX_gc[iw], 0,0, \
  AX_size[iw].width, AX_size[iw].height, 0,0); AXSYNC(iw); }
/* then every video mode asks the server to truly finish the job. */

/* For AX_VIDEOMODE_NO_SHM and AX_VIDEOMODE_SHM_IMAGE video modes, */
/* #define AX_capture(iw) { AXSYNC(iw); \ */
/* if (AX_videomode == AX_VIDEOMODE_NO_SHM) \ */
/* XGetSubImage ( AX_display[iw], AX_drawable[iw], 0,0, AX_size[iw].width, \ */
/* AX_size[iw].height, ~0L, ZPixmap, AX_img[iw], 0,0 ); \ */
/* else if (AX_videomode == AX_VIDEOMODE_SHM_IMAGE) \ */
/* XShmGetImage ( AX_display[iw], AX_drawable[iw], AX_img[iw], 0,0, ~0L ); \ */
/* AXSYNC(iw); } */
#define AX_capture(iw) {}
/* capture what is on screen to client datamemory via AX_img[].    */

#define AXDisplayWidth(iw) DisplayWidth(AX_display[iw],AX_screen[iw])
#define AXDisplayHeight(iw) DisplayHeight(AX_display[iw],AX_screen[iw])

/* XGetGeometry() via a type AXSize "newsize" */
#define AXGetGeometry(iw,newsize) XGetGeometry(AX_display[iw], \
  AX_win[iw], &AX_root[iw], &AX_hint[iw].x, &AX_hint[iw].y, \
  (unsigned int *) &newsize.width, (unsigned int *) &newsize.height, \
  (unsigned int *) &AX_borderwidth[iw], &AX_depth[iw])

#define AXResizeWindow(iw) XResizeWindow(AX_display[iw], AX_win[iw], \
  AX_size[iw].width, AX_size[iw].height)

#define AXDestroyWindow(iw) XDestroyWindow(AX_display[iw],AX_win[iw])

#define AXCloseDisplay(iw) XCloseDisplay(AX_display[iw])

#define AXSetForeground(iw,pixel) \
  XSetForeground(AX_display[iw],AX_gc[iw],pixel)

#define AXSetNamedForeground(iw,namedpixel_idx) \
  AXSetForeground(iw,AX_namedpixel[namedpixel_idx])

#define AXSetBackground(iw,pixel) \
  XSetBackground(AX_display[iw],AX_gc[iw],pixel)

#define AXSetNamedBackground(iw,namedpixel_idx) \
  AXSetBackground(iw,AX_namedpixel[namedpixel_idx])

#define AXSetLineAttributes(iw,line_width,line_style,cap_style,join_style) \
  XSetLineAttributes(AX_display[iw],AX_gc[iw],\
  line_width,line_style,cap_style,join_style)

#define AXSetDashes(iw,dash_offset,dash_list,n) \
  XSetDashes(AX_display[iw],AX_gc[iw],dash_offset,dash_list,n)

#define AXSetArcMode(iw,arc_mode) \
  XSetArcMode(AX_display[iw],AX_gc[iw],arc_mode)

#define AXSetFillStyle(iw,fill_style) \
  XSetFillStyle(AX_display[iw],AX_gc[iw],fill_style)

#define AXSetFillRule(iw,fill_rule) \
  XSetFillRule(AX_display[iw],AX_gc[iw],fill_rule)

#define AXSetFont(iw,font)  XSetFont(AX_display[iw],AX_gc[iw],font)

#define AXDrawPoint(iw,x,y) \
  XDrawPoint(AX_display[iw],AX_drawable[iw],AX_gc[iw],x,y)

#define AXDrawPoints(iw,points,npoints,mode) \
  XDrawPoints(AX_display[iw],AX_drawable[iw],AX_gc[iw],points,npoints,mode)

#define AXDrawLine(iw,x1,y1,x2,y2) \
  XDrawLine(AX_display[iw],AX_drawable[iw],AX_gc[iw],x1,y1,x2,y2)

#define AXDrawLines(iw,points,npoints,mode) \
  XDrawLines(AX_display[iw],AX_drawable[iw],AX_gc[iw],points,npoints,mode)

#define AXDrawPolygon(iw,points,npoints) \
  ( AXDrawLines(iw,points,npoints,CoordModeOrigin), AXDrawLine \
  (iw,points[npoints-1].x,points[npoints-1].y,points[0].x,points[0].y) )

#define AXDrawSegments(iw,segments,nsegments) \
  XDrawSegments(AX_display[iw],AX_drawable[iw],AX_gc[iw],segments,nsegments)

#define AXDrawRectangle(iw,x,y,width,height) \
  XDrawRectangle(AX_display[iw],AX_drawable[iw],AX_gc[iw],x,y,width,height)

#define AXDrawRectangles(iw,rectangles,nrectangles) \
  XDrawRectangles(AX_display[iw],AX_drawable[iw],AX_gc[iw],\
  rectangles,nrectangles)

#define AXFillRectangle(iw,x,y,width,height) \
  XFillRectangle(AX_display[iw],AX_drawable[iw],AX_gc[iw],x,y,width,height)

#define AXFillRectangles(iw,rectangles,nrectangles) \
  XFillRectangles(AX_display[iw],AX_drawable[iw],AX_gc[iw],rectangles,\
  nrectangles)

#define AXFillPolygon(iw,points,npoints,shape,mode) \
  XFillPolygon(AX_display[iw],AX_drawable[iw],AX_gc[iw],points,npoints,\
  shape,mode)

#define AXDrawArc(iw,x,y,width,height,angle1,angle2) \
  XDrawArc(AX_display[iw],AX_drawable[iw],AX_gc[iw],x,y,width,height,\
  angle1,angle2)

#define AXDrawArcs(iw,arcs,narcs) \
  XDrawArcs(AX_display[iw],AX_drawable[iw],AX_gc[iw],arcs,narcs)

#define AXFillArc(iw,x,y,width,height,angle1,angle2) \
  XFillArc(AX_display[iw],AX_drawable[iw],AX_gc[iw],x,y,width,height,\
  angle1,angle2)

#define AXFillArcs(iw,arcs,narcs) \
  XFillArcs(AX_display[iw],AX_drawable[iw],AX_gc[iw],arcs,narcs)

#define AXDrawString(iw,x,y,string,length) \
  XDrawString(AX_display[iw],AX_drawable[iw],AX_gc[iw],x,y,string,length)

#define AXDrawString16(iw,x,y,string,length) \
  XDrawString16(AX_display[iw],AX_drawable[iw],AX_gc[iw],x,y,string,length)

#define AXFlush(iw) XFlush(AX_display[iw])

#define AXSync(iw,discard) XSync(AX_display[iw],discard)

/* No, keep the events please. */
#define AXSYNC(iw) AXSync(iw,False)

#define AXEventsQueued(iw,mode) XEventsQueued(AX_display[iw],mode)

#define AXPending(iw) XPending(AX_display[iw])

#define AXNextEvent(iw) XNextEvent(AX_display[iw],&AX_event[iw])

#define AXPressKey(iw) while(1) { AXNextEvent(iw); \
  if (AX_event[iw].type == KeyPress) break; }

/* with modifiers */
#define AXKeysym(iw) XKeycodeToKeysym\
  (AX_display[iw], AX_event[iw].xkey.keycode, AX_event[iw].xkey.state)

/* naked key: "suppose without modifiers" */
#define AXKeysym0(iw) XKeycodeToKeysym\
  (AX_display[iw], AX_event[iw].xkey.keycode, 0)

/* http://www.mit.edu/afs/net/tools/diff.src/sun4src/play/play.c */
#define AXKeySYM(iw) ( AXKeysym(iw) ? AXKeysym(iw) : AXKeysym0(iw) )

#define AXCTRL(iw) (((XKeyEvent *)&AX_event[iw])->state & ControlMask)
#define AXCtrl(iw) ONE_OR_ZERO( AXCTRL(iw) )
#define AXSHFT(iw) (((XKeyEvent *)&AX_event[iw])->state & ShiftMask)
#define AXShft(iw) ONE_OR_ZERO( AXSHFT(iw) )
#define AXLOCK(iw) (((XKeyEvent *)&AX_event[iw])->state & LockMask)
#define AXLock(iw) ONE_OR_ZERO( AXLOCK(iw) )

#ifdef _SunOS
#define AXMETA(iw) (((XKeyEvent *)&AX_event[iw])->state & \
  (Mod1Mask | Mod2Mask | Mod3Mask | Mod4Mask | Mod5Mask))
#else
#define AXMETA(iw) (((XKeyEvent *)&AX_event[iw])->state & \
  (Mod1Mask))
#endif

#define AXMeta(iw) ONE_OR_ZERO( AXMETA(iw) )
#define AXPRESSBTN(iw) (AX_event[iw].xbutton.button)
#define AXBTNTIME(iw) (AX_event[iw].xbutton.time)
/* MotionNotify events */
#define AXMOTIONBTN1(iw) (AX_event[iw].xmotion.state & Button1Mask)
#define AXMOTIONBTN2(iw) (AX_event[iw].xmotion.state & Button2Mask)
#define AXMOTIONBTN3(iw) (AX_event[iw].xmotion.state & Button3Mask)

#define AXPeekEvent(iw) XPeekEvent(AX_display[iw],&AX_event[iw])

#define AXWindowEvent(iw,event_mask) \
  XWindowEvent(AX_display[iw],AX_win[iw],event_mask,&AX_event[iw])

#define AXCheckWindowEvent(iw,event_mask) \
  XCheckWindowEvent(AX_display[iw],AX_win[iw],event_mask,&AX_event[iw])

#define AXMaskEvent(iw,event_mask) \
  XMaskEvent(AX_display[iw],event_mask,&AX_event[iw])

#define AXCheckMaskEvent(iw,event_mask) \
  XCheckMaskEvent(AX_display[iw],event_mask,&AX_event[iw])

#define AXCheckTypedEvent(iw,event_type) \
  XCheckTypedEvent(AX_display[iw],event_type,&AX_event[iw])

#define AXCheckTypedWindowEvent(iw,event_type) \
  XCheckTypedWindowEvent(AX_display[iw],AX_win[iw],event_type,&AX_event[iw])

#define AXIfEvent(iw,predicate,arg) \
  XIfEvent(AX_display[iw],&AX_event[iw],predicate,arg)

#define AXCheckIfEvent(iw,predicate,arg) \
  XCheckIfEvent(AX_display[iw],&AX_event[iw],predicate,arg)

#define AXPeekIfEvent(iw,predicate,arg) \
  XPeekIfEvent(AX_display[iw],&AX_event[iw],predicate,arg)

#define AXPutBackEvent(iw) XPutBackEvent(AX_display[iw],&AX_event[iw])

#define AXKeycodeToKeysym(iw,keycode,index) \
  XKeycodeToKeysym(AX_display[iw],keycode,index)

#define AXQueryPointer(iw) \
 XQueryPointer(AX_display[iw],AX_win[iw],&AX_root[iw],&AX_child[iw],\
 &AX_root_x[iw],&AX_root_y[iw],&AX_win_x[iw],&AX_win_y[iw],&AX_pointermask[iw])

#define AXQueryPointerInWindow(iw) \
  ( AXQueryPointer(iw) && AX_IJINWIN(iw,AX_win_x[iw],AX_win_y[iw]) )

#define AXWarpPointer(iw,x,y) \
  XWarpPointer(AX_display[iw],None,AX_root[iw],0,0,0,0,x,y)

/* window title bar name */
#define AXStoreName(iw) XStoreName(AX_display[iw],AX_win[iw],AX_title[iw])
/* window icon name */
#define AXSetIconName(iw) XSetIconName(AX_display[iw],AX_win[iw],AX_title[iw])
/* both */
#define AXSetName(iw) (AXStoreName(iw),AXSetIconName(iw))

#define AXLowerWindow(iw) XLowerWindow(AX_display[iw],AX_win[iw])
#define AXLOWERWindow(iw) (AXLowerWindow(iw), AXSYNC(iw))
#define AXRaiseWindow(iw) XRaiseWindow(AX_display[iw],AX_win[iw])
#define AXRAISEWindow(iw) (AXRaiseWindow(iw), AXSYNC(iw))

#define AXQLength(iw) XQLength(AX_display[iw])

/** Scan.c: **/

/************************************************************/
/* Scan-convert 2D primitives with floating-point precision */
/************************************************************/

/* physical limits */
#if defined(_SunOS)
#define AX_MAXWIDTH   1536
#define AX_MAXHEIGHT  AX_MAXWIDTH
#elif defined(_alpha)
#define AX_MAXWIDTH   1200
#define AX_MAXHEIGHT  AX_MAXWIDTH
#else
/* #define AX_MAXWIDTH   16384 */
/* #define AX_MAXWIDTH   8192 */
/* #define AX_MAXWIDTH   4096 */
/* #define AX_MAXWIDTH   3072 */
#define AX_MAXWIDTH   2560
#define AX_MAXHEIGHT  AX_MAXWIDTH
#endif
/* floating-point drawing board:  [0,width] x [0,height] */
#ifdef AX_MEMORY_EFFICIENT
typedef float AX_Float;
#define AX_INFINITY   SINGLE_PRECISION_INFINITY
#define AX_TINY      (3e-6)
#else
#define AX_INFINITY   DOUBLE_PRECISION_INFINITY
#define AX_TINY      (1e-12)
typedef double AX_Float;
#endif
#define AX_FLOAT(f) ((AX_Float)(f))
/* seven effective digits is way enough for raster devices */
#define AX_XINWIN(iw,x)      XINW(x,AX_size[iw].width)
#define AX_XOUWIN(iw,x)      XOUW(x,AX_size[iw].width)
#define AX_YINWIN(iw,y)      XINW(y,AX_size[iw].height)
#define AX_YOUWIN(iw,y)      XOUW(y,AX_size[iw].height)
#define AX_XYINWIN(iw,x,y)  (AX_XINWIN(iw,x) && AX_YINWIN(iw,y))
#define AX_XYOUWIN(iw,x,y)  (AX_XOUWIN(iw,x) || AX_YOUWIN(iw,y))

/* labelling pixels of the drawing board: [0,width-1] x [0,height-1] */
#ifdef AX_MEMORY_EFFICIENT
typedef short int AX_IJ;
#else
typedef short int AX_IJ;
#endif
/* raster coordination won't use up 30K */
#define AX_IINWIN(iw,i)      INW(i,AX_size[iw].width)
#define AX_IOUWIN(iw,i)      OUW(i,AX_size[iw].width)
#define AX_JINWIN(iw,j)      INW(j,AX_size[iw].height)
#define AX_JOUWIN(iw,j)      OUW(j,AX_size[iw].height)
#define AX_IJINWIN(iw,i,j)  (AX_IINWIN(iw,i) && AX_JINWIN(iw,j))
#define AX_IJOUWIN(iw,i,j)  (AX_IOUWIN(iw,i) || AX_JOUWIN(iw,j))

/* The xy<->ij conversion should be context-dependent: */
/* here I just give one that does not cause seg-fault. */
#define AX_x_to_i(iw,x,i)  { (i)=(x); if ((i)==AX_size[iw].width)  (i)--; }
#define AX_y_to_j(iw,y,j)  { (j)=(y); if ((j)==AX_size[iw].height) (j)--; }
#define AX_xy_to_ij(iw,x,y,i,j)  { AX_x_to_i(iw,x,i); AX_y_to_j(iw,y,j) }

/* pixel is represented by its center */
#define AX_i_to_x(i,x)  ( (x) = (i) + 0.5 )
#define AX_j_to_y(j,y)  ( (y) = (j) + 0.5 )
#define AX_ij_to_xy(i,j,x,y)  ( AX_i_to_x(i,x), AX_j_to_y(j,y) )

/***************************/
/* Scan Conversion Manager */
/***************************/

typedef struct
{
    AX_IJ      i;    /* [0,width-1]  */
    AX_IJ      j;    /* [0,height-1] */
} AX_PAIR;

typedef union
{
    AX_PAIR    p;
    AX_Offset  offset;
} AX_BP;  /* aligned */

typedef union
{
    AX_Float   area;
    AX_Offset  offset;
} AX_Content;  /* aligned */

typedef struct
{
    AX_BP       b;
    AX_Content  c;
} AX_AP;  /* aligned */

/* pixel-based alias operation stack */
extern AX_AP *AX_aop [AX_MAXWIN];
/* alias stack memory limit */
extern int AX_maxaop [AX_MAXWIN];

/* maximal number of alias stack objects */
#define AX_AOP_MAXOBJ       1
/* alias stack against the amount of window pixels */
#define AX_AOP_MAXRATIO    (0.125)
/* margin to catch alias stack overflow */
#define AX_AOP_WARNINGZONE  128
#define AX_AOP_top(iw)      AX_aop[iw]
#define AX_AOP_TOP(iw,i)    AX_aop[iw][i]

/* local aop operations */
#define AX_clearaop()  ( aop[0].b.offset = 1 )
#define AX_Printaop(aop,m) { printf ("AOP list:\n"); \
  for (m=1; m<(aop)[0].b.offset; m++) printf("j=%d i=%d area=%f\n", \
  (aop)[m].b.p.j, (aop)[m].b.p.i, (aop)[m].c.area); }
#define AX_printaop(m) AX_Printaop(aop,m)
#define AX_addaop(ia,ja,AREA) \
  ( aop[aop[0].b.offset].b.p.i = (ia), aop[aop[0].b.offset].b.p.j = (ja), \
  aop[aop[0].b.offset++].c.area = (AREA) )

/* pixel-based block operation stack */
extern AX_BP *AX_bop [AX_MAXWIN];
/* block stack space limit */
extern int AX_maxbop [AX_MAXWIN];

/* maximal number of block stack objects */
#define AX_BOP_MAXOBJ       AX_AOP_MAXOBJ
/* block stack against the amount of window pixels */
#define AX_BOP_MAXRATIO    (1.0)
/* margin to catch block stack overflow */
#define AX_BOP_WARNINGZONE  0
#define AX_BOP_top(iw)      AX_bop[iw]
#define AX_BOP_TOP(iw,i)    AX_bop[iw][i]

/* local bop operations */
#define AX_clearbop()  (bop[0].offset = 1)
#define AX_Printbop(bop,m) { printf ("BOP list:\n"); \
  for (m=1; m<(bop)[0].offset; m++) \
  printf("j=%d i=%d\n", (bop)[m].p.j, (bop)[m].p.i); }
#define AX_printbop(m) AX_Printbop(bop,m)
#define AX_addbop(ia,ja) \
  ( bop[bop[0].offset].p.i = (ia), bop[bop[0].offset++].p.j = (ja) )

/* Allocate scan conversion stacks AX_aop[iw] and AX_bop[iw], */
/* sizes of which are based on an estimate using AX_size[iw]. */
int AX_plugin_Scan_module (int iw);

/* resize scan conversion stacks AX_aop[iw] and AX_bop[iw] */
/* according to current AX_size[iw] using realloc().       */
int AX_resize_Scan_module (int iw);

/* Free scan conversion scratch stack AX_aop[iw] and AX_bop[iw] */
int AX_plugout_Scan_module (int iw);


/**********************/
/* Clipping Functions */
/**********************/

/* Test if line segment (x1,y1)-(x2,y2) intersects with */
/* (x3,y3)-(x4,y4); if it does, save intersection x, y. */
#define AX_LINES_PARALLEL                    -1
#define AX_LINES_INTERSECT_BUT_NOT_SEGMENTS   0
#define AX_SEGMENTS_INTERSECT                 1
int AX_SegmentsIntersect ( AX_Float x1, AX_Float y1,
                           AX_Float x2, AX_Float y2,
                           AX_Float x3, AX_Float y3,
                           AX_Float x4, AX_Float y4,
                           AX_Float *x, AX_Float *y );

/* The Cohen-Sutherland line clipping algorithm */
bool AX_CS_LineClipWindow ( AX_Float *x0, AX_Float *y0,
                            AX_Float *x1, AX_Float *y1,
                            AX_Float width, AX_Float height );
bool AX_CS_LineClipRectangle ( AX_Float *x0, AX_Float *y0,
                               AX_Float *x1, AX_Float *y1,
                               AX_Float xmin, AX_Float ymin,
                               AX_Float xmax, AX_Float ymax );

/* Mark S. Sobkow, Paul Pospisil and Yee-Hong Yang, Computers & */
/* Graphics Vol. 11, No. 4, pp. 459-467, 1987, adopted by John  */
/* Schultz in Notes/clip2d.readme, clip2d.mod.                  */
bool AX_SPY_LineClipWindow ( AX_Float *x0, AX_Float *y0,
                             AX_Float *x1, AX_Float *y1,
                             AX_Float width, AX_Float height );
bool AX_SPY_LineClipRectangle ( AX_Float *x0, AX_Float *y0,
                                AX_Float *x1, AX_Float *y1,
                                AX_Float xmin, AX_Float ymin,
                                AX_Float xmax, AX_Float ymax );

#define AX_LineClipWindow(x0,y0,x1,y1,width,height) \
  AX_SPY_LineClipWindow(x0,y0,x1,y1,width,height)
#define AX_LineClipRectangle(x0,y0,x1,y1,xmin,ymin,xmax,ymax) \
  AX_SPY_LineClipRectangle(x0,y0,x1,y1,xmin,ymin,xmax,ymax)

typedef struct
{
    AX_Float x;
    AX_Float y;
} AX_Vertex;

#define AX_VertexAssign(v,X,Y)  ( (v).x=(X), (v).y=(Y) )
#define AX_RandomVertex(v,width,height) \
  ( (v).x=Frandom()*(width), (v).y=Frandom()*(height) )
#define AX_SameVertex(v0,v1) ( EQ((v0).x,(v1).x) && EQ((v0).y,(v1).y) )
#define AX_VertexShift(v,dx,dy) ( (v).x+=(dx), (v).y+=(dy) )
#define AX_CopyVertex(v,vc)  ( (vc).x = (v).x, (vc).y = (v).y )

typedef struct
{
    AX_Float outNormal_x;  /* unnormalized normal direction */
    AX_Float outNormal_y;
    AX_Float outNormal_self;
} AX_ClipBoundary;
/* outward "distance" (unnormalized) from ClipBoundary */
#define AX_VertexOutwardDistance(Vertex,ClipBoundary) \
  ((Vertex).x*(ClipBoundary).outNormal_x+(Vertex).y*(ClipBoundary).outNormal_y)

/* determine if Vertex is inside (left) of an infinite ClipBoundary */
#define AX_VertexInsideClipBoundary(Vertex,ClipBoundary) \
 (AX_VertexOutwardDistance(Vertex,ClipBoundary)<=(ClipBoundary).outNormal_self)
/* (AX_VertexOutwardDistance(Vertex,ClipBoundary) < \ */
/* (ClipBoundary).outNormal_self + AX_TINY) */

/* Convert directed edge v0->v1 to an infinite clipping boundary */
void AX_EdgeToClipBoundary (AX_Vertex *v0, AX_Vertex *v1, AX_ClipBoundary *c);

/* Calculate the intersection of two clip boundaries */
void AX_ClipIntersectClip
(AX_ClipBoundary *c0, AX_ClipBoundary *c1, AX_Vertex *intersection);

/* Given that "v0" is inside ClipBoundary "c" with "v1" */
/* outside, or vice versa, get the intersection point.  */
void AX_IntersectClipBoundary
(AX_Vertex *v0, AX_Vertex *v1, AX_ClipBoundary *c, AX_Vertex *intersection);

#define AX_PolygonMaxVertex 8
typedef struct
{
    int nVertex;
    AX_Vertex Vertex [AX_PolygonMaxVertex];
} AX_Polygon;

/* Allocate memory for a polygon (should be freed */
/* later!) and initialize its vertex coordinates. */
AX_Polygon *AX_NewPolygon (int nVertex, ...);

/* (Re)assign polygon vertices */
AX_Polygon *AX_PolygonAssign (AX_Polygon *p, int nVertex, ...);
AX_Polygon *AX_PolygonASSIGN (AX_Polygon *p, int nVertex, XPoint *V);
#define AX_PolygonToXPoints(p,V,i)  for (i=0; i<(p).nVertex; i++) { \
  V[i].x = (p).Vertex[i].x; V[i].y = (p).Vertex[i].y; }

/* Add a vertex to polygon */
#define AX_PolygonAddVertex(v,p) ((p).Vertex[(p).nVertex].x = (v).x, \
  (p).Vertex[(p).nVertex].y = (v).y, (p).nVertex++)

/* Copy polygon p to pc */
#define AX_PolygonCopy(p,pc,i) CAN ( (pc).nVertex = (p).nVertex; \
  for (i=0; i<(p).nVertex; i++) (pc).Vertex[i] = (p).Vertex[i]; )

/* Canonical definition of polygon in Linear Programming: Ax >= b */
typedef struct
{
    int nClipBoundary;
    AX_ClipBoundary ClipBoundary [AX_PolygonMaxVertex];
} AX_Clipper;

/* from polygon vertices, get the set of clip boundaries */
void AX_PolygonToClipper (AX_Polygon *p, AX_Clipper *c);

/* Vertex indices (0,1) corresponds to edge 0 */
#define AX_VertexSourceToEdge(Vertex_i,Edge_j) ( (Vertex_i) == (Edge_j) )
#define AX_VertexSinkToEdge(Vertex_i,Edge_j,nVertex) \
  ( (Vertex_i) == ((Edge_j) + 1) % (nVertex) )
#define AX_VertexIncidentToEdge(Vertex_i,Edge_j,nVertex) \
  ( AX_VertexSourceToEdge(Vertex_i,Edge_j) || \
  AX_VertexSinkToEdge(Vertex_i,Edge_j,nVertex) )

/* from a set of clip boundaries, get polygon vertices */
void AX_ClipperToPolygon (AX_Clipper *c, AX_Polygon *p);

/* Allocate memory for a Clipper (should be freed */
/* later!) and initialize by vertex coordinates.  */
AX_Clipper *AX_NewClipper (int nVertex, ...);

/* (Re)assign Clipper by vertices */
AX_Clipper *AX_ClipperAssign (AX_Clipper *c, int nVertex, ...);
AX_Clipper *AX_ClipperASSIGN (AX_Clipper *c, int nVertex, XPoint *V);

/* Test whether a vertex is inside Clipper */
bool AX_InsideClipper (AX_Vertex *v, AX_Clipper *c);
bool AX_INSIDEClipper (XPoint *V, AX_Clipper *c);

/* Check if a polygon Clipper region is "feasible" */
bool AX_feasibleClipper (AX_Clipper *c);

/* Sutherland-Hodgman polygon-polygon clipping algorithm */
void AX_SH_PolygonPolygonClip
(AX_Polygon *p, AX_Clipper *c, AX_Polygon *result);
void AX_SH_PolygonWindowClip
(AX_Polygon *p, AX_Float width, AX_Float height, AX_Polygon *result);

/* Check if segment (xmin,y)-(xmax,y) has finite */
/* length intersection with circle (x0,y0,r).    */
bool AX_HorizontalSegmentIntersectCircle
(AX_Float y, AX_Float xmin, AX_Float xmax,
 AX_Float x0, AX_Float y0, AX_Float radius);

/* Check if segment (x,ymin)-(x,ymax) has finite */
/* length intersection with circle (x0,y0,r).    */
bool AX_VerticalSegmentIntersectCircle
(AX_Float x, AX_Float ymin, AX_Float ymax,
 AX_Float x0, AX_Float y0, AX_Float radius);

#define AX_CircleIntersectRectangleBorder(x0,y0,radius,xmin,ymin,xmax,ymax) \
  AX_HorizontalSegmentIntersectCircle(ymin,xmin,xmax,x0,y0,radius) || \
  AX_HorizontalSegmentIntersectCircle(ymax,xmin,xmax,x0,y0,radius) || \
  AX_VerticalSegmentIntersectCircle(xmin,ymin,ymax,x0,y0,radius)   || \
  AX_VerticalSegmentIntersectCircle(xmax,ymin,ymax,x0,y0,radius)
/* Determine if circle (x0,y0,r) has finite area */
/* in rectangle (xmin,ymin) - (xmax,ymax).       */
#define AX_CircleIntersectRectangle(x0,y0,radius,xmin,ymin,xmax,ymax) \
  ( ( XIN(x0,xmin,xmax) && XIN(y0,ymin,ymax) ) || \
  AX_CircleIntersectRectangleBorder(x0,y0,radius,xmin,ymin,xmax,ymax) )
#define AX_CircleIntersectWindowBorder(x0,y0,radius,width,height) \
  AX_CircleIntersectRectangleBorder(x0,y0,radius,0,0,width,height)
#define AX_CircleIntersectWindow(x0,y0,radius,width,height) \
  AX_CircleIntersectRectangle(x0,y0,radius,0,0,width,height)

/* Calculate length of line segment (xmin,y)-(xmax,y) in circle (x0,y0,r) */
AX_Float AX_HorizontalSegmentInCircle
(AX_Float y, AX_Float xmin, AX_Float xmax,
 AX_Float x0, AX_Float y0, AX_Float radius2);

/* Calculate length of line segment (x,ymin)-(x,ymax) in circle (x0,y0,r) */
AX_Float AX_VerticalSegmentInCircle
(AX_Float x, AX_Float ymin, AX_Float ymax,
 AX_Float x0, AX_Float y0, AX_Float radius2);

/* If segment (*xmin,y)-(*xmax,y) intersects circle (x0,y0,r), */
/* revise *xmin,*xmax and return TRUE; otherwise return FALSE. */
bool AX_HorizontalSegmentClipCircle
(AX_Float y, AX_Float *xmin, AX_Float *xmax,
 AX_Float x0, AX_Float y0, AX_Float radius2);

/* If segment (x,*ymin)-(x,*ymax) intersects circle (x0,y0,r), */
/* revise *ymin,*ymax and return TRUE; otherwise return FALSE. */
bool AX_VerticalSegmentClipCircle
(AX_Float x, AX_Float *ymin, AX_Float *ymax,
 AX_Float x0, AX_Float y0, AX_Float radius2);

typedef struct
{
    AX_Float x0, y0;  /* ellipse center */
    AX_Float ax, ay;  /* normalized major axis direction */
    AX_Float a, a2;  /* length and squared length of major axis */
    AX_Float b, b2;  /* length and squared length of minor axis */
    AX_Float Dxx,Dxy,Dyy,D0;  /* Dxx x^2 + 2Dxy xy + Dyy y^2 - D0 = 0 */
    AX_Float ymin,ymax,ymin_x,ymax_x;  /* y-container */
} AX_Ellipse; /* aligned */

/* Specify an ellipse by major axis displacements and minor axis length */
AX_Ellipse *AX_EllipseAssign (AX_Float x0, AX_Float y0, AX_Float ax,
                              AX_Float ay, AX_Float b, AX_Ellipse *e);

/* Specify an ellipse by major axis angle (in radians) */
/* with x-axis and major, minor axis lengths.          */
AX_Ellipse *AX_EllipseASSIGN (AX_Float x0, AX_Float y0, AX_Float angle,
                              AX_Float a, AX_Float b, AX_Ellipse *e);
/* uses degree instead of radian */
#define AX_ELLIPSEASSIGN(x0,y0,angle,a,b,e) \
  AX_EllipseASSIGN(x0,y0,DEGREE_TO_RADIAN(angle),a,b,e)

/* Specify an ellipse by A*(x-x0)^2+2*B*(x-x0)*(y-y0)+C*(y-y0)^2=1 */
AX_Ellipse *AX_Ellipseassign (AX_Float x0, AX_Float y0, AX_Float A,
                              AX_Float B, AX_Float C, AX_Ellipse *e);

/* Check whether segment (xmin,y)-(xmax,y) has */
/* finite length intersection with ellipse.    */
bool AX_HorizontalSegmentIntersectEllipse
(AX_Float y, AX_Float xmin, AX_Float xmax, AX_Ellipse *e);

/* Calculate length of line segment (xmin,y)-(xmax,y) in ellipse */
AX_Float AX_HorizontalSegmentInEllipse
(AX_Float y, AX_Float xmin, AX_Float xmax, AX_Ellipse *e);

/* If segment (*xmin,y)-(*xmax,y) intersects ellipse, revise */
/* *xmin, *xmax and return TRUE; otherwise return FALSE.     */
bool AX_HorizontalSegmentClipEllipse
(AX_Float y, AX_Float *xmin, AX_Float *xmax, AX_Ellipse *e);

/* Check whether segment (x,ymin)-(x,ymax) has */
/* finite length intersection with ellipse.    */
bool AX_VerticalSegmentIntersectEllipse
(AX_Float x, AX_Float ymin, AX_Float ymax, AX_Ellipse *e);

/* Calculate length of line segment (x,ymin)-(x,ymax) in ellipse */
AX_Float AX_VerticalSegmentInEllipse
(AX_Float x, AX_Float ymin, AX_Float ymax, AX_Ellipse *e);

/* If segment (x,*ymin)-(x,*ymax) intersects ellipse, revise */
/* *ymin, *ymax and return TRUE; otherwise return FALSE.     */
bool AX_VerticalSegmentClipEllipse
(AX_Float x, AX_Float *ymin, AX_Float *ymax, AX_Ellipse *e);

/* determine if ellipse has finite area in rectangle (xmin,ymin)-(xmax,ymax) */
#define AX_EllipseIntersectRectangle(e,xmin,ymin,xmax,ymax) \
  ( ( XIN((e)->x0,xmin,xmax) && XIN((e)->y0,ymin,ymax) ) || \
  AX_HorizontalSegmentIntersectEllipse(ymin,xmin,xmax,e) || \
  AX_HorizontalSegmentIntersectEllipse(ymax,xmin,xmax,e) || \
  AX_VerticalSegmentIntersectEllipse(xmin,ymin,ymax,e)   || \
  AX_VerticalSegmentIntersectEllipse(xmax,ymin,ymax,e) )

#define AX_EllipseIntersectWindow(e,width,height) \
  AX_EllipseIntersectRectangle(e,0,0,width,height)


/***************/
/* Color Field */
/***************/

/* Color field with an origin (& its color) and linear gradients */
typedef struct
{
    AX_Float r0; AX_Float rx; AX_Float ry;
    AX_Float g0; AX_Float gx; AX_Float gy;
    AX_Float b0; AX_Float bx; AX_Float by;
} AX_GradientColor;

/* Get r,g,b brightness in floating-point [0,1] */
#define AX_r(c) (((AX_Float)((c)&AX_rmask))/AX_rmask)
#define AX_g(c) (((AX_Float)((c)&AX_gmask))/AX_gmask)
#define AX_b(c) (((AX_Float)((c)&AX_bmask))/AX_bmask)

#define AX_extrapolate(a0,ax,ay,dx,dy) ((a0)+(ax)*(dx)+(ay)*(dy))

#define AX_colorfield(r0,rx,ry,g0,gx,gy,b0,bx,by,dx,dy) AX_Colorcarrier \
  ( AX_extrapolate(r0,rx,ry,dx,dy), AX_extrapolate(g0,gx,gy,dx,dy), \
  AX_extrapolate(b0,bx,by,dx,dy) )

#define AX_gradientcolor(c,x,y) AX_colorfield \
  ( (c).r0,(c).rx,(c).ry, (c).g0,(c).gx,(c).gy, (c).b0,(c).bx,(c).by, x,y )

/* Assert [0,1] limit on color field brightnesses */
AX_Carrier AX_safe_gradientcolor (AX_GradientColor *c, AX_Float x, AX_Float y);

/* initialize a color field with no gradients */
#define AX_ZeroGradientColor(c,r,g,b) \
 ((c).r0=(r),(c).g0=(g),(c).b0=(b),(c).rx=(c).ry=(c).gx=(c).gy=(c).bx=(c).by=0)

/* Deduce color field given two points: it assumes */
/* there is no variation in the normal direction.  */
void AX_LinearGradientColor
(AX_Float x0, AX_Float y0, AX_Float r0, AX_Float g0, AX_Float b0,
 AX_Float x1, AX_Float y1, AX_Float r1, AX_Float g1, AX_Float b1,
 AX_GradientColor *c);

/* deduce color field given three points */
void AX_TriangulateGradientColor
(AX_Float x0, AX_Float y0, AX_Float r0, AX_Float g0, AX_Float b0,
 AX_Float x1, AX_Float y1, AX_Float r1, AX_Float g1, AX_Float b1,
 AX_Float x2, AX_Float y2, AX_Float r2, AX_Float g2, AX_Float b2,
 AX_GradientColor *c);


/**********************/
/* Various Primitives */
/**********************/

/* read ScanModuleAPI.txt */

/* Interest Begins: interest begins at this pixel */
#define AX_IB(l,IB) ( IB=(l) ) /* IB,LP,RP,IE are all inclusive */
/* Left Peace: left foot does not disturb this pixel */
#define AX_LP(l,LP) { LP=(l)+1; if ((l)+1==LP) LP--; }
/* Right Peace: right foot does not disturb this pixel */
#define AX_RP(r,RP) ( RP=((AX_IJ)(r))-1 )
/* Interest Ends: interest ends with this pixel */
#define AX_IE(r,IE) { IE=(r); if ((r)==IE) IE--; }


/*****************/
/* straight line */
/*****************/

/* Simple line demo of the ScanModule API. The valid domain */
/* of x,y for AX_SimpleLineScan() is [0,height) x [0,width) */
void AX_SimpleLineScan
(AX_Float x0, AX_Float y0, AX_Float x1, AX_Float y1, AX_BP *bop);

/* Because line is a 1D instead of 2D object, its scan conversion is  */
/* not perfectly defined. The drawing board is [0,width) x [0,height) */
/* for AX_SimpleLine(), which is also a superset of XDrawLine().      */
void AX_SimpleLine (int iw, AX_Float x0, AX_Float y0,
                    AX_Float x1, AX_Float y1, AX_Carrier c);
void AX_SimpleLineGradientColor
(int iw, AX_Float x0, AX_Float y0, AX_Carrier c0,
 AX_Float x1, AX_Float y1, AX_Carrier c1);

/* Gupta-Sproull anti-aliased scan conversion of straight line */

/* From signed distance to anti-alias intensity: */
/* sharp cone filter for line of thickness 1.    */
/* #define AX_GS_intensity(d)  (1-ABS(d)/2.35) */
/* #define AX_GS_intensity(d)  (1-ABS(d)/2) */
#define AX_GS_intensity(d)  (1-ABS(d)/1.5)

/* Because it involves three pixels, the valid domain of */
/* x,y for AX_GSLineScan() is [1,height-1) x [1,width-1) */
void AX_GSLineScan(AX_Float x0,AX_Float y0,AX_Float x1,AX_Float y1,AX_AP *aop);

/* Because line is a 1D instead of 2D object, its scan conversion is  */
/* not perfectly defined. The drawing board is [1,width) x [1,height) */
void AX_GSLine (int iw, AX_Float x0, AX_Float y0,
                AX_Float x1, AX_Float y1, AX_Carrier c);
void AX_GSLineGradientColor
(int iw, AX_Float x0, AX_Float y0, AX_Carrier c0,
 AX_Float x1, AX_Float y1, AX_Carrier c1);


/***********/
/* polygon */
/***********/

/* Anti-alias intensity with area weight saturation */
#define AX_SAT(area) (((area)+0.65)/1.65)
/* #define AX_SAT(area) (area) */
#define AX_MIX(c0,c1,a1) ((c0==0)?AX_carrier(AX_SAT(a1)*((c1)&AX_rmask), \
  AX_SAT(a1)*((c1)&AX_gmask), AX_SAT(a1)*((c1)&AX_bmask)):AX_mix(c0,c1,a1))

/* increase area[] by those covered by trapezoidal slice 0 < Dy < 1 */
void AX_sliceConversion ( AX_Float Dy, AX_Float l, AX_Float lslope,
                          AX_Float r, AX_Float rslope, AX_Float *area );

/* process a thick trapezoidal slice that lies between two basecuts */
void AX_basecutsConversion
(int j0, int Dj, AX_Float l0, AX_Float lslope, AX_Float r0, AX_Float rslope,
 AX_AP *aop, AX_BP *bop);

typedef struct
{
    AX_Float y;
    AX_Float l;
    AX_Float r;
    int j;
    AX_Float lslope;
    AX_Float rslope;
    int is_basecut;
    int no_basecut; 
} AX_Cut;

#define AX_CIRCLE_TANDEM 3
#define AX_MAXCUT (AX_MAXHEIGHT+1+AX_CIRCLE_TANDEM)
typedef struct
{
    int ncut;
    int IB;
    int IE;  /* [IB,IE] delimits the width of the entire shape */
    AX_Cut cut [AX_MAXCUT];
} AX_Tandem;

/* Scan conversion of shapes parsed into tandem horizontal cuts */
void AX_TandemConversion (AX_Tandem *t0, AX_AP *aop, AX_BP *bop);

/* calculate the leftmost and rightmost intersections of a */
/* horizontal cut with the polygon, if they do intersect.  */
bool AX_PolygonHorizontalCut (AX_Polygon *p, AX_Float y, AX_Cut *c);

/* Scan convert polygons with exact anti-aliasing */
void AX_PolygonScan (AX_Polygon *p, AX_AP *aop, AX_BP *bop);

/* Draw filled polygon in single color "c" with exact anti-aliasing */
void AX_POLYGON (int iw, AX_Polygon *p, AX_Carrier c);

/* Draw filled triangle in single color "c" with exact anti-aliasing */
void AX_Triangle (int iw, AX_Float x0, AX_Float y0, AX_Float x1,
                  AX_Float y1, AX_Float x2, AX_Float y2, AX_Carrier c);

/* Draw rectangular bar with centerline (x0,y0)-(x1,y1) and thickness */
/* +-"barwidth"/2 in a single color "c" with exact anti-aliasing.     */
void AX_Bar (int iw, AX_Float x0, AX_Float y0, AX_Float x1,
             AX_Float y1, AX_Float barwidth, AX_Carrier c);

/* Draw gradient colorbar with exact anti-aliasing */
void AX_BarGradientColor
(int iw, AX_Float x0, AX_Float y0, AX_Float r0, AX_Float g0, AX_Float b0,
 AX_Float x1, AX_Float y1, AX_Float r1, AX_Float g1, AX_Float b1,
 AX_Float barwidth);

/* Draw filled color-triangle with exact anti-aliasing */
void AX_TriangleGradientColor
(int iw, AX_Float x0, AX_Float y0, AX_Float r0, AX_Float g0, AX_Float b0,
 AX_Float x1, AX_Float y1, AX_Float r1, AX_Float g1, AX_Float b1,
 AX_Float x2, AX_Float y2, AX_Float r2, AX_Float g2, AX_Float b2);

/**********/
/* circle */
/**********/

#define AX_CAPHEIGHT  2

/***********************************************************/
/* Scan convert circle clipped by window with almost exact */
/* anti-aliasing. The circle must have finite area in the  */
/* window (use AX_CircleIntersectWindow() beforehand).     */
/***********************************************************/
void AX_CircleWindowScan
(AX_Float x0, AX_Float y0, AX_Float radius, AX_Float width, AX_Float height,
 AX_AP *aop, AX_BP *bop);

/* Draw filled circle in single color "c" */
void AX_Circle
(int iw, AX_Float x0, AX_Float y0, AX_Float radius, AX_Carrier c);

/***********/
/* ellipse */
/***********/

/************************************************************/
/* Scan convert ellipse clipped by window with almost exact */
/* anti-aliasing. The ellipse must have finite area in the  */
/* window (use AX_EllipseIntersectWindow() beforehand).     */
/************************************************************/
void AX_EllipseWindowScan
(AX_Ellipse *e, AX_Float width, AX_Float height, AX_AP *aop, AX_BP *bop);

/* Draw filled ellipse in single color "c" */
void AX_ELLIPSE (int iw, AX_Ellipse *e, AX_Carrier c);


/* Sorter.c: */

/********************************************/
/* Fast sorting routines for AX_Float array */
/********************************************/

/* AX memory management */
#define AXFmem(N) ((AX_Float *)malloc((N)*sizeof(AX_Float)))
#define AXImem(N) ((int *)malloc((N)*sizeof(int)))
#define AXBmem(size_in_bits) \
  ((Bmap *)AXImem(SEIL(BITS_TO_BYTES(size_in_bits),sizeof(int))))
#define AXFREE(x) free((void *)(x))

/* Assign index array idx[] := a..b (1..10, 10..1); return idx[] */
int *AX_sequentially_index (int idx[], int a, int b);

/* assign index array idx[] := 0..N-1 */
void AX_Sequentially_index (int N, int idx[]);

/* randomly permute index array a[min..max] (e.g. a[0..n-1]) */
int *AX_r_permute (int *a, int min, int max);

/* determine if values in idx[] is a permutation of min..max */
int AX_is_permutation (int idx[], int min, int max);

/* determine whether x[idx[i]], i=0..N-1, is non-decreasing */
int AX_is_nondecreasing (int N, AX_Float x[], int idx[]);

/* rearrange array of AX_Floats: y[i] := x[idx[i]], i=0..N-1 */
void AX_Vrearrange (int N, AX_Float x[], int idx[], AX_Float y[]);

/* rearrange array of AX_Floats: x[] := x[idx[]], then set idx[] := 0..N-1 */
void AX_VRearrange (int N, AX_Float x[], int idx[]);

#define AX_USE_OLD_IDX  0  /* incoming index permutation but not sequential */
#define AX_USE_NEW_IDX  1  /* incoming index undefined: must initialize */
#define AX_USE_SEQ_IDX  2  /* incoming index sequential */

/* non-recursive quicksort category:  */
void AX_qsort_glibc (int N,  AX_Float x[], int idx[], int option);
void AX_qsort_numerical_recipes (int N,  AX_Float x[], int idx[], int option);

/* recursive mergesort category: */
void AX_msort_Lee (int N,  AX_Float x[], int idx[], int option);
/* non-recursive mergesort category: */
void AX_msort_Ju (int N,  AX_Float x[], int idx[], int option);


/* 3D.c: */

/**************/
/* 3D Support */
/**************/

#define AX_V3ZERO(a)  ((a)[0] = (a)[1] = (a)[2] = 0)
/* generate three independent random components on (-0.5,0.5) */
#define AX_V3FRANDOM(a) ((a)[0]=FRANDOM(),(a)[1]=FRANDOM(),(a)[2]=FRANDOM())

#define AX_V3DOT(a,b)  ((a)[0]*(b)[0] + (a)[1]*(b)[1] + (a)[2]*(b)[2])
#define AX_V3LENGTH2(a)  AX_V3DOT(a,a)
#define AX_V3LENGTH(a)  sqrt((double)AX_V3LENGTH2(a))
#define AX_V3MUL(multiplier,a)  \
  ((a)[0]*=(multiplier),(a)[1]*=(multiplier),(a)[2]*=(multiplier))
#define AX_V3DIV(a,divisor)  \
  ((a)[0]/=(divisor),(a)[1]/=(divisor),(a)[2]/=(divisor))
#define AX_V3div(a,divisor,b)  \
  ((b)[0]=(a)[0]/(divisor),(b)[1]=(a)[1]/(divisor),(b)[2]=(a)[2]/(divisor))
#define AX_V3NORMALIZE(a,r) (r=AX_V3LENGTH(a),AX_V3DIV(a,r))
#define AX_V3normalize(a,b,r) (r=AX_V3LENGTH(a),AX_V3div(a,r,b))

#define AX_V3EQV(a,b) ((b)[0]=(a)[0],(b)[1]=(a)[1],(b)[2]=(a)[2])
#define AX_V3NEG(a,b) ((b)[0]=-(a)[0],(b)[1]=-(a)[1],(b)[2]=-(a)[2])

#define AX_V3EQZERO(a) ( (0 == (a)[0]) && (0 == (a)[1]) && (0 == (a)[2]) )
#define AX_V3NEZERO(a) ( (0 != (a)[0]) || (0 != (a)[1]) || (0 != (a)[2]) )

#define AX_V3TRIM(a,b) \
  ((b)[0]=TRIM((a)[0]),(b)[1]=TRIM((a)[1]),(b)[2]=TRIM((a)[2]))
#define AX_V3IMAGE(a,b) \
  ((b)[0]=IMAGE((a)[0]),(b)[1]=IMAGE((a)[1]),(b)[2]=IMAGE((a)[2]))

#define AX_V3add(a,b,c) ((c)[0] = (a)[0] + (b)[0], \
  (c)[1] = (a)[1] + (b)[1], (c)[2] = (a)[2] + (b)[2])
#define AX_V3ADD(a,b) ((a)[0]+=(b)[0], (a)[1]+=(b)[1], (a)[2]+=(b)[2])

#define AX_V3sub(a,b,c) ((c)[0] = (a)[0] - (b)[0], \
  (c)[1] = (a)[1] - (b)[1], (c)[2] = (a)[2] - (b)[2])
#define AX_V3SUB(a,b) ((a)[0]-=(b)[0], (a)[1]-=(b)[1], (a)[2]-=(b)[2])

#define AX_V3CROSS(a,b,c)  ((c)[0] = (a)[1]*(b)[2] - (a)[2]*(b)[1], \
  (c)[1] = (a)[2]*(b)[0]-(a)[0]*(b)[2], (c)[2] = (a)[0]*(b)[1]-(a)[1]*(b)[0])

/* c[] := a[] - multiplier * b[] */
#define AX_V3submul(a,multiplier,b,c) ( (c)[0]=(a)[0]-(multiplier)*(b)[0], \
  (c)[1] = (a)[1]-(multiplier)*(b)[1], (c)[2] = (a)[2]-(multiplier)*(b)[2] )
/* a[] := a[] - multiplier * b[] */
#define AX_V3SUBMUL(a,multiplier,b) ( (a)[0] -= (multiplier) * (b)[0], \
  (a)[1] -= (multiplier) * (b)[1], (a)[2] -= (multiplier) * (b)[2] )

/* c[] := a[] + multiplier * b[] */
#define AX_V3addmul(a,multiplier,b,c) ( (c)[0]=(a)[0]+(multiplier)*(b)[0], \
  (c)[1] = (a)[1]+(multiplier)*(b)[1], (c)[2] = (a)[2]+(multiplier)*(b)[2] )
/* a[] := a[] + multiplier * b[] */
#define AX_V3ADDMUL(a,multiplier,b) ( (a)[0] += (multiplier) * (b)[0], \
  (a)[1] += (multiplier) * (b)[1], (a)[2] += (multiplier) * (b)[2] )

/* c[] := a[] + b[] / divisor */
#define AX_V3adddiv(a,b,divisor,c) ( (c)[0] = (a)[0]+(b)[0]/(divisor), \
  (c)[1] = (a)[1]+(b)[1]/(divisor), (c)[2] = (a)[2]+(b)[2]/(divisor) )
/* a[] := a[] + b[] / divisor */
#define AX_V3ADDDIV(a,b,divisor) ( (a)[0] += (b)[0] / (divisor), \
  (a)[1] += (b)[1] / (divisor), (a)[2] += (b)[2] / (divisor) )

#define AX_M3EQV(A,B) (B[0][0] = A[0][0], B[0][1] = A[0][1], \
  B[0][2] = A[0][2], B[1][0] = A[1][0], B[1][1] = A[1][1], \
  B[1][2] = A[1][2], B[2][0] = A[2][0], B[2][1] = A[2][1], B[2][2] = A[2][2])

#define AX_M3ZERO(A) (A[0][0]=0, A[0][1]=0, A[0][2]=0, \
  A[1][0]=0, A[1][1]=0, A[1][2]=0, A[2][0]=0, A[2][1]=0, A[2][2]=0)

#define AX_M3IDENTITY(A) ( A[0][0]=1, A[0][1]=0, A[0][2]=0, \
  A[1][0]=0, A[1][1]=1, A[1][2]=0, A[2][0]=0, A[2][1]=0, A[2][2]=1 )

#define AX_M3DIAGONAL(a0,a1,a2,A) (A[0][0]=(a0), A[0][1]=0, A[0][2]=0, \
  A[1][0]=0, A[1][1]=(a1), A[1][2]=0, A[2][0]=0, A[2][1]=0, A[2][2]=(a2))

#define AX_M3FRANDOM(A) (A[0][0] = FRANDOM(), \
  A[0][1] = FRANDOM(), A[0][2] = FRANDOM(), A[1][0] = FRANDOM(), \
  A[1][1] = FRANDOM(), A[1][2] = FRANDOM(), A[2][0] = FRANDOM(), \
  A[2][1] = FRANDOM(), A[2][2] = FRANDOM())

#define AX_M3NEG(A,B) (B[0][0] = -A[0][0], \
  B[0][1] = -A[0][1], B[0][2] = -A[0][2], B[1][0] = -A[1][0], \
  B[1][1] = -A[1][1], B[1][2] = -A[1][2], B[2][0] = -A[2][0], \
  B[2][1] = -A[2][1], B[2][2] = -A[2][2])

#define AX_M3TRANSPOSE(A,B) ( \
  B[0][0] = A[0][0], B[0][1] = A[1][0], B[0][2] = A[2][0], \
  B[1][0] = A[0][1], B[1][1] = A[1][1], B[1][2] = A[2][1], \
  B[2][0] = A[0][2], B[2][1] = A[1][2], B[2][2] = A[2][2] )

#define AX_M3SYMMETRIZE(A,B) ( B[0][0] = A[0][0], \
  B[1][1] = A[1][1], B[2][2] = A[2][2], \
  B[0][1] = B[1][0] = (A[1][0]+A[0][1])/2., \
  B[0][2] = B[2][0] = (A[2][0]+A[0][2])/2., \
  B[1][2] = B[2][1] = (A[2][1]+A[1][2])/2. )
#define AX_M3Symmetrize(A) ( A[0][1] = A[1][0] = (A[1][0]+A[0][1])/2., \
  A[0][2] = A[2][0] = (A[2][0]+A[0][2])/2., \
  A[1][2] = A[2][1] = (A[2][1]+A[1][2])/2.)

#define AX_M3TR(A) (A[0][0]+A[1][1]+A[2][2])

#define AX_M3DETERMINANT(A) ( \
  A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1]) + \
  A[0][1]*(A[1][2]*A[2][0]-A[1][0]*A[2][2]) + \
  A[0][2]*(A[1][0]*A[2][1]-A[2][0]*A[1][1]) )

/* B[][] := A[][]^-1; determinant := det(A) */
#define AX_M3inv(A,B,determinant) ( \
  B[0][0] = A[1][1]*A[2][2]-A[1][2]*A[2][1], \
  B[1][1] = A[2][2]*A[0][0]-A[2][0]*A[0][2], \
  B[2][2] = A[0][0]*A[1][1]-A[0][1]*A[1][0], \
  B[1][0] = A[1][2]*A[2][0]-A[1][0]*A[2][2], \
  B[2][1] = A[2][0]*A[0][1]-A[2][1]*A[0][0], \
  B[0][2] = A[0][1]*A[1][2]-A[0][2]*A[1][1], \
  B[2][0] = A[1][0]*A[2][1]-A[2][0]*A[1][1], \
  B[0][1] = A[2][1]*A[0][2]-A[0][1]*A[2][2], \
  B[1][2] = A[0][2]*A[1][0]-A[1][2]*A[0][0], \
  determinant = A[0][0]*B[0][0] + A[0][1]*B[1][0] + A[0][2]*B[2][0], \
  B[0][0] /= determinant, B[1][1] /= determinant, B[2][2] /= determinant, \
  B[1][0] /= determinant, B[2][1] /= determinant, B[0][2] /= determinant, \
  B[2][0] /= determinant, B[0][1] /= determinant, B[1][2] /= determinant )

/* B[][] := A[][]^-T; determinant := det(A) */
#define AX_M3invtranspose(A,B,determinant) ( \
  B[0][0] = A[1][1]*A[2][2]-A[1][2]*A[2][1], \
  B[1][1] = A[2][2]*A[0][0]-A[2][0]*A[0][2], \
  B[2][2] = A[0][0]*A[1][1]-A[0][1]*A[1][0], \
  B[0][1] = A[1][2]*A[2][0]-A[1][0]*A[2][2], \
  B[1][2] = A[2][0]*A[0][1]-A[2][1]*A[0][0], \
  B[2][0] = A[0][1]*A[1][2]-A[0][2]*A[1][1], \
  B[0][2] = A[1][0]*A[2][1]-A[2][0]*A[1][1], \
  B[1][0] = A[2][1]*A[0][2]-A[0][1]*A[2][2], \
  B[2][1] = A[0][2]*A[1][0]-A[1][2]*A[0][0], \
  determinant = A[0][0]*B[0][0] + A[0][1]*B[0][1] + A[0][2]*B[0][2], \
  B[0][0] /= determinant, B[1][1] /= determinant, B[2][2] /= determinant, \
  B[1][0] /= determinant, B[2][1] /= determinant, B[0][2] /= determinant, \
  B[2][0] /= determinant, B[0][1] /= determinant, B[1][2] /= determinant )

#define AX_M3MULTIPLY(multiplier,A,B) ( B[0][0] = (multiplier) * A[0][0], \
  B[0][1] = (multiplier) * A[0][1], B[0][2] = (multiplier) * A[0][2], \
  B[1][0] = (multiplier) * A[1][0], B[1][1] = (multiplier) * A[1][1], \
  B[1][2] = (multiplier) * A[1][2], B[2][0] = (multiplier) * A[2][0], \
  B[2][1] = (multiplier) * A[2][1], B[2][2] = (multiplier) * A[2][2] )

#define AX_M3DIVIDE(A,divisor,B) ( B[0][0] = A[0][0] / (divisor), \
  B[0][1] = A[0][1] / (divisor), B[0][2] = A[0][2] / (divisor), \
  B[1][0] = A[1][0] / (divisor), B[1][1] = A[1][1] / (divisor), \
  B[1][2] = A[1][2] / (divisor), B[2][0] = A[2][0] / (divisor), \
  B[2][1] = A[2][1] / (divisor), B[2][2] = A[2][2] / (divisor) )

#define AX_M3SUBDIAG(A,a,B) ( B[0][0] = A[0][0] - (a), \
  B[0][1] = A[0][1], B[0][2] = A[0][2], B[1][0] = A[1][0], \
  B[1][1] = A[1][1] - (a), B[1][2] = A[1][2], B[2][0] = A[2][0], \
  B[2][1] = A[2][1], B[2][2] = A[2][2] - (a) )

/* c[] := A[][] * b[] */
#define AX_M3mulV3(A,b,c) ( (c)[0] = AX_V3DOT(A[0], b), \
  (c)[1] = AX_V3DOT(A[1], b), (c)[2] = AX_V3DOT(A[2], b) )

/* b[] := A[][] * b[] */
#define AX_M3MULV3(A,b,tmp) ( AX_M3mulV3(A,b,tmp), AX_V3EQV(tmp,b) )

#define AX_M3ADD(A,B,C) ( C[0][0] = A[0][0] + B[0][0], \
  C[0][1] = A[0][1] + B[0][1], C[0][2] = A[0][2] + B[0][2], \
  C[1][0] = A[1][0] + B[1][0], C[1][1] = A[1][1] + B[1][1], \
  C[1][2] = A[1][2] + B[1][2], C[2][0] = A[2][0] + B[2][0], \
  C[2][1] = A[2][1] + B[2][1], C[2][2] = A[2][2] + B[2][2] )

#define AX_M3SUB(A,B,C) ( C[0][0] = A[0][0] - B[0][0], \
  C[0][1] = A[0][1] - B[0][1], C[0][2] = A[0][2] - B[0][2], \
  C[1][0] = A[1][0] - B[1][0], C[1][1] = A[1][1] - B[1][1], \
  C[1][2] = A[1][2] - B[1][2], C[2][0] = A[2][0] - B[2][0], \
  C[2][1] = A[2][1] - B[2][1], C[2][2] = A[2][2] - B[2][2] )

#define AX_M3MUL(A,B,C) ( \
  C[0][0] = A[0][0]*B[0][0] + A[0][1]*B[1][0] + A[0][2]*B[2][0], \
  C[0][1] = A[0][0]*B[0][1] + A[0][1]*B[1][1] + A[0][2]*B[2][1], \
  C[0][2] = A[0][0]*B[0][2] + A[0][1]*B[1][2] + A[0][2]*B[2][2], \
  C[1][0] = A[1][0]*B[0][0] + A[1][1]*B[1][0] + A[1][2]*B[2][0], \
  C[1][1] = A[1][0]*B[0][1] + A[1][1]*B[1][1] + A[1][2]*B[2][1], \
  C[1][2] = A[1][0]*B[0][2] + A[1][1]*B[1][2] + A[1][2]*B[2][2], \
  C[2][0] = A[2][0]*B[0][0] + A[2][1]*B[1][0] + A[2][2]*B[2][0], \
  C[2][1] = A[2][0]*B[0][1] + A[2][1]*B[1][1] + A[2][2]*B[2][1], \
  C[2][2] = A[2][0]*B[0][2] + A[2][1]*B[1][2] + A[2][2]*B[2][2] )

/* this is cubic random, NOT spherical random */
#define AX_M3RANDOMROTATION(R,tmp) ( AX_V3FRANDOM(R[0]), \
  AX_V3NORMALIZE(R[0],tmp), AX_V3FRANDOM(R[1]), AX_V3CROSS(R[0],R[1],R[2]), \
  AX_V3NORMALIZE(R[2],tmp), AX_V3CROSS(R[2],R[0],R[1]) )

/* white color spot light */
#define AX_3D_LIGHT_R (1.)
#define AX_3D_LIGHT_G (1.)
#define AX_3D_LIGHT_B (1.)
/* z_angle = 20, x_angle = 55 degrees */
#define AX_3D_LIGHT_X (0.53898554469576)
#define AX_3D_LIGHT_Y (-0.76975113132006)
#define AX_3D_LIGHT_Z (-0.34202014332567)
/* light influence tries to complement the color */
/* difference between light itself and object.   */
#define AX_3D_LRGB(r,g,b,lr,lg,lb) \
  (lr=AX_3D_LIGHT_R-(r), lg=AX_3D_LIGHT_G-(g), lb=AX_3D_LIGHT_B-(b))
/* U is the outgoing surface normal */
#define AX_3D_LDOT(U) \
  (AX_3D_LIGHT_X*U[0] + AX_3D_LIGHT_Y*U[1] + AX_3D_LIGHT_Z*U[2])

/***************/
/* 3D Renderer */
/***************/

typedef struct
{
    double d0;      /* abscissa */
    double dx[3];   /* normal */
} AX_3D_Filter_Plane;

#define AX_3D_MAX_FILTER_PLANE 16

typedef struct
{
    double x[3];  /* coordinates of the viewpoint */
    double k;  /* conversion factor from radian to window pixels */
    double V[3][3];  /* V[i][0-2] is the normalized ith axis of viewport */
    double zcut;   /* (0,zcut] of the view frustum */
    double wx,wy;  /* viewpoint coordinates in window frame (pixels) */
    AX_3D_Filter_Plane fp[AX_3D_MAX_FILTER_PLANE]; /* filter planes */
} AX_3D_Viewport;  

extern AX_3D_Viewport AX_3D [AX_MAXWIN];

#define AX_3D_DEF_VIEW_ANGLE_IN_DEGREE (60)

/* Assign default observer (viewport) AX_3D[iw] */
void AX_3D_DEFViewport (int iw);

/* floating point z-buffer */
extern AX_Float *AX_zmap[AX_MAXWIN];

#ifdef AX_MEMORY_EFFICIENT
#define AX_clearzmap(iw) \
  FLOATset(AX_zmap[iw], AX_INFINITY, AX_size[iw].width*AX_size[iw].height)
#else
#define AX_clearzmap(iw) \
  DOUBLEset(AX_zmap[iw], AX_INFINITY, AX_size[iw].width*AX_size[iw].height)
#endif

/* Allocate AX_zmap[iw] according to AX_size[iw] and set viewport defaults */
int AX_plugin_3D_module (int iw);
#define AX_PLUGIN_3D_module(iw) (AX_plugin_3D_module(iw), AX_clearzmap(iw))

/* Free z-buffer memory allocation AX_zmap[iw] */
int AX_plugout_3D_module (int iw);

/* resize AX_zmap[iw] according to current AX_size[iw] by */
/* realloc(), and reassign viewport parameters sensibly.  */
int AX_resize_3D_module (int iw);

typedef union
{
    AX_Float   area;
    AX_Offset  carrier;
} AX_content;  /* aligned */

/* ready-to-deploy pixel */
typedef struct
{
    AX_Offset   offset;
    AX_content  c;
    AX_Float    z;
} AX_3D_Pixel;

#define AX_3D_AssignRGB(P,R,G,B) ((P).r=(R),(P).g=(G),(P).b=(B))

/****************/
/* Line Support */
/****************/

typedef struct
{
    AX_Float x0[3], x1[3];  /* 3D coordinates */
    AX_Float X0[2], X1[2];  /* Viewport coordinates */
    AX_Float r, g, b;       /* line color */
} AX_3D_Line;

typedef struct
{
    int n_lines;
    AX_3D_Line *LINE;
} AX_3D_Lines;

/* Realloc n_lines; if used as malloc, assign L->LINE=NULL first */
void AX_3D_Lines_Realloc (AX_3D_Lines *L, int n_lines);

/* Realloc n_lines followed by x0,y0,z0,x1,y1,z1,r,g,b,... interface */
void AX_3D_Lines_Recreate (AX_3D_Lines *L, int n_lines, ...);

/* free allocated memory and set NULL */
void AX_3D_Lines_Free (AX_3D_Lines *L);

/* Generate line scan and draw to pixmap alias pixels according to z-buffer */
void AX_3D_Lines_Zprint (int iw, AX_3D_Lines *L);

/* define stack memory for a Parallelepiped Wireframe */
#define AX_3D_Define_PW(L) \
  AX_3D_Line AX_3D_PW_LINE[12]; AX_3D_Lines L[1]; \
  L[0].n_lines=12; L[0].LINE=AX_3D_PW_LINE;

/* Assign parallelepiped wireframe H[][] originated at x0,y0,z0 */
AX_3D_Lines *AX_3D_PW_Assign
(AX_Float x0, AX_Float y0, AX_Float z0, AX_Float H[3][3],
 AX_Float r, AX_Float g, AX_Float b, AX_3D_Lines *L);


/*******************/
/* Polygon Support */
/*******************/

/***********/
/* Shading */
/***********/
#define AX_PS_EDGE  0.65
#define AX_PS_TOPI  0.35
#define AX_PS_TOP  (AX_PS_EDGE+AX_PS_TOPI)
#define AX_PS_LIGHT 0.85
/* back surface */
#define AX_SP_EDGE  0.
#define AX_SP_TOPI  0.4
#define AX_SP_TOP  (AX_SP_EDGE+AX_SP_TOPI)
#define AX_SP_LIGHT 0.1

typedef union
{
    int offset;
    struct AX_3D_PIXEL *p;
} AX_OffsetP; /* aligned */

struct AX_3D_PIXEL
{
    AX_OffsetP  o;
    AX_Float    area;
    AX_Offset   carrier;
    AX_Float    z;
};

/* this represents a plane of the same material (r,g,b) but */
/* different shading when looking at the front surface (PS) */
/* (vertices 0,1,2..n-1 forming a right-handed system with  */
/* respect to the front normal) than the back surface (SP). */
#define AX_3D_PolygonMaxVertex 5

typedef struct
{ /* cut in two corners by two parallel planes */
    int nVertex;
    AX_Float x[3*(AX_3D_PolygonMaxVertex+2)];
} AX_3D_polygon;

typedef struct
{
    int nVertex;
    AX_Float x[3*AX_3D_PolygonMaxVertex];
    AX_Float front_normal[3];
    AX_3D_polygon q[1];  /* polygon in viewport frame */
    AX_Polygon p[1];  /* projected polygon on screen */
    AX_Float r,g,b;
    struct AX_3D_PIXEL *IMG;  /* for the alias pixels */
    AX_3D_Pixel *img;  /* for the block pixels */
} AX_3D_Polygon;  /* aligned */

/* (Re)assign 3D polygon by x0,y0,z0, x1,y1,z1, ... interface */
AX_3D_Polygon *AX_3D_PolygonAssign
(AX_3D_Polygon *p, int nVertex, AX_Float r, AX_Float g, AX_Float b, ...);

/* (Re)assign 3D polygon by *x0, *x1, ... interface */
AX_3D_Polygon *AX_3D_PolygonASSIGN
(AX_3D_Polygon *p, int nVertex, AX_Float r, AX_Float g, AX_Float b, ...);

#define AX_3D_PolygonCopy(p,pc,i) { (pc)->nVertex = (p)->nVertex; \
  for (i=0; i<3*(p)->nVertex; i++) (pc)->x[i] = (p)->x[i]; }

#define AX_3D_Polygontranslate(p,dx,dy,dz,pc,i) { \
  for (i=0; i<(p).nVertex; i++) { (pc).x[3*i]=(p).x[3*i]+dx; \
  (pc).x[3*i+1]=(p).x[3*i+1]+dy; (pc).x[3*i+2]=(p).x[3*i+2]+dz; } }

#define AX_3D_PolygonTranslate(p,dx,pc,i) \
  AX_3D_Polygontranslate(p,(dx)[0],(dx)[1],(dx)[2],pc,i)

#define AX_3D_PolygonAddVertex(v,P) \
  ( AX_V3EQV(v,&((P).x[3*(P).nVertex])), (P).nVertex++ )

/* Sutherland-Hodgman 3D polygon/plane (n*x>=n[3]) clipping algorithm */
void AX_3D_SH_PolygonClipPlane
(AX_3D_polygon *p, AX_Float n[4], AX_3D_polygon *result);

typedef struct
{
    int n_polygons;
    AX_3D_Polygon *POLYGON;
    int n_shown;
    /* z-order to preempt polygons that are "visible" but covered */
    AX_Float *vz;
    int *idx;
} AX_3D_Polygons;

/* Realloc n_polygons; if used as malloc, let P->POLYGON=P->idx=NULL first */
void AX_3D_Polygons_Realloc (AX_3D_Polygons *P, int n_polygons);

/* free allocated memory and set NULL */
void AX_3D_Polygons_Free (AX_3D_Polygons *P);

/* generate polygon scans and register them in the z-buffer */
void AX_3D_Polygons_Zpaper(int iw, AX_3D_Polygons *P);

/* draw the polygons to pixmap (block pixels only) according to z-buffer */
void AX_3D_Polygons_ZBprint(int iw, AX_3D_Polygons *P);

/* draw the polygons to pixmap (alias pixels only) according to z-buffer */
void AX_3D_Polygons_ZAprint(int iw, AX_3D_Polygons *P);

/* Use z-buffer morphed linklist to correctly handle */
/* multiple aliases from polygons WITHIN this group. */
void AX_3D_Polygons_Zaprint(int iw, AX_3D_Polygons *P);

/* Assign parallelepiped (6 polygons) for H[][] originated at x0,y0,z0 */
void AX_3D_ParallelepipedAssign
(AX_Float x0, AX_Float y0, AX_Float z0, AX_Float H[3][3],
 AX_Float r, AX_Float g, AX_Float b, AX_3D_Polygons *P);

/****************/
/* Ball Shading */
/****************/
#define AX_BS_EDGE  0.4  /* color intensity at ball edge */
#define AX_BS_TOPI  0.6  /* color intensity increment if at top */
#define AX_BS_TOP  (AX_BS_EDGE+AX_BS_TOPI)

/****************/
/* Ball Caching */
/****************/
/* tabulation with respect to the circle radius */
#define AX_BC_RMAX   128
#define AX_BC_RMESH  1    /* total of AX_BC_RMAX*AX_BC_RMESH balls cached */
#define AX_FREEZE_WIDTH (2*AX_BC_RMAX+1)
#define AX_BC_BP_RATIO  2.1

typedef struct
{
    AX_Offset offset;  /* in current window configuration */
    AX_Float  c1;   /* covered area for AOP, delta z for BOP */
} AX_BC_Unit;  /* aligned */

typedef struct
{
    AX_Float radius;
    int nc;  /* max row - min row */
    int *c;  /* for each row, min col to max col + 1 */
    int nb;  /* number of block pixels */
    AX_BC_Unit *b;  /* block pixel basic data */
    AX_PAIR *B;
    AX_Float *light_influence;  /* block pixel light intensity */
    int na;  /* number of alias pixels */
    AX_BC_Unit *a;  /* alias pixel basic data */
    AX_PAIR *A;
} AX_BC; /* aligned */

/* Install ball cache for this window */
int AX_plugin_BC_module (int iw);
#define AX_PLUGIN_BC_module(iw) \
  (AX_plugin_BC_module(iw), printf ("ball caching finished\n"))

/* Update ball cache information according to new AX_size[iw] */
int AX_resize_BC_module (int iw);

/* Free ball cache allocations */
int AX_plugout_BC_module (int iw);

/****************/
/* Ball Support */
/****************/

/* ready-to-deploy pixel */
typedef struct
{
    AX_Offset   offset;
    AX_Carrier  carrier;
    AX_Float    c1;  /* covered area for AOP, z for BOP */
} AX_3D_pixel;

/* ball image pointer */
typedef union
{
    AX_BC *bc;        /* if by ball cache */
    AX_3D_pixel *bp;  /* if scanned on the fly */
} AX_3D_BallScan;

/* object status */
#define AX_3D_BALL_COMPLETE  1   /* if not, then clipped */
#define AX_3D_BALL_CACHED    2   /* if not, then scan-on-the-fly */

typedef struct
{
    AX_Float x[3];  /* ball center in real space */
    AX_Float radius;
    AX_Float r,g,b;
    AX_Float vx[3];
    AX_IJ i0,j0;  /* image data */
    AX_3D_BallScan vs;
} AX_3D_Ball;

typedef struct
{
    int n_balls;
    AX_3D_Ball *BALL;
    int n_shown;
    /* z-order to preempt balls that are "visible" but covered */
    int *idx;
    unsigned char *status;  /* in viewport frame */
    /* occupied "on-the-fly" stack space */
    int bp_occupied;
} AX_3D_Balls;

/* Realloc n_balls */
void AX_3D_Balls_Realloc (AX_3D_Balls *B, int n_balls);

/* Realloc n_balls followed by x0,y0,z0,radius,r,g,b,... interface */
void AX_3D_Balls_Recreate (AX_3D_Balls *B, int n_balls, ...);

/* free allocated memory and set NULL */
void AX_3D_Balls_Free (AX_3D_Balls *B);

/* Calculate coordinates in viewport frame and determine visibility */
void AX_3D_Balls_Zdet (int iw, AX_3D_Balls *B);

/* register the balls in z-buffer */
void AX_3D_Balls_Zpaper (int iw, AX_3D_Balls *B);

/* Draw the balls to pixmap (block pixels only) according to z-buffer */
void AX_3D_Balls_ZBprint (int iw, AX_3D_Balls *B);

/* draw the balls to pixmap (alias pixels only) according to z-buffer */
void AX_3D_Balls_ZAprint (int iw, AX_3D_Balls *B);

/* Find out the index of the topmost ball at a given pixel */
/* after Zpaper/ZBprint(); return -1 if there is no match. */
int AX_3D_Balls_Zgrab (int iw, AX_3D_Balls *B, int gi, int gj);


/*********************/
/* Ellipsoid Support */
/*********************/

/***********/
/* Shading */
/***********/
#define AX_ES_EDGE  0.4  /* color intensity at edge */
#define AX_ES_TOPI  0.6  /* color intensity increment if at top */
#define AX_ES_TOP  (AX_ES_EDGE+AX_ES_TOPI)

typedef struct
{
    AX_Float x[3]; /* real coordinates of center */
    AX_Float G[3][3]; /* (x-x0)' * G * (x-x0) = 1 */
    AX_Float r,g,b;
    AX_Float vx[3]; /* coordinates in viewport frame */
    AX_Float a,B,c,d,e,f,D; /* ellipsoid determinant */
    AX_Ellipse p[1]; /* projected ellipse on screen */
    AX_3D_Pixel *img;
} AX_3D_Ellipsoid; /* aligned */

typedef struct
{
    int n_ellipsoids;
    AX_3D_Ellipsoid *ELLIPSOID;
    int n_shown;
    /* z-order to preempt ellipsoids that are "visible" but covered */
    int *idx;
} AX_3D_Ellipsoids;

/* Realloc n_ellipsoids */
void AX_3D_Ellipsoids_Realloc (AX_3D_Ellipsoids *E, int n_ellipsoids);

/* Realloc n_ellipsoids followed by x0,y0,z0,G,r,g,b,... interface */
void AX_3D_Ellipsoids_Recreate (AX_3D_Ellipsoids *E, int n_ellipsoids, ...);

/* free allocated memory and set NULL */
void AX_3D_Ellipsoids_Free (AX_3D_Ellipsoids *E);

/* Calculate coordinates in viewport frame and determine visibility */
void AX_3D_Ellipsoids_Zdet (int iw, AX_3D_Ellipsoids *E);

/* Generate ellipsoid scans and register them in the z-buffer */
void AX_3D_Ellipsoids_Zpaper (int iw, AX_3D_Ellipsoids *E);

/* draw the ellipsoids to pixmap (block pixels only) according to z-buffer */
void AX_3D_Ellipsoids_ZBprint (int iw, AX_3D_Ellipsoids *E);

/* draw the ellipsoids to pixmap (alias pixels only) according to z-buffer */
void AX_3D_Ellipsoids_ZAprint (int iw, AX_3D_Ellipsoids *E);


/********************/
/* Cylinder Support */
/********************/

/***********/
/* Shading */
/***********/
#define AX_CS_EDGE  0.4  /* color intensity at edge */
#define AX_CS_TOPI  0.6  /* color intensity increment if at top */
#define AX_CS_TOP  (AX_CS_EDGE+AX_CS_TOPI)

typedef struct
{
    AX_Float x0[3];
    AX_Float axis[4];  /* normalized x1-x0 and its length */
    AX_Float radius;  /* radius <= 0 is invisibility flag */
    AX_Float r,g,b;
    AX_3D_Pixel *img;
} AX_3D_Cylinder; /* aligned */

#define AX_3D_CYLINDER_ARM_AXIS(C,x1) ( \
  AX_V3sub(x1,(C).x0,(C).axis),AX_V3NORMALIZE((C).axis,(C).axis[3]) )

typedef struct
{
    double dx[3], y0[3], U[4];
    AX_Polygon p[1];
} AX_3D_Cylinder_Power; /* aligned */

typedef struct
{
    int n_cylinders;
    AX_3D_Cylinder *CYLINDER;
    int n_power;
    AX_3D_Cylinder_Power *power;
    AX_3D_Cylinder **cylinder;
    /* z-order to preempt cylinders that are "visible" but covered */
    AX_Float *vz;
    int *idx;
    AX_3D_Pixel *stack;
} AX_3D_Cylinders; /* aligned */

/* Realloc n_cylinders */
void AX_3D_Cylinders_Realloc (AX_3D_Cylinders *C, int n_cylinders);

/* Realloc n_cylinders followed by x0,y0,z0,x1,y1,z1,radius,r,g,b,... */
void AX_3D_Cylinders_Recreate (AX_3D_Cylinders *C, int n_cylinders, ...);

/* Realloc n_cylinders followed by x0[],x1[],radius,r,g,b,... interface */
void AX_3D_Cylinders_RECREATE (AX_3D_Cylinders *C, int n_cylinders, ...);

/* free allocated memory and set NULL */
void AX_3D_Cylinders_Free (AX_3D_Cylinders *C);

/* Calculate coordinates in viewport frame and determine visibility */
void AX_3D_Cylinders_Zdet (int iw, AX_3D_Cylinders *C);

/* generate cylinder scans and register them in the z-buffer */
void AX_3D_Cylinders_Zpaper(int iw, AX_3D_Cylinders *C);

/* draw the cylinders to pixmap (block pixels only) according to z-buffer */
void AX_3D_Cylinders_ZBprint (int iw, AX_3D_Cylinders *C);

/* draw the cylinders to pixmap (alias pixels only) according to z-buffer */
void AX_3D_Cylinders_ZAprint (int iw, AX_3D_Cylinders *C);

/* Find out the index of the topmost cylinder at a given pixel */
/* after Zpaper/Zprint(); return -1 if there is no match.      */
int AX_3D_Cylinders_Zgrab (int iw, AX_3D_Cylinders *C, int gi, int gj);

#endif  /* _AX_h */
