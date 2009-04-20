/* 2004-2005 Futoshi Shimizu */
#ifndef CUI_H
#define CUI_H

#include "A.h"

#define CUI_GREETINGS       "Threaded atomistic configuration viewer V3.0"
#define CUI_N_FILE          5
#define CUI_TIMEOUT         {0, 10000}
#define CUI_LINEMAX         2048
#define CUI_SCRIPT_DEFAULT  "scr_cmds"
#define CUI_SCRIPT_STARTUP  ".A"
#define CUI_PORT_FROM       21040
#define CUI_PORT_TO         21049
#define CUI_PROMPT          "> "
#define CUI_PROTOCOL_OK     "ok"
#define CUI_PROTOCOL_NG     "error"
#define CUI_PROTOCOL_QUIT   "quit"
#define CUI_ARG_REQUIRED    "arg_required"
#define CUI_EPS_RESOLUTION_DEFAULT  2560
#define CUI_XTERM_DEFAULT \
                "xterm -bg gray40 -fg white -sb -sl 3000 -cr yellow -fn 7x13"



#define IS_MANAGER          1

#ifdef  CUI_GLOBAL
    int cui_startup_width = DEFAULT_WIDTH, cui_startup_height = DEFAULT_HEIGHT;
    void cui_AX_dump(int iw) { if (AX_display[iw]) AX_dump(iw); }
#else
#   define  CUI_GLOBAL extern
    CUI_GLOBAL int cui_startup_width, cui_startup_height;
    void cui_AX_dump(int iw);
#endif

CUI_GLOBAL int cui_enabled;
CUI_GLOBAL int cui_listen;
CUI_GLOBAL int cui_connect;
CUI_GLOBAL int cui_stdin;
CUI_GLOBAL int cui_startup_width;
CUI_GLOBAL int cui_startup_height;
CUI_GLOBAL int cui_argc;
CUI_GLOBAL int cui_iw;
CUI_GLOBAL char *cui_argv[128];
CUI_GLOBAL char *cui_hostport;
CUI_GLOBAL char *cui_scrfname;
CUI_GLOBAL char *cui_geometry;
CUI_GLOBAL char cui_title[CUI_LINEMAX];
CUI_GLOBAL bool cui_diligence;
CUI_GLOBAL double cui_time;
CUI_GLOBAL Window cui_xterm_win;

/* cui.c */
double cui_wtime(void);
#ifdef ATOMEYE_LIB
int cui_init(int *argc, char ***argv,  void (*on_click)(int atom), void (*on_close)(), void (*on_advance)(char *instr));
#else
int cui_init(int *argc, char ***argv);
#endif
bool cui_treatevent(int iw);
bool gui_treatevent(int iw);
int cui_config_load_A3(Alib_Declare_Config);

/* cui_AX.c */
int cui_AX_openwindow
(int cid, char *title, unsigned int width, unsigned int height);
void cui_AX_resizewindow(int iw, Bool do_window_too);
void cui_AX_closewindow(int iw);

/* A */

#undef  DEFAULT_WIDTH
#define DEFAULT_WIDTH   cui_startup_width
#undef  DEFAULT_HEIGHT
#define DEFAULT_HEIGHT  cui_startup_height

/* Atoms */

#undef  CONFIG_LOAD
#define CONFIG_LOAD(fname,Config_Alib_to_Alib)\
                ((strcmp(fname,"/dev/null") == 0)\
                    ?       cui_config_load_A3(Config_Aapp_to_Alib)\
                    : Config_Load(fname,stdout,Config_Alib_to_Alib))

/* X & AX */
#define XStoreName(display, w, window_name)\
            ((display&&w!=xterm_win) ? XStoreName(  display,w,window_name) : 0)
#define XSetIconName(display, w, window_name)\
            ((display&&w!=xterm_win) ? XSetIconName(display,w,window_name) : 0)
#define XMapWindow(display, w)\
            ((display&&w) ? XMapWindow(display, w) : 0)
#define XIconifyWindow(display, w, screen_number)\
            (display ? XIconifyWindow(display, w,screen_number) : (Status)0)

#define CUI_AX(IW, type, value) (AX_display[IW] ? (value) : (type)0)
#define CUI_TITLE(iw) ((!cui_enabled||iw!=cui_iw) ? AX_title[iw] :\
            strncpy(cui_title+1, AX_title[iw], sizeof(cui_title))-1)

#undef  AXQueryPointer
#define AXQueryPointer(iw) CUI_AX(iw, Bool,\
         XQueryPointer(AX_display[iw],AX_win[iw],&AX_root[iw],&AX_child[iw],\
        &AX_root_x[iw],&AX_root_y[iw],&AX_win_x[iw],&AX_win_y[iw],\
        &AX_pointermask[iw]))
#undef  AXSetName
#define AXSetName(iw)   CUI_AX(iw, int, (\
            XStoreName(AX_display[iw],AX_win[iw],CUI_TITLE(iw)),\
            XSetIconName(AX_display[iw],AX_win[iw],CUI_TITLE(iw))))
#undef  AXIdentifyWIN
#define AXIdentifyWIN(iw,identifier) CUI_AX(iw, Window,\
        AXIdentifyWin(iw, AX_root[iw], identifier))
#undef  AXSETICON
#define AXSETICON(iw)   CUI_AX(iw, int, AXSetICON(iw,icon))
#undef  AXRaiseWindow
#define AXRaiseWindow(iw) CUI_AX(iw,int,XRaiseWindow(AX_display[iw],AX_win[iw]))

#undef  AX_show
#define AX_show(iw) if (AX_display[iw])\
                    { if (AX_videomode == AX_VIDEOMODE_SHM_PIXMAP) \
  XCopyArea(AX_display[iw], AX_pixmap[iw], AX_win[iw], AX_gc[iw], 0,0, \
  AX_size[iw].width, AX_size[iw].height, 0,0); AXSYNC(iw); }

#undef  AX_dump
#define AX_dump(iw)                 cui_AX_dump(iw)
#define AX_openwindow(cid, title, width, height) \
                    ((cui_enabled >= 0)\
                            ? AX_openwindow(    cid,  title, width, height)\
                        : cui_AX_openwindow(-1*(cid), title, width, height))
#define AX_resizewindow(iw, dwt)\
                    switch ((AX_cid[iw] >= 0)) {\
                    case 0: cui_AX_resizewindow(iw, dwt); break;\
                    default:    AX_resizewindow(iw, dwt); break;\
                    }
#define AX_closewindow(iw)\
                    switch ((AX_cid[iw] >= 0)) {\
                    case 0: cui_AX_closewindow(iw); break;\
                    default:    AX_closewindow(iw); break;\
                    }
#endif
