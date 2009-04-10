/* 2004-2005 Futoshi Shimizu */
#include "cui.h"
#undef  AX_openwindow
#undef  AX_resizewindow
#undef  AX_closewindow
#define nonvacant(iw) (AX_cid[iw]!=0)

static int cui_AX_plugin_ShmGC_module (int iw)
{
    AX_AssertWinID ("cui_AX_plugin_ShmGC_module", iw);
    AX_Check_Module_Dependence ("cui_AX_plugin_ShmGC_module", iw,MODULE_ShmGC);
    if ( AX_Module_Is_PluggedIn [iw] [MODULE_ShmGC] )
        pe ("cui_AX_plugin_ShmGC_module:\n"
            "you are crazy, this module is already plugged in\n");

    /*switch (AX_depth[iw]) {
    case 8:  AX_bytes = 1; break;
    case 16: AX_bytes = 2; break;
    case 24: case 32: AX_bytes = 4; break;
    default: AX_bytes = 4;
    }*/
    AX_bytes = 4;

    /* non-shared client data-memory */
    MALLOC ( AX_plugin_ShmGC_module, AX_mem[iw].uc,
             AX_size[iw].width * AX_size[iw].height * AX_bytes,
             unsigned char );

/*  AX_bytes = AX_img[iw] -> bytes_per_line / AX_size[iw].width; */
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

    AX_namedbg(iw,AX_DEF_BACKGROUND);
    AX_Module_Is_PluggedIn [iw] [MODULE_ShmGC] = 1;
    return (1);
} /* end cui_AX_plugin_ShmGC_module() */


static int cui_AX_plugout_ShmGC_module (int iw)
{
    AX_AssertWinID ("cui_AX_plugout_ShmGC_module", iw);
    if ( ! AX_Module_Is_PluggedIn [iw] [MODULE_ShmGC] )
        pe ("cui_AX_plugout_ShmGC_module:\n"
            "you are crazy, this module is not plugged in\n");

    Free(AX_mem[iw].uc);
    AX_Module_Is_PluggedIn [iw] [MODULE_ShmGC] = 0;
    return (1);
} /* end cui_AX_plugout_ShmGC_module() */


static int cui_AX_resize_ShmGC_module (int iw)
{
    REALLOC( AX_plugin_ShmGC_module, AX_mem[iw].uc,
             AX_size[iw].width * AX_size[iw].height * AX_bytes,
             unsigned char );
    return (1);
} /* end cui_AX_resize_ShmGC_module() */


static struct {
    int (*plugin)  (int iw);
    int (*plugout) (int iw);
    int (*resize)  (int iw);
} cui_AX_modules [MODULE_MAX] = {
    {&cui_AX_plugin_ShmGC_module,
                        &cui_AX_plugout_ShmGC_module,
                                                 &cui_AX_resize_ShmGC_module},
    {&AX_plugin_Scan_module, &AX_plugout_Scan_module, &AX_resize_Scan_module},
    {&AX_plugin_3D_module,   &AX_plugout_3D_module,   &AX_resize_3D_module},
    {&AX_plugin_BC_module,   &AX_plugout_BC_module,   &AX_resize_BC_module}
};


int cui_AX_openwindow
(int cid, char *title, unsigned int width, unsigned int height)
{
    int iw, minor;

    for (iw = 0; (iw < AX_MAXWIN) && nonvacant(iw); iw++);
    if (iw == AX_MAXWIN) {
        fprintf(stderr,"cui_AX_openwindow: AX_MAXWIN = %d reached\n",AX_MAXWIN);
        return -1;
    }
    
    AX_display[iw]          = NULL;
    AX_cid[iw]              = cid;
    AX_size[iw].width       = width;
    AX_size[iw].height      = height;
    AX_borderwidth[iw]      = 0;
    AX_depth[iw]            = 32;
    AX_rmask                = 0xff0000;
    AX_gmask                = 0x00ff00;
    AX_bmask                = 0x0000ff;
    AX_namedpixel[AX_RED]   = AX_Colorpixel(1,0,0);
    AX_namedpixel[AX_GREEN] = AX_Colorpixel(0,1,0);
    AX_namedpixel[AX_BLUE]  = AX_Colorpixel(0,0,1);
    AX_noneedcolormap       = 1;
    AX_videomode            = AX_VIDEOMODE_NO_SHM;

    for (minor=0; minor<=AX_DEF_Module_Support_Level; minor++)
        cui_AX_modules[minor].plugin(iw);

    return iw;
}

/* Resize all modules that are already plugged in, to AX_size[iw] */
void cui_AX_resizewindow (int iw, Bool do_window_too)
{
    int module;
    AX_AssertWinID ("AX_resizewindow", iw);
    for (module=0; module<MODULE_MAX; module++)
        if (AX_Module_Is_PluggedIn[iw][module])
            cui_AX_modules[module].resize(iw);
    return;
} /* end AX_resizewindow() */


/* close the window and free all resources */
void cui_AX_closewindow (int iw)
{
    int jw;
    AX_AssertWinID ("AX_closewindow", iw);
    for (jw=0; jw<MODULE_MAX; jw++)
        if (AX_Module_Is_PluggedIn [iw] [jw])
            cui_AX_modules[jw].plugout(iw);
    /* close the window */
    if (AX_display[iw]) {
        AXDestroyWindow(iw);
        AXSYNC(iw);
    }
    for (jw = 0; jw < AX_MAXWIN; jw++)
        if ( (jw != iw) && (AX_cid[jw] == AX_cid[iw]) )
            break;
    /* if no officemate, also pull the phone line */
    if (jw == AX_MAXWIN && AX_display[iw]) {
        XFreeColormap (AX_display[iw], AX_colormap[iw]);
        AXCloseDisplay(iw);
    }
    /* unregister this care provider */
    AX_cid [iw] = 0;
    return;
} /* end AX_closewindow() */
