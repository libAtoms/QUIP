/* 2004-2005 Futoshi Shimizu */
#ifndef P3DPATCH_H
#define P3DPATCH_H

#ifndef USE_CUI
#   error
#else
#   include "cui.h"
#endif

#include "P3DExt.h"
#include "P3D.h"

#include "A.h"
#include <setjmp.h>

#define P3DP_NEIGHBORLIST_BINATOM_RATIO 4.0 /* 3.0 */

#define Neighborlist_Recreate(Config_Alib_to_Alib,info,ct,tp,N) do {\
            if (p3dp_enabled)\
                p3dp_Neighborlist_Recreate(Config_Alib_to_Alib,info,ct,tp,N);\
            else     Neighborlist_Recreate(Config_Alib_to_Alib,info,ct,tp,N);\
        } while (0)

#ifndef P3DP_GLOBAL

#   define  CalculateSimpleStatistics(N, series, bytes_separation, valtype, s)\
       p3dp_CalculateSimpleStatistics(N, series, bytes_separation, valtype, s)

#   undef   AXNextEvent
#   define  AXNextEvent(iw)             p3dp_AXNextEvent(iw)
#   undef   AXQLength
#   define  AXQLength(iw)               p3dp_AXQLength(iw)
#   undef   AXPending
#   define  AXPending(iw)               p3dp_AXPending(iw)
#   undef   AXGetGeometry
#   define  AXGetGeometry(iw,newsize)   p3dp_AXGetGeometry(iw,&(newsize))
#   undef   AXQueryPointer
#   define  AXQueryPointer(iw)          p3dp_AXQueryPointer(iw)
#   undef   AXKeysym0
#   define  AXKeysym0(iw)               p3dp_AXKeysym0(iw)

#   undef   AX_openwindow
#   define  AX_openwindow p3dp_AX_openwindow
/*
#   define  AX_openwindow(cid, title, width, height) \
            ((cui_enabled >= 0 && (p3dp_enabled < 1 || !p3d_rank(p3dp_cell)))\
                            ? AX_openwindow(    cid,  title, width, height)\
                        : cui_AX_openwindow(-1*(cid), title, width, height))
*/
#   undef   AX_resizewindow
#   define  AX_resizewindow p3dp_AX_resizewindow
#   undef   AX_closewindow
#   define  AX_closewindow p3dp_AX_closewindow
int p3dp_AX_openwindow
(int cid, char *title, unsigned int width, unsigned int height);
void p3dp_AX_resizewindow (int iw, Bool do_window_too);
void p3dp_AX_closewindow (int iw);

#   define  AX_3D_Balls_Zdet(iw, B)\
            switch (p3dp_enabled) {\
            case 0:       AX_3D_Balls_Zdet(iw, B); break;\
            default: p3dp_AX_3D_Balls_Zdet(iw, B); break;}

#   define  AX_3D_Cylinders_Zdet(iw, C)\
            switch (p3dp_enabled) {\
            case 0:       AX_3D_Cylinders_Zdet(iw, C); break;\
            default: p3dp_AX_3D_Cylinders_Zdet(iw, C); break;}

#   define  P3DP_GLOBAL extern
#endif

P3DP_GLOBAL Cell p3dp_cell;
P3DP_GLOBAL jmp_buf quit_env;
P3DP_GLOBAL int p3dp_enabled;
P3DP_GLOBAL char *p3dp_decomp;
P3DP_GLOBAL char **p3dp_mpiexec_argv;
P3DP_GLOBAL int p3dp_mpiexec_argc;
P3DP_GLOBAL int p3dp_rank_grab;
P3DP_GLOBAL char p3dp_spec[(MENDELEYEV_MAX+1)*SYMBOL_CHAR+1];

typedef struct {
    int rank_stack[ATOM_STACK_SIZE];
} P3DP_Navigator;
P3DP_GLOBAL P3DP_Navigator p3dp_n[AX_MAXWIN];

#undef  IS_MANAGER
#define IS_MANAGER  (!p3dp_enabled || !p3d_rank(p3dp_cell))

#undef  CONFIG_LOAD
#define CONFIG_LOAD(fname,Config_Alib_to_Alib) \
                p3dp_Config_Load(fname,stdout,Config_Alib_to_Alib)

#undef  rebind_CT
#define rebind_CT(Config_Alib_to_Alib, specification, ct, tp) (!p3dp_enabled) ?\
                rebind_ct(Config_Alib_to_Alib, specification, ct, tp, stdout) :\
                rebind_ct(Config_Alib_to_Alib, p3dp_spec,     ct, tp,\
                                                (IS_MANAGER) ? stdout : NULL)

#define p3dp_Bcast(buf, count, type) do { if (p3dp_enabled) \
    p3d_bcast(p3dp_cell, buf, count, type, 0); } while (0)


int p3dp_init(int *argc, char ***argv);
int p3dp_finalize(void);
int p3dp_Config_Load(char *fname, FILE *info, Alib_Declare_Config);
void p3dp_config_scatter(Alib_Declare_Config);
void p3dp_zmap_reduce(int iw);
void p3dp_atom_stack_insert(int iw, int atom, int rank);
void p3dp_CalculateSimpleStatistics
(int N, char *series, int bytes_separation, int valtype, SimpleStatistics *s);
void p3dp_bcast_int2(int *val0, int *val1);

int p3dp_AXNextEvent(int iw);
int p3dp_AXQLength(int iw);
int p3dp_AXPending(int iw);
int p3dp_AX_plugin_3D_module(int iw);
Status p3dp_AXGetGeometry(int iw, AXSize *newsize);
Bool p3dp_AXQueryPointer(int iw);
KeySym p3dp_AXKeysym0(int iw);

void p3dp_Neighborlist_Recreate
(Alib_Declare_Config, FILE *info, Chemtab *ct, Tp **tp, Neighborlist *N);

/* p3dp_info.c */
void p3dp_s(double *si, int i, int irank);
void p3dp_atom_pair_s(double *sj, double *si, double dxji[4]);
double p3dp_atom_triplet_s
(double *sk, double *sj, double *si, double dxkj[4], double dxij[4]);
void p3dp_print_atom_pair_info(int iw, int j, int i, int jrank, int irank);
void p3dp_print_atom_triplet_info
(int iw, int k, int j, int i, int krank, int jrank, int irank);
void p3dp_print_atom_quartet_info
(int iw, int l, int k, int j, int i, int lrank, int krank, int jrank,int irank);

/* p3dp_AX.c */
void p3dp_AX_3D_Balls_Zdet(int iw, AX_3D_Balls *B);
int p3dp_AX_3D_Balls_Zgrab(int iw, AX_3D_Balls *B, int gi, int gj, int *rank);
void p3dp_AX_3D_Cylinders_Zdet(int iw, AX_3D_Cylinders *C);
int p3dp_AX_3D_Cylinders_Zgrab(int iw,AX_3D_Cylinders*C,int gi,int gj,int*rank);
#endif
