/* $Id: P3D.h,v 1.10 2005/07/23 17:36:36 shimizu Exp $
 * 2004-2005 Futoshi SHIMIZU
 */
#ifndef P3D_H
#define P3D_H

#include "mpi.h"
#include "VecMat3.h"
#include "P3DExt.h"

#define P3D_INFINITY 1.e+308
#define P3D_LINEMAX 1024

#ifndef CONFIG_MAX_AUXILIARY
#define CONFIG_MAX_AUXILIARY 32
#endif


typedef struct {
#ifdef P3DAtomPRE
    P3DAtomPRE
#endif
    V3 r, v;                /* position, velocity */
    double m;               /* mass */
    char sym[sizeof(int)];  /* symbol */
    int iw;                 /* work (for P3D) */
    double aux[CONFIG_MAX_AUXILIARY];
#ifdef P3DAtomMD
    V3 a, f;        /* acceleration, force */
    double p;       /* potential */
    M3 vir;         /* virial */
#endif
#ifdef P3DAtomPOST
    P3DAtomPOST
#endif
} P3DAtom;

typedef struct Auxiliary_tag {
    int n_aux;
    char name[CONFIG_MAX_AUXILIARY][P3D_LINEMAX];
    char unit[CONFIG_MAX_AUXILIARY][P3D_LINEMAX];
} *Auxiliary;


typedef struct Cell_tag {
    M3 h, virial;
    int idx_trans;
    int pbc[3];
    int reduced_coordinates;
    Auxiliary auxiliary;
    struct Cell_private_tag *private;
} *Cell;

typedef struct Pair_tag {
    P3DAtom **begin, **end;
} *Pair;

typedef struct List_tag *List;



enum p3d_index_flags {
    I_ENT, I_ONT, I_WNT,
    I_EOT, I_OOT, I_WOT,
    I_EST, I_OST, I_WST,

    I_ENO, I_ONO, I_WNO,
    I_EOO, I_OOO, I_WOO,
    I_ESO, I_OSO, I_WSO,

    I_ENB, I_ONB, I_WNB,
    I_EOB, I_OOB, I_WOB,
    I_ESB, I_OSB, I_WSB,
    I_XYZ
};

enum p3d_flags_index {
    F_OOO = 0x00,
    F_EOO = 0x20,
    F_WOO = 0x10,
    F_ONO = 0x08,
    F_ENO =     F_EOO | F_ONO,
    F_WNO =     F_WOO | F_ONO,
    F_OSO = 0x04,
    F_ESO =     F_EOO | F_OSO,
    F_WSO =     F_WOO | F_OSO,
    F_OOT = 0x02,
    F_EOT =     F_EOO |         F_OOT,
    F_WOT =     F_WOO |         F_OOT,
    F_ONT =             F_ONO | F_OOT,
    F_ENT =     F_EOO | F_ONO | F_OOT,
    F_WNT =     F_WOO | F_ONO | F_OOT,
    F_OST =             F_OSO | F_OOT,
    F_EST =     F_EOO | F_OSO | F_OOT,
    F_WST =     F_WOO | F_OSO | F_OOT,
    F_OOB = 0x01,
    F_EOB =     F_EOO |         F_OOB,
    F_WOB =     F_WOO |         F_OOB,
    F_ONB =             F_ONO | F_OOB,
    F_ENB =     F_EOO | F_ONO | F_OOB,
    F_WNB =     F_WOO | F_ONO | F_OOB,
    F_OSB =             F_OSO | F_OOB,
    F_ESB =     F_EOO | F_OSO | F_OOB,
    F_WSB =     F_WOO | F_OSO | F_OOB
};



#ifndef P3D_GLOBAL
#   define P3D_GLOBAL extern
#endif
P3D_GLOBAL MPI_Comm p3d_global_comm;
P3D_GLOBAL int p3d_global_nprocs, p3d_global_rank;
P3D_GLOBAL int p3d_global_io_any, p3d_global_rank_fprintf;
P3D_GLOBAL char *p3d_global_decomp_string;


/*#undef P3D_QUIET*/

/**
#define p3d_abort(errorcode) (fprintf(stderr, "[%d] Abort, %s %u(%s): %d\n",\
        p3d_global_rank, __FILE__, __LINE__, __func__, errorcode),\
        fflush(stderr), MPI_Abort(MPI_COMM_WORLD, errorcode))
**/
/**/
#define p3d_abort(errorcode) (fprintf(stderr, "[%d] Abort, %s %u: %d\n",\
        p3d_global_rank, __FILE__, __LINE__, errorcode),\
        fflush(stderr), MPI_Abort(MPI_COMM_WORLD, errorcode))
/**/


#define P3D_LOOP(cell,pnt) \
    for (pnt = p3d_atom_begin(cell); pnt < p3d_atom_end(cell); pnt++)

#define P3D_LOOP_PAIR_HEAD(pair,pnt)\
        {P3DAtom **pp; for (pnt = *(pp=pair->begin); pp < pair->end; pnt=*++pp){

#define P3D_LOOP_PAIR_TAIL(pair) }} pair++

#define p3d_zero(cell,member) \
    do {P3DAtom *pnt; P3D_LOOP(cell,pnt) pnt->member=0; } while(0)

#define p3d_V3zero(cell,member) \
    do {P3DAtom *pnt; P3D_LOOP(cell,pnt) V3ZERO(pnt->member); } while(0)

#define p3d_M3zero(cell,member) \
    do {P3DAtom *pnt; P3D_LOOP(cell,pnt) M3ZERO(pnt->member); } while(0)

#define DWIDTH "24"
#define DDIGIT "16"
#define DFORMAT "%" DWIDTH "." DDIGIT "e"
#define VFORMAT DFORMAT DFORMAT DFORMAT
#define MFORMAT VFORMAT VFORMAT VFORMAT
#define FFORMAT " %." DWIDTH "f"
#define FFORMAT3 FFORMAT FFORMAT FFORMAT


/* Core.c */

int p3d_init(int *argc,char ***argv, MPI_Comm comm);
int p3d_finalize(void);

int p3d_fflush(FILE *fp);
int p3d_fprintf(FILE *fp, const char *format, ...);
int p3d_message(const char *format, ...);

MPI_Comm p3d_comm(Cell cell);
int p3d_nprocs(Cell cell);
int p3d_rank(Cell cell);
int p3d_n_atoms(Cell cell);
int p3d_n_images(Cell cell);
P3DAtom *p3d_atom_begin(Cell cell);
P3DAtom *p3d_atom_end(Cell cell);
P3DAtom *p3d_image_begin(Cell cell);
P3DAtom *p3d_image_end(Cell cell);

Cell p3d_new_cell(char *decomp, MPI_Datatype m_image, MPI_Comm comm);
int p3d_set_atoms(Cell cell, P3DAtom *atoms, int n);
int p3d_cat_atoms(Cell cell, P3DAtom *atoms, int n);
int p3d_remove_atom(Cell cell, P3DAtom *atom);
int p3d_reset_cell_c(Cell cell, double crust);
#define p3d_reset_cell(cell) p3d_reset_cell_c(cell, 0.0)
int p3d_update_cell(Cell cell);

int p3d_bcast(Cell cell, void *buff, int count,
                        MPI_Datatype datatype, int root);
int p3d_reduce(Cell cell, void *sbuf, void *rbuf, int count,
                        MPI_Datatype datatype, MPI_Op op);

int p3d_read_config(Cell cell, char *path);
int p3d_write_config(Cell cell, char *path,
        char *s_formats, char*velocity_formats, int n_aux, ...);

void p3d_coord_reduced2real(Cell cell);


/* P3DNeighborList.c */

List p3d_new_list(Cell cell);
void p3d_realloc_list(List list, int n_max, int p_max, int i_max);
int p3d_reset_list(List list, Cell cell, double cutoff);
Pair p3d_pair_atom_begin(List list);
Pair p3d_pair_image_begin(List list);

#define p3d_reset_both(cell, list, rc_list) \
    do {\
        p3d_reset_cell_c(cell, rc_list); p3d_reset_list(list, cell, rc_list);\
    } while (0)


/* P3DUtils.c */

#define p3d_set_a_from_f(cell) \
    do { P3DAtom *pnt; P3D_LOOP(cell,pnt)\
        V3MUL(1.0/pnt->m, pnt->f, pnt->a);\
    } while (0)


enum {
    F_FORCE     = 0x01,
    F_POTENTIAL = 0x02,
    F_VIRIAL    = 0x04,
    F_ALL       = F_FORCE | F_POTENTIAL | F_VIRIAL
};

void calc_interaction(Cell cell, List list, int flags);

void p3d_update_r_v_verlet(Cell cell, double dt);
void p3d_update_va_v_verlet(Cell cell, double dt);
void p3d_update_gear3_pred(Cell cell, double dt);
void p3d_update_gear3_corr(Cell cell, double dt);

double p3d_step_sd(Cell cell, double alpha, double drmax);

double p3d_sum_kinetic_energy(Cell cell);
double p3d_sum_potential_energy(Cell cell);
double p3d_sum_momentum(Cell cell, V3 momentum);
void p3d_sum_stress(Cell cell, M3 stress);
double p3d_calc_temperature(Cell cell, double k_blt);
double p3d_calc_temperature2(Cell cell, double k_blt);
void p3d_velocity_scaling(Cell cell, double temperature, double k_blt);
void p3d_velocity_scaling2(Cell cell, double temperature, double k_blt);
void p3d_set_temperature(Cell cell, double temp, double k_blt);
void p3d_cm_vzero(Cell cell);
double p3d_sum_msd(Cell cell, int offset_0, int offset_s, int reset);

void p3d_set_number(Cell cell, int offset);
void p3d_scale_cell(Cell cell, double target_pressure, double pk, int ip[3]);

/*
*/

#endif
