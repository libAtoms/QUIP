/* 2004-2005 Futoshi Shimizu */
#include "P3DExt.h"
#include "P3D_p.h"
#define  P3DP_GLOBAL
#include "p3dp.h"
#include "Atoms.h"
#include <stdlib.h>


#ifdef AX_MEMORY_EFFICIENT
#   define MPI_AX_FLOAT        MPI_FLOAT
#   define MPI_AX_FLOAT_INT    MPI_FLOAT_INT
#else
#   define MPI_AX_FLOAT        MPI_DOUBLE
#   define MPI_AX_FLOAT_INT    MPI_DOUBLE_INT
#endif

static MPI_Datatype m_sym;
static MPI_Op m_minsym;
static void minsym_func(void *in, void *inout, int *len, MPI_Datatype *datatype)
{
    int i, j;
    char *c_in = in, *c_inout = inout;

    for (j = 0; j < *len; j++) {
        if (strncmp(c_in, c_inout, SYMBOL_CHAR) < 0)
            for (i = 0; i < SYMBOL_CHAR; i++)
                c_inout[i] = c_in[i];
        c_in    += SYMBOL_CHAR;
        c_inout += SYMBOL_CHAR;
    }
}

int p3dp_init(int *argc, char ***argv)
{
    int i;
    static char buf[CUI_LINEMAX] = "";

    p3dp_mpiexec_argc = 0;
    p3dp_enabled = 0;
    p3dp_decomp = "";

    if (argv && *argv) {
        if ((i = strlen((*argv)[0])) > 0 && (*argv)[0][i-1] == 'p') {
            p3dp_enabled = cui_enabled = 1;
            MPI_Init(argc, argv);
        }

	for (i = 1; i < *argc; i++) {
            if (strncmp((*argv)[i], "-mpi=", strlen("-mpi=")) == 0) {
                strncat(buf, " ", sizeof(buf));
                strncat(buf, (*argv)[i]+strlen("-mpi="), sizeof(buf));
                (*argv)[i] = 0;
            }
            else if (strncmp((*argv)[i], "-decomp=", strlen("-decomp=")) == 0 ||
                     strncmp((*argv)[i], "-DECOMP=", strlen("-DECOMP=")) == 0) {
                p3dp_decomp = (*argv)[i];
                (*argv)[i] = 0;
            }
        }
    }
    {
        int j = i = 0;
        while (j < *argc) {
            while ((*argv)[j] == 0 && j < *argc) j++;
            if (j < *argc) (*argv)[i++] = (*argv)[j++];
        }
        if (!(*argv)[i-1]) i--;
        *argc = i;
    }

    if (!p3dp_enabled && buf[0]) { /* mpiexec */
        int len = strlen(buf);
        char **p;

        p3dp_enabled = cui_enabled =  1;
        if (!*p3dp_decomp)
            p3dp_decomp = "-decomp=0";

        i = 0;
        while (i < len) {
            while (isspace(buf[i]) && i < len) i++;
            if (i < len) {
                p3dp_mpiexec_argc++;
                while (!isspace(buf[i]) && i < len) i++;
            }
        }
        if (p3dp_mpiexec_argc == 0) {
            strcpy(buf, " mpiexec -n 1");
            p3dp_mpiexec_argc = 3;
        }
        else if (p3dp_mpiexec_argc == 1 && (i = atoi(buf)) > 0) {
            sprintf(buf, " mpiexec -n %d", i);
            p3dp_mpiexec_argc = 3;
        }

        fprintf(stderr, "[MPI process startup command: \"%s\"]\n", buf+1);
        len = strlen(buf);
        for (i = 0; i < len; i++)
            if (isspace(buf[i]))
                buf[i] = 0;

        MALLOC(p3dp_init, p3dp_mpiexec_argv, p3dp_mpiexec_argc+1, char*);
        p = p3dp_mpiexec_argv;
        i = 0;
        while (i < len) {
            while ((!buf[i] || isspace(buf[i])) && i < len) i++;
            if (i < len) {
                *p++ = buf+i;
                while (buf[i] && !isspace(buf[i]) && i < len) i++;
            }
        }
        *p = NULL;
        return 0;
    }

    if (*p3dp_decomp)
        p3dp_enabled = (p3dp_decomp[1] == 'D') ? -1: 1;

    if (p3dp_enabled) { /*p3dp_init0(argc, argv);*/
        MPI_Datatype m_image, types[2] = {MPI_CHAR, MPI_UB};
        MPI_Aint displacements[2] = {0, sizeof(P3DAtom)};
        int blocklengths[2] = {(size_t)&((P3DAtom*)0)->iw, 1};

        MPI_Initialized(&i);
        if (!i)
            MPI_Init(argc, argv);
        p3d_init(argc, argv, MPI_COMM_NULL);

        MPI_Type_contiguous(SYMBOL_CHAR, MPI_CHAR, &m_sym);
        MPI_Type_commit(&m_sym);
        MPI_Op_create(minsym_func, True, &m_minsym);

        MPI_Type_struct(2, blocklengths, displacements, types,&m_image);
        MPI_Type_commit(&m_image);

        if (!*p3dp_decomp)
            p3dp_decomp = "-decomp=0";
        p3dp_cell = p3d_new_cell(p3dp_decomp+strlen("-decomp="),
                                                m_image, MPI_COMM_WORLD);
        p3dp_cell->reduced_coordinates = 1;
    }

    return 0;
}


int p3dp_finalize(void)
{
    if (p3dp_enabled) {
        p3d_finalize();
        MPI_Finalize();
    }
    return 0;
}


static int compsym(const void *s1, const void *s2)
{
    return strncmp((char*) s1, (char*) s2, SYMBOL_CHAR);
}


static void p3dp_config_add_image(Alib_Declare_Config)
{
    Cell_private cp = p3dp_cell->private;
    P3DAtom *p;
    V3 lb, ub;
    enum p3d_flags_index f_xyz[3] = {F_EOO|F_WOO, F_ONO|F_OSO, F_OOT|F_OOB};
    int i, j, k, n, nspec_l, nspec_g, mask = 0;
    double rcut, rcut_max = 0.0;
    char spec[(MENDELEYEV_MAX+1)*SYMBOL_CHAR+1] = "";
    char sym_i[SYMBOL_SIZE], sym_j[SYMBOL_SIZE], sym_end[SYMBOL_SIZE];

    for (i = 0; i < SYMBOL_CHAR; i++)
        sym_end[i] = 0x7f;
    sym_end[SYMBOL_CHAR] = 0;

    nspec_l = 0;
    for (p = p3d_atom_begin(p3dp_cell); p < p3d_atom_end(p3dp_cell); p++) {
        for (i = 0; i < nspec_l; i++)
            if (strncmp(p->sym, &spec[i*SYMBOL_CHAR], SYMBOL_CHAR) == 0)
                break;
        if (i == nspec_l) {
            for (j = 0; j < SYMBOL_CHAR; j++)
                spec[nspec_l*SYMBOL_CHAR+j] = p->sym[j];
            nspec_l++;
        }
    }
    qsort(spec, nspec_l, SYMBOL_CHAR, compsym);
    spec[(nspec_l+1)*SYMBOL_CHAR] = 0;
    strncat(spec, sym_end, SYMBOL_CHAR);
    i = 0;
    for (nspec_g = 0; nspec_g < MENDELEYEV_MAX; nspec_g++) {
        p3d_reduce(p3dp_cell, &spec[i*SYMBOL_CHAR],
                          &p3dp_spec[nspec_g*SYMBOL_CHAR], 1, m_sym, m_minsym);
        if (strncmp(&p3dp_spec[nspec_g*SYMBOL_CHAR], sym_end, SYMBOL_CHAR) == 0)
            break;
        if (i < nspec_l && strncmp(&p3dp_spec[nspec_g*SYMBOL_CHAR],
                                  &spec[i*SYMBOL_CHAR], SYMBOL_CHAR) == 0)
            i++;
    };
    p3dp_spec[(nspec_g)*SYMBOL_CHAR] = 0;
    
    sym_i[SYMBOL_CHAR] = sym_j[SYMBOL_CHAR] = 0;
    for (i = 0; i < nspec_g*SYMBOL_CHAR; i += SYMBOL_CHAR) {
        for (k = 0; k < SYMBOL_CHAR; k++)
            sym_i[k] = p3dp_spec[i + k];
        for (j = i; j < nspec_g*SYMBOL_CHAR; j += SYMBOL_CHAR) {
            for (k = 0; k < SYMBOL_CHAR; k++)
                sym_j[k] = p3dp_spec[j + k];
            rcut = NEIGHBORLIST_RCUT_RATIO *
                               (ATOM_RADIUS_IN_A(search_atom_by_symbol(sym_i))+
                                ATOM_RADIUS_IN_A(search_atom_by_symbol(sym_j)));
            if (rcut > rcut_max) rcut_max = rcut;
        }
    }
    rcut_max += 0.1;

    p3d_reset_cell_c(p3dp_cell, rcut_max);

    V3EQV(cp->lb, lb);
    V3EQV(cp->ub, ub);
    for (i = 0; i < 3; i++) {
        if (cp->dims[i] == 1)
            mask |= f_xyz[i];
        if (lb[i] < 0)
            lb[i] = 0;
        if (ub[i] >= 1.0)
            ub[i]  = 1.0;
    }

    n = 0;
    for (p = p3d_image_begin(p3dp_cell); p < p3d_image_end(p3dp_cell); p++)
        if (!(p->iw = mask & p3d_flags_boundary(p->r, lb, ub))) n++;

    *np = p3d_n_atoms(p3dp_cell) + n;
    Config_realloc (Config_Alib_to_Alib);
    for (k = 0; k < CONFIG_num_auxiliary; k++)
        REALLOC(p3dp_config_add_image, CONFIG_auxiliary[k], *np, double);

    n = 0;
    for (p = p3d_atom_begin(p3dp_cell); p < p3d_atom_end(p3dp_cell); p++) {
        (*mass)[n] = p->m;
        strncpy(SYMBOL(n), p->sym, 3);
        (*s )[n * DIMENSION    ] = p->r[0];
        (*s )[n * DIMENSION + 1] = p->r[1];
        (*s )[n * DIMENSION + 2] = p->r[2];
        for (k = 0; k < CONFIG_num_auxiliary; k++)
            CONFIG_auxiliary[k][n]= p->aux[k];
        n++;
    }

    for (p = p3d_image_begin(p3dp_cell); p < p3d_image_end(p3dp_cell); p++) {
        if (p->iw) continue;
        (*mass)[n] = p->m;
        strncpy(SYMBOL(n), p->sym, 3);
        (*s )[n * DIMENSION    ] = p->r[0];
        (*s )[n * DIMENSION + 1] = p->r[1];
        (*s )[n * DIMENSION + 2] = p->r[2];
        for (k = 0; k < CONFIG_num_auxiliary; k++)
            CONFIG_auxiliary[k][n]= p->aux[k];
        n++;
    }
}

static int p3dp_Config_load(char *fname, FILE *info, Alib_Declare_Config)
{
    Auxiliary auxiliary = p3dp_cell->auxiliary;
    int k;
    extern void Config_analyze_H (Alib_Declare_Config, FILE *out);

    if (info) Fprintf(info, "Loading configuration from \"%s\":\n", fname);

    p3d_read_config(p3dp_cell, fname);
    M3EQV(p3dp_cell->h, H);

    *np = p3d_n_atoms(p3dp_cell);
    if (info) Fprintf(info, "%d atoms found.. reallocating memory\n\n", *np);
    if ((CONFIG_num_auxiliary = auxiliary->n_aux)) {
        for (k = 0; k < CONFIG_num_auxiliary; k++) {
            strncpy(CONFIG_auxiliary_name[k], auxiliary->name[k], TERMSIZE);
            strncpy(CONFIG_auxiliary_unit[k], auxiliary->unit[k], TERMSIZE);
        }
    }

    p3dp_config_add_image(Config_Alib_to_Alib);

    if (info) {
        Config_analyze_H       (Config_Alib_to_Alib, info); fcr(info);
        Config_analyze_density (Config_Alib_to_Alib, info); fcr(info);
        Config_analyze_species (Config_Alib_to_Alib, info); fcr(info);
    }

    return(CONFIG_CFG_LOADED);
}

void p3dp_config_scatter(Alib_Declare_Config)
{
    int mg, ml, k;
    Auxiliary auxiliary;

    if (!p3dp_enabled)
        return;

    p3d_set_atoms(p3dp_cell, p3d_atom_begin(p3dp_cell), 0);
    auxiliary = p3dp_cell->auxiliary;
    if (!p3d_rank(p3dp_cell)) {
        int n;
        P3DAtom *atoms = (P3DAtom*) malloc(*np * sizeof(P3DAtom)), *p = atoms;

        if (!atoms) p3d_abort(-1);
        M3EQV(H, p3dp_cell->h);
        for (n = 0; n < *np; n++,p++) {
            p->m = (*mass)[n];
            strncpy(p->sym, SYMBOL(n), 3);
            p->r[0] = (*s )[n * DIMENSION    ];
            p->r[1] = (*s )[n * DIMENSION + 1];
            p->r[2] = (*s )[n * DIMENSION + 2];
            for (k = 0; k < CONFIG_num_auxiliary; k++)
                p->aux[k] = CONFIG_auxiliary[k][n];
        }
        p3d_set_atoms(p3dp_cell, atoms, *np);
        free(atoms);

        if ((auxiliary->n_aux = CONFIG_num_auxiliary)) {
            for (k = 0; k < CONFIG_num_auxiliary; k++) {
               strncpy(auxiliary->name[k],CONFIG_auxiliary_name[k],P3D_LINEMAX);
               strncpy(auxiliary->unit[k],CONFIG_auxiliary_unit[k],P3D_LINEMAX);
            }
        }
    }

    p3d_bcast(p3dp_cell, p3dp_cell->h, 9, MPI_DOUBLE, 0);
    do {
        ml = p3d_reset_cell(p3dp_cell);
        p3d_reduce(p3dp_cell, &ml, &mg, 1, MPI_INT, MPI_SUM);
    } while (mg);
    M3EQV(p3dp_cell->h, H);

    p3d_bcast(p3dp_cell, &CONFIG_num_auxiliary, 1, MPI_INT, 0);
    if (CONFIG_num_auxiliary) {
        p3d_bcast(p3dp_cell, auxiliary,sizeof(struct Auxiliary_tag),MPI_CHAR,0);
        if (p3d_rank(p3dp_cell)) {
            for (k = 0; k < CONFIG_num_auxiliary; k++) {
                strncpy(CONFIG_auxiliary_name[k], auxiliary->name[k], TERMSIZE);
                strncpy(CONFIG_auxiliary_unit[k], auxiliary->unit[k], TERMSIZE);
            }
        }
    }

    p3dp_config_add_image(Config_Alib_to_Alib);
}


int p3dp_Config_Load(char *fname, FILE *info, Alib_Declare_Config)
{
    char *p, *q;

    if (strcmp(config_fname, "/dev/null") == 0) {
        if (IS_MANAGER)
            cui_config_load_A3(Config_Alib_to_Alib);
        if (p3dp_enabled)
            p3dp_config_scatter(Config_Alib_to_Alib);
        return CONFIG_CFG_LOADED;
    }
    else if (!p3dp_enabled)
        return Config_Load(fname, info, Config_Alib_to_Alib);

    if (!IS_MANAGER)
        info = NULL;

    Fprintf(info, "Guessing file format of \"%s\":\n", fname);
    p = strrchr(fname, '/');
    if (p == NULL) p = fname;
    q = strchr(p, '.');
    if (q != NULL) {
        p = strrstr(q, "cfg");
        if (p != NULL) {
            if (strstr(p, "pdb") || strstr(p, "PDB"))
                goto PDB;
        }
        else if (strstr(q, "pdb") || strstr(q, "PDB"))
            goto PDB;
    }
    Fprintf(info, "should be Ju Li's CFG format.\n\n");
    p3dp_Config_load(fname, info, Config_Alib_to_Alib);
    return CONFIG_CFG_LOADED;
  PDB:
    Fprintf(info, "should be Protein Data Bank format.\n\n");

    if (IS_MANAGER) {
        Config_load_from_pdb(fname, info, Config_Alib_to_Alib);
        Config_to_bounding_box_config(Config_Alib_to_Alib,info);
    }
    p3dp_config_scatter(Config_Alib_to_Alib);

    return CONFIG_CFG_LOADED; /*return CONFIG_PDB_LOADED;*/
} /* end Config_Load() */



void p3dp_CalculateSimpleStatistics
(int N, char *series, int bytes_separation, int valtype, SimpleStatistics *s)
{
    int i, nn;
    char *p;
    double val;

    if (!p3dp_enabled) {
        CalculateSimpleStatistics(N, series, bytes_separation, valtype, s);
        return;
    }

    nn = p3d_n_atoms(p3dp_cell);
    p3d_reduce(p3dp_cell, &nn, &N, 1, MPI_INT, MPI_SUM);

    if (N < 1)
    {
        s->N = 0;
        s->idx_min = -1;
        s->min = 0;
        s->idx_max = -1;
        s->max = 0;
        s->average = 0;
        s->variance = 0;
        s->standard_deviation = 0;
        return;
    }

    s->N = N;
    s->idx_min = 0;
    s->min = IOVAL(series,valtype);
    s->idx_max = 0;
    s->max = IOVAL(series,valtype);
    s->average = 0;
    s->variance = 0;

    for (i=0,p=series; i<nn; i++,p+=bytes_separation)
    {
        val = IOVAL(p,valtype);
        if (val < s->min)
        {
            s->idx_min = i;
            s->min = val;
        }
        if (val > s->max)
        {
            s->idx_max = i;
            s->max = val;
        }
        s->average += val;
        s->variance += val*val;
    }

    {
        double og[2], ol[2], min_l, max_l;
        ol[0] = s->average;
        ol[1] = s->variance;
        p3d_reduce(p3dp_cell, ol, og, 2, MPI_DOUBLE, MPI_SUM);
        s->average  = og[0];
        s->variance = og[1];
        min_l = s->min;
        p3d_reduce(p3dp_cell, &min_l, &s->min, 1, MPI_DOUBLE, MPI_MIN);
        if (s->min != min_l)
            s->idx_min = -1;
        max_l = s->max;
        p3d_reduce(p3dp_cell, &max_l, &s->max, 1, MPI_DOUBLE, MPI_MAX);
        if (s->max != max_l)
            s->idx_max = -1;
    }

    s->average /= N;
    s->variance -= N*s->average*s->average;
    if (N>1) s->variance /= (N-1);
    else s->variance = 0;
    s->standard_deviation = sqrt(s->variance);
}


void p3dp_bcast_int2(int *val0, int *val1)
{
    if (p3dp_enabled) {
        int val[2];
        val[0] = *val0;
        val[1] = *val1;
        p3d_bcast(p3dp_cell, val, 2, MPI_INT, 0);
        *val0 = val[0];
        *val1 = val[1];
    }
}


void p3dp_zmap_reduce(int iw)
{
    if (p3dp_enabled) {
        static int size_zmap_mem_uc[AX_MAXWIN] = {0};
        static struct ax_float_int {
            AX_Float zmap;
            int mem_uc;
        } *zmap_mem_uc_l[AX_MAXWIN] = {0}, *zmap_mem_uc_g[AX_MAXWIN] = {0};
        int i, n = AX_size[iw].width*AX_size[iw].height;
        /*if (AX_bytes > sizeof(int)) p3d_abort(-1);*/

        if (n > size_zmap_mem_uc[iw]) {
            int size = n * sizeof(struct ax_float_int);
            if (!(zmap_mem_uc_l[iw] =
                    (struct ax_float_int*) realloc(zmap_mem_uc_l[iw], size)) ||
                !(zmap_mem_uc_g[iw] =
                    (struct ax_float_int*) realloc(zmap_mem_uc_g[iw], size)))
                p3d_abort(-2);
            size_zmap_mem_uc[iw] = n;
        }

        for(i = 0; i < n; i++ ) {
            zmap_mem_uc_l[iw][i].zmap = AX_zmap[iw][i];
            zmap_mem_uc_l[iw][i].mem_uc = *(int *)&AX_mem[iw].uc[AX_bytes*i];
        }
        p3d_reduce(p3dp_cell, zmap_mem_uc_l[iw], zmap_mem_uc_g[iw],
                                        n, MPI_AX_FLOAT_INT, MPI_MINLOC);
        if (!p3d_rank(p3dp_cell)) {
            for (i = 0; i < n; i++)
                *(int *)&AX_mem[iw].uc[AX_bytes*i]= zmap_mem_uc_g[iw][i].mem_uc;
        }
        for(i = 0; i < n; i++)
            AX_zmap[iw][i] = zmap_mem_uc_g[iw][i].zmap;
    }
}



void p3dp_atom_stack_insert(int iw, int atom, int rank)
{
    int i;
    for (i=ATOM_STACK_SIZE-1; i>0; i--)
        n[iw].atom_stack[i] = n[iw].atom_stack[i-1];
    n[iw].atom_stack[0] = atom;
    if (p3dp_enabled) {
        for (i=ATOM_STACK_SIZE-1; i>0; i--) {
            p3dp_n[iw].rank_stack[i] = p3dp_n[iw].rank_stack[i-1];
        }
        p3dp_n[iw].rank_stack[0] = rank;
    }
    return;
}


/* AX */

int p3dp_AXNextEvent(int iw)
{
    if (p3dp_enabled) {
        int value[AX_MAXWIN];
        if (!p3d_rank(p3dp_cell))
            value[iw] = AXNextEvent(iw);
        p3d_bcast(p3dp_cell, &AX_event[iw], sizeof(XEvent), MPI_BYTE, 0);
        p3d_bcast(p3dp_cell, &value[iw], 1, MPI_INT, 0);
        return value[iw];
    }
    else
        return AXNextEvent(iw);
}


int p3dp_AXQLength(int iw)
{
    if (p3dp_enabled) {
        int value[AX_MAXWIN];
        if (!p3d_rank(p3dp_cell))
            value[iw] = AXQLength(iw);
        p3d_bcast(p3dp_cell, &value[iw], 1, MPI_INT, 0);
        return value[iw];
    }
    else
        return AXQLength(iw);
}

int p3dp_AXPending(int iw)
{
    if (p3dp_enabled) {
        int value[AX_MAXWIN];
        if (!p3d_rank(p3dp_cell))
            value[iw] = AXPending(iw);
        p3d_bcast(p3dp_cell, &value[iw], 1, MPI_INT, 0);
        return value[iw];
    }
    else
        return AXPending(iw);
}


Status p3dp_AXGetGeometry(int iw, AXSize *newsize)
{
    if (p3dp_enabled) {
        Status value[AX_MAXWIN];
        if (!p3d_rank(p3dp_cell))
            value[iw] = AXGetGeometry(iw, (*newsize));
        p3d_bcast(p3dp_cell, newsize,    sizeof(AXSize), MPI_BYTE, 0);
        p3d_bcast(p3dp_cell, &value[iw], sizeof(Status), MPI_BYTE, 0);
        return value[iw];
    }
    else
        return AXGetGeometry(iw, (*newsize));
}


Bool p3dp_AXQueryPointer(int iw)
{
    if (p3dp_enabled) {
        Bool value[AX_MAXWIN];
        if (!p3d_rank(p3dp_cell))
            value[iw] = AXQueryPointer(iw);
        p3d_bcast(p3dp_cell, &value[iw], sizeof(Bool), MPI_CHAR, 0);
        return value[iw];
    }
    else
        return AXQueryPointer(iw);
}


KeySym p3dp_AXKeysym0(int iw)
{
    if (p3dp_enabled) {
        KeySym value[AX_MAXWIN];
        if (!p3d_rank(p3dp_cell))
            value[iw] = AXKeysym0(iw);
        p3d_bcast(p3dp_cell, &value[iw], sizeof(KeySym), MPI_CHAR, 0);
        return value[iw];
    }
    else
        return AXKeysym0(iw);
}


/* Atoms/Neighborlist.c */

#undef  NEIGHBORLIST_BINATOM_RATIO
#define NEIGHBORLIST_BINATOM_RATIO P3DP_NEIGHBORLIST_BINATOM_RATIO
#undef  NEIGHBORLIST_BINDEXS
#define NEIGHBORLIST_BINDEXS(N,s,i) p3dp_neighborlist_bindexs(N,s,i)
#define LIST_DISPOSABLE(N) \
  (((N)->min_shrinkage == 1) && ((N)->max_atom_displacement == 0))
static int nbin_g[3], nbin_o[3];
static int p3dp_neighborlist_bindexs(Neighborlist *N, double *s, int i)
{
    int j, ii[3];
    for (j = 0; j < DIMENSION; j++) {
        ii[j] = INT(s[DIMENSION*i+j]*nbin_g[j]) - nbin_o[j]; 
        if      (ii[j] < 0)              ii[j] = N->nbin[j] - 1;
        else if (ii[j] > N->nbin[j] - 1) ii[j] = 0;
    }
    return NEIGHBORLIST_BINDEX(N, ii[0], ii[1], ii[2]);
}

void p3dp_Neighborlist_Recreate
(Alib_Declare_Config, FILE *info, Chemtab *ct, Tp **tp, Neighborlist *N)
{
    register int i, m;
    int j, k, di, dj, dk, ii, jj, kk, l, n;
    double rmax, thickness[DIMENSION], ds[DIMENSION], dx[DIMENSION];

    /* Verify the only assumption of this program */
    for (i=0; i<*np; i++)
    {
        /* if (IN( (*s)[DIMENSION*i],   -EPS, 0 )) (*s)[DIMENSION*i]   = 0; */
        /* if (IN( (*s)[DIMENSION*i+1], -EPS, 0 )) (*s)[DIMENSION*i+1] = 0; */
        /* if (IN( (*s)[DIMENSION*i+2], -EPS, 0 )) (*s)[DIMENSION*i+2] = 0; */
        if ( V3NEED_TRIM(&((*s)[DIMENSION*i])) )
        {
            Fprintf(info, "** %s atom (1-%d) has\n""** s = (%g %g %g)\n",
                    word_for_order(i+1), *np, (*s)[DIMENSION*i],
                    (*s)[DIMENSION*i+1], (*s)[DIMENSION*i+2]);
            V3mM3( &((*s)[DIMENSION*i]), H, dx );
            V3MuL( ulength_IN_A, dx );
            Fprintf(info, "** or x = (%g %g %g) A.\n\n", dx[0], dx[1], dx[2]);
            Fprintf(info, "This is impermissible for constructing "
                    "neighborlist,\n");
            if (N->s_overflow_err_handler ==
                NEIGHBORLIST_S_OVERFLOW_ERR_HANDLER_QUIT)
                pe ("Neighborlist_Recreate() exiting...\n");
            else if (N->s_overflow_err_handler ==
                     NEIGHBORLIST_S_OVERFLOW_ERR_HANDLER_FOLD_INTO_PBC)
            { /* even though the user may not mean PBC */
                Fprintf(info, "so we have to fold all atoms into [0,1).\n\n");
                Config_fold_into_PBC (Config_Alib_to_Alib);
                break;
            }
            else  /* NEIGHBORLIST_S_OVERFLOW_ERR_HANDLER_BOUNDING_BOX */
            { /* even though the user may mean PBC */
                Fprintf(info, "so we use BOUNDING BOX as H[][] instead:\n\n");
                Config_to_bounding_box_config(Config_Alib_to_Alib,info);
                Fcr(info);
                break;
            }
        }
    }

    /* Also check RCUT symmetry in the case of pairwise interaction */
    if (N->pairwise)
        for (i=0; i<ct->t; i++)
            for (j=i+1; j<ct->t; j++)
                if (NEIGHBOR_TABLE(N->rcut,ct,i,j) !=
                    NEIGHBOR_TABLE(N->rcut,ct,j,i))
                {
                    Fprintf(info, "Warning: rcut[%d][%d] != "
                            "rcut[%d][%d]\n", i,j, j,i);
                    Fprintf(info, "-> %f will be used for both.\n\n",
                            NEIGHBOR_TABLE(N->rcut,ct,i,j)=
                            NEIGHBOR_TABLE(N->rcut,ct,j,i)=
                            MAX(NEIGHBOR_TABLE(N->rcut,ct,i,j),
                                NEIGHBOR_TABLE(N->rcut,ct,j,i)));
                }

    /* Print out the form */
    Fprintf(info, "Create neighborlist with the following parameters:\n");
    Fprintf(info, "pairwise saving = %s,  track_dx = %s\n",
            BWORD(N->pairwise), BWORD(N->track_dx));
    for (i=0; i<DIMENSION; i++)
        Fprintf(info, "pbc[%d] = %s,  ", i, BWORD(N->pbc[i]));
    Fcr(info);
    Fprintf(info, "Strain Session min_shrinkage = %g,\n"
            "All Atom Tether max_atom_displacement = %g (reduced),\n",
            N->min_shrinkage, N->max_atom_displacement);
    if ( (N->min_shrinkage > 1) || (N->min_shrinkage <= 0) )
        pe ("Neighborlist_Recreate: min_shrinkage must be (0,1].\n\n");
    if ( N->max_atom_displacement < 0 )
        pe ("Neighborlist_Recreate: max_atom_displacement must be >= 0.\n\n");
        
    REALLOC( Neighborlist_Recreate_Form, N->rlist, SQUARE(ct->t), double );
    if (LIST_DISPOSABLE(N))
    {
        Fprintf(info, "-> COMPRESSED disposable atom-atom list\n"
                "with Rlist_ij = Rcut_ij;\n\n");
        VEQV( SQUARE(ct->t), N->rcut, N->rlist );
    }
    else
    {
        Fprintf(info, "-> UNCOMPRESSED reusable atom-atom list\n"
                "with %g x Rlist_ij = Rcut_ij + 2 x %g;\n\n",
                N->min_shrinkage, N->max_atom_displacement);
        if (info) Sfpr(info, "Rcut  = %M (reduced)\n ", N->rcut, ct->t,ct->t);
        for (i=0; i<ct->t; i++)
            for (j=0; j<ct->t; j++)
                NEIGHBOR_TABLE(N->rlist,ct,i,j) =
                    (NEIGHBOR_TABLE(N->rcut,ct,i,j) +
                     2 * N->max_atom_displacement) / N->min_shrinkage;
    }

    for (i=0; i<ct->t; i++)
        for (j=0; j<ct->t; j++)
            if (NEIGHBOR_TABLE(N->rcut,ct,i,j) <= 0)
                NEIGHBOR_TABLE(N->rlist,ct,i,j) = 0;

    REALLOC( Neighborlist_Recreate_Form, N->rlist2, SQUARE(ct->t), double );
    for (i=0; i<ct->t; i++)
        for (j=0; j<ct->t; j++)
            NEIGHBOR_TABLE(N->rlist2,ct,i,j) =
                SQUARE(NEIGHBOR_TABLE(N->rlist,ct,i,j));

    if (info) Sfpr(info, "rlist = %M (reduced)\n ", N->rlist, ct->t,ct->t);
    rmax = VMAX (SQUARE(ct->t), N->rlist);
    Fprintf(info, "MAX(rlist) = %g  ->  bin thicknesses should >= %g\n",
            rmax, rmax);

  restart:
    M3rowthicknesses (H, thickness);
    if (info) Vfpr(info, "since H thicknesses = %M", thickness, DIMENSION);

    if ( N->small_cell_err_handler ==
         NEIGHBORLIST_SMALL_CELL_ERR_HANDLER_MULTIPLY )
    {
        ii = ceil(2*rmax / thickness[0]);
        jj = ceil(2*rmax / thickness[1]);
        kk = ceil(2*rmax / thickness[2]);
        if ( (ii > 1) || (jj > 1) || (kk > 1) )
        {
            Config_multiply ( ii, jj, kk, Config_Alib_to_Alib );
            rebind_CT (Config_Alib_to_Alib, "", ct, tp);
            goto restart;
        }
    }

    N->nbin[DIMENSION] = 1;
    for (i=0; i<DIMENSION; i++)
    {
        if ( N->pbc[i] && (thickness[i] < 2*rmax) &&
             (N->small_cell_err_handler !=
              NEIGHBORLIST_SMALL_CELL_ERR_HANDLER_NOCHECK) )
            pe ("Neighborlist_Recreate: losing r<rmax images may happen\n"
                "in direction %d, please increase H[%d][] "
                "(rmax=%g,thickness=%g).\n", i, i, rmax,thickness[i]);
        N->nbin[i] = thickness[i] / rmax / p3dp_cell->private->dims[i];
        if (N->nbin[i] < 1) N->nbin[i] = 1;
        if (N->nbin[i] > NEIGHBORLIST_MAXBIN) N->nbin[i]=NEIGHBORLIST_MAXBIN;
        nbin_g[i] = N->nbin[i] * p3dp_cell->private->dims[i];
        nbin_o[i] = N->nbin[i] * p3dp_cell->private->coords[i];
        N->nbin[DIMENSION] *= N->nbin[i];
    }
    Fprintf(info, "-> %d x ", N->nbin[0]);
    for (i=1; i<DIMENSION-1; i++)
        Fprintf(info, "%d x ", N->nbin[i]);
    Fprintf(info, "%d = %d bins.\n", N->nbin[DIMENSION-1],
            N->nbin[DIMENSION]);
    
    /* bin-bin: prepare for maximal bin connectivity */
    N->binbinlist = Irecreatelist
        (&N->binbindx, N->nbin[DIMENSION], DIMENSION3);
    for (i=0; i<N->nbin[0]; i++)
        for (j=0; j<N->nbin[1]; j++)
            for (k=0; k<N->nbin[2]; k++)
            {
                l = NEIGHBORLIST_BINDEX(N,i,j,k);
                for (di=-1; di<=1; di++)
                    for (dj=-1; dj<=1; dj++)
                        for (dk=-1; dk<=1; dk++)
                        {
                            ii = i + di;
                            jj = j + dj;
                            kk = k + dk;
                            if (ii >= N->nbin[0])
                            {
                                if (N->pbc[0]) ii-=N->nbin[0]; else continue;
                            }
                            if (ii < 0)
                            {
                                if (N->pbc[0]) ii+=N->nbin[0]; else continue;
                            }
                            if (jj >= N->nbin[1])
                            {
                                if (N->pbc[1]) jj-=N->nbin[1]; else continue;
                            }
                            if (jj < 0)
                            {
                                if (N->pbc[1]) jj+=N->nbin[1]; else continue;
                            }
                            if (kk >= N->nbin[2])
                            {
                                if (N->pbc[2]) kk-=N->nbin[2]; else continue;
                            }
                            if (kk < 0)
                            {
                                if (N->pbc[2]) kk+=N->nbin[2]; else continue;
                            }
                            m = NEIGHBORLIST_BINDEX(N,ii,jj,kk);
                            /* make sure it is a new record */
                            for (n=N->binbindx[2*l]; n<N->binbindx[2*l+1]; n++)
                                if (m == N->binbinlist[n]) break;
                            if (n == N->binbindx[2*l+1])
                                N->binbinlist[N->binbindx[2*l+1]++] = m;
                        }
            }
    Fprintf(info, "bin-bin list created.\n\n");

    /* it is reasonable to assume that bin-connectivity never changes */
    N->binbinlist = Icompresslist (&N->binbindx, N->nbin[DIMENSION]);
    Ilist_statistics ("bin-bin", N->binbindx, N->nbin[DIMENSION], TRUE, info);
    Fcr(info);

    Fprintf(info, "On average, there are %.2f atoms per bin,\n"
            "but we will offer space of %d for each.\n",
            DOUBLE(*np) / N->nbin[DIMENSION], n = ceil
            (DOUBLE(*np) / N->nbin[DIMENSION] * NEIGHBORLIST_BINATOM_RATIO));
    N->binlist = Irecreatelist (&N->bindx, N->nbin[DIMENSION], n);
    REALLOC( Neighborlist_Recreate, N->mybin, *np, int );

    /* bin-atom: fill atoms into bins */
    for (i=0; i<*np; i++)
        Iappend(N->bindx, N->binlist, i, N->mybin[i]=
                NEIGHBORLIST_BINDEXS(N,*s,i), 0, N->nbin[DIMENSION]);
    Fprintf(info, "bin-atom list created.\n\n");

    Ilist_statistics("bin-atom", N->bindx, N->nbin[DIMENSION], FALSE, info);
    Fcr(info);

    /* atom-atom */
    for (n=i=0; i<ct->t; i++)
    {
        if ( N->pairwise && (!N->keepcounter) )
            N->neighbormax[i] = N->neighbormax[i] / 2;
        n += N->neighbormax[i] * ct->count[i];
    }
    if (N->pairwise) Fprintf(info, "After pairwise/own_pair saving,\n");
    if (info)
        ifpr(info, "likely MAX(neighbor) = %M atoms,", N->neighbormax, ct->t);
    REALLOC( Neighborlist_Recreate, N->idx, 2*(*np)+1+n, int );
    N->idx[0] = N->idx[1] = 0;
    for (i=1; i<*np; i++)
        N->idx[2*i] = N->idx[2*i+1] =
            N->idx[2*i-1]+N->neighbormax[(int)(*tp)[i-1]];
    if ( (N->idx[2*(*np)]=n) != N->idx[2*i-1]+N->neighbormax[(int)(*tp)[i-1]] )
        pe ("Neighborlist_Recreate: checksum failed.\n");
    N->list = &N->idx[2*(*np)] + 1;
    for (i=0; i<*np; i++)
    { /* already compressed bin-bin list */
        ii = N->binbindx[N->mybin[i]];
        jj = N->binbindx[N->mybin[i]+1];
        for (j=ii; j<jj; j++)
        { /* neighboring bins */
            l = N->binbinlist[j];
            for (k=N->bindx[2*l]; k<N->bindx[2*l+1]; k++)
            { /* particles in neighboring bins */
                m = N->binlist[k];
                if (N->pairwise && (!own_pair(i,m))) continue;
                if (i == m) continue;
                V3SUB ( &((*s)[DIMENSION*m]), &((*s)[DIMENSION*i]), ds);
                if (N->pbc[0]) Image(ds[0]);
                if (N->pbc[1]) Image(ds[1]);
                if (N->pbc[2]) Image(ds[2]);
                V3mM3 (ds, H, dx);
                if (V3LENGTH2(dx) <=
                    NEIGHBOR_TABLE(N->rlist2,ct,(*tp)[i],(*tp)[m]))
                    Iappend (N->idx, N->list, m, i, 0, *np);
            }
        }
    }
    Fprintf(info, "atom-atom list created.\n\n");
    Ilist_statistics ("atom-atom", N->idx, *np, FALSE, info);
    Fcr(info);

    M3EQV (H, N->H0);
    if (LIST_DISPOSABLE(N))
    {
        N->list = Icompresslist (&N->idx, *np);
        Ilist_statistics ("atom-atom", N->idx, *np, TRUE, info);
        Fcr(info);
        Ifreelist (&N->binbindx, &N->binbinlist);
        Ifreelist (&N->bindx, &N->binlist);
        Free (N->mybin);
        Fprintf(info, "All bin-related allocations freed.\n");
    }
    else  /* LIST_REUSABLE(N)) */
    {
        if (!N->keepcounter) N->remesh_counter = 0;
        REALLOC (Neighborlist_Recreate, N->sa, DIMENSION*(*np), double);
        VEQV (DIMENSION*(*np), (*s), N->sa);
        Fprintf (info, "Particle anchors established.\n");
        REALLOC (Neighborlist_Recreate, N->maxtether2, ct->t, double);
        VMaX (ct->t, N->neighbormax, i, N->stack_max);
        N->stack_max *= NEIGHBORLIST_STACK_MAX_RATIO;
        REALLOC (Neighborlist_Recreate, N->stack, N->stack_max, int);
        if (!N->keepcounter) N->reanchor_counter = 0;
        if (N->track_dx)
        {
            REALLOC (Neighborlist_Recreate, N->dx, 2*DIMENSION*(*np), double);
            N->dxa = N->dx + DIMENSION*(*np);
            if (!N->keepcounter)
            {
                VZERO (2*DIMENSION*(*np), N->dx);
                Fprintf (info, "dx tracking initiated, counters cleared.\n");
            }
        }
        if (!N->keepcounter) N->maintenance_counter = 0;
    }
    N->keepcounter = TRUE;
    return;
}
