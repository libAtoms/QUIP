/* $Id: P3DCore.c,v 1.10 2005/08/19 02:42:36 shimizu Exp $
 * 2004-2005 Futoshi SHIMIZU
 */

#define P3D_GLOBAL
#include "P3D_p.h"
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>
#include <stdarg.h>

enum p3d_index_flags
    direction_flags[] = {F_EOO, F_WOO, F_ONO, F_OSO, F_OOT, F_OOB};

#define P3D_DECOMP_KEY "-decomp="

static int initialized;

int p3d_init(int *argc, char ***argv, MPI_Comm comm)
{
    int flag, *attr_val; 

    p3d_global_rank_fprintf = 0;

    MPI_Initialized(&initialized);
    if (!initialized) MPI_Init(argc, argv);

    if (comm == MPI_COMM_NULL) comm = MPI_COMM_WORLD;
    MPI_Comm_dup(comm, &p3d_global_comm);
    MPI_Comm_size(p3d_global_comm, &p3d_global_nprocs);
    MPI_Comm_rank(p3d_global_comm, &p3d_global_rank);

    MPI_Attr_get(MPI_COMM_WORLD, MPI_IO, &attr_val, &flag);
    p3d_global_io_any = (flag && *attr_val == MPI_ANY_SOURCE);

    if (!(p3d_global_decomp_string = (char*) malloc(1))) p3d_abort(-1);
    *p3d_global_decomp_string = '\0';

    if (argv && *argv) {
        int i, len, keylen = strlen(P3D_DECOMP_KEY);
	for (i = 1; i < *argc; i++) {
            if (strncmp((*argv)[i], P3D_DECOMP_KEY, keylen) == 0) {
                if ((len = strlen((*argv)[i])) > keylen) {
                    p3d_global_decomp_string =
                        (char*)realloc(p3d_global_decomp_string, len-keylen+1);
                    if (p3d_global_decomp_string)
                        strcpy(p3d_global_decomp_string, (*argv)[i] + keylen);
                }
                (*argv)[i] = 0;
            }
        }
        { /* mpich: MPID_ArgSqueeze */
            int j = i = 0;
            while (j < *argc) {
	        while ((*argv)[j] == 0 && j < *argc) j++;
	        if (j < *argc) (*argv)[i++] = (*argv)[j++];
            }
            if (!(*argv)[i-1]) i--;
            *argc = i;
        }
    }

    return 0;
}


int p3d_finalize(void)
{
    MPI_Comm_free(&p3d_global_comm);
    if (!initialized) MPI_Finalize();
    return 0;
}


MPI_Comm p3d_comm(Cell cell) { return cell->private->comm; }
int p3d_nprocs(Cell cell) { return cell->private->nprocs; }
int p3d_rank(Cell cell) { return cell->private->rank; }
int p3d_n_atoms(Cell cell) { return cell->private->point_end
                                - cell->private->point_begin; }
int p3d_n_images(Cell cell) { return cell->private->image_end
                                - cell->private->point_end; }
P3DAtom *p3d_atom_begin(Cell cell) { return cell->private->point_begin; }
P3DAtom *p3d_atom_end(Cell cell) { return cell->private->point_end; } 
P3DAtom *p3d_image_begin(Cell cell) { return cell->private->point_end; } 
P3DAtom *p3d_image_end(Cell cell) { return cell->private->image_end; }


static void p3d_realloc_cell(Cell cell, int n_max, int p_max)
{
    Cell_private p = cell->private;
    /*int i;*/

    if (n_max > p->n_max) {
        P3DAtom *point_begin_org = p->point_begin;

        /*p3d_fprintf(stderr, "n_max: %d -> %d\n", p->n_max, n_max);*/
        p->point_begin=(P3DAtom*)realloc(p->point_begin, n_max*sizeof(P3DAtom));
        if (p->point_begin == NULL) p3d_abort(1);
        p->n_max = n_max;

        if ((char*)p->point_begin != (char*)point_begin_org) {
            P3DAtom **pp;
            p->point_end = p->point_begin + (p->point_end - point_begin_org);
            p->image_end = p->point_begin + (p->image_end - point_begin_org);
            for (pp = p->p_point_begin; pp < p->p_point_end; pp++)
                *pp =      p->point_begin + (*pp          - point_begin_org);
        }

        p->f_realloc = 1;
    }
    if (p_max > p->p_max) {
        P3DAtom **p_point_begin_org = p->p_point_begin;

        /*p3d_fprintf(stderr, "p_max: %d -> %d\n", p->p_max, p_max);*/
        p->p_point_begin =
            (P3DAtom**) realloc(p->p_point_begin, p_max*sizeof(P3DAtom*));
        if (p->point_begin == NULL) p3d_abort(1);
        p->p_max = p_max;

        if ((char*)p->p_point_begin != (char*)p_point_begin_org) {
            p->p_point_end=p->p_point_begin+(p->p_point_end-p_point_begin_org);
            /*
            for (i = 0; i <= 6; i++)
                p->update[i]=p->p_point_begin+(p->update[i]-p_point_begin_org);
            */
        }

        p->f_realloc = 1;
    }
}

static V3 v_eps = {1.e-7, 1.e-7, 1.e-7};

static int p3d_decompose(Cell cell, char *decomp, MPI_Comm comm)
{
    Cell_private p = cell->private;
    char spec[3][512] = {"0", "0", "0"};
    int i;

    
    for (i = 0; i < 6; i++) {
        if (p->req_send[i]!=MPI_REQUEST_NULL) MPI_Request_free(&p->req_send[i]);
        if (p->req_recv[i]!=MPI_REQUEST_NULL) MPI_Request_free(&p->req_recv[i]);
    }

    if (!decomp || !strlen(decomp)) decomp = p3d_global_decomp_string;
    {
        char *s = strstr(decomp, "io_any");
        if (s) {
            s += strlen("io_any");
            if (*s == '=')
                sscanf(s+1, "%d", &p3d_global_io_any);
            else
                p3d_global_io_any = 1;
        }
    }

    {
        int j, k;
        sscanf(decomp, "%[^,],%[^,],%[^,]", spec[0], spec[1], spec[2]);
        for (i = 0; i < 3; i++) {
            k = 0;
            for (j = 0; j < strlen(spec[i]); j++) {
                switch (spec[i][j]) {
                case '(':
                case '/':
                case ')':
                    k++;
                    break;
                default:
                    break;
                }
            }
            if (k != 0) p->dims[i] = k;
            else        p->dims[i] = atoi(spec[i]);
        }
    }
    {
        int nprocs, free, c[3], periods[3] = {1, 1, 1};
        int *rank_send = p->rank_send, *rank_recv = p->rank_recv;

        if (comm == MPI_COMM_NULL) {
            if (p->comm != MPI_COMM_NULL) {
                comm = p->comm;
                free = 1;
            }
            else
                comm = p3d_global_comm;
        }
        free = (p->comm != MPI_COMM_NULL);

        MPI_Comm_size(comm, &nprocs);
        MPI_Dims_create(nprocs, 3, p->dims);
        MPI_Cart_create(comm, 3, p->dims, periods, 1, &p->comm);
        if (free) MPI_Comm_free(&comm);

        MPI_Comm_size(p->comm, &p->nprocs);
        MPI_Comm_rank(p->comm, &p->rank);
        MPI_Cart_coords(p->comm, p->rank, 3, p->coords);

        p->trans_mask = F_OOO;
        if (p->coords[0] == p->dims[0] - 1) p->trans_mask |= F_EOO;
        if (p->coords[0] ==              0) p->trans_mask |= F_WOO;
        if (p->coords[1] == p->dims[1] - 1) p->trans_mask |= F_ONO;
        if (p->coords[1] ==              0) p->trans_mask |= F_OSO;
        if (p->coords[2] == p->dims[2] - 1) p->trans_mask |= F_OOT;
        if (p->coords[2] ==              0) p->trans_mask |= F_OOB;

        c[1] = p->coords[1];        c[2] = p->coords[2];
        c[0] = p->coords[0] + 1;    MPI_Cart_rank(p->comm, c, &rank_send[0]);
        c[0] = p->coords[0] - 1;    MPI_Cart_rank(p->comm, c, &rank_send[1]);
        c[0] = p->coords[0];        c[2] = p->coords[2];
        c[1] = p->coords[1] + 1;    MPI_Cart_rank(p->comm, c, &rank_send[2]);
        c[1] = p->coords[1] - 1;    MPI_Cart_rank(p->comm, c, &rank_send[3]);
        c[0] = p->coords[0];        c[1] = p->coords[1];
        c[2] = p->coords[2] + 1;    MPI_Cart_rank(p->comm, c, &rank_send[4]);
        c[2] = p->coords[2] - 1;    MPI_Cart_rank(p->comm, c, &rank_send[5]);
        rank_recv[0] = rank_send[1]; rank_recv[1] = rank_send[0];
        rank_recv[2] = rank_send[3]; rank_recv[3] = rank_send[2];
        rank_recv[4] = rank_send[5]; rank_recv[5] = rank_send[4];
    }
    {
        double *lb = (double*) &p->lb, *ub = (double*) &p->ub;
        char *s;
        int j;
        for (i = 0; i < 3; i++) {
            s = strchr(spec[i], '(');
            if (s != NULL) {
                if (p->coords[i] == 0) {
                    lb[i] = 0.0;
                    sscanf(++s, "%lf", &ub[i]);
                }
                else {
                    for (j = 1; j < p->coords[i]; j++) s = strchr(++s, '/');
                    sscanf(++s, "%lf", &lb[i]);
                    if (p->coords[i] == p->dims[i] - 1)
                        ub[i] = 1.0;
                    else {
                        s = strchr(++s, '/');
                        sscanf(++s, "%lf", &ub[i]);
                    }
                }
            }
            else {
                double d = 1.0 / p->dims[i];
                lb[i] = d *  p->coords[i];
                ub[i] = d * (p->coords[i] + 1);
            }
        }
#ifndef P3D_QUIET
        {
            char line[P3D_LINEMAX];
            snprintf(line, sizeof(line),
                    "# %d: %d[%g,%g] %d[%g,%g] %d[%g,%g]\n",
                    p->rank,
                    p->coords[0], p->lb[0], p->ub[0],
                    p->coords[1], p->lb[1], p->ub[1],
                    p->coords[2], p->lb[2], p->ub[2]);
            if (p->rank)
                MPI_Send(line, sizeof(line), MPI_CHAR, 0, 0, p->comm);
            else {
                MPI_Status stat;
                int rank;
                fprintf(stderr, "# \"%s\" -> %d procs (%dx%dx%d)\n%s",
                    decomp, p->nprocs, p->dims[0], p->dims[1], p->dims[2],line);
                for (rank = 1; rank < p->nprocs; rank++) {
                    MPI_Recv(line, sizeof(line), MPI_CHAR,rank,0,p->comm,&stat);
                    fprintf(stderr, line);
                }
            }
        }
#endif
        V3SuB(lb, v_eps);
        V3AdD(v_eps, ub);
    }

    return 0;
}


Cell p3d_new_cell(char *decomp, MPI_Datatype m_image, MPI_Comm comm)
{
    Cell cell = (Cell) malloc(sizeof(struct Cell_tag));
    Cell_private p = (Cell_private) malloc(sizeof(struct Cell_private_tag));

    if (cell == NULL) p3d_abort(1);
    if (p == NULL) p3d_abort(1);
    cell->private = p;

    M3ASSIGN(1., 0., 0., 0., 1., 0., 0., 0., 1., cell->h);

    p->crust = 0.0;
    p->n_max = 0;
    p->point_begin = p->point_end = p->image_end = NULL;

    p->p_max = 0;
    p->p_point_begin = p->p_point_end = NULL;

    p->f_realloc = 1;

    {
        int i;
        for (i = 0; i <= 6; i++) p->update[i] = NULL;
    }

    MPI_Type_contiguous(sizeof(P3DAtom), MPI_CHAR, &p->M_POINT);
    MPI_Type_commit(&p->M_POINT);

    if (m_image != MPI_DATATYPE_NULL) {
        p->M_IMAGE = m_image;
    }
    else {
        int blocklengths[2] = {offsetof(P3DAtom, r) + sizeof(V3), 1};
        MPI_Datatype types[2] = {MPI_CHAR, MPI_UB};
        MPI_Aint displacements[2] = {0, sizeof(P3DAtom)};
        MPI_Type_struct(2, blocklengths, displacements, types, &p->M_IMAGE);
        MPI_Type_commit(&p->M_IMAGE);
    }

    {
        int i;
        for (i = 0; i < 6; i++) p->req_send[i]=p->req_recv[i]=MPI_REQUEST_NULL;
    }

    p->comm = MPI_COMM_NULL;
    p3d_decompose(cell, decomp, comm);

    cell->reduced_coordinates = 0;

    cell->auxiliary = (Auxiliary) malloc(sizeof(struct Auxiliary_tag));
    cell->auxiliary->n_aux = 0;

    return cell;
}


int p3d_set_atoms(Cell cell, P3DAtom *point, int n)
{
    Cell_private p = cell->private;

    if (p->point_begin != point) {
        p3d_realloc_cell(cell, n, 0);
        memcpy(p->point_begin, point, sizeof(P3DAtom) * n);
    }
    p->image_end = p->point_end = p->point_begin + n;

    return (p->point_begin == point);
}


int p3d_cat_atoms(Cell cell, P3DAtom *point, int n)
{
    Cell_private p = cell->private;
    int n_max = p->point_end - p->point_begin + n;

    p3d_realloc_cell(cell, n_max, 0);

    memcpy(p->point_end, point, sizeof(P3DAtom) * n);
    p->point_end += n;
    p->image_end = p->point_end;

    return n_max;
}


int p3d_remove_atom(Cell cell, P3DAtom *atom)
{
    Cell_private p = cell->private;

    if (atom - p->point_begin < 0 || atom - p->point_end >= 0)
        return -1;
    *atom = *(--p->point_end);
    p->image_end = p->point_end;

    return p->point_end - p->point_begin;
}


int p3d_flags_boundary(V3 s, V3 lb, V3 ub)
{
    int flags = F_OOO;
    if (s[0] <   lb[0]) flags |= F_WOO;
    if (s[0] >=  ub[0]) flags |= F_EOO;
    if (s[1] <   lb[1]) flags |= F_OSO;
    if (s[1] >=  ub[1]) flags |= F_ONO;
    if (s[2] <   lb[2]) flags |= F_OOB;
    if (s[2] >=  ub[2]) flags |= F_OOT;
    return flags;
}


int p3d_flags_boundary2(V3 s, V3 lb, V3 ub)
{
    int flags = F_OOO;

    if      (s[0] <  0.0) s[0] += 1.0;
    else if (s[0] >= 1.0) s[0] -= 1.0;
    if      (s[1] <  0.0) s[1] += 1.0;
    else if (s[1] >= 1.0) s[1] -= 1.0;
    if      (s[2] <  0.0) s[2] += 1.0;
    else if (s[2] >= 1.0) s[2] -= 1.0;

    if (s[0] <   lb[0]) flags |= F_WOO;
    if (s[0] >=  ub[0]) flags |= F_EOO;
    if (s[1] <   lb[1]) flags |= F_OSO;
    if (s[1] >=  ub[1]) flags |= F_ONO;
    if (s[2] <   lb[2]) flags |= F_OOB;
    if (s[2] >=  ub[2]) flags |= F_OOT;
    return flags;
}


int p3d_reset_cell_c(Cell cell, double crust)
{
    Cell_private p = cell->private;
    P3DAtom *pnt;
    M3 h_inv;
    int dir0, n_max, p_max, n_send_total = 0;

    p->p_point_end = p->p_point_begin;
    for (dir0 = 0; dir0 < 6; dir0++) {
        if (p->req_recv[dir0] != MPI_REQUEST_NULL)
            MPI_Request_free(&p->req_recv[dir0]);
        if (p->req_send[dir0] != MPI_REQUEST_NULL)
            MPI_Request_free(&p->req_send[dir0]);
    }

    M3inv(cell->h, h_inv);

    if (cell->reduced_coordinates) {
        for (pnt = p->point_begin; pnt < p->point_end; pnt++) {
            pnt->iw = p3d_flags_boundary(pnt->r, p->lb, p->ub);
        }
    }
    else {
        for (pnt = p->point_begin; pnt < p->point_end; pnt++) {
            V3 s;
            V3mM3(pnt->r, h_inv, s);
            pnt->iw = p3d_flags_boundary(s, p->lb, p->ub);
        }
    }

    for (dir0 = 0; dir0 < 6; dir0++) {
        int n_send, n_recv, dir_flags;
        P3DAtom *point_send;
        MPI_Status status;
        P3DAtom temp;

        dir_flags = direction_flags[dir0];

        n_send = 0;
        for (pnt = p->point_begin; pnt < p->point_end; pnt++)
            if (pnt->iw & dir_flags) n_send++;
        n_send_total += n_send;
        MPI_Sendrecv(   &n_send, 1, MPI_INT, p->rank_send[dir0], dir0,
                        &n_recv, 1, MPI_INT, p->rank_recv[dir0], dir0,
                        p->comm, &status);

        n_max = (p->point_end - p->point_begin) + n_recv;
        p3d_realloc_cell(cell, n_max, 0);

        point_send = p->point_begin + p->n_max;
        for (pnt = p->point_begin; pnt < p->point_end; pnt++) {
            if (pnt->iw & dir_flags) {
                temp = *pnt;
                *pnt-- = *--p->point_end;
                *--point_send = temp;
            }
        }
        if (dir_flags & p->trans_mask) { /* PBC */
            int idx = dir0 / 2;
            V3 trans;
            V3ZERO(trans);
            trans[idx] = 1.0;
            if (!cell->reduced_coordinates) {
                V3 tmp;
                V3EQV(trans, tmp);
                V3mM3(tmp, cell->h, trans);
            }
            if (dir0 % 2 == 0) {
                for (pnt = point_send; pnt < p->point_begin + p->n_max; pnt++)
                    V3SuB(pnt->r, trans);
            }
            else {
                for (pnt = point_send; pnt < p->point_begin + p->n_max; pnt++)
                    V3AdD(trans, pnt->r);
            }
        }
        MPI_Sendrecv(point_send, n_send, p->M_POINT, p->rank_send[dir0], dir0,
                   p->point_end, n_recv, p->M_POINT, p->rank_recv[dir0], dir0,
                   p->comm, &status);

        pnt = p->point_end;
        p->point_end += n_recv;
        if (cell->reduced_coordinates) {
            while (pnt < p->point_end) {
                pnt->iw = p3d_flags_boundary(pnt->r, p->lb, p->ub);
                pnt++;
            }
        }
        else {
            while (pnt < p->point_end) {
                V3 s;
                V3mM3(pnt->r, h_inv, s);
                pnt->iw = p3d_flags_boundary(s, p->lb, p->ub);
                pnt++;
            }
        }
    }

    p->image_end = p->point_end;
    if ((p->crust = crust)) {
        V3 lb_c, ub_c, c;
        P3DAtom *image_end, **p_point;
        int n_send[6], n_recv[6], max_nn_send = 0;

        V3ASSIGN(h_inv[0][0], h_inv[1][1], h_inv[2][2], c);
        V3MuL(crust, c);
        V3ADD(p->lb, c, lb_c); V3AdD(v_eps, lb_c); V3AdD(v_eps, lb_c);
        V3SUB(p->ub, c, ub_c); V3SuB(ub_c, v_eps); V3SuB(ub_c, v_eps);

        if (cell->reduced_coordinates) {
            for (pnt = p->point_begin; pnt < p->point_end; pnt++) {
                pnt->iw = p3d_flags_boundary(pnt->r, lb_c, ub_c);
            }
        }
        else {
            for (pnt = p->point_begin; pnt < p->point_end; pnt++) {
                V3 s;
                V3mM3(pnt->r, h_inv, s);
                pnt->iw = p3d_flags_boundary(s, lb_c, ub_c);
            }
        }

        for (dir0 = 0; dir0 < 6; dir0 += 2) {
            P3DAtom *point_send0, *point_send1;
            V3 trans;
            MPI_Request req_send[2], req_recv[2];
            MPI_Status status[2];
            int nn_send, nn_recv;
            int dir_flags0, dir_flags1, dir1 = dir0 + 1, idx = dir0 / 2;
        
        /* count */
            dir_flags0 = direction_flags[dir0];
            dir_flags1 = direction_flags[dir1];
            n_send[dir0] = n_send[dir1] = 0;
            for (pnt = p->point_begin; pnt < p->image_end; pnt++) {
                if (pnt->iw & dir_flags0) n_send[dir0]++;
                if (pnt->iw & dir_flags1) n_send[dir1]++;
            }
            nn_send = n_send[dir0] + n_send[dir1];
            if (nn_send > max_nn_send) max_nn_send = nn_send;
            MPI_Sendrecv(   &n_send[dir0], 1, MPI_INT, p->rank_send[dir0], dir0,
                            &n_recv[dir0], 1, MPI_INT, p->rank_recv[dir0], dir0,
                            p->comm, status);
            MPI_Sendrecv(   &n_send[dir1], 1, MPI_INT, p->rank_send[dir1], dir1,
                            &n_recv[dir1], 1, MPI_INT, p->rank_recv[dir1], dir1,
                            p->comm, status);
            nn_recv = n_recv[dir0] + n_recv[dir1];
            n_max = p->image_end - p->point_begin + max_nn_send + nn_recv;
            p_max = p->p_point_end - p->p_point_begin + nn_send;
            p3d_realloc_cell(cell, n_max, p_max);

        /* sendrecv */
            image_end = p->image_end;
            if (n_recv[dir0] > 0) {
                MPI_Irecv(image_end, n_recv[dir0], p->M_POINT,
                        p->rank_recv[dir0], dir0, p->comm, &req_recv[0]);
                image_end += n_recv[dir0];
            }
            else
                req_recv[0] = MPI_REQUEST_NULL;
            if (n_recv[dir1] > 0) {
                MPI_Irecv(image_end, n_recv[dir1], p->M_POINT,
                        p->rank_recv[dir1], dir1, p->comm, &req_recv[1]);
                image_end += n_recv[dir1];
            }
            else
                req_recv[1] = MPI_REQUEST_NULL;

            V3ZERO(trans);
            trans[idx] = 1.0;
            if (!cell->reduced_coordinates) {
                V3 tmp;
                V3EQV(trans, tmp);
                V3mM3(tmp, cell->h, trans);
            }

            p_point = p->p_point_end;
            p->p_point_end += nn_send;

            point_send0 = p->point_begin + p->n_max;
            if (n_send[dir0] > 0) {
                for (pnt = p->point_begin; pnt < p->image_end; pnt++) {
                    if (pnt->iw & dir_flags0) {
                        *p_point++ = pnt;
                        *--point_send0 = *pnt;
                    }
                }
                if (dir_flags0 & p->trans_mask) {
                    for (pnt=point_send0; pnt<p->point_begin+p->n_max; pnt++)
                        V3SuB(pnt->r, trans);
                }
                MPI_Isend(point_send0, n_send[dir0], p->M_POINT,
                        p->rank_send[dir0], dir0, p->comm, &req_send[0]);
            }
            else
                req_send[0] = MPI_REQUEST_NULL;
        
            point_send1 = point_send0;
            if (n_send[dir1] > 0) {
                for (pnt = p->point_begin; pnt < p->image_end; pnt++) {
                    if (pnt->iw & dir_flags1) {
                        *p_point++ = pnt;
                        *--point_send1 = *pnt;
                    }
                }
                if (dir_flags1 & p->trans_mask) {
                    for (pnt = point_send1; pnt < point_send0; pnt++) 
                        V3AdD(trans, pnt->r);
                }
                MPI_Isend(point_send1, n_send[dir1], p->M_POINT,
                        p->rank_send[dir1], dir1, p->comm, &req_send[1]);
            }
            else
                req_send[1] = MPI_REQUEST_NULL;

            MPI_Waitall(2, req_recv, status);
            if (cell->reduced_coordinates) {
                for (pnt = p->image_end; pnt < image_end; pnt++) {
                    pnt->iw = p3d_flags_boundary(pnt->r, lb_c, ub_c);
                }
            }
            else {
                for (pnt = p->image_end; pnt < image_end; pnt++) {
                    V3 s;
                    V3mM3(pnt->r, h_inv, s);
                    pnt->iw = p3d_flags_boundary(s, lb_c, ub_c);
                }
            }
            p->image_end = image_end;
            MPI_Waitall(2, req_send, status);
        }

        /* init */
        image_end = p->point_end;
        p_point = p->p_point_begin;
        for (dir0 = 0; dir0 < 6; dir0 += 2) {
            P3DAtom *point_send = p->point_begin + p->n_max;
            int dir1 = dir0 + 1;
        
            if (n_recv[dir0] > 0) {
                MPI_Recv_init(image_end, n_recv[dir0], p->M_IMAGE,
                        p->rank_recv[dir0], dir0, p->comm, &p->req_recv[dir0]);
                image_end += n_recv[dir0];
            }
            if (n_recv[dir1] > 0) {
                MPI_Recv_init(image_end, n_recv[dir1], p->M_IMAGE,
                        p->rank_recv[dir1], dir1, p->comm, &p->req_recv[dir1]);
                image_end += n_recv[dir1];
            }

            p->update[dir0] = p_point;
            if (n_send[dir0] > 0) {
                p_point    += n_send[dir0];
                point_send -= n_send[dir0];
                MPI_Send_init(point_send, n_send[dir0], p->M_IMAGE,
                        p->rank_send[dir0], dir0, p->comm, &p->req_send[dir0]);
            }
        
            p->update[dir1] = p_point;
            if (n_send[dir1] > 0) {
                p_point    += n_send[dir1];
                point_send -= n_send[dir1];
                MPI_Send_init(point_send, n_send[dir1], p->M_IMAGE,
                        p->rank_send[dir1], dir1, p->comm, &p->req_send[dir1]);
            }

        }
        p->update[6] = p_point;
    }

    p->f_realloc = 0;

    return n_send_total;
}


int p3d_update_cell(Cell cell)
{
    Cell_private p = cell->private;
    P3DAtom *pnt, *point_send0, *point_send1, **p_point;
    V3 trans;
    MPI_Status status[2];
    int dir0, dir1, dir2;

    if (p->f_realloc) p3d_abort(-1);

    for (dir0 = 0; dir0 < 6; dir0 += 2) {
        int idx = dir0 / 2;
        dir1 = dir0 + 1;
    
        if (p->req_recv[dir0]!=MPI_REQUEST_NULL) MPI_Start(&p->req_recv[dir0]);
        if (p->req_recv[dir1]!=MPI_REQUEST_NULL) MPI_Start(&p->req_recv[dir1]);

        V3ZERO(trans);
        trans[idx] = 1.0;
        if (!cell->reduced_coordinates) {
            V3 tmp;
            V3EQV(trans, tmp);
            V3mM3(tmp, cell->h, trans);
        }
    
        point_send0 = p->point_begin + p->n_max;
        if (p->req_send[dir0] != MPI_REQUEST_NULL) {
            for (p_point=p->update[dir0]; p_point<p->update[dir1]; p_point++)
                *--point_send0 = **p_point;
            if (direction_flags[dir0] & p->trans_mask) {
                for (pnt = point_send0; pnt < p->point_begin+p->n_max; pnt++)
                    V3SuB(pnt->r, trans);
            }
            MPI_Start(&p->req_send[dir0]);
        }
    
        point_send1 = point_send0;
        if (p->req_send[dir1] != MPI_REQUEST_NULL) {
            dir2 = dir1 + 1;
            for (p_point=p->update[dir1]; p_point<p->update[dir2]; p_point++)
                *--point_send1 = **p_point;
            if (direction_flags[dir1] & p->trans_mask) {
                for (pnt = point_send1; pnt < point_send0; pnt++)
                    V3AdD(trans, pnt->r);
            }
            MPI_Start(&p->req_send[dir1]);
        }
    
        MPI_Waitall(2, &p->req_recv[dir0], status);
        MPI_Waitall(2, &p->req_send[dir0], status);
    }

    return 0;
}


/* Comm */
int p3d_bcast(Cell cell, void *buff, int count, MPI_Datatype datatype, int root)
{
    return MPI_Bcast(buff, count, datatype, root, cell->private->comm);
}

int p3d_reduce(Cell cell, void *sbuf, void *rbuf, int count,
                        MPI_Datatype datatype, MPI_Op op)
{
    return MPI_Allreduce(sbuf, rbuf, count, datatype, op, cell->private->comm);
}


/* I/O */

int p3d_fflush(FILE *fp)
{
    int ret = 0;

    if ((p3d_global_rank_fprintf==2)                       ||
        (p3d_global_rank_fprintf==0 && p3d_global_rank==0) ||
        (p3d_global_rank_fprintf==1 && p3d_global_rank==p3d_global_nprocs-1)) {
        ret = fflush(fp);
    }

    return ret;
}


int p3d_fprintf(FILE *fp, const char *format, ...)
{
    int ret;

    ret=(p3d_global_rank_fprintf==2)? fprintf(fp, "[%d] ", p3d_global_rank) : 0;
    if ((p3d_global_rank_fprintf==2)                       ||
        (p3d_global_rank_fprintf==0 && p3d_global_rank==0) ||
        (p3d_global_rank_fprintf==1 && p3d_global_rank==p3d_global_nprocs-1)) {
        char *s;
        for (s = (char*)format; s; s++) {
            if (*s == '%' && *(s+1) != '%')
                break;
        }
        if (*s) {
            va_list va;
            va_start(va, format);
            ret += vfprintf(fp, format, va);
            va_end(va);
        }
        else
            ret += fprintf(fp, format);
    }

    return ret;
}


static char *fgets2(char *s, int size, FILE *stream)
{
    do {
        if (NULL == fgets(s, size, stream)) return NULL;
    } while ('#' == *s || '\n' == *s);
    return s;
}


static int GZIP_MAGIC[]  = {0x1f, 0x8b};
static int BZIP2_MAGIC[] = {'B','Z','h'};
static int PMDS_MAGIC[]  = {'<','M','D','S','t','e','n','c','i','l','A'};

struct filter_t {
    int num;
    int *magic;
    char *command;
    char *option;
};

static struct filter_t filters[] = {
    { sizeof(GZIP_MAGIC) /sizeof(int), GZIP_MAGIC,  "gzip",  "-cd" },
    { sizeof(BZIP2_MAGIC)/sizeof(int), BZIP2_MAGIC, "bzip2", "-cd" },
    { sizeof(PMDS_MAGIC) /sizeof(int), PMDS_MAGIC,  "pmds_cat", "-v" },
    { 0, NULL, NULL, NULL }
};


static FILE *fopen2(const char *path, const char *mode)
{
    FILE *fp = fopen(path, mode);

    if (fp && *mode == 'r') {
        struct filter_t *f;

        for (f = filters; f->num; f++) {
            int i;

            for (i = 0; i < f->num; i++) {
                if (f->magic[i] != fgetc(fp))
                    break;
            }
            if (i == f->num) {
                int pipefds[2];
                pid_t pid;

                fclose(fp);
                if (pipe(pipefds) == -1 || (pid = fork()) == -1)
                    return NULL;

                if (pid == 0) {
                    if (fork() == 0) {
                        char *argv[4];
                        close(1);
                        dup(pipefds[1]);
                        close(pipefds[1]);
                        close(pipefds[0]);
                        argv[0] = f->command;
                        argv[1] = f->option;
                        argv[2] = (char*)path;
                        argv[3] = NULL;
                        execvp(f->command, argv);
                    }
                    _exit(0);
                }
                else {
                    wait(NULL);
                    close(pipefds[1]);
                    fp = fdopen(pipefds[0], mode);
                    return fp;
                }
            }
            else
                rewind(fp);
        }
    }

    return fp;
}

#ifndef USE_ZLIB

    typedef FILE *P3DFileP;
#   define P3DFopen(path, mode) ((P3DFileP)fopen2(path, mode))
#   define P3DFclose(fp)       fclose(fp)       
#   define P3DFputs(s, fp)     fputs(s, fp)
#   define P3DFgets(s, n, fp)  fgets(s, n, fp)
#   define P3DFgets2(s, n, fp) fgets2(s, n, fp)

#else

#   include "zlib.h"
    typedef struct File_tag {
        void *fp;
        int compress;
    } *P3DFileP;
#   define P3DFputs(s, fp)     (((fp)->compress) ? gzputs((gzFile)(fp)->fp, s)\
                                              : fputs(s, (FILE*)(fp)->fp))
#   define P3DFgets(s, n, fp)  (((fp)->compress) ? gzgets((gzFile)(fp)->fp,s,n)\
                                              : fgets(s, n, (FILE*)(fp)->fp))

    static int chknamegz(const char *path)
    {
        int len = strlen(path);
        return (len > 2 && path[len-3] == '.'
                        && path[len-2] == 'g' && path[len-1] == 'z');
    }
    
    static P3DFileP P3DFopen(const char *path, const char *mode)
    {
        P3DFileP fp = (P3DFileP) malloc(sizeof(struct File_tag));

        if (fp) {
            if (*mode == 'r') {
                fp->compress = 0;
                fp->fp = (P3DFileP) fopen2(path, mode);
            }
            else if (chknamegz(path)) {
                fp->compress = 1;
                fp->fp = (P3DFileP) gzopen(path, mode);
            }
            else {
                fp->compress = 0;
                fp->fp = (P3DFileP)  fopen(path, mode);
            }
        }
        return fp;
    }
    
    int P3DFclose(P3DFileP fp)
    {
        int value = EOF;

        if (fp) {
            value = (fp->compress) ?
                    gzclose((gzFile)(fp->fp)) : fclose((FILE*)(fp->fp));
            free(fp);
        }

        return value;
    }
    
    static char *P3DFgets2(char *s, int n, P3DFileP fp)
    {
        if (fp->compress) {
            do {
                if (NULL == gzgets((gzFile)fp->fp, s, n)) return NULL;
            } while ('#' == *s || '\n' == *s);
            return s;
        }
        else
            return fgets2(s, n, (FILE*)fp->fp);
    
    }
#endif

static int chkpathname(char *path)
{
   int i, divided = 0;

   for (i=0; i<strlen(path)-1; i++) {
       if ( path[i] == '%' && path[++i] != '%') {
           divided = 1;
           break;
       }
   }

   return divided;
}


static char *p3d_filename(char *path, Cell cell)
{
    int i;
    static char newpath[FILENAME_MAX];

    newpath[0] = '\0';
    for (i=0; i<strlen(path)-1; i++) {
        if ( path[i] == '%' && path[++i] != '%') {
            sprintf(newpath, path, cell->private->rank);
            break;
        }
    }

    return (newpath[0] != '\0') ? newpath : path;
}




enum p3d_cfg_tags {
    CHG_TAG_NOMATCH,
    CHG_TAG_NUM, CHG_TAG_A, CHG_TAG_R,
    CHG_TAG_H0_11, CHG_TAG_H0_12, CHG_TAG_H0_13,
    CHG_TAG_H0_21, CHG_TAG_H0_22, CHG_TAG_H0_23,
    CHG_TAG_H0_31, CHG_TAG_H0_32, CHG_TAG_H0_33,
    CHG_TAG_Tr_11, CHG_TAG_Tr_12, CHG_TAG_Tr_13,
    CHG_TAG_Tr_21, CHG_TAG_Tr_22, CHG_TAG_Tr_23,
    CHG_TAG_Tr_31, CHG_TAG_Tr_32, CHG_TAG_Tr_33,
    CHG_TAG_Et_11, CHG_TAG_Et_12, CHG_TAG_Et_13,
    CHG_TAG_Et_21, CHG_TAG_Et_22, CHG_TAG_Et_23,
    CHG_TAG_Et_31, CHG_TAG_Et_32, CHG_TAG_Et_33,
    CHG_TAG_NO_V, CHG_TAG_ENTRY, CHG_TAG_AUX
};

typedef struct {
    enum p3d_cfg_tags id;
    char *name;
    char *comment;
} ChgTag;

static ChgTag CHG_TAGS[] = {
    {CHG_TAG_NUM,   "Number of particles =",
                    "(required) this must be the first line"},
    {CHG_TAG_A,     "A =",
                    "(optional) basic length-scale: default A = 1.0"},
    {CHG_TAG_R,     "R =",
                    "(optional) basic rate-scale: default R = 1.0 [ns^-1]"},
    {CHG_TAG_H0_11, "H0(1,1) =", NULL},
    {CHG_TAG_H0_12, "H0(1,2) =", NULL},
    {CHG_TAG_H0_13, "H0(1,3) =",
                    "(required) this is the supercell's 1st edge, in A"},
    {CHG_TAG_H0_21, "H0(2,1) =", NULL},
    {CHG_TAG_H0_22, "H0(2,2) =", NULL},
    {CHG_TAG_H0_23, "H0(2,3) =",
                    "(required) this is the supercell's 2nd edge, in A"},
    {CHG_TAG_H0_31, "H0(3,1) =", NULL},
    {CHG_TAG_H0_32, "H0(3,2) =", NULL},
    {CHG_TAG_H0_33, "H0(3,3) =",
                    "(required) this is the supercell's 3rd edge, in A"},
    {CHG_TAG_Tr_11, "Transform(1,1) =", NULL},
    {CHG_TAG_Tr_12, "Transform(1,2) =", NULL},
    {CHG_TAG_Tr_13, "Transform(1,3) =", NULL},
    {CHG_TAG_Tr_21, "Transform(2,1) =", NULL},
    {CHG_TAG_Tr_22, "Transform(2,2) =", NULL},
    {CHG_TAG_Tr_23, "Transform(2,3) =", NULL},
    {CHG_TAG_Tr_31, "Transform(3,1) =", NULL},
    {CHG_TAG_Tr_32, "Transform(3,2) =", NULL},
    {CHG_TAG_Tr_33, "Transform(3,3) =",
                    "(optional) apply additional transformation on H0:"
                    "  H = H0 * Transform;\n# default = Identity matrix."},
    {CHG_TAG_Et_11, "eta(1,1) =", NULL},
    {CHG_TAG_Et_12, "eta(1,2) =", NULL},
    {CHG_TAG_Et_13, "eta(1,3) =", NULL},
    {CHG_TAG_Et_21, "eta(2,1) =", NULL},
    {CHG_TAG_Et_22, "eta(2,2) =", NULL},
    {CHG_TAG_Et_23, "eta(2,3) =", NULL},
    {CHG_TAG_Et_31, "eta(3,1) =", NULL},
    {CHG_TAG_Et_32, "eta(3,2) =", NULL},
    {CHG_TAG_Et_33, "eta(3,3) =",
                    "(optional) apply additional Lagrangian strain on H0:"
                    "\n# H = H0 * sqrt(Identity_matrix + 2 * eta);"
                    "\n# default = zero matrix."},
    {CHG_TAG_NO_V,  ".NO_VELOCITY.", NULL},
    {CHG_TAG_ENTRY, "entry_count =", NULL},
    {CHG_TAG_AUX,   "auxiliary[", NULL},
};
#define N_CHG_TAGS sizeof(CHG_TAGS)/sizeof(ChgTag)


static int chk_tag(char *line, char **s)
{
    ChgTag *ct;
    for (ct = CHG_TAGS; ct < CHG_TAGS + N_CHG_TAGS; ct++) {
        if ((*s = strstr(line, ct->name))) {
            *s += strlen(ct->name);
            return ct->id;
        }
    }
    return CHG_TAG_NOMATCH;
}

#define SCANFMT3    "%lf%lf%lf"
#define AUXSCANFMT  "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf"
#define AUX16(pnt)  (pnt)->aux,   (pnt)->aux+1, (pnt)->aux+2, (pnt)->aux+3,\
                    (pnt)->aux+4, (pnt)->aux+5, (pnt)->aux+6, (pnt)->aux+7,\
                    (pnt)->aux+8, (pnt)->aux+9, (pnt)->aux+10,(pnt)->aux+11,\
                    (pnt)->aux+12,(pnt)->aux+13,(pnt)->aux+14,(pnt)->aux+15\

#define ulength_IN_A 1.0
#define utime_IN_NS 1.0
#define umass_IN_AMU 1.0


static int read_config_header(Cell cell, P3DFileP fp, char *line, int size,
                            double *r, int *entry_count, int *no_velocity)
{
    int n;
    double a = 1.0;

    *r = 1.0;
    *entry_count = *no_velocity = 0;

    while (P3DFgets2(line, size, fp)) {
        char *s;
        int k;
        chk_tag(line, &s);
        switch (chk_tag(line, &s)) {
        case CHG_TAG_NOMATCH:
            if (a != 1.0) {
                a /= ulength_IN_A;
                M3MultiplY(a, cell->h);
            }
            if (*r != 1.0) {
                *r *= utime_IN_NS;
            }
            return n;
            break;
        case CHG_TAG_NUM:
            sscanf(s, "%d", &n);
            break;
        case CHG_TAG_A:
            sscanf(s, "%lf", &a);
            break;
        case CHG_TAG_R:
            sscanf(s, "%lf", r);
            break;
        case CHG_TAG_H0_11: sscanf(s, "%lf", &cell->h[0][0]);   break;
        case CHG_TAG_H0_12: sscanf(s, "%lf", &cell->h[0][1]);   break;
        case CHG_TAG_H0_13: sscanf(s, "%lf", &cell->h[0][2]);   break;
        case CHG_TAG_H0_21: sscanf(s, "%lf", &cell->h[1][0]);   break;
        case CHG_TAG_H0_22: sscanf(s, "%lf", &cell->h[1][1]);   break;
        case CHG_TAG_H0_23: sscanf(s, "%lf", &cell->h[1][2]);   break;
        case CHG_TAG_H0_31: sscanf(s, "%lf", &cell->h[2][0]);   break;
        case CHG_TAG_H0_32: sscanf(s, "%lf", &cell->h[2][1]);   break;
        case CHG_TAG_H0_33: sscanf(s, "%lf", &cell->h[2][2]);   break;
        case CHG_TAG_NO_V:
            *no_velocity = 1;
            break;
        case CHG_TAG_ENTRY:
            sscanf(s, "%d" , entry_count);
            cell->auxiliary->n_aux = *entry_count;
            if (*no_velocity) cell->auxiliary->n_aux -= 3;
            break;
        case CHG_TAG_AUX:
            if (sscanf(s, "%d", &k) == 1 && k < CONFIG_MAX_AUXILIARY) {
                char *name, *unit;
                if (k >= cell->auxiliary->n_aux)
                         cell->auxiliary->n_aux = k + 1;
                (name = cell->auxiliary->name[k])[0] = 0;
                (unit = cell->auxiliary->unit[k])[0] = 0;
                if ((s = strchr(s, '='))) {
                    s++;
                    while (*s && isspace(*s)) s++;
                    if (*s) {
                        int i;
                        sscanf(s, "%[^[]", name);
                        for (i = strlen(name) - 1; i >= 0; i--) {
                            if (!isspace(name[i]))
                                break;
                            else
                                name[i] = 0;
                        }
                        if ((s = strchr(s, '['))) {
                            strncpy(unit, s, P3D_LINEMAX);
                            for (i = strlen(unit) - 1; i >= 0; i--) {
                                if (!isspace(unit[i]))
                                    break;
                                else
                                    unit[i] = 0;
                            }
                        }
                    }
                }
            }
            break;
        default:
            break;
        }
    }

    return -1;
}


static int read_config_body(Cell cell, P3DFileP fp, char *line, int size,
                P3DAtom *point, int max, int *remain, double r,
                int entry_count, int no_velocity)
{
    static char symbol[sizeof(int)] = " ";
    P3DAtom *point_end;
    int i, n = (*remain < max) ? *remain : max;
    *remain -= n;

    if (!entry_count) { /* Standard cfg */
        P3DAtom *pnt = point;
        sscanf(line, "%lf%s%lf%lf%lf%lf%lf%lf",
                &pnt->m, &symbol[1], V3e(pnt->r), V3e(pnt->v));
        strncpy(pnt->sym, &symbol[strlen(symbol)-2], sizeof(int));
        pnt++;

        for (i = 1; i < n; i++) {
            P3DFgets2(line, size, fp);
            sscanf(line, "%lf%s%lf%lf%lf%lf%lf%lf",
                    &pnt->m, &symbol[1], V3e(pnt->r), V3e(pnt->v));
            strncpy(pnt->sym, &symbol[strlen(symbol)-2],sizeof(int));
            pnt++;
        }
        point_end = pnt;
    }
    else { /* Extended cfg */
        static char *sym;
        static double mass;
        P3DAtom *pnt = point;
        if (no_velocity) {
            if (sscanf(line,SCANFMT3 AUXSCANFMT,V3e(pnt->r),AUX16(pnt))==1) {
                mass = pnt->r[0];
                P3DFgets2(line, size, fp);
                sscanf(line, "%s", &symbol[1]);
                sym = &symbol[strlen(symbol) - 2];
                P3DFgets2(line, size, fp);
                sscanf(line,SCANFMT3 AUXSCANFMT, V3e(pnt->r), AUX16(pnt));
            }
            pnt->m = mass;
            strncpy(pnt->sym, sym, sizeof(int));
            V3ZERO(pnt->v);
            pnt++;

            for (i = 1; i < n; i++) {
                P3DFgets2(line, size, fp);
                if (sscanf(line,SCANFMT3 AUXSCANFMT,V3e(pnt->r),AUX16(pnt))==1){
                    mass = pnt->r[0];
                    P3DFgets2(line, size, fp);
                    sscanf(line, "%s", &symbol[1]);
                    sym = &symbol[strlen(symbol) - 2];
                    P3DFgets2(line, size, fp);
                    sscanf(line, SCANFMT3 AUXSCANFMT, V3e(pnt->r), AUX16(pnt));
                }
                pnt->m = mass;
                strncpy(pnt->sym, sym, sizeof(int));
                V3ZERO(pnt->v);
                pnt++;
            }

        }
        else { /* with velocity */
            if (sscanf(line, "%lf%lf%lf%lf%lf%lf"AUXSCANFMT,
                                V3e(pnt->r), V3e(pnt->v), AUX16(pnt)) == 1) {
                mass = pnt->r[0];
                P3DFgets2(line, size, fp);
                sscanf(line, "%s", &symbol[1]);
                sym = &symbol[strlen(symbol) - 2];
                P3DFgets2(line, size, fp);
                sscanf(line, "%lf%lf%lf%lf%lf%lf"AUXSCANFMT,
                                V3e(pnt->r), V3e(pnt->v), AUX16(pnt));
            }
            pnt->m = mass;
            strncpy(pnt->sym, sym, sizeof(int));
            pnt++;

            for (i = 1; i < n; i++) {
                P3DFgets2(line, size, fp);
                if (sscanf(line, "%lf%lf%lf%lf%lf%lf"AUXSCANFMT,
                                    V3e(pnt->r), V3e(pnt->v), AUX16(pnt))==1) {
                    mass = pnt->r[0];
                    P3DFgets2(line, size, fp);
                    sscanf(line, "%s", &symbol[1]);
                    sym = &symbol[strlen(symbol) - 2];
                    P3DFgets2(line, size, fp);
                    sscanf(line, "%lf%lf%lf%lf%lf%lf"AUXSCANFMT,
                                    V3e(pnt->r), V3e(pnt->v), AUX16(pnt));
                }
                pnt->m = mass;
                strncpy(pnt->sym, sym, sizeof(int));
                pnt++;
            }
        }
        point_end = pnt;
    }


    if (!entry_count || !no_velocity) {
        P3DAtom *pnt;
        if (r != 1.0) {
            for (pnt = point; pnt < point_end; pnt++) {
                V3MuL(r, pnt->v);
            }
        }
        if (!cell->reduced_coordinates) {
            V3 sv;
            for (pnt = point; pnt < point_end; pnt++) {
                V3EQV(pnt->v, sv);
                V3mM3(sv, cell->h, pnt->v);
            }
        }
    }

    if (!cell->reduced_coordinates) {
        P3DAtom *pnt;
        V3 sr;
        for (pnt = point; pnt < point_end; pnt++) {
            V3EQV(pnt->r, sr);
            V3mM3(sr, cell->h, pnt->r);
        }
    }

    if (umass_IN_AMU != 1.0) {
        P3DAtom *pnt;
        for (pnt = point; pnt < point_end; pnt++) {
            pnt->m /= umass_IN_AMU;
        }
    }

    if (*remain) P3DFgets2(line, size, fp);

    return point_end - point;
}


static int read_config_body_filter(Cell cell, P3DFileP fp, char *line, int size,
                P3DAtom *point, int max, int *remain, double r,
                int entry_count, int no_velocity)
{
    static char symbol[sizeof(int)] = " ";
    P3DAtom *start = point, *end = start + max, *pnt = start;
    P3DAtom dummy, *dp = &dummy;
    V3 lb, ub;
    int i;

    for (i = 0; i < 3; i++) {
        lb[i] = ((double)cell->private->coords[i]    )/ cell->private->dims[i];
        ub[i] = ((double)cell->private->coords[i] + 1)/ cell->private->dims[i];
    }
        
    if (!entry_count) { /* Standard cfg */
        sscanf(line, "%lf%s%lf%lf%lf%lf%lf%lf", 
                        &dp->m, &symbol[1], V3e(dp->r), V3e(dp->v));
        if (!p3d_flags_boundary2(dp->r, lb, ub)) {
            if (pnt == end)
                goto scale;
            strncpy(dp->sym, &symbol[strlen(symbol)-2], sizeof(int));
            *pnt++ = *dp;
        }

        for (i = 1; i < *remain; i++) {
            P3DFgets2(line, size, fp);
            sscanf(line, "%lf%s%lf%lf%lf%lf%lf%lf",
                            &dp->m, &symbol[1], V3e(dp->r), V3e(dp->v));
            if (!p3d_flags_boundary2(dp->r, lb, ub)) {
                if (pnt == end) {
                    *remain -= i;
                    goto scale;
                }
                strncpy(dp->sym, &symbol[strlen(symbol)-2], sizeof(int));
                *pnt++ = *dp;
            }
        }
    }
    else { /* Extended cfg */
        static char *sym;
        static double mass;
        if (no_velocity) {
            V3ZERO(dp->v);
            if (sscanf(line,"%lf%lf%lf"AUXSCANFMT,V3e(dp->r),AUX16(dp)) == 1) {
                mass = dp->r[0];
                P3DFgets2(line, size, fp);
                sscanf(line, "%s", &symbol[1]);
                sym = &symbol[strlen(symbol) - 2];
                P3DFgets2(line, size, fp);
                sscanf(line,"%lf%lf%lf"AUXSCANFMT,V3e(dp->r),AUX16(dp));
            }
            if (!p3d_flags_boundary2(dp->r, lb, ub)) {
                if (pnt == end)
                    goto scale;
                dp->m = mass;
                strncpy(dp->sym, sym, sizeof(int));
                *pnt++ = *dp;
            }

            for (i = 1; i < *remain; i++) {
                P3DFgets2(line, size, fp);
                if (sscanf(line,"%lf%lf%lf"AUXSCANFMT,V3e(dp->r),AUX16(dp))==1){
                    mass = dp->r[0];
                    P3DFgets2(line, size, fp);
                    sscanf(line, "%s", &symbol[1]);
                    sym = &symbol[strlen(symbol) - 2];
                    P3DFgets2(line, size, fp);
                    sscanf(line, "%lf%lf%lf"AUXSCANFMT,V3e(dp->r),AUX16(dp));
                }
                if (!p3d_flags_boundary2(dp->r, lb, ub)) {
                    if (pnt == end) {
                        *remain -= i;
                        goto scale;
                    }
                    dp->m = mass;
                    strncpy(dp->sym, sym, sizeof(int));
                    *pnt++ = *dp;
                }
            }
        }
        else { /* with velocity */
            if (sscanf(line, "%lf%lf%lf%lf%lf%lf"AUXSCANFMT,
                                V3e(dp->r), V3e(dp->v), AUX16(dp)) == 1) {
                mass = dp->r[0];
                P3DFgets2(line, size, fp);
                sscanf(line, "%s", &symbol[1]);
                sym = &symbol[strlen(symbol) - 2];
                P3DFgets2(line, size, fp);
                sscanf(line, "%lf%lf%lf%lf%lf%lf"AUXSCANFMT,
                                V3e(dp->r), V3e(dp->v), AUX16(dp));
            }
            if (!p3d_flags_boundary2(dp->r, lb, ub)) {
                if (pnt == end)
                    goto scale;
                dp->m = mass;
                strncpy(dp->sym, sym, sizeof(int));
                *pnt++ = *dp;
            }

            for (i = 1; i < *remain; i++) {
                P3DFgets2(line, size, fp);
                if (sscanf(line, "%lf%lf%lf%lf%lf%lf"AUXSCANFMT,
                                    V3e(dp->r), V3e(dp->v), AUX16(dp))==1) {
                    mass = dp->r[0];
                    P3DFgets2(line, size, fp);
                    sscanf(line, "%s", &symbol[1]);
                    sym = &symbol[strlen(symbol) - 2];
                    P3DFgets2(line, size, fp);
                    sscanf(line, "%lf%lf%lf%lf%lf%lf"AUXSCANFMT,
                                    V3e(dp->r), V3e(dp->v), AUX16(dp));
                }
                if (!p3d_flags_boundary2(dp->r, lb, ub)) {
                    if (pnt == end) {
                        *remain -= i;
                        goto scale;
                    }
                    dp->m = mass;
                    strncpy(dp->sym, sym, sizeof(int));
                    *pnt++ = *dp;
                }
            }
        }
    }
    *remain = 0;
    end = pnt;

scale:
    if (!entry_count || !no_velocity) {
        if (r != 1.0) {
            for (pnt = start; pnt < end; pnt++) {
                V3MuL(r, pnt->v);
            }
        }
        if (!cell->reduced_coordinates) {
            V3 sv;
            for (pnt = start; pnt < end; pnt++) {
                V3EQV(pnt->v, sv);
                V3mM3(sv, cell->h, pnt->v);
            }
        }
    }

    if (!cell->reduced_coordinates) {
        V3 sr;
        for (pnt = start; pnt < end; pnt++) {
            V3EQV(pnt->r, sr);
            V3mM3(sr, cell->h, pnt->r);
        }
    }

    if (umass_IN_AMU != 1.0) {
        for (pnt = start; pnt < end; pnt++) {
            pnt->m /= umass_IN_AMU;
        }
    }

    return end - start;
}


int p3d_read_config(Cell cell, char *path)
{
    Cell_private p = cell->private;
    int n, ml, mg, nl, ndiv;
    int entry_count, no_velocity, remain = 0, ntmp = 0;
    int multi, k = 0;
    P3DAtom *ptmp = NULL;
    P3DFileP fp = NULL;
    double r;
    char line[P3D_LINEMAX];

    /*p3d_global_io_any = 0;*/
    multi = chkpathname(path);

    p3d_set_atoms(cell, p3d_atom_begin(cell), 0);
    cell->auxiliary->n_aux = 0;
    ndiv = (multi) ? 1 : p->nprocs;

    if (multi && !p3d_global_io_any && p->nprocs > 1) {
        fprintf(stderr, "open error\n");
        p3d_abort(-1);
    }

    if (p3d_global_io_any || p->rank == 0) {
        if ((fp = P3DFopen(p3d_filename(path, cell), "r")) == NULL)
            p3d_abort(-1);
    }

    if (p3d_global_io_any) {
        int (*read_body)(Cell,P3DFileP,char*,int,P3DAtom*,int,int*,
            double,int,int)
                        = (multi) ? read_config_body : read_config_body_filter;
        for (k = 0; ; k++) {
            if (remain == 0)
                remain = read_config_header(cell, fp, line, sizeof(line),
                                            &r, &entry_count, &no_velocity);
            if (remain < 0)
                break;
            else if (remain > 0) {
                if (k == 0) {
                    n = remain / ndiv + 1;
                    if (p->n_max < n)
                        p3d_realloc_cell(cell, n, 0);
                    nl = read_body(cell, fp, line, sizeof(line), p->point_begin,
                                p->n_max, &remain, r, entry_count, no_velocity);
                    p3d_set_atoms(cell, p->point_begin, nl);
                }
                else {
                    n = remain;
                    if(ntmp<n&&!(ptmp=(P3DAtom*)malloc(sizeof(P3DAtom)*(ntmp=n))))
                        p3d_abort(-2);
                    nl = read_body(cell, fp, line,sizeof(line), ptmp,
                                ntmp, &remain, r, entry_count, no_velocity);
                    p3d_cat_atoms(cell, ptmp, nl);
                }
            }
        }
        if (multi) {
            do {
                ml = p3d_reset_cell(cell);
                p3d_reduce(cell, &ml, &mg, 1, MPI_INT, MPI_SUM);
            } while (mg);
        }
    }
    else {
        for (k = 0; ; k++) {
            if (remain == 0 && p->rank == 0)
                remain = read_config_header(cell, fp, line, sizeof(line),
                                            &r, &entry_count, &no_velocity);
            p3d_bcast(cell, &remain, 1, MPI_INT, 0);
            if (remain < 0)
                break;
            else if (remain > 0) {
                p3d_bcast(cell, &cell->h[0][0], 9, MPI_DOUBLE, 0);
                p3d_bcast(cell, &cell->auxiliary->n_aux, 1, MPI_INT, 0);
                if (cell->auxiliary->n_aux) {
                    p3d_bcast(cell, cell->auxiliary,
                                sizeof(struct Auxiliary_tag), MPI_CHAR,0);
                }

                if (p->rank == 0) {
                    if (k == 0) {
                        n = remain / ndiv;
                        if (p->n_max < n)
                            p3d_realloc_cell(cell, n, 0);
                        nl = read_config_body(cell,fp,line,sizeof(line),
                                    p->point_begin, p->n_max, &remain,
                                    r, entry_count, no_velocity);
                        p3d_set_atoms(cell, p->point_begin, nl);
                    }
                    else {
                        n = remain;
                        if (ntmp < nl && !(ptmp =
                                    (P3DAtom*)malloc(sizeof(P3DAtom)*(ntmp=n))))
                            p3d_abort(-2);
                        nl = read_config_body(cell, fp, line, sizeof(line),ptmp,
                                    ntmp, &remain, r, entry_count, no_velocity);
                        p3d_cat_atoms(cell, ptmp, nl);
                    }
                }
                do {
                    ml = p3d_reset_cell(cell);
                    p3d_reduce(cell, &ml, &mg, 1, MPI_INT, MPI_SUM);
                } while (mg);
            }
        }
    }

    if (ptmp) free(ptmp);

    return p3d_n_atoms(cell);
}


#define TAG_WRITE(fp, buf, fmt, num, val) \
        do {\
            ChgTag *ct;\
            for (ct = CHG_TAGS; ct < CHG_TAGS + N_CHG_TAGS; ct++) {\
                if (ct->id == num) {\
                    sprintf(buf, "%s" fmt "\n", ct->name, val);\
                    P3DFputs(line, fp);\
                    if (ct->comment) {\
                        sprintf(buf, "# %s\n\n", ct->comment);\
                        P3DFputs(line, fp);\
                    }\
                    break;\
                }\
            }\
        } while (0)

int p3d_write_config(Cell cell, char *path,
        char *s_formats, char*velocity_formats, int n_aux, ...)
{
    Cell_private p = cell->private;
    V3 sr, sv;
    M3 h_inv;
    P3DAtom *pnt;
    P3DFileP fp = NULL;
    char line[P3D_LINEMAX], line2[P3D_LINEMAX], oldsymbol[sizeof(int)];
    double cv_inv = 1.0 / utime_IN_NS, oldmass;
    int n = p3d_n_atoms(cell), write_each = chkpathname(path);
    int entry_count = 3;
    struct auxiliary {
        char *description;
        char *format;
        void *pointer;
        int bytes_separation;
    } aux[CONFIG_MAX_AUXILIARY];

    if (!s_formats) /* Standerd */
        n_aux = 0;
    else {          /* Extended */
        va_list ap;
        int i;
        if (n_aux > CONFIG_MAX_AUXILIARY || n_aux < 0) p3d_abort(-1);
        va_start(ap, n_aux);
        for (i = 0; i < n_aux; i++) {
            aux[i].description = va_arg(ap, char *);
            aux[i].format = va_arg(ap, char *);
            aux[i].pointer = va_arg(ap, char *);
            aux[i].bytes_separation = va_arg(ap, int);
        }
        va_end(ap);
    }

    if (!write_each) {
        int m = n;
        MPI_Reduce(&m, &n, 1, MPI_INT, MPI_SUM, 0, p->comm);
    }

    if (write_each || p->rank == 0) {
        if ((fp = P3DFopen(p3d_filename(path,cell),"w")) == NULL) p3d_abort(-1);

        if (p->nprocs > 1 && write_each ) {
            sprintf(line2, "%d \t\t# %%%d (/%d)", n, p->rank, p->nprocs);
            TAG_WRITE(fp, line, " %s", CHG_TAG_NUM, line2);
        }
        else
            TAG_WRITE(fp, line, " %d", CHG_TAG_NUM, n);

        TAG_WRITE(fp, line,"%s",CHG_TAG_A," 1.0 Angstrom (basic length-scale)");

        TAG_WRITE(fp, line, DFORMAT" A", CHG_TAG_H0_11, cell->h[0][0]);
        TAG_WRITE(fp, line, DFORMAT" A", CHG_TAG_H0_12, cell->h[0][1]);
        TAG_WRITE(fp, line, DFORMAT" A", CHG_TAG_H0_13, cell->h[0][2]);
        TAG_WRITE(fp, line, DFORMAT" A", CHG_TAG_H0_21, cell->h[1][0]);
        TAG_WRITE(fp, line, DFORMAT" A", CHG_TAG_H0_22, cell->h[1][1]);
        TAG_WRITE(fp, line, DFORMAT" A", CHG_TAG_H0_23, cell->h[1][2]);
        TAG_WRITE(fp, line, DFORMAT" A", CHG_TAG_H0_31, cell->h[2][0]);
        TAG_WRITE(fp, line, DFORMAT" A", CHG_TAG_H0_32, cell->h[2][1]);
        TAG_WRITE(fp, line, DFORMAT" A", CHG_TAG_H0_33, cell->h[2][2]);
        P3DFputs("# H = H0\n\n", fp);

        if (!s_formats) { /* Standard cfg */
            P3DFputs(
                    "# ENSUING ARE THE ATOMS, EACH ATOM DESCRIBED BY A ROW\n"
                    "# 1st entry is atomic mass in a.m.u.\n"
                    "# 2nd entry is the chemical symbol (max 2 chars)\n"
                    "\n"
                    "# 3rd entry is reduced coordinate s1 (dimensionless)\n"
                    "# 4th entry is reduced coordinate s2 (dimensionless)\n"
                    "# 5th entry is reduced coordinate s3 (dimensionless)\n"
                    "# real coordinates x = s * H,  x, s are 1x3 row vectors\n"
                    "\n"
                    "# 6th entry is d(s1)/dt in basic rate-scale R\n"
                    "# 7th entry is d(s2)/dt in basic rate-scale R\n"
                    "# 8th entry is d(s3)/dt in basic rate-scale R\n"
                    , fp);
            TAG_WRITE(fp, line, "%s", CHG_TAG_R, " 1.0 [ns^-1]");
        }
        else { /* Extended cfg */
            int i;
            P3DFputs(
                    "# Each atom is described by a row:\n"
                    "\n"
                    "# 1st entry is reduced coordinate s1 (dimensionless)\n"
                    "# 2nd entry is reduced coordinate s1 (dimensionless)\n"
                    "# 3rd entry is reduced coordinate s1 (dimensionless)\n"
                    "# real coordinates x = s * H,  x, s are 1x3 row vectors\n"
                    "\n"
                    , fp);
            if (velocity_formats) {
                P3DFputs(
                    "# 4th entry is d(s1)/dt in basic rate-scale R\n"
                    "# 5th entry is d(s2)/dt in basic rate-scale R\n"
                    "# 6th entry is d(s3)/dt in basic rate-scale R\n"
                    , fp);
                TAG_WRITE(fp, line, "%s", CHG_TAG_R, " 1.0 [ns^-1]");
                entry_count += 3;
            }
            else {
                TAG_WRITE(fp, line, "%s", CHG_TAG_NO_V, "");
                P3DFputs(
                    "# Atom velocities are deemed irrelevant for\n"
                    "# this configuration so they are not stored.\n"
                    "\n"
                    , fp);
            }

            entry_count += n_aux;
            TAG_WRITE(fp, line, " %d\n", CHG_TAG_ENTRY, entry_count);

            for (i = 0; i < n_aux; i++) {
                sprintf(line2, "%d] = %s", i, aux[i].description);
                TAG_WRITE(fp, line, "%s", CHG_TAG_AUX, line2);
            }

            pnt = p3d_atom_begin(cell);
            oldmass = pnt->m;
            strncpy(oldsymbol, pnt->sym, sizeof(oldsymbol));
            sprintf(line,
                "#\n"
                "# These properties are piece-wise uniform:\n"
                "%.15g\n"
                "# (required) atomic mass in a.m.u.\n"
                "%2s\n"
                "# (required) chemical symbol (max 2 chars)\n\n",
                oldmass * umass_IN_AMU, oldsymbol);
            P3DFputs(line, fp);
        }

        M3inv(cell->h, h_inv);

        if (!s_formats) {
            if (cell->reduced_coordinates) {
                for (pnt=p3d_atom_begin(cell);pnt<p3d_atom_end(cell);pnt++) {
                    V3MUL(cv_inv, pnt->v, sv);
                    sprintf(line, "%f %s"FFORMAT3 FFORMAT3"\n",
                            pnt->m*umass_IN_AMU, pnt->sym, V3E(pnt->r),V3E(sv));
                    P3DFputs(line, fp);
                }
            }
            else {
                for (pnt=p3d_atom_begin(cell);pnt<p3d_atom_end(cell);pnt++) {
                    V3mM3(pnt->r, h_inv, sr);
                    V3mM3(pnt->v, h_inv, sv);
                    V3MuL(cv_inv, sv);
                    sprintf(line, "%f %s"FFORMAT3 FFORMAT3"\n",
                            pnt->m*umass_IN_AMU, pnt->sym, V3E(sr), V3E(sv));
                    P3DFputs(line, fp);
                }
            }
        }
        else {
            char *pointer[CONFIG_MAX_AUXILIARY];
            int i;
            for (i = 0; i < n_aux; i++) pointer[i] = aux[i].pointer;
            for (pnt=p3d_atom_begin(cell); pnt<p3d_atom_end(cell); pnt++) {
                if (oldmass != pnt->m || strncmp(oldsymbol, pnt->sym, 2)) {
                    oldmass =  pnt->m;
                    strncpy(oldsymbol, pnt->sym, sizeof(oldsymbol));
                    sprintf(line,"%.15g\n%2s\n",oldmass*umass_IN_AMU,oldsymbol);
                    P3DFputs(line, fp);
                }
                if (cell->reduced_coordinates)
                    sprintf(line, s_formats, V3E(pnt->r));
                else {
                    V3mM3(pnt->r, h_inv, sr);
                    sprintf(line, s_formats, V3E(sr));
                }
                P3DFputs(line, fp);
                if (velocity_formats) {
                    if (cell->reduced_coordinates)
                        V3MUL(cv_inv, pnt->v, sv);
                    else {
                        V3mM3(pnt->v, h_inv, sv);
                        V3MuL(cv_inv, sv);
                    }
                    sprintf(line, velocity_formats, V3E(sv));
                    P3DFputs(line, fp);
                }
                for (i = 0; i < n_aux; i++) {
                    sprintf(line, aux[i].format, *((double*)pointer[i]));
                    P3DFputs(line, fp);
                    pointer[i] += aux[i].bytes_separation;
                }
                P3DFputs("\n", fp);
            }
        }
    }

    if (!write_each) {
        static int ntmp = 0;
        static double *atmp = NULL;
        double *pp;
        MPI_Datatype m_point = cell->private->M_POINT;
        int tag = 0;
        if (p->rank) {
            MPI_Send(&n, 1, MPI_INT, 0, tag, p->comm);
            MPI_Send(p3d_atom_begin(cell), n, m_point, 0, tag, p->comm);
            if (n_aux) {
                char *pointer[CONFIG_MAX_AUXILIARY];
                int i, j;
                if (n > ntmp) {
                    ntmp = n;
                    atmp = (double*) realloc(atmp, sizeof(double)*ntmp*n_aux);
                    if (atmp == NULL) p3d_abort(-2);
                }
                for (i = 0; i < n_aux; i++) pointer[i] = aux[i].pointer;
                pp = atmp;
                for (j = 0; j < n; j++) {
                    for (i = 0; i < n_aux; i++) {
                        *pp++ = *((double*)pointer[i]);
                        pointer[i] += aux[i].bytes_separation;
                    }
                }
                MPI_Send(atmp, n*n_aux, MPI_DOUBLE, 0, tag, p->comm);
            }
        }
        else {
            MPI_Status stat;
            int rank, nn = n;
            static P3DAtom *ptmp = NULL;
            for (rank = 1; rank < p->nprocs; rank++) {
                MPI_Recv(&n, 1, MPI_INT, rank, tag, p->comm, &stat);
                if (n > ntmp) {
                    ntmp = n;
                    ptmp = (P3DAtom*) realloc(ptmp, sizeof(P3DAtom) * ntmp);
                    if (ptmp == NULL) p3d_abort(-2);
                    if (n_aux) {
                        atmp = (double*)realloc(atmp,sizeof(double)*ntmp*n_aux);
                        if (atmp == NULL) p3d_abort(-2);
                    }
                }
                MPI_Recv(ptmp, n, m_point, rank, tag, p->comm, &stat);
                if (n_aux) {
                    MPI_Recv(atmp, n*n_aux, MPI_DOUBLE, rank,tag,p->comm,&stat);
                    pp = atmp;
                }
                if (!s_formats) {
                    if (cell->reduced_coordinates) {
                        for (pnt = ptmp; pnt < ptmp + n; pnt++) {
                            V3MUL(cv_inv, pnt->v, sv);
                            sprintf(line, "%f %s"FFORMAT3 FFORMAT3"\n",
                                        pnt->m, pnt->sym, V3E(pnt->r),V3E(sv));
                            P3DFputs(line, fp);
                        }
                    }
                    else {
                        for (pnt = ptmp; pnt < ptmp + n; pnt++) {
                            V3mM3(pnt->r, h_inv, sr);
                            V3mM3(pnt->v, h_inv, sv);
                            V3MuL(cv_inv, sv);
                            sprintf(line, "%f %s"FFORMAT3 FFORMAT3"\n",
                                        pnt->m, pnt->sym, V3E(sr), V3E(sv));
                            P3DFputs(line, fp);
                        }
                    }
                } else {
                    for (pnt = ptmp; pnt < ptmp + n; pnt++) {
                        int i;
                        if (oldmass != pnt->m||strncmp(oldsymbol,pnt->sym,2)) {
                            oldmass =  pnt->m;
                            strncpy(oldsymbol, pnt->sym, sizeof(oldsymbol));
                            sprintf(line, "%.15g\n%2s\n",
                                        oldmass*umass_IN_AMU, oldsymbol);
                            P3DFputs(line, fp);
                        }
                        if (cell->reduced_coordinates)
                            sprintf(line, s_formats, V3E(pnt->r));
                        else {
                            V3mM3(pnt->r, h_inv, sr);
                            sprintf(line, s_formats, V3E(sr));
                        }
                        P3DFputs(line, fp);
                        if (velocity_formats) {
                            if (cell->reduced_coordinates)
                                V3MUL(cv_inv, pnt->v, sv);
                            else {
                                V3mM3(pnt->v, h_inv, sv);
                                V3MuL(cv_inv, sv);
                            }
                            sprintf(line, velocity_formats, V3E(sv));
                            P3DFputs(line, fp);
                        }
                        for (i = 0; i < n_aux; i++) {
                            sprintf(line, aux[i].format, *pp++);
                            P3DFputs(line, fp);
                        }
                        P3DFputs("\n", fp);
                    }
                }
            }
            n = nn;
        }
    }

    if (fp) P3DFclose(fp);

    return p3d_n_atoms(cell);
}


void p3d_coord_reduced2real(Cell cell)
{
    P3DAtom *pnt;
    V3 s;
    P3D_LOOP(cell, pnt) {
        V3EQV(pnt->r, s);
        V3mM3(s, cell->h, pnt->r);
        V3EQV(pnt->v, s);
        V3mM3(s, cell->h, pnt->v);
    }
    cell->reduced_coordinates = 0;
}
