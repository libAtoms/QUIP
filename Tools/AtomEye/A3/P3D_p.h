/* $Id: P3D_p.h,v 1.1.1.1 2005/03/20 20:14:24 shimizu Exp $
 * 2004-2005 Futoshi SHIMIZU
 */
#ifndef P3D_P_H
#define P3D_P_H

#include "P3D.h"

struct Cell_private_tag {
    P3DAtom *point_begin, *point_end, *image_end;
    P3DAtom **p_point_begin, **p_point_end, **update[7];
    int n_max, p_max, f_realloc;
    V3 lb, ub;
    double crust;
    int nprocs, rank;
    int dims[3], coords[3], trans_mask;
    int rank_send[6], rank_recv[6];
    MPI_Comm comm;
    MPI_Request req_send[6], req_recv[6];
    MPI_Datatype M_POINT, M_IMAGE;
};

typedef struct Cell_private_tag *Cell_private;

struct List_tag {
    int *num;
    Pair pair_point, pair_image;
    P3DAtom ***p3d_point, ***p3d_image, **p_point_begin;
    int n_max, p_max, i_max;
};

/* P3DCore.c */

int p3d_flags_boundary(V3 s, V3 lb, V3 ub);

#endif
