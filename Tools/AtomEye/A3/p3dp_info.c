#include "A.h"
#include "p3dp.h"

void p3dp_s(double *si, int i, int irank)
{
    if (p3d_rank(p3dp_cell) == irank)
        V3EQV(&s[DIMENSION*i], si);
    p3d_bcast(p3dp_cell, si, 3, MPI_DOUBLE, irank);
}

void p3dp_atom_pair_s(double *sj, double *si, double dxji[4])
{
    double ds[3];
    V3SUB (sj, si, ds);
    V3ImagE (ds);
    V3M3LENGTH2 (ds, H, dxji);
    return;
}


static void print_atom_pair_info_s
(int iw, int j, int i, int jrank, int irank, double *sj, double *si)
{
    double dxji[4];
    p3dp_atom_pair_s(sj, si, dxji);
    if (IS_MANAGER)
        printf ("x[%d(%d)]-x[%d(%d)] = (%g %g %g), distance = %g A;\n",
                j, jrank, i, irank, V3E(dxji), sqrt(dxji[3]));
    return;
}

void p3dp_print_atom_pair_info(int iw, int j, int i, int jrank, int irank)
{
    double sj[3], si[3];
    p3dp_s(sj, j, jrank);
    p3dp_s(si, i, irank);
    print_atom_pair_info_s(iw, j, i, jrank, irank, sj, si);
    return;
}


double p3dp_atom_triplet_s
(double *sk, double *sj, double *si, double dxkj[4], double dxij[4])
{
    p3dp_atom_pair_s(sk, sj, dxkj);
    p3dp_atom_pair_s(si, sj, dxij);
    if ( (dxkj[3]>0) && (dxij[3]>0) )
        return( acos( V3DOT(dxkj, dxij) / sqrt(dxkj[3]) / sqrt(dxij[3]) ) );
    else return (0);
}


void p3dp_print_atom_triplet_info
(int iw, int k, int j, int i, int krank, int jrank, int irank)
{
    double dxkj[4], dxij[4], angle;
    double sk[3], sj[3], si[3];
    p3dp_s(sk, k, krank);
    p3dp_s(sj, j, jrank);
    p3dp_s(si, i, irank);
    angle = p3dp_atom_triplet_s(sk, sj, si, dxkj, dxij);
    print_atom_pair_info_s(iw, i, j, irank, jrank, si, sj);
    if (IS_MANAGER)
        printf ("bond angle = %g degrees.\n", RADIAN_TO_DEGREE(angle));
    print_atom_pair_info_s(iw, k, j, krank, jrank, sk, sj);
    return;
}


void p3dp_print_atom_quartet_info
(int iw, int l, int k, int j, int i, int lrank, int krank, int jrank, int irank)
{
    double dxkj[4], dxij[4], angle, normal[4], dxlk[4], dxjk[4], dihedral;
    double sl[3], sk[3], sj[3], si[3];
    p3dp_s(sl, l, lrank);
    p3dp_s(sk, k, krank);
    p3dp_s(sj, j, jrank);
    p3dp_s(si, i, irank);
    angle = p3dp_atom_triplet_s(sk, sj, si, dxkj, dxij);
    print_atom_pair_info_s(iw, i, j, irank, jrank, si, sj);
    if (IS_MANAGER)
        printf ("bond angle = %g degrees.\n", RADIAN_TO_DEGREE(angle));
    print_atom_pair_info_s(iw, k, j, krank, jrank, sk, sj);
    V3CROSS (dxkj, dxij, normal);
    normal[3] = V3LENGTH2 (normal);
    angle = p3dp_atom_triplet_s(sl, sk, sj, dxlk, dxjk);
    if (IS_MANAGER)
        printf ("bond angle = %g degrees.\n", RADIAN_TO_DEGREE(angle));
    print_atom_pair_info_s(iw, l, k, lrank, krank, sl, sk);
    /* right-handed helix gives positive dihedral angle */
    if ( (normal[3]>0) && (dxlk[3]>0) )
        dihedral = acos( V3DOT(normal,dxlk)/sqrt(normal[3])/sqrt(dxlk[3]) ) -
            PI / 2;
    else dihedral = 0;
    if (IS_MANAGER)
        printf ("dihedral angle = %g degrees.\n", RADIAN_TO_DEGREE(dihedral));
    return;
}
