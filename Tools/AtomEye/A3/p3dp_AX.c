#include "AX/3D.c"
#include "p3dp.h"

typedef struct {
    double d;
    int i;
} double_int;

/* Calculate coordinates in viewport frame and determine visibility */
void p3dp_AX_3D_Balls_Zdet (int iw, AX_3D_Balls *B)
{
    register AX_3D_Ball *ball=NULL;
    register int i, m;
    register AX_Float tmp=0;
    int j, fp_activated, activated[AX_3D_MAX_FILTER_PLANE], i_g, i_l;
    AX_Float dx[3];
    for (fp_activated=0,j=0; j<AX_3D_MAX_FILTER_PLANE; j++)
    {
        activated[j] = AX_V3NEZERO(AX_3D[iw].fp[j].dx);
        fp_activated = fp_activated || activated[j];
    }
    for (m=0;;m++)
    { /* determine which balls are visible */
        B->n_shown = 0;
        for (i=B->n_balls; i--;)
        {
            ball = B->BALL + i;
            if (ball->r < 0) continue; /* invisibility flag */
            if (fp_activated)
            {
                for (j=0; j<AX_3D_MAX_FILTER_PLANE; j++)
                    if ( activated[j] )
                        if ( AX_V3DOT(AX_3D[iw].fp[j].dx, ball->x) <
                             AX_3D[iw].fp[j].d0-AX_TINY )
                            break;
                if (j<AX_3D_MAX_FILTER_PLANE) continue; /* invisibility flag */
            }
            AX_V3sub(ball->x, AX_3D[iw].x, dx);
            /* if viewpoint inside the ball, initiate ejection */
            tmp = AX_V3LENGTH2(dx);
            if (tmp <= SQUARE(ball->radius)) break;
            /* convert to viewframe coordinates */
            ball->vx[2] = AX_V3DOT(AX_3D[iw].V[2], dx);
            /* must be in (0,zcut) view frustum */
            if ((ball->vx[2] <= 0) ||
                (ball->vx[2] > AX_3D[iw].zcut)) continue;
            /* now ball is at front and not enclosing the viewpoint */
            tmp = AX_3D[iw].k / ball->vx[2];
            ball->vx[0] = AX_V3DOT(AX_3D[iw].V[0], dx) * tmp + AX_3D[iw].wx;
            ball->vx[1] = AX_V3DOT(AX_3D[iw].V[1], dx) * tmp + AX_3D[iw].wy;
            tmp *= ball->radius;
            if ( OU(ball->vx[0],-tmp,AX_size[iw].width+tmp) ||
                 OU(ball->vx[1],-tmp,AX_size[iw].height+tmp) ) continue;
            B->status[i] = 0;
            if (tmp < AX_BC_rmax[iw])
            {
                B->status[i] |= AX_3D_BALL_CACHED;
                ball->vs.bc = &AX_bc[iw][INT(tmp*AX_BC_RMESH)];
                /* circle center quantization (periodic) */
                FLOOR(ball->vx[0], ball->i0);
                FLOOR(ball->vx[1], ball->j0);
                ball->vx[0] = ball->i0 + 0.5;
                ball->vx[1] = ball->j0 + 0.5;
                tmp = ball->vs.bc->radius;
            }
            if ( AX_CircleIntersectWindowBorder
                 ( ball->vx[0], ball->vx[1], tmp,
                   AX_size[iw].width, AX_size[iw].height ) )
                B->idx[B->n_shown++] = i;
            else if ( AX_XYINWIN(iw, ball->vx[0], ball->vx[1]) )
            {
                B->idx[B->n_shown++] = i;
                B->status[i] |= AX_3D_BALL_COMPLETE;
            }
        }

        i_l = i;
        p3d_reduce(p3dp_cell, &i_l, &i_g, 1, MPI_INT, MPI_MAX);

        if (i_g >= 0) /* viewpoint ejection */
        {
            double_int di_g, di_l;
            di_l.d = tmp;
            di_l.i = p3d_rank(p3dp_cell);
            p3d_reduce(p3dp_cell, &di_l, &di_g,1,MPI_DOUBLE_INT,MPI_MINLOC);
            if (di_g.i == p3d_rank(p3dp_cell)) {
                if (tmp == 0)
                    AX_V3ADDMUL( AX_3D[iw].x, ball->radius * (1+AX_TINY),
                                 AX_3D[iw].V[2] );
                else
                {
                    tmp = ball->radius * (1+AX_TINY) / sqrt(tmp) - 1;
                    AX_V3SUBMUL( AX_3D[iw].x, tmp, dx );
                }
                if (m >= 3)
                {
                    AX_V3sub (AX_3D[iw].x, ball->x, dx);
                    AX_3D[iw].x[0] += Frandom() * dx[0];
                    AX_3D[iw].x[1] += Frandom() * dx[1];
                    AX_3D[iw].x[2] += Frandom() * dx[2];
                }
            }
            if (m >= 1000) pe ("AX_3D_Balls_Zdet: viewpoint ejection failed\n");
            p3d_bcast(p3dp_cell, AX_3D[iw].x, 3, MPI_DOUBLE, di_g.i);
        }
        else break;
    }
    AX_sort_3D_Balls(B);
    return;
} /* end AX_3D_Balls_Zdet() */

/* Find out the index of the topmost ball at a given pixel */
/* after Zpaper/ZBprint(); return -1 if there is no match. */
int p3dp_AX_3D_Balls_Zgrab (int iw, AX_3D_Balls *B, int gi, int gj, int *rank)
{
    register AX_3D_Ball *ball;
    register int m,goffset,k;
    register AX_Float *base, vradius, x2, y2, r2;
    register AX_BC_Unit *b;
    int i;
    double_int di_g, di_l = {AX_INFINITY, 0};
    if (AX_IJOUWIN(iw,gi,gj)) goto not_found;
    goffset = AX_OFF(iw,gi,gj);
    for (i=0; i<B->n_shown; i++)
    {
        ball = B->BALL + B->idx[i];
        vradius = AX_3D[iw].k / ball->vx[2] * ball->radius;
        x2 = SQUARE(gi+0.5-ball->vx[0]);
        y2 = SQUARE(gj+0.5-ball->vx[1]);
        r2 = SQUARE(vradius);
        if (x2 + y2 < r2)
        {
            if (B->status[B->idx[i]] & AX_3D_BALL_CACHED)
            {
                k = AX_OFF(iw, ball->i0, ball->j0);
                base = AX_zmap[iw] + k;
                k = goffset - k;
                b = ball->vs.bc->b;
                for (m=ball->vs.bc->nb; m--;)
                    if ((b[m].offset == k) &&
                        (ball->vx[2] - b[m].c1 * ball->radius <
                         base[b[m].offset] * (1+AX_TINY)))
                        goto found;
                b = ball->vs.bc->a;
                for (m=ball->vs.bc->na; m--;)
                    if ((b[m].offset == k) &&
                        (ball->vx[2] < base[b[m].offset]))
                        goto found;
            }
            else  /* scan it on the fly */
            {
                AX_CircleWindowScan
                    ( ball->vx[0], ball->vx[1], vradius,
                      AX_size[iw].width, AX_size[iw].height,
                      AX_AOP_top(iw), AX_BOP_top(iw) );
                for (m=1; m<AX_BOP_TOP(iw,0).offset; m++)
                {
                    if ( AX_OFF(iw, AX_BOP_TOP(iw,m).p.i,
                                AX_BOP_TOP(iw,m).p.j) != goffset )
                        continue;
                    if ( ball->vx[2] - sqrt(1.-(x2+y2)/r2) *
                         ball->radius <= AX_zmap[iw][goffset] * (1+AX_TINY) )
                        goto found;
                }
                for (m=1; m<AX_AOP_TOP(iw,0).b.offset; m++)
                {
                    if ( AX_OFF(iw, AX_AOP_TOP(iw,m).b.p.i,
                                AX_AOP_TOP(iw,m).b.p.j) != goffset )
                        continue;
                    if ( ball->vx[2] <= AX_zmap[iw][goffset] )
                        goto found;
                }
            }
        }
    }
    goto not_found;

found:
    di_l.d = ball->vx[2];
not_found:
    di_l.i = p3d_rank(p3dp_cell);
    p3d_reduce(p3dp_cell, &di_l, &di_g, 1, MPI_DOUBLE_INT, MPI_MINLOC);
    if (di_g.d == AX_INFINITY)
        return -1;
    else {
        int idx;
        *rank = di_g.i;
        if (p3d_rank(p3dp_cell) == *rank)
            idx = B->idx[i];
        p3d_bcast(p3dp_cell, &idx, 1, MPI_INT, *rank);
        return idx;
    }
} /* end AX_3D_Balls_Zgrab() */



void p3dp_AX_3D_Cylinders_Zdet (int iw, AX_3D_Cylinders *C)
{
    register AX_3D_Cylinder *cylinder;
    register AX_3D_Cylinder_Power *power;
    register double tmp=0;
    int i, m, i_g, i_l;
    double x1[4], radius2, UU[3],VV[3], xx[3], X0[4],X1[4], A,B,c;
    int j, fp_activated, activated[AX_3D_MAX_FILTER_PLANE];
    MALLOC (AX_3D_Cylinders_Zdet, C->power, C->n_cylinders,
            AX_3D_Cylinder_Power );
    MALLOC (AX_3D_Cylinders_Zdet, C->vz, C->n_cylinders, AX_Float);
    MALLOC (AX_3D_Cylinders_Zdet, C->idx, C->n_cylinders, int);
    REALLOC (AX_3D_Cylinders_Zdet, C->cylinder, C->n_cylinders,
             AX_3D_Cylinder *);
    for (fp_activated=0,j=0; j<AX_3D_MAX_FILTER_PLANE; j++)
    {
        activated[j] = AX_V3NEZERO(AX_3D[iw].fp[j].dx);
        fp_activated = fp_activated || activated[j];
    }
    for (power=C->power,m=0; ; power=C->power,m++)
    {
        for (i=0; i<C->n_cylinders; i++)
        {
            cylinder = C->CYLINDER + i;
            /* invisibility flag */
            if (cylinder->radius <= 0) continue;
            if (fp_activated)
            {
                AX_V3addmul (cylinder->x0, cylinder->axis[3],
                             cylinder->axis, x1);
                for (j=0; j<AX_3D_MAX_FILTER_PLANE; j++)
                    if ( activated[j] )
                        if ( ( AX_V3DOT(AX_3D[iw].fp[j].dx, cylinder->x0) <
                               AX_3D[iw].fp[j].d0-AX_TINY ) ||
                             ( AX_V3DOT(AX_3D[iw].fp[j].dx, x1) <
                               AX_3D[iw].fp[j].d0-AX_TINY ) )
                            break;
                if (j<AX_3D_MAX_FILTER_PLANE) continue; /* invisibility flag */
            }
            AX_V3sub (cylinder->x0, AX_3D[iw].x, power->dx);
            A = AX_V3DOT (power->dx, cylinder->axis);
            AX_V3submul (power->dx, A, cylinder->axis, power->y0);
            tmp = AX_V3LENGTH2 (power->y0);
            radius2 = SQUARE (cylinder->radius);
            if (tmp <= radius2 * (1+AX_3D_CYLINDER_REIL))
            {
                if ((A<=0) && (A>=-cylinder->axis[3])) break;
                else continue; /* invisible */
            }
            /* Convert to viewport frame */
            AX_M3mulV3 (AX_3D[iw].V, power->dx, power->y0);
            AX_M3mulV3 (AX_3D[iw].V, cylinder->axis, power->U);
            AX_V3addmul (power->y0, cylinder->axis[3], power->U, x1);
            if ( (power->y0[2] <= 0) || (power->y0[2] > AX_3D[iw].zcut) ||
                 (x1[2] <= 0) || (x1[2] > AX_3D[iw].zcut) )
            { /* intercepted by the z=0 plane */
                if ((power->y0[2] <= 0) && (x1[2] <= 0)) continue;
                else if ( (power->y0[2] <= 0) && (x1[2] > 0) )
                {
                    tmp = AX_3D_CYLINDER_VEIL * x1[2] - power->y0[2];
                    power->y0[0] += tmp * (x1[0] - power->y0[0]) /
                        (x1[2] - power->y0[2]);
                    power->y0[1] += tmp * (x1[1] - power->y0[1]) /
                        (x1[2] - power->y0[2]);
                    power->y0[2] = AX_3D_CYLINDER_VEIL * x1[2];
                }
                else if ( (power->y0[2] > 0) && (x1[2] <= 0) )
                {
                    tmp = AX_3D_CYLINDER_VEIL * power->y0[2] - x1[2];
                    x1[0] += tmp * (power->y0[0] - x1[0]) /
                        (power->y0[2] - x1[2]);
                    x1[1] += tmp * (power->y0[1] - x1[1]) /
                        (power->y0[2] - x1[2]);
                    x1[2] = AX_3D_CYLINDER_VEIL * power->y0[2];
                }
                /* intercepted by the z=AX_3D[iw].zcut plane */
                if ( (power->y0[2] > AX_3D[iw].zcut) &&
                     (x1[2] > AX_3D[iw].zcut) ) continue;
                else if ( (power->y0[2] > AX_3D[iw].zcut) &&
                          (x1[2] <= AX_3D[iw].zcut) )
                {
                    tmp = AX_3D[iw].zcut - power->y0[2];
                    power->y0[0] += tmp * (x1[0] - power->y0[0]) /
                        (x1[2] - power->y0[2]);
                    power->y0[1] += tmp * (x1[1] - power->y0[1]) /
                        (x1[2] - power->y0[2]);
                    power->y0[2] = AX_3D[iw].zcut;
                }
                else if ( (power->y0[2] <= AX_3D[iw].zcut) &&
                          (x1[2] > AX_3D[iw].zcut) )
                {
                    tmp = AX_3D[iw].zcut - x1[2];
                    x1[0] += tmp * (power->y0[0] - x1[0]) /
                        (power->y0[2] - x1[2]);
                    x1[1] += tmp * (power->y0[1] - x1[1]) /
                        (power->y0[2] - x1[2]);
                    x1[2] = AX_3D[iw].zcut;
                }
                power->U[3] = DISTANCE( x1[0] - power->y0[0],
                                        x1[1] - power->y0[1],
                                        x1[2] - power->y0[2] );
            }
            else power->U[3] = cylinder->axis[3];
            /* both centers are at front and we are not in the cylinder */
            X0[0] = AX_3D[iw].k * power->y0[0] / power->y0[2];
            X0[1] = AX_3D[iw].k * power->y0[1] / power->y0[2];
            X1[0] = AX_3D[iw].k * x1[0] / x1[2];
            X1[1] = AX_3D[iw].k * x1[1] / x1[2];
            if ( radius2 * AX_3D[iw].k/power->y0[2] * AX_3D[iw].k/x1[2] *
                 (SQUARE(X1[0]-X0[0])+SQUARE(X1[1]-X0[1])) <
                 AX_3D_CYLINDERS_MINAREA2 ) continue;
            tmp = AX_V3LENGTH2(power->y0);
            power->dx[0] = -AX_V3DOT(power->U, power->y0);
            power->dx[1] = tmp - power->dx[0]*power->dx[0] - radius2;
            power->dx[2] = AX_3D[iw].k * AX_3D[iw].k;
            /* generate perpendicular wheels */
            AX_V3CROSS (power->U, x1, UU);
            tmp = cylinder->radius / AX_V3LENGTH(UU);
            AX_V3MUL (tmp, UU);
            AX_V3CROSS (power->U, UU, VV);
            /* scan sequence */
            ENSURE (X0[0], X1[0], tmp);
            ENSURE (X0[1], X1[1], tmp);
            if ( X1[0]-X0[0] >= X1[1]-X0[1] )
            {
                for (m=0; m<AX_3D_CYLINDER_MESH; m++)
                {
                    xx[0] = power->y0[0] + (table[2*m]*UU[0] +
                                            table[2*m+1]*VV[0]);
                    xx[1] = power->y0[1] + (table[2*m]*UU[1] +
                                            table[2*m+1]*VV[1]);
                    xx[2] = power->y0[2] + (table[2*m]*UU[2] +
                                            table[2*m+1]*VV[2]);
                    if (xx[2] > 0)
                    {
                        tmp = AX_3D[iw].k * xx[0] / xx[2];
                        if (tmp > X1[0]) X1[0] = tmp;
                        else if (tmp < X0[0]) X0[0] = tmp;
                    }
                    xx[0] = x1[0] + (table[2*m]*UU[0] + table[2*m+1]*VV[0]);
                    xx[1] = x1[1] + (table[2*m]*UU[1] + table[2*m+1]*VV[1]);
                    xx[2] = x1[2] + (table[2*m]*UU[2] + table[2*m+1]*VV[2]);
                    if (xx[2] > 0)
                    {
                        tmp = AX_3D[iw].k * xx[0] / xx[2];
                        if (tmp > X1[0]) X1[0] = tmp;
                        else if (tmp < X0[0]) X0[0] = tmp;
                    }
                }
                if (X0[0] < -AX_3D[iw].wx) X0[0] = -AX_3D[iw].wx;
                if (X1[0] > AX_size[iw].width-AX_3D[iw].wx)
                    X1[0] = AX_size[iw].width-AX_3D[iw].wx;
                /* to get the y-limits */
                X0[2] = power->U[0]*X0[0]/AX_3D[iw].k+power->U[2];
                X0[3] = -power->y0[0]*X0[0]/AX_3D[iw].k-
                    power->y0[2]-power->dx[0]*X0[2];
                A = (SQUARE(power->y0[1]+power->dx[0]*power->U[1])-
                     (1-power->U[1]*power->U[1])*power->dx[1])
                    / power->dx[2];
                B = (X0[3]*(-power->y0[1]-power->dx[0]*power->U[1])+
                     X0[2]*power->U[1]*power->dx[1]) / AX_3D[iw].k;
                c = X0[3]*X0[3] - (X0[0]*X0[0]/power->dx[2]+1-X0[2]*X0[2])
                    *(power->dx[1]);
                AX_SOLVE_QUADRATIC(A,B,c,tmp,X0[2],X0[3]);
                X1[2] = power->U[0]*X1[0]/AX_3D[iw].k+power->U[2];
                X1[3] = -power->y0[0]*X1[0]/AX_3D[iw].k-power->y0[2]-
                    power->dx[0]*X1[2];
                A = (SQUARE(power->y0[1]+power->dx[0]*power->U[1])-
                     (1-power->U[1]*power->U[1])*power->dx[1])
                    / power->dx[2];
                B = (X1[3]*(-power->y0[1]-power->dx[0]*power->U[1])+
                     X1[2]*power->U[1]*power->dx[1]) / AX_3D[iw].k;
                c = X1[3]*X1[3] - (X1[0]*X1[0]/power->dx[2]+1-X1[2]*X1[2])
                    *(power->dx[1]);
                AX_SOLVE_QUADRATIC(A,B,c,tmp,X1[2],X1[3]);
                power->p->nVertex = 4;
                power->p->Vertex[0].x = power->p->Vertex[1].x =
                    AX_3D[iw].wx + X0[0];
                power->p->Vertex[0].y = AX_3D[iw].wy + X0[2];
                power->p->Vertex[1].y = AX_3D[iw].wy + X0[3]; 
                power->p->Vertex[2].x = power->p->Vertex[3].x =
                    AX_3D[iw].wx + X1[0];
                power->p->Vertex[2].y = AX_3D[iw].wy + X1[3];
                power->p->Vertex[3].y = AX_3D[iw].wy + X1[2];
            }
            else  /* X1[0]-X0[0] < X1[1]-X0[1] */
            {
                for (m=0; m<AX_3D_CYLINDER_MESH; m++)
                {
                    xx[0] = power->y0[0] + (table[2*m]*UU[0] +
                                            table[2*m+1]*VV[0]);
                    xx[1] = power->y0[1] + (table[2*m]*UU[1] +
                                            table[2*m+1]*VV[1]);
                    xx[2] = power->y0[2] + (table[2*m]*UU[2] +
                                            table[2*m+1]*VV[2]);
                    if (xx[2] > 0)
                    {
                        tmp = AX_3D[iw].k * xx[1] / xx[2];
                        if (tmp > X1[1]) X1[1] = tmp;
                        else if (tmp < X0[1]) X0[1] = tmp;
                    }
                    xx[0] = x1[0] + (table[2*m]*UU[0] + table[2*m+1]*VV[0]);
                    xx[1] = x1[1] + (table[2*m]*UU[1] + table[2*m+1]*VV[1]);
                    xx[2] = x1[2] + (table[2*m]*UU[2] + table[2*m+1]*VV[2]);
                    if (xx[2] > 0)
                    {
                        tmp = AX_3D[iw].k * xx[1] / xx[2];
                        if (tmp > X1[1]) X1[1] = tmp;
                        else if (tmp < X0[1]) X0[1] = tmp;
                    }
                }
                if (X0[1] < -AX_3D[iw].wy) X0[1] = -AX_3D[iw].wy;
                if (X1[1] > AX_size[iw].height-AX_3D[iw].wy)
                    X1[1] = AX_size[iw].height-AX_3D[iw].wy;
                /* to get the y-limits */
                X0[2] = power->U[1]*X0[1]/AX_3D[iw].k+power->U[2];
                X0[3] = -power->y0[1]*X0[1]/AX_3D[iw].k-
                    power->y0[2]-power->dx[0]*X0[2];
                A = (SQUARE(power->y0[0]+power->dx[0]*power->U[0])-
                     (1-power->U[0]*power->U[0])*power->dx[1])
                    / power->dx[2];
                B = (X0[3]*(-power->y0[0]-power->dx[0]*power->U[0])+
                     X0[2]*power->U[0]*power->dx[1]) / AX_3D[iw].k;
                c = X0[3]*X0[3] - (X0[1]*X0[1]/power->dx[2]+1-X0[2]*X0[2])
                    *(power->dx[1]);
                AX_SOLVE_QUADRATIC(A,B,c,tmp,X0[2],X0[3]);
                X1[2] = power->U[1]*X1[1]/AX_3D[iw].k+power->U[2];
                X1[3] = -power->y0[1]*X1[1]/AX_3D[iw].k-
                    power->y0[2]-power->dx[0]*X1[2];
                A = (SQUARE(power->y0[0]+power->dx[0]*power->U[0])-
                     (1-power->U[0]*power->U[0])*power->dx[1]) /
                    power->dx[2];
                B = (X1[3]*(-power->y0[0]-power->dx[0]*power->U[0])+
                     X1[2]*power->U[0]*power->dx[1]) / AX_3D[iw].k;
                c = X1[3]*X1[3] - (X1[1]*X1[1]/power->dx[2]+1-X1[2]*X1[2])
                    *(power->dx[1]);
                AX_SOLVE_QUADRATIC(A,B,c,tmp,X1[2],X1[3]);
                power->p->nVertex = 4;
                power->p->Vertex[1].y = power->p->Vertex[0].y =
                    AX_3D[iw].wy + X0[1];
                power->p->Vertex[0].x = AX_3D[iw].wx + X0[2];
                power->p->Vertex[1].x = AX_3D[iw].wx + X0[3]; 
                power->p->Vertex[2].y = power->p->Vertex[3].y =
                    AX_3D[iw].wy + X1[1];
                power->p->Vertex[2].x = AX_3D[iw].wx + X1[3];
                power->p->Vertex[3].x = AX_3D[iw].wx + X1[2];
            }
            AX_SH_PolygonWindowClip (power->p, AX_size[iw].width,
                                     AX_size[iw].height, power->p);
            if (power->p->nVertex >= 3)
            {
                C->cylinder[power-C->power] = cylinder;
                C->vz[power-C->power] = power->y0[2] + x1[2];
                power++;
            }
        }
        /*if (i < C->n_cylinders)*/
        i_l = (i < C->n_cylinders);
        p3d_reduce(p3dp_cell, &i_l, &i_g, 1, MPI_INT, MPI_MAX);
        if (i_g)
        {  /* viewpoint ejection */
            double_int di_g, di_l;
            di_l.d = 0.0;
            di_l.i = p3d_rank(p3dp_cell);
            if (i_l) di_l.d = V3LENGTH2(power->y0);
            p3d_reduce(p3dp_cell, &di_l, &di_g, 1, MPI_DOUBLE_INT, MPI_MAXLOC);
            if (di_g.i == p3d_rank(p3dp_cell)) {
                AX_V3ADD(AX_3D[iw].x, power->y0);
                while (tmp == 0)
                { /* get a random orientation */
                    AX_V3FRANDOM(power->dx);
                    tmp = AX_V3DOT(power->dx, cylinder->axis);
                    AX_V3submul(power->dx, tmp, cylinder->axis, power->y0);
                    tmp = AX_V3LENGTH2(power->y0);
                }
                tmp = (1+AX_3D_CYLINDER_REIL)*cylinder->radius / sqrt(tmp);
                AX_V3SUBMUL(AX_3D[iw].x, tmp, power->y0);
            }
            if (m>=1000)
                pe ("AX_3D_Cylinders_Zdet: viewpoint ejection failed\n");
            p3d_bcast(p3dp_cell, AX_3D[iw].x, 3, MPI_DOUBLE, di_g.i);
        }
        else break;
    }
    C->n_power = power - C->power;
    REALLOC(AX_3D_Cylinders_Zdet, C->power, C->n_power, AX_3D_Cylinder_Power);
    REALLOC(AX_3D_Cylinders_Zdet, C->cylinder, C->n_power, AX_3D_Cylinder *);
    REALLOC(AX_3D_Cylinders_Zdet, C->idx, C->n_power, int);
    AX_qsort_numerical_recipes(C->n_power, C->vz, C->idx, AX_USE_NEW_IDX);
    Free (C->vz);
    return;
} /* end AX_3D_Cylinders_Zdet() */


int p3dp_AX_3D_Cylinders_Zgrab(int iw,AX_3D_Cylinders*C,int gi,int gj,int*rank)
{
    register AX_3D_Pixel *img;
    register int m,goffset;
    int i;
    double_int di_g, di_l = {AX_INFINITY, 0};
    if (AX_IJOUWIN(iw,gi,gj)) goto not_found;
    goffset = AX_OFF(iw,gi,gj);
    for (i=0; i<C->n_power; i++)
    {
        img = C->cylinder[i]->img;
        for (m=1; m<img[0].offset; m++)
            if (goffset == img[m].offset)
                if (img[m].z <= AX_zmap[iw][goffset])
                    goto found;
        for (m=img[0].offset+1;
             m<img[0].offset+img[img[0].offset].offset; m++)
            if (goffset == img[m].offset)
                if (img[m].z <= AX_zmap[iw][goffset])
                    goto found;
    }
    goto not_found;

found:
    di_l.d = img[m].z;
not_found:
    di_l.i = p3d_rank(p3dp_cell);
    p3d_reduce(p3dp_cell, &di_l, &di_g, 1, MPI_DOUBLE_INT, MPI_MINLOC);
    if (di_g.d == AX_INFINITY)
        return -1;
    else {
        int idx;
        *rank = di_g.i;
        if (p3d_rank(p3dp_cell) == *rank)
            idx = C->cylinder[i]-C->CYLINDER;
        p3d_bcast(p3dp_cell, &idx, 1, MPI_INT, *rank);
        return idx;
    }
} /* end AX_3D_Cylinders_Zgrab() */


#undef  AX_openwindow
#undef  AX_resizewindow
#undef  AX_closewindow
#define nonvacant(iw) (AX_cid[iw]!=0)

static int p3dp_AX_plugin_ShmGC_module (int iw)
{
    AX_AssertWinID ("p3dp_AX_plugin_ShmGC_module", iw);
    AX_Check_Module_Dependence ("p3dp_AX_plugin_ShmGC_module", iw,MODULE_ShmGC);
    if ( AX_Module_Is_PluggedIn [iw] [MODULE_ShmGC] )
        pe ("p3dp_AX_plugin_ShmGC_module:\n"
            "you are crazy, this module is already plugged in\n");

    switch (AX_depth[iw]) {
    case 8:  AX_bytes = 1; break;
    case 16: AX_bytes = 2; break;
    case 24: case 32: AX_bytes = 4; break;
    default: AX_bytes = 4;
    }

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
} /* end p3dp_AX_plugin_ShmGC_module() */


static int p3dp_AX_plugout_ShmGC_module (int iw)
{
    AX_AssertWinID ("p3dp_AX_plugout_ShmGC_module", iw);
    if ( ! AX_Module_Is_PluggedIn [iw] [MODULE_ShmGC] )
        pe ("p3dp_AX_plugout_ShmGC_module:\n"
            "you are crazy, this module is not plugged in\n");

    Free(AX_mem[iw].uc);
    AX_Module_Is_PluggedIn [iw] [MODULE_ShmGC] = 0;
    return (1);
} /* end p3dp_AX_plugout_ShmGC_module() */


static int p3dp_AX_resize_ShmGC_module (int iw)
{
    REALLOC( AX_plugin_ShmGC_module, AX_mem[iw].uc,
             AX_size[iw].width * AX_size[iw].height * AX_bytes,
             unsigned char );
    return (1);
} /* end p3dp_AX_resize_ShmGC_module() */


static struct {
    int (*plugin)  (int iw);
    int (*plugout) (int iw);
    int (*resize)  (int iw);
} p3dp_AX_modules [MODULE_MAX] = {
    {&p3dp_AX_plugin_ShmGC_module,
                        &p3dp_AX_plugout_ShmGC_module,
                                                 &p3dp_AX_resize_ShmGC_module},
    {&AX_plugin_Scan_module, &AX_plugout_Scan_module, &AX_resize_Scan_module},
    {&AX_plugin_3D_module,   &AX_plugout_3D_module,   &AX_resize_3D_module},
    {&AX_plugin_BC_module,   &AX_plugout_BC_module,   &AX_resize_BC_module}
};


int p3dp_AX_openwindow
(int cid, char *title, unsigned int width, unsigned int height)
{
    int iw, minor;
    AX_Carrier mask[3];
    AX_Pixel   pixel[3];

    if (cui_enabled < 0)
        return cui_AX_openwindow(-1*cid, title, width, height);
    else if (p3dp_enabled <= 0)
        return AX_openwindow(cid, title, width, height);

    if (p3d_rank(p3dp_cell) == 0) {
        iw = AX_openwindow(cid, title, width, height);
        mask[0] = AX_rmask;
        mask[1] = AX_gmask;
        mask[2] = AX_bmask;
        pixel[0] = AX_namedpixel[AX_RED];
        pixel[1] = AX_namedpixel[AX_GREEN];
        pixel[2] = AX_namedpixel[AX_BLUE];
        p3d_bcast(p3dp_cell, &AX_depth[iw], 1, MPI_INT, 0);
        p3d_bcast(p3dp_cell, mask,  3*sizeof(AX_Carrier), MPI_CHAR, 0);
        p3d_bcast(p3dp_cell, pixel, 3*sizeof(AX_Pixel),   MPI_CHAR, 0);
    }
    else {
        for (iw = 0; (iw < AX_MAXWIN) && nonvacant(iw); iw++);
        if (iw == AX_MAXWIN) {
            fprintf(stderr,"cui_AX_openwindow: AX_MAXWIN = %d reached\n",
                    AX_MAXWIN);
            return -1;
        }
        
        AX_display[iw]          = NULL;
        AX_cid[iw]              = cid;
        AX_size[iw].width       = width;
        AX_size[iw].height      = height;
        AX_borderwidth[iw]      = 0;
        AX_noneedcolormap       = 1;
        AX_videomode            = AX_VIDEOMODE_NO_SHM;

        p3d_bcast(p3dp_cell, &AX_depth[iw], 1, MPI_INT, 0);
        p3d_bcast(p3dp_cell, mask,  3*sizeof(AX_Carrier), MPI_CHAR, 0);
        p3d_bcast(p3dp_cell, pixel, 3*sizeof(AX_Pixel),   MPI_CHAR, 0);

        AX_rmask = mask[0];
        AX_gmask = mask[1];
        AX_bmask = mask[2];
        AX_namedpixel[AX_RED]   = pixel[0];
        AX_namedpixel[AX_GREEN] = pixel[1];
        AX_namedpixel[AX_BLUE]  = pixel[2];

        for (minor=0; minor<=AX_DEF_Module_Support_Level; minor++)
            p3dp_AX_modules[minor].plugin(iw);
    }

    return iw;
}

/* Resize all modules that are already plugged in, to AX_size[iw] */
void p3dp_AX_resizewindow (int iw, Bool do_window_too)
{
    if (cui_enabled < 0)
        cui_AX_resizewindow(iw, do_window_too);
    else if (p3dp_enabled <= 0 || p3d_rank(p3dp_cell) == 0)
        AX_resizewindow(iw, do_window_too);
    else {
        int module;
        AX_AssertWinID ("AX_resizewindow", iw);
        for (module=0; module<MODULE_MAX; module++)
            if (AX_Module_Is_PluggedIn[iw][module])
                p3dp_AX_modules[module].resize(iw);
    }
    return;
} /* end AX_resizewindow() */


/* close the window and free all resources */
void p3dp_AX_closewindow (int iw)
{
    if (cui_enabled < 0)
        cui_AX_closewindow(iw);
    else if (p3dp_enabled <= 0 || p3d_rank(p3dp_cell) == 0)
        AX_closewindow(iw);
    else {
        int jw;
        AX_AssertWinID ("AX_closewindow", iw);
        for (jw=0; jw<MODULE_MAX; jw++)
            if (AX_Module_Is_PluggedIn [iw] [jw])
                p3dp_AX_modules[jw].plugout(iw);
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
    }
    return;
} /* end AX_closewindow() */
