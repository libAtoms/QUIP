/********************************************/
/* libAX: -lX11 -lXext -lpng -lz -ljpeg -lm */
/*        -lScalar -lIO -lTimer             */
/*                                          */
/* Accelerated low-level graphics library   */
/* with X Window Shared Memory Extension.   */
/*                                          */
/* Jan.13, 2000 Ju Li <liju99@mit.edu>      */
/********************************************/

#include "AX.h"

/***************/
/* 3D Renderer */
/***************/

AX_3D_Viewport AX_3D [AX_MAXWIN];

/* Assign default observer (viewport) AX_3D[iw] */
void AX_3D_DEFViewport (int iw)
{
    int j;
    AX_V3ZERO(AX_3D[iw].x);
    AX_3D[iw].wx = AX_size[iw].width / 2.;
    AX_3D[iw].wy = AX_size[iw].height / 2.;
    AX_3D[iw].k = sqrt( AX_3D[iw].wx * AX_3D[iw].wy ) /
        tan( DEGREE_TO_RADIAN(AX_3D_DEF_VIEW_ANGLE_IN_DEGREE/2.) );
    AX_M3IDENTITY(AX_3D[iw].V);
    AX_3D[iw].zcut = AX_INFINITY;
    for (j=0; j<AX_3D_MAX_FILTER_PLANE; j++)
        AX_V3ZERO(AX_3D[iw].fp[j].dx);
    return;
} /* end AX_3D_DEFViewport() */

/* floating point z-buffer */
AX_Float *AX_zmap[AX_MAXWIN];

/* Allocate AX_zmap[iw] according to AX_size[iw] and set viewport defaults */
int AX_plugin_3D_module (int iw)
{
    AX_AssertWinID ("AX_plugin_3D_module", iw);
    AX_Check_Module_Dependence ("AX_plugin_3D_module", iw, MODULE_3D);
    if ( AX_Module_Is_PluggedIn [iw] [MODULE_3D] )
        pe ("AX_plugin_3D_module:\n"
            "you are crazy, this module is already plugged in\n");
    MALLOC ( AX_plugin_3D_module, AX_zmap[iw],
             AX_size[iw].width * AX_size[iw].height, AX_Float );
    AX_3D_DEFViewport (iw);
    AX_Module_Is_PluggedIn [iw] [MODULE_3D] = 1;
    return(1);
} /* end AX_plugin_3D_module() */


/* Free z-buffer memory allocation AX_zmap[iw] */
int AX_plugout_3D_module (int iw)
{
    AX_AssertWinID ("AX_plugout_3D_module", iw);
    if ( !AX_Module_Is_PluggedIn [iw] [MODULE_3D] )
        pe ("AX_plugout_3D_module:\n"
            "you are crazy, this module is not plugged in\n");
    Free(AX_zmap[iw]);
    AX_Module_Is_PluggedIn [iw] [MODULE_3D] = 0;
    return(1);
} /* end AX_plugout_3D_module() */


/* Resize AX_zmap[iw] according to current AX_size[iw] by */
/* realloc(), and reassign viewport parameters sensibly.  */
int AX_resize_3D_module (int iw)
{
    AX_AssertWinID ("AX_resize_3D_module", iw);
    if ( ! AX_Module_Is_PluggedIn [iw] [MODULE_3D] )
        pe ("AX_resize_3D_module:\n"
            "you are crazy, this module is not plugged in\n");
    AX_zmap[iw] = realloc ( AX_zmap[iw], AX_size[iw].width *
                            AX_size[iw].height * sizeof(AX_Float) );
    if (AX_zmap[iw] == NULL)
        pe ("AX_resize_3D_module: realloc failed\n");
    AX_3D[iw].k *= sqrt ( AX_size[iw].width / 2. *
                          AX_size[iw].height / 2. /
                          AX_3D[iw].wx / AX_3D[iw].wy );
    AX_3D[iw].wx = AX_size[iw].width  / 2.;
    AX_3D[iw].wy = AX_size[iw].height / 2.;
    return(1);
} /* end AX_resize_3D_module() */


/****************/
/* Line Support */
/****************/


/* Realloc n_lines; if used as malloc, assign L->LINE=NULL first */
void AX_3D_Lines_Realloc (AX_3D_Lines *L, int n_lines)
{
    L->n_lines = n_lines;
    REALLOC (AX_3D_Lines_Realloc, L->LINE, n_lines, AX_3D_Line);
    return;
} /* end AX_3D_Lines_Realloc() */


/* Realloc n_lines followed by x0,y0,z0,x1,y1,z1,r,g,b, ... interface */
void AX_3D_Lines_Recreate (AX_3D_Lines *L, int n_lines, ...)
{
    register int i;
    va_list ap;
    AX_3D_Lines_Realloc (L, n_lines);
    va_start(ap, n_lines);
    for (i=0; i<n_lines; i++)
    {
        L->LINE[i].x0[0] = va_arg(ap, double);
        L->LINE[i].x0[1] = va_arg(ap, double);
        L->LINE[i].x0[2] = va_arg(ap, double);
        L->LINE[i].x1[0] = va_arg(ap, double);
        L->LINE[i].x1[1] = va_arg(ap, double);
        L->LINE[i].x1[2] = va_arg(ap, double);
        AX_3D_AssignRGB (L->LINE[i], va_arg(ap, double), va_arg(ap, double),
                         va_arg(ap, double));
    }
    va_end(ap);
    return;
} /* end AX_3D_Lines_Recreate() */


/* free allocated memory and set NULL */
void AX_3D_Lines_Free (AX_3D_Lines *L)
{
    Free (L->LINE);
    return;
} /* end AX_3D_Lines_Free() */


#define AX_3D_LINE_VEIL 0.001
/* Generate line scan and draw to pixmap alias pixels according to z-buffer */
void AX_3D_Lines_Zprint (int iw, AX_3D_Lines *L)
{
    int i;
    AX_3D_Line *line;
    AX_Float x0[3], x1[3], dx[3];
    register int m;
    register AX_Float tmp;
    register AX_Offset offset;
    register AX_Carrier carrier;
    for (i=0; i<L->n_lines; i++)
    {
        line = L->LINE + i;
        if (line->r < 0) continue;  /* invisibility flag */
        AX_V3sub(line->x0, AX_3D[iw].x, dx);
        AX_M3mulV3 (AX_3D[iw].V, dx, x0);
        AX_V3sub(line->x1, AX_3D[iw].x, dx);
        AX_M3mulV3 (AX_3D[iw].V, dx, x1);
        AX_V3sub(x1, x0, dx);
        if ((x0[2] <= 0) && (x1[2] <= 0)) continue;
        else if ( (x0[2] <= 0) && (x1[2] > 0) )
        {
            tmp = AX_3D_LINE_VEIL * x1[2] - x0[2];
            x0[0] += tmp * dx[0] / dx[2];
            x0[1] += tmp * dx[1] / dx[2];
            x0[2] = AX_3D_LINE_VEIL * x1[2];
        }
        else if ( (x0[2] > 0) && (x1[2] <= 0) )
        {
            tmp = AX_3D_LINE_VEIL * x0[2] - x1[2];
            x1[0] += tmp * dx[0] / dx[2];
            x1[1] += tmp * dx[1] / dx[2];
            x1[2] = AX_3D_LINE_VEIL * x0[2];
        }
        if ( (x0[2] > AX_3D[iw].zcut) && (x1[2] > AX_3D[iw].zcut) )
            continue;
        else if ( (x0[2] > AX_3D[iw].zcut) && (x1[2] <= AX_3D[iw].zcut) )
        {
            tmp = AX_3D[iw].zcut - x0[2];
            x0[0] += tmp * dx[0] / dx[2];
            x0[1] += tmp * dx[1] / dx[2];
            x0[2] = AX_3D[iw].zcut;
        }
        else if ( (x0[2] <= AX_3D[iw].zcut) && (x1[2] > AX_3D[iw].zcut) )
        {
            tmp = AX_3D[iw].zcut - x1[2];
            x1[0] += tmp * dx[0] / dx[2];
            x1[1] += tmp * dx[1] / dx[2];
            x1[2] = AX_3D[iw].zcut;
        }
        line->X0[0] = AX_3D[iw].k * x0[0] / x0[2] + AX_3D[iw].wx;
        line->X0[1] = AX_3D[iw].k * x0[1] / x0[2] + AX_3D[iw].wy;
        line->X1[0] = AX_3D[iw].k * x1[0] / x1[2] + AX_3D[iw].wx;
        line->X1[1] = AX_3D[iw].k * x1[1] / x1[2] + AX_3D[iw].wy;
        if ( (line->X0[0]==line->X1[0]) && (line->X0[1]==line->X1[1]) )
            continue;
        if ( ! AX_LineClipRectangle
             ( &line->X0[0], &line->X0[1], &line->X1[0], &line->X1[1],
               1, 1, (AX_size[iw].width-1)*(1-AX_TINY),
               (AX_size[iw].height-1)*(1-AX_TINY) ) ) continue;
        AX_GSLineScan (line->X0[0], line->X0[1], line->X1[0], line->X1[1],
                       AX_AOP_top(iw));
        /* normalized line direction */
        AX_V3NORMALIZE(dx, tmp);
        carrier = AX_Colorcarrier(line->r, line->g, line->b);

        AX_C(

            for (m=1; m<AX_AOP_TOP(iw,0).b.offset; m++)
            {
                offset = AX_OFF
                    (iw, AX_AOP_TOP(iw,m).b.p.i, AX_AOP_TOP(iw,m).b.p.j);
                x1[0] = (AX_AOP_TOP(iw,m).b.p.i + 0.5 - AX_3D[iw].wx)
                    / AX_3D[iw].k;
                x1[1] = (AX_AOP_TOP(iw,m).b.p.j + 0.5 - AX_3D[iw].wy)
                    / AX_3D[iw].k;
                x1[2] = 1;
                tmp = AX_V3DOT(x1, dx);
                AX_V3SUBMUL(x1, tmp, dx);
                tmp = AX_V3DOT(x1,x0) / AX_V3LENGTH2(x1);
                if (tmp < AX_zmap[iw][offset])
                    AX_mem[iw].uc[offset] = AX_mix
                        ( AX_mem[iw].uc[offset],
                          carrier, AX_AOP_TOP(iw,m).c.area );
            },

            for (m=1; m<AX_AOP_TOP(iw,0).b.offset; m++)
            {
                offset = AX_OFF
                    (iw, AX_AOP_TOP(iw,m).b.p.i, AX_AOP_TOP(iw,m).b.p.j);
                x1[0] = (AX_AOP_TOP(iw,m).b.p.i + 0.5 - AX_3D[iw].wx)
                    / AX_3D[iw].k;
                x1[1] = (AX_AOP_TOP(iw,m).b.p.j + 0.5 - AX_3D[iw].wy)
                    / AX_3D[iw].k;
                x1[2] = 1;
                tmp = AX_V3DOT(x1,dx);
                AX_V3SUBMUL(x1,tmp,dx);
                tmp = AX_V3DOT(x1,x0) / AX_V3LENGTH2(x1);
                if (tmp < AX_zmap[iw][offset])
                    AX_mem[iw].i2[offset] = AX_mix
                        ( AX_mem[iw].i2[offset],
                          carrier, AX_AOP_TOP(iw,m).c.area );
            },

            for (m=1; m<AX_AOP_TOP(iw,0).b.offset; m++)
            {
                offset = AX_OFF
                    (iw, AX_AOP_TOP(iw,m).b.p.i, AX_AOP_TOP(iw,m).b.p.j);
                x1[0] = (AX_AOP_TOP(iw,m).b.p.i + 0.5 - AX_3D[iw].wx)
                    / AX_3D[iw].k;
                x1[1] = (AX_AOP_TOP(iw,m).b.p.j + 0.5 - AX_3D[iw].wy)
                    / AX_3D[iw].k;
                x1[2] = 1;
                tmp = AX_V3DOT(x1, dx);
                AX_V3SUBMUL(x1, tmp, dx);
                tmp = AX_V3DOT(x1,x0) / AX_V3LENGTH2(x1);
                if (tmp < AX_zmap[iw][offset])
                    AX_mem[iw].i4[offset] = AX_mix
                        ( AX_mem[iw].i4[offset],
                          carrier, AX_AOP_TOP(iw,m).c.area );
            }

            );
    }
    return;
} /* end AX_3D_Lines_Zprint() */


/* Assign parallelepiped wireframe H[][] originated at x0,y0,z0 */
AX_3D_Lines *AX_3D_PW_Assign
(AX_Float x0, AX_Float y0, AX_Float z0, AX_Float H[3][3], AX_Float r,
 AX_Float g, AX_Float b, AX_3D_Lines *L)
{
    register int i;
    for (i=0; i<12; i++)
        AX_3D_AssignRGB (L->LINE[i], r,g,b);
    L->LINE[0].x0[0] = x0;
    L->LINE[0].x0[1] = y0;
    L->LINE[0].x0[2] = z0;
    L->LINE[1].x0[0] = L->LINE[0].x1[0] = x0 + H[0][0];
    L->LINE[1].x0[1] = L->LINE[0].x1[1] = y0 + H[0][1];
    L->LINE[1].x0[2] = L->LINE[0].x1[2] = z0 + H[0][2];
    L->LINE[2].x0[0] = L->LINE[1].x1[0] = L->LINE[0].x1[0] + H[1][0];
    L->LINE[2].x0[1] = L->LINE[1].x1[1] = L->LINE[0].x1[1] + H[1][1];
    L->LINE[2].x0[2] = L->LINE[1].x1[2] = L->LINE[0].x1[2] + H[1][2];
    L->LINE[3].x0[0] = L->LINE[2].x1[0] = x0 + H[1][0];
    L->LINE[3].x0[1] = L->LINE[2].x1[1] = y0 + H[1][1];
    L->LINE[3].x0[2] = L->LINE[2].x1[2] = z0 + H[1][2];
    L->LINE[8].x0[0] = L->LINE[3].x1[0] = x0;
    L->LINE[8].x0[1] = L->LINE[3].x1[1] = y0;
    L->LINE[8].x0[2] = L->LINE[3].x1[2] = z0;
    L->LINE[8].x1[0] = x0 + H[2][0];
    L->LINE[8].x1[1] = y0 + H[2][1];
    L->LINE[8].x1[2] = z0 + H[2][2];
    for (i=0; i<4; i++)
    {
        AX_V3add(L->LINE[i].x0, H[2], L->LINE[4+i].x0);
        AX_V3add(L->LINE[i].x1, H[2], L->LINE[4+i].x1);
    }
    AX_V3add(L->LINE[8].x0, H[0], L->LINE[9].x0);
    AX_V3add(L->LINE[8].x1, H[0], L->LINE[9].x1);
    AX_V3add(L->LINE[9].x0, H[1], L->LINE[10].x0);
    AX_V3add(L->LINE[9].x1, H[1], L->LINE[10].x1);
    AX_V3add(L->LINE[8].x0, H[1], L->LINE[11].x0);
    AX_V3add(L->LINE[8].x1, H[1], L->LINE[11].x1);
    return (L);
} /* end AX_3D_PW_Assign() */


#ifdef _3D_Lines_TEST
#define Z0 4
#define boxwidth  2
#define boxheight 2
#define boxthickness 2
#define PNG_FILENAME "/tmp/lines.png"
#define EPS_FILENAME "/tmp/lines.eps"
int main (int argc, char *argv[])
{
    AX_Float H[3][3],R[3][3],A[3][3];
    AX_3D_Define_PW(L);
    AX_M3RANDOMROTATION(R,H[0][0]);
    AX_M3TRANSPOSE(R,H);
    AX_M3MUL(R,H,A);
    AX_M3DIAGONAL(boxwidth,boxheight,boxthickness,A);
    AX_M3MUL(A,R,H);
    /* printf ("%f %f %f\n", H[0][0], H[1][1], H[0][1]); */
    AX_openwindow (getpid(), NULL, AX_DEFWIDTH, AX_DEFHEIGHT);
    AX_PLUGIN_3D_module(0);
    AX_namedmop(0,AX_WHITE);
    AX_3D_Lines_Zprint (0, AX_3D_PW_Assign
                        (-boxwidth/2, -boxheight/2, Z0, H, 1, 0, 0, L));
    AX_dump(0); AX_show(0);
    AXPressKey(0);
    AX_save_pixmap_as_PNG(0, PNG_FILENAME);
    printf ("Image saved on \"%s\".\n", PNG_FILENAME);
    AX_save_pixmap_as_EPS(0, EPS_FILENAME);
    printf ("Image saved on \"%s\".\n", EPS_FILENAME);
    AX_closewindow(0);
    Press_return();
    return(0);
}
#endif /* _3D_Lines_TEST */


/*******************/
/* Polygon Support */
/*******************/

/* (Re)assign 3D polygon by x0,y0,z0, x1,y1,z1, ... interface */
AX_3D_Polygon *AX_3D_PolygonAssign
(AX_3D_Polygon *p, int nVertex, AX_Float r, AX_Float g, AX_Float b, ...)
{
    register int i;
    AX_Float x0[3], x1[3], tmp=0;
    va_list ap;
    if (nVertex > AX_3D_PolygonMaxVertex)
        pe ("AX_3D_PolygonAssign: nVertex = %d exceeds "
            "AX_3D_PolygonMaxVertex = %d\n", nVertex, AX_3D_PolygonMaxVertex);
    p->nVertex = nVertex;
    AX_3D_AssignRGB (*p, r, g, b);
    va_start(ap, b);
    for (i=0; i<nVertex; i++)
    {
        p->x[3*i]   = va_arg(ap, double);
        p->x[3*i+1] = va_arg(ap, double);
        p->x[3*i+2] = va_arg(ap, double);
    }
    for (i=0; i<nVertex; i++)
    { /* work out the surface normal by triangulating 0,1,n-1 vertices */
        AX_V3sub ( &(p->x[3*((i+1)%nVertex)]), &(p->x[3*i]), x0 );
        AX_V3sub ( &(p->x[3*((i+nVertex-1)%nVertex)]), &(p->x[3*i]), x1 );
        AX_V3CROSS (x0, x1, p->front_normal);
        tmp = AX_V3LENGTH2(p->front_normal);
        if (tmp != 0) break;
    }
    if ( i == nVertex )
        pe ("AX_3D_PolygonAssign: polygon has zero area.\n");
    else
    {
        tmp = sqrt(tmp);
        AX_V3DIV (p->front_normal, tmp);
    }
    va_end (ap);
    return (p);
} /* end AX_3D_PolygonAssign() */


/* (Re)assign 3D polygon by *x0, *x1, ... interface */
AX_3D_Polygon *AX_3D_PolygonASSIGN
(AX_3D_Polygon *p, int nVertex, AX_Float r, AX_Float g, AX_Float b, ...)
{
    register int i;
    AX_Float *xp, x0[3], x1[3], tmp=0;
    va_list ap;
    if (nVertex > AX_3D_PolygonMaxVertex)
        pe ("AX_3D_PolygonASSIGN: nVertex = %d exceeds "
            "AX_3D_PolygonMaxVertex = %d\n", nVertex, AX_3D_PolygonMaxVertex);
    p->nVertex = nVertex;
    AX_3D_AssignRGB (*p, r, g, b);
    va_start(ap, b);
    for (i=0; i<nVertex; i++)
    {
        xp = va_arg(ap,AX_Float *);
        AX_V3EQV(xp, &(p->x[3*i]));
    }
    for (i=0; i<nVertex; i++)
    { /* work out the surface normal by triangulating 0,1,n-1 vertices */
        AX_V3sub ( &(p->x[3*((i+1)%nVertex)]), &(p->x[3*i]), x0 );
        AX_V3sub ( &(p->x[3*((i+nVertex-1)%nVertex)]), &(p->x[3*i]), x1 );
        AX_V3CROSS (x0, x1, p->front_normal);
        tmp = AX_V3LENGTH2(p->front_normal);
        if (tmp != 0) break;
    }
    if ( i == nVertex )
        pe ("AX_3D_PolygonASSIGN: polygon has zero area.\n");
    else
    {
        tmp = sqrt(tmp);
        AX_V3DIV (p->front_normal,tmp);
    }
    va_end (ap);
    return (p);
} /* end AX_3D_PolygonASSIGN() */


/* Sutherland-Hodgman 3D polygon/plane (n*x>=n[3]) clipping algorithm */
void AX_3D_SH_PolygonClipPlane
(AX_3D_polygon *p, AX_Float n[4], AX_3D_polygon *result)
{
    register int j;
    AX_3D_Polygon from[1];
    AX_Float *v0,*v1,d0,d1,dx[3];
    AX_3D_PolygonCopy(p,from,j);
    result->nVertex = 0;
    /* start with the last vertex of the polygon */
    v0 = from->x + 3*(from->nVertex-1);
    for (j=0; j<from->nVertex; j++)
    {
        v1 = from->x + 3*j;
        if  ( (d1=AX_V3DOT(v1,n)) >= n[3] )
        { /* case 1 and 4 */
            if ( (d0=AX_V3DOT(v0,n)) >= n[3] )
                AX_3D_PolygonAddVertex (v1, *result);  /* case 1 */
            else
            { /* case 4 */
                AX_V3sub(v1,v0,dx);
                d0 = (n[3]-d0)/(d1-d0);
                AX_V3addmul(v0,d0,dx,&(result->x[3*result->nVertex]));
                result->nVertex++;
                AX_3D_PolygonAddVertex(v1,*result);
            }
        }
        else /* case 2 and 3 */
            if ( (d0=AX_V3DOT(v0,n)) >= n[3] )
            { /* case 2 */
                AX_V3sub(v1,v0,dx);
                d0 = (n[3]-d0)/(d1-d0);
                AX_V3addmul(v0,d0,dx,&(result->x[3*result->nVertex]));
                result->nVertex++;
            }  /* no action for case 3 */
        v0 = v1;
    } /* for j */
    return;
} /* end AX_3D_SH_PolygonClipPlane() */


/* Realloc n_polygons; if used as malloc, let P->POLYGON=P->idx=NULL first */
void AX_3D_Polygons_Realloc (AX_3D_Polygons *P, int n_polygons)
{
    int i;
    P->n_polygons = P->n_shown = n_polygons;
    REALLOC (AX_3D_Polygons_Realloc, P->POLYGON, n_polygons, AX_3D_Polygon);
    REALLOC (AX_3D_Polygons_Realloc, P->vz, n_polygons, AX_Float);
    REALLOC (AX_3D_Polygons_Realloc, P->idx, n_polygons, int);
    for (i=0; i<n_polygons; i++) P->idx[i] = i;
    return;
} /* end AX_3D_Polygons_Realloc() */


/* free allocated memory and set NULL */
void AX_3D_Polygons_Free (AX_3D_Polygons *P)
{
    Free (P->POLYGON);
    Free (P->vz);
    Free (P->idx);
    return;
} /* end AX_3D_Polygons_Free() */


#define AX_3D_POLYGON_VEIL                  0.001
#define AX_3D_POLYGON_EDGE_INFLUENCE_VEIL   0.01
#define AX_3D_POLYGON_BACKFACE_RETARDATION  0.25

/* Generate polygon scans and register them in the z-buffer */
void AX_3D_Polygons_Zpaper (int iw, AX_3D_Polygons *P)
{
    register AX_3D_Polygon *polygon;
    register int j,m;
    register AX_Float tmp;
    int i, n_shown;
    AX_Float dx[4], U[3];
    AX_Float edge_influence, lr,lg,lb, light_influence;
    for (P->n_shown=i=0; i<P->n_polygons; i++)
    {
        polygon = P->POLYGON + i;
        m = 0;    /* number of vertices with z <= 0 */
        tmp = 0;  /* max z of those vertices with z > 0 */
        polygon->q->nVertex = 0;
        /* number of vertices with z > AX_3D[iw].zcut */
        for (j=0; j<polygon->nVertex; j++)
        {
            AX_V3sub (&polygon->x[3*j], AX_3D[iw].x, dx);
            AX_M3mulV3 (AX_3D[iw].V, dx, &polygon->q->x[3*j]);
            if (polygon->q->x[3*j+2] <= 0) m++;
            else
            {
                if (polygon->q->x[3*j+2] > tmp)
                    tmp = polygon->q->x[3*j+2];
                if (polygon->q->x[3*j+2] > AX_3D[iw].zcut)
                    polygon->q->nVertex++;
            }
        }
        if ( (m == polygon->nVertex) ||
             (polygon->q->nVertex == polygon->nVertex) ) continue;
        if (polygon->q->nVertex)
        {
            polygon->q->nVertex = polygon->nVertex;
            dx[0] = dx[1] = 0;
            dx[2] = -1;
            dx[3] = -AX_3D[iw].zcut;
            AX_3D_SH_PolygonClipPlane (polygon->q, dx, polygon->q);
        }
        else polygon->q->nVertex = polygon->nVertex;
        if ( m )
        {
            dx[0] = dx[1] = 0;
            dx[2] = 1;
            dx[3] = tmp * AX_3D_POLYGON_VEIL;
            AX_3D_SH_PolygonClipPlane (polygon->q, dx, polygon->q);
        }
        /* then clip against the window */
        polygon->p->nVertex = polygon->q->nVertex;
        for (j=0; j<polygon->q->nVertex; j++)
        {
            polygon->p->Vertex[j].x = AX_3D[iw].wx + AX_3D[iw].k *
                polygon->q->x[3*j] / polygon->q->x[3*j+2];
            polygon->p->Vertex[j].y = AX_3D[iw].wy + AX_3D[iw].k *
                polygon->q->x[3*j+1] / polygon->q->x[3*j+2];
        }
        AX_SH_PolygonWindowClip
            (polygon->p, (AX_Float)AX_size[iw].width,
             (AX_Float)AX_size[iw].height, polygon->p);
        if (polygon->p->nVertex >= 3)
        {
            for (P->vz[i]=j=0; j<polygon->q->nVertex; j++)
                P->vz[i] += polygon->q->x[3*j+2];
            P->vz[i] /= polygon->q->nVertex;
            P->idx[P->n_shown++] = i;
        }
    }
    AX_qsort_numerical_recipes
        (P->n_shown, P->vz, P->idx, AX_USE_OLD_IDX);
    n_shown = P->n_shown;
    for (P->n_shown=i=0; i<n_shown; i++)
    {
        polygon = P->POLYGON + P->idx[i];
        /* front surface normal */
        AX_M3mulV3 (AX_3D[iw].V, polygon->front_normal, dx);
        dx[3] = AX_V3DOT(dx, polygon->q->x); /* plane determinant */
        AX_PolygonScan (polygon->p, AX_AOP_top(iw), AX_BOP_top(iw));
        polygon->IMG = (struct AX_3D_PIXEL *) malloc
            ( AX_AOP_TOP(iw,0).b.offset * sizeof(struct AX_3D_PIXEL) +
              AX_BOP_TOP(iw,0).offset * sizeof(AX_3D_Pixel) );
        polygon->img = (AX_3D_Pixel *)
            ( ((char *)polygon->IMG) +
              AX_AOP_TOP(iw,0).b.offset * sizeof(struct AX_3D_PIXEL) );
        if (dx[2] < 0)  /* front surface */
            AX_3D_LRGB (AX_PS_TOP*polygon->r, AX_PS_TOP*polygon->g,
                        AX_PS_TOP*polygon->b, lr, lg, lb);
        else  /* back surface */
            AX_3D_LRGB (AX_SP_TOP*polygon->r, AX_SP_TOP*polygon->g,
                        AX_SP_TOP*polygon->b, lr, lg, lb);
        U[2] = 1; /* direction of view-sight */
        for (j=m=1; m<AX_BOP_TOP(iw,0).offset; m++)
        {
            polygon->img[j].offset = AX_OFF
                (iw, AX_BOP_TOP(iw,m).p.i, AX_BOP_TOP(iw,m).p.j);
            U[0] = (AX_BOP_TOP(iw,m).p.i + 0.5 - AX_3D[iw].wx)
                / AX_3D[iw].k;
            U[1] = (AX_BOP_TOP(iw,m).p.j + 0.5 - AX_3D[iw].wy)
                / AX_3D[iw].k;
            edge_influence = AX_V3DOT(U, dx);
            polygon->img[j].z = dx[3] / edge_influence;
            if (polygon->img[j].z < AX_zmap[iw][polygon->img[j].offset])
            {
                edge_influence /= -AX_V3LENGTH(U);
                if (edge_influence > 0)
                {
                    light_influence = AX_PS_LIGHT * edge_influence *
                        AX_3D_LDOT(dx);
                    edge_influence = AX_PS_EDGE + AX_PS_TOPI *
                        edge_influence;
                }
                else
                {
                    light_influence = AX_SP_LIGHT * edge_influence *
                        AX_3D_LDOT(dx);
                    edge_influence = AX_SP_EDGE - AX_SP_TOPI *
                        edge_influence;
                }
                if (light_influence < 0) light_influence = 0;
                polygon->img[j].c.carrier = AX_Colorcarrier
                    ( edge_influence * polygon->r + light_influence * lr,
                      edge_influence * polygon->g + light_influence * lg,
                      edge_influence * polygon->b + light_influence * lb );
                AX_zmap[iw][polygon->img[j].offset] = polygon->img[j].z;
                j++;
            }
        }
        if (j-1 <= -1)
        {
            free (polygon->IMG);
            continue;
        }
        else
        {
            polygon->img[0].offset = j;
            P->idx[P->n_shown++] = P->idx[i];
        }
        for (j=m=1; m<AX_AOP_TOP(iw,0).b.offset; m++)
        {
            polygon->IMG[j].o.offset = AX_OFF
                (iw, AX_AOP_TOP(iw,m).b.p.i, AX_AOP_TOP(iw,m).b.p.j);
            U[0] = (AX_AOP_TOP(iw,m).b.p.i+0.5 - AX_3D[iw].wx)
                / AX_3D[iw].k;
            U[1] = (AX_AOP_TOP(iw,m).b.p.j+0.5 - AX_3D[iw].wy)
                / AX_3D[iw].k;
            edge_influence = AX_V3DOT(U,dx);
            polygon->IMG[j].z = dx[3] / edge_influence;
            if (polygon->IMG[j].z < AX_zmap[iw][polygon->IMG[j].o.offset])
            {
                edge_influence /= -AX_V3LENGTH(U);
                if (edge_influence > 0)
                { /* screen surface pixels that's too oblique with viewsight */
                    if (edge_influence < AX_3D_POLYGON_EDGE_INFLUENCE_VEIL)
                        continue;
                    light_influence = AX_PS_LIGHT * edge_influence *
                        AX_3D_LDOT(dx);
                    edge_influence = AX_PS_EDGE + AX_PS_TOPI * edge_influence;
                }
                else /* (edge_influence <= 0) */
                { /* screen surface pixels that's too oblique with viewsight */
                    if (edge_influence > -AX_3D_POLYGON_EDGE_INFLUENCE_VEIL)
                        continue;
                    /* prevent backface from overcoming frontface at border */
                    polygon->IMG[j].z *=
                        (1 - edge_influence *
                         AX_3D_POLYGON_BACKFACE_RETARDATION);
                    light_influence = AX_SP_LIGHT * edge_influence *
                        AX_3D_LDOT(dx);
                    edge_influence = AX_SP_EDGE - AX_SP_TOPI * edge_influence;
                }
                if (light_influence < 0) light_influence = 0;
                polygon->IMG[j].carrier = AX_Colorcarrier
                    ( edge_influence * polygon->r + light_influence * lr,
                      edge_influence * polygon->g + light_influence * lg,
                      edge_influence * polygon->b + light_influence * lb );
                polygon->IMG[j].area = AX_AOP_TOP(iw,m).c.area;
                j++;
            }
        }
        polygon->IMG[0].o.offset = j;
    }
    return;
} /* end AX_3D_Polygons_Zpaper() */


/* Draw polygons to pixmap (block pixels only) according to z-buffer */
void AX_3D_Polygons_ZBprint(int iw, AX_3D_Polygons *P)
{
    int i;
    register int j;
    register AX_3D_Polygon *polygon;

    AX_C(

        for (i=0; i<P->n_shown; i++)
        {
            polygon = P->POLYGON + P->idx[i];
            for (j=1; j<polygon->img[0].offset; j++)
                if (polygon->img[j].z <=
                    AX_zmap[iw][polygon->img[j].offset])
                    AX_mem[iw].uc[polygon->img[j].offset] =
                        polygon->img[j].c.carrier;
        },

        for (i=0; i<P->n_shown; i++)
        {
            polygon = P->POLYGON + P->idx[i];
            for (j=1; j<polygon->img[0].offset; j++)
                if (polygon->img[j].z <=
                    AX_zmap[iw][polygon->img[j].offset])
                    AX_mem[iw].i2[polygon->img[j].offset] =
                        polygon->img[j].c.carrier;
        },

        for (i=0; i<P->n_shown; i++)
        {
            polygon = P->POLYGON + P->idx[i];
            for (j=1; j<polygon->img[0].offset; j++)
                if (polygon->img[j].z <=
                    AX_zmap[iw][polygon->img[j].offset])
                    AX_mem[iw].i4[polygon->img[j].offset] =
                        polygon->img[j].c.carrier;
        }

        );
    return;
} /* end AX_3D_Polygons_ZBprint() */


/* Draw polygons to pixmap (alias pixels only) according to z-buffer */
void AX_3D_Polygons_ZAprint(int iw, AX_3D_Polygons *P)
{
    int i;
    register int j;
    register AX_3D_Polygon *polygon;

    AX_C(

        for (i=P->n_shown-1; i>=0; i--)
        {
            polygon = P->POLYGON + P->idx[i];
            for (j=1; j<polygon->IMG[0].o.offset; j++)
                if (polygon->IMG[j].z <=
                    AX_zmap[iw][polygon->IMG[j].o.offset])
                    AX_mem[iw].uc[polygon->IMG[j].o.offset] = AX_mix
                        ( AX_mem[iw].uc[polygon->IMG[j].o.offset],
                          polygon->IMG[j].carrier, polygon->IMG[j].area );
            free(polygon->IMG);
        },

        for (i=P->n_shown-1; i>=0; i--)
        {
            polygon = P->POLYGON + P->idx[i];
            for (j=1; j<polygon->IMG[0].o.offset; j++)
                if (polygon->IMG[j].z <=
                    AX_zmap[iw][polygon->IMG[j].o.offset])
                    AX_mem[iw].i2[polygon->IMG[j].o.offset] = AX_mix
                        ( AX_mem[iw].i2[polygon->IMG[j].o.offset],
                          polygon->IMG[j].carrier, polygon->IMG[j].area );
            free(polygon->IMG);
        },

        for (i=P->n_shown-1; i>=0; i--)
        {
            polygon = P->POLYGON + P->idx[i];
            for (j=1; j<polygon->IMG[0].o.offset; j++)
                if (polygon->IMG[j].z <=
                    AX_zmap[iw][polygon->IMG[j].o.offset])
                    AX_mem[iw].i4[polygon->IMG[j].o.offset] = AX_mix
                        ( AX_mem[iw].i4[polygon->IMG[j].o.offset],
                          polygon->IMG[j].carrier, polygon->IMG[j].area );
            free(polygon->IMG);
        }

        );
    return;
} /* end AX_3D_Polygons_ZAprint() */


#define AX_3D_SAMEZ(z0,z1,dz) \
  (dz=(z1)-(z0), ABS(dz)<(0.2*AX_3D_POLYGON_BACKFACE_RETARDATION)*(z0))
#define AX_3D_BLOCKFORMED(area) ((area)>=0.99)

/* Use z-buffer morphed linklist to correctly handle */
/* multiple aliases from polygons WITHIN this group. */
void AX_3D_Polygons_Zaprint(int iw, AX_3D_Polygons *P)
{
    register AX_3D_Polygon *polygon;
    register int j,m;
    register struct AX_3D_PIXEL *p,*q=NULL;
    register AX_Float tmp,zmin;
    int i;
    /* initial screening of covered alias pixels and polygons */
    for (i=P->n_shown-1; i>=0; i--)
        if ((polygon=P->POLYGON+P->idx[i])->IMG)
        {
            for (m=0,j=1; j<polygon->IMG[0].o.offset; j++)
                if (polygon->IMG[j].z > AX_zmap[iw][polygon->IMG[j].o.offset])
                    polygon->IMG[j].o.offset = -1;
                else
                { /* codebook change but no loss of information */
                    polygon->IMG[j].o.offset = - 2 - polygon->IMG[j].o.offset;
                    m++;
                }
            if (m == 0)
            { /* early retirement of polygon */
                free((void *)polygon->IMG);
                polygon->IMG = NULL;
            }
        }
    /* link from zmap entry to any pixel that points to it: lose original z */
    for (i=P->n_shown-1; i>=0; i--)
        if ((polygon=P->POLYGON+P->idx[i])->IMG)
            for (j=1; j<polygon->IMG[0].o.offset; j++)
                if (polygon->IMG[j].o.offset < -1)
                    *((struct AX_3D_PIXEL **)
                      (&AX_zmap[iw][-polygon->IMG[j].o.offset-2])) =
                        &(polygon->IMG[j]);
    /* Since any needed zmap entry now points to a certain pixel, any other */
    /* pixel can append behind, or coalesce with members of the linklist if */
    /* they have the same z. The last member always keeps the zmap offset.  */
    for (i=P->n_shown-1; i>=0; i--)
        if ((polygon=P->POLYGON+P->idx[i])->IMG)
            for (j=1; j<polygon->IMG[0].o.offset; j++)
                if (polygon->IMG[j].o.offset < -1) /* code still sound */
                    if ( (p=*((struct AX_3D_PIXEL **)
                              (&AX_zmap[iw][-polygon->IMG[j].o.offset-2]))) !=
                         &(polygon->IMG[j]) )
                        while (1)
                        {
                            if ( AX_3D_SAMEZ (polygon->IMG[j].z, p->z, tmp) )
                            { /* merge the two alias pixels */
                                p->area += polygon->IMG[j].area;
                                tmp = polygon->IMG[j].area / p->area;
                                p->carrier = AX_mix
                                    (p->carrier, polygon->IMG[j].carrier, tmp);
                                p->z = (1-tmp) * p->z + tmp *
                                    polygon->IMG[j].z;
                                /* release the current pixel */
                                polygon->IMG[j].o.offset = -1;
                                break;
                            }
                            if ( p->o.offset >= 0 ) p = p->o.p;
                            else
                            { /* create a new leaf */
                                p->o.p = &polygon->IMG[j];
                                break;
                            }
                        }

    AX_C( 
    
        for (i=P->n_shown-1; i>=0; i--)
            if ((polygon=P->POLYGON+P->idx[i])->IMG)
                for (j=1; j<polygon->IMG[0].o.offset; j++)
                    if (polygon->IMG[j].o.offset < -1)  /* leaf node found */
                    { /* leaf node has zmap offset */
                        p = *((struct AX_3D_PIXEL **)
                              (&AX_zmap[iw][-polygon->IMG[j].o.offset-2]));
                        zmin = tmp = AX_INFINITY;
                        while(1)
                        { /* find the closest block pixel, if any */
                            if (AX_3D_BLOCKFORMED(p->area))
                            { /* both correctly and incorrectly merged */
                                if (p->z < tmp)
                                {
                                    tmp = p->z;
                                    AX_zmap[iw][-polygon->IMG[j].o.offset-2]
                                        = tmp;
                                    AX_mem[iw].uc[-polygon->IMG[j].o.offset-2]
                                        = p->carrier;
                                }
                            }
                            else if (p->z < zmin)
                            { /* find the closest "true" alias pixel, if any */
                                zmin = p->z;
                                q = p;
                            }
                            if (p->o.offset >= 0) p = p->o.p;
                            else break;
                        }
                        if (zmin < tmp)
                        { /* "true" alias pixel would appear */
                            AX_mem[iw].uc[-polygon->IMG[j].o.offset-2] = AX_MIX
                                ( AX_mem[iw].uc[-polygon->IMG[j].o.offset-2],
                                  q->carrier, q->area );
                            /* somewhat restore zmap */
                            if (tmp == AX_INFINITY)
                                AX_zmap[iw][-polygon->IMG[j].o.offset-2] =
                                    zmin;
                        }
                    },

        for (i=P->n_shown-1; i>=0; i--)
            if ((polygon=P->POLYGON+P->idx[i])->IMG)
                for (j=1; j<polygon->IMG[0].o.offset; j++)
                    if (polygon->IMG[j].o.offset < -1)  /* leaf node found */
                    { /* leaf node has zmap offset */
                        p = *((struct AX_3D_PIXEL **)
                              (&AX_zmap[iw][-polygon->IMG[j].o.offset-2]));
                        zmin = tmp = AX_INFINITY;
                        while(1)
                        { /* find the closest block pixel, if any */
                            if (AX_3D_BLOCKFORMED(p->area))
                            { /* both correctly and incorrectly merged */
                                if (p->z < tmp)
                                {
                                    tmp = p->z;
                                    AX_zmap[iw][-polygon->IMG[j].o.offset-2]
                                        = tmp;
                                    AX_mem[iw].i2[-polygon->IMG[j].o.offset-2]
                                        = p->carrier;
                                }
                            }
                            else if (p->z < zmin)
                            { /* find the closest "true" alias pixel, if any */
                                zmin = p->z;
                                q = p;
                            }
                            if (p->o.offset >= 0) p = p->o.p;
                            else break;
                        }
                        if (zmin < tmp)
                        { /* "true" alias pixel would appear */
                            AX_mem[iw].i2[-polygon->IMG[j].o.offset-2] = AX_MIX
                                ( AX_mem[iw].i2[-polygon->IMG[j].o.offset-2],
                                  q->carrier, q->area );
                            /* somewhat restore zmap */
                            if (tmp == AX_INFINITY)
                                AX_zmap[iw][-polygon->IMG[j].o.offset-2]
                                    = zmin;
                        }
                    },

        for (i=P->n_shown-1; i>=0; i--)
            if ((polygon=P->POLYGON+P->idx[i])->IMG)
                for (j=1; j<polygon->IMG[0].o.offset; j++)
                    if (polygon->IMG[j].o.offset < -1)  /* leaf node found */
                    { /* leaf node has zmap offset */
                        p = *((struct AX_3D_PIXEL **)
                              (&AX_zmap[iw][-polygon->IMG[j].o.offset-2]));
                        zmin = tmp = AX_INFINITY;
                        while(1)
                        { /* find the closest block pixel, if any */
                            if (AX_3D_BLOCKFORMED(p->area))
                            { /* both correctly and incorrectly merged */
                                if (p->z < tmp)
                                {
                                    tmp = p->z;
                                    AX_zmap[iw][-polygon->IMG[j].o.offset-2]
                                        = tmp;
                                    AX_mem[iw].i4[-polygon->IMG[j].o.offset-2]
                                        = p->carrier;
                                }
                            }
                            else if (p->z < zmin)
                            { /* find the closest "true" alias pixel, if any */
                                zmin = p->z;
                                q = p;
                            }
                            if (p->o.offset >= 0) p = p->o.p;
                            else break;
                        }
                        if (zmin < tmp)
                        { /* "true" alias pixel would appear */
                            AX_mem[iw].i4[-polygon->IMG[j].o.offset-2] = AX_MIX
                                ( AX_mem[iw].i4[-polygon->IMG[j].o.offset-2],
                                  q->carrier, q->area );
                            /* somewhat restore zmap */
                            if (tmp == AX_INFINITY)
                                AX_zmap[iw][-polygon->IMG[j].o.offset-2] =
                                    zmin;
                        }
                    }
        
        );
    for (i=P->n_shown-1; i>=0; i--)
        if ((polygon=P->POLYGON+P->idx[i])->IMG)
            free((void *)polygon->IMG);
    return;
} /* end AX_3D_Polygons_Zaprint() */


/* Assign parallelepiped (6 polygons) for H[][] originated at x0,y0,z0 */
void AX_3D_ParallelepipedAssign
(AX_Float x0, AX_Float y0, AX_Float z0, AX_Float H[3][3],
 AX_Float r, AX_Float g, AX_Float b, AX_3D_Polygons *P)
{
    register int i;
    register AX_Float tmp;
    AX_Float x[3], G[3][3];
    x[0] = x0;
    x[1] = y0;
    x[2] = z0;
    for (i=0; i<6; i++)
    {
        AX_3D_AssignRGB(P->POLYGON[i],r,g,b);
        P->POLYGON[i].nVertex = 4;
    }
    AX_M3invtranspose (H, G, tmp);
    AX_V3EQV (x, P->POLYGON[0].x);
    AX_V3add (x, H[2], &(P->POLYGON[0].x[3]));
    AX_V3add (&(P->POLYGON[0].x[3]), H[1], &(P->POLYGON[0].x[6]));
    AX_V3add (x, H[1], &(P->POLYGON[0].x[9]));
    AX_3D_PolygonTranslate (P->POLYGON[0], H[0], P->POLYGON[1], i);
    AX_V3normalize (G[0], P->POLYGON[1].front_normal, tmp);
    AX_V3NEG (P->POLYGON[1].front_normal, P->POLYGON[0].front_normal);
    AX_V3EQV (x, P->POLYGON[2].x);
    AX_V3add (x, H[0], &(P->POLYGON[2].x[3]));
    AX_V3add (&(P->POLYGON[2].x[3]), H[2], &(P->POLYGON[2].x[6]));
    AX_V3add (x, H[2], &(P->POLYGON[2].x[9]));
    AX_3D_PolygonTranslate (P->POLYGON[2], H[1], P->POLYGON[3], i);
    AX_V3normalize (G[1], P->POLYGON[3].front_normal, tmp);
    AX_V3NEG (P->POLYGON[3].front_normal, P->POLYGON[2].front_normal);
    AX_V3EQV (x, P->POLYGON[4].x);
    AX_V3add (x, H[1], &(P->POLYGON[4].x[3]));
    AX_V3add (&(P->POLYGON[4].x[3]), H[0], &(P->POLYGON[4].x[6]));
    AX_V3add (x, H[0], &(P->POLYGON[4].x[9]));
    AX_3D_PolygonTranslate (P->POLYGON[4], H[2], P->POLYGON[5], i);
    AX_V3normalize (G[2], P->POLYGON[5].front_normal, tmp);
    AX_V3NEG (P->POLYGON[5].front_normal, P->POLYGON[4].front_normal);
    return;
} /* end AX_3D_ParallelepipedAssign() */


#ifdef _3D_Parallelepiped_TEST
#define Z0           4.
#define boxwidth     2.
#define boxheight    2.
#define boxthickness 2.
int main (int argc, char *argv[])
{
    AX_3D_Polygons P[1] = {0};
    AX_Float H[3][3],R[3][3],A[3][3];
    AX_openwindow (getpid(), NULL, AX_DEFWIDTH, AX_DEFHEIGHT);
    AX_plugin_3D_module(0);
    AX_3D_Polygons_Realloc(P,6);
    TimeRandomize();
    while(1)
    {
        AX_namedmop(0,AX_WHITE);
        AX_clearzmap(0);
        AX_M3RANDOMROTATION(R,H[0][0]);
        AX_M3DIAGONAL(boxwidth,boxheight,boxthickness,A);
        AX_M3MUL(A,R,H);
        AX_3D_ParallelepipedAssign
            (-boxwidth/2, -boxheight/2, Z0, H, 0.5,0.5,0.5, P);
        AX_3D_Polygons_Zpaper (0, P);
        printf ("%d\n", P->n_shown);
        AX_3D_Polygons_ZBprint (0, P);
        AX_dump(0); AX_show(0);
        AXPressKey(0);
        AX_3D_Polygons_Zaprint (0, P);
        AX_dump(0); AX_show(0);
        AXPressKey(0);
    }
    AX_3D_Polygons_Free(P);
    Press_return();
    return (0);
}
#endif /* _3D_Parallelepiped_TEST */


#ifdef _3D_Polygons_TEST
#define Z0   10.
#define mesh 30
#define n_polygons (mesh+1)
#define cylheight 4.
#define cylradius 2.
#define platewidth 3.
#define PNG_FILENAME "/tmp/polygons.png"
#define JPG_FILENAME "/tmp/polygons.jpg"
#define EPS_FILENAME "/tmp/polygons.eps"
int main (int argc, char *argv[])
{
    register int i,j;
    AX_3D_Polygons P[1] = {0};
    AX_Float R[3][3], x0,y0,x1,y1, tmp[3];
    AX_3D_Polygons_Realloc(P,n_polygons);
    for (i=0; i<mesh; i++)
    {
        x0 = cylradius * cos(i*2*PI/mesh);
        y0 = cylradius * sin(i*2*PI/mesh);
        x1 = cylradius * cos((i+1)*2*PI/mesh);
        y1 = cylradius * sin((i+1)*2*PI/mesh);
        AX_3D_PolygonAssign (P->POLYGON+i, 4, 1,0,0,
                             x0,y0,-cylheight,
                             x1,y1,-cylheight,
                             x1,y1,cylheight,
                             x0,y0,cylheight);
    }
    AX_3D_PolygonAssign (P->POLYGON+i, 4, 0,0,1,
                         -platewidth,-platewidth,0.,
                         platewidth,-platewidth,0.,
                         platewidth,platewidth,0.,
                         -platewidth,platewidth,0.);
    AX_M3RANDOMROTATION(R,tmp[0]);
    for (i=0; i<n_polygons; i++)
    {
        for (j=0; j<P->POLYGON[i].nVertex; j++)
        {
            AX_M3MULV3(R,&(P->POLYGON[i].x[3*j]),tmp);
            P->POLYGON[i].x[3*j+2] += Z0;
        }
        AX_M3MULV3(R,P->POLYGON[i].front_normal,tmp);
    }
    AX_openwindow (getpid(), "Zaprint", AX_DEFWIDTH, AX_DEFHEIGHT);
    AX_PLUGIN_3D_module(0);
    AX_namedmop(0, AX_WHITE);
    AX_3D_Polygons_Zpaper (0, P);
    AX_3D_Polygons_ZBprint (0, P);
    AX_3D_Polygons_Zaprint (0, P);
    AX_dump(0); AX_show(0);
    AX_openwindow (getpid(), "ZAprint", AX_DEFWIDTH, AX_DEFHEIGHT);
    AX_namedmop(1, AX_WHITE);
    AX_PLUGIN_3D_module(1);
    AX_3D_Polygons_Zpaper (1, P);
    AX_3D_Polygons_ZBprint (1, P);
    AX_3D_Polygons_ZAprint (1, P);
    AX_dump(1); AX_show(1);
    Press_return();
    AX_save_pixmap_as_PNG(0, PNG_FILENAME);
    printf ("Image saved on \"%s\".\n", PNG_FILENAME);
    AX_save_pixmap_as_JPG(0, JPG_FILENAME);
    printf ("Image saved on \"%s\".\n", JPG_FILENAME);
    AX_save_pixmap_as_EPS(0, EPS_FILENAME);
    printf ("Image saved on \"%s\".\n", EPS_FILENAME);
    Press_return();
    AX_closewindow(0);
    Press_return();
    AX_closewindow(1);
    Press_return();
    return (0);
}
#endif /* _3D_Polygons_TEST */


/****************/
/* Ball Shading */
/****************/

static AX_Float AX_B_z (AX_Float rx, AX_Float ry, AX_Float radius)
{
    register AX_Float r2;
    rx /= radius;
    ry /= radius;
    r2 = rx*rx + ry*ry;
    return (sqrt(1-r2));
} /* end AX_B_z() */

/* highlight spot (1.2 x 50 degrees) and its normalization */
#define HX   0.68829172362126
#define HY  -0.98298245314679
#define HM2  4.84
/* color influence of light on ball and normalized z */
AX_Float AX_B_light_influence_and_z
(AX_Float rx, AX_Float ry, AX_Float radius, AX_Float *z)
{
    register AX_Float r2, light_influence;
    rx /= radius;
    ry /= radius;
    r2 = rx*rx + ry*ry;
    *z = sqrt(1-r2);
    rx -= HX;
    ry -= HY;
    light_influence = 1 - (rx*rx + ry*ry) / HM2;
    light_influence *= light_influence;
    light_influence *= light_influence * (1 - r2*r2);
    return (light_influence);
} /* end AX_B_light_influence_and_z() */
#undef HM2
#undef HY
#undef HX


/****************/
/* Ball Caching */
/****************/
static AX_BC AX_bc [AX_MAXWIN] [AX_BC_RMAX * AX_BC_RMESH];
/* used because the window width might be smaller than 2*AX_BC_RMAX+1 */
static int AX_BC_rmax [AX_MAXWIN];
/* offset shift AX_BC_rmax[iw] x AX_BC_rmax[iw] to decode */
/* ball scan by putting it in the AX_size[iw].width mesh. */
static int AX_BC_offshift [AX_MAXWIN];
/* this may sound ironic, but BC is going to provide */
/* scratch space for on-the-fly scan conversions.    */
static int AX_BC_bp_top [AX_MAXWIN];
static AX_3D_pixel *AX_BC_bp [AX_MAXWIN];

/* Install ball cache for this window */
int AX_plugin_BC_module (int iw)
{ /* for memory testing purpose when AX_BC_RMAX is changed */
    /* AX_Float c0_ratio = 2; */
    /* AX_Float b0_ratio = 2; */
    /* AX_Float a0_ratio = 2; */
    AX_Float c0_ratio = 0.98450;
    AX_Float b0_ratio = 0.98387;
    AX_Float a0_ratio = 1.26331;
    int c0, b0, a0, c, b, a, k, m, width, jw;
    AX_IJ i0=AX_BC_RMAX, j0=AX_BC_RMAX;
    AX_Float x0=i0+0.5, y0=j0+0.5;
    AX_AssertWinID ("AX_plugin_BC_module", iw);
    AX_Check_Module_Dependence ("AX_plugin_BC_module", iw, MODULE_BC);
    if ( AX_Module_Is_PluggedIn [iw] [MODULE_BC] )
        pe ("AX_plugin_BC_module:\n"
            "you are crazy, this module is already plugged in\n");
    c0 = 2*AX_BC_RMAX * AX_BC_RMAX*AX_BC_RMESH * c0_ratio;
    b0 = PI*SQUARE(AX_BC_RMAX) * AX_BC_RMAX*AX_BC_RMESH / 3 * b0_ratio;
    a0 = PI*AX_BC_RMAX * AX_BC_RMAX*AX_BC_RMESH * a0_ratio;
    /* try to copy other's ball cache */
    for ( jw = 0; jw < AX_MAXWIN; jw++ )
        if ( AX_Module_Is_PluggedIn [jw] [MODULE_BC] )
        {
            for (k=0; k<AX_BC_RMAX*AX_BC_RMESH; k++)
                AX_bc[iw][k] = AX_bc[jw][k];
            MALLOC ( AX_plugin_BC_module, AX_bc[iw][0].a, a0, AX_BC_Unit );
            MALLOC ( AX_plugin_BC_module, AX_bc[iw][0].b, b0, AX_BC_Unit );
            MALLOC ( AX_plugin_BC_module, AX_bc[iw][0].c, c0, int );
            memcpy ( AX_bc[iw][0].a, AX_bc[jw][0].a, a0 * sizeof(AX_BC_Unit) );
            memcpy ( AX_bc[iw][0].b, AX_bc[jw][0].b, b0 * sizeof(AX_BC_Unit) );
            memcpy ( AX_bc[iw][0].c, AX_bc[jw][0].c, c0 * sizeof(int) );
            for (k=1; k<AX_BC_RMAX*AX_BC_RMESH; k++)
            {
                AX_bc[iw][k].a = AX_bc[iw][k-1].a + AX_bc[iw][k-1].na;
                AX_bc[iw][k].b = AX_bc[iw][k-1].b + AX_bc[iw][k-1].nb;
                AX_bc[iw][k].c = AX_bc[iw][k-1].c + 2 * AX_bc[iw][k-1].nc;
            }
            AX_BC_rmax [iw] = AX_BC_rmax [jw];
            AX_BC_offshift [iw] = AX_BC_offshift [jw];
            AX_Module_Is_PluggedIn [iw] [MODULE_BC] = 1;
            AX_resize_BC_module (iw);
            goto exit;
        }
    /* live flat plate preconditioner: private */
    MALLOC ( AX_plugin_BC_module, AX_bc[iw][0].c, c0, int );
    /* live-offset/dead-z pixels: private */
    MALLOC ( AX_plugin_BC_module, AX_bc[iw][0].b, b0, AX_BC_Unit );
    /* dead i,j pixels: shared */
    MALLOC ( AX_plugin_BC_module, AX_bc[iw][0].B, b0, AX_PAIR );
    /* dead light_influence pixels: shared */
    MALLOC ( AX_plugin_BC_module, AX_bc[iw][0].light_influence, b0, AX_Float );
    /* live-offset/dead-z pixels: private */
    MALLOC ( AX_plugin_BC_module, AX_bc[iw][0].a, a0, AX_BC_Unit );
    /* dead i,j pixels: shared */
    MALLOC ( AX_plugin_BC_module, AX_bc[iw][0].A, a0, AX_PAIR );
    AX_BC_rmax [iw] = (AX_size[iw].width - 1) / 2;
    if (AX_BC_rmax [iw] >= AX_BC_RMAX) AX_BC_rmax [iw] = AX_BC_RMAX;
    AX_BC_offshift[iw] = (AX_size[iw].width+1) * AX_BC_rmax [iw];
    /* descending radius for cache coherency under z-order */
    for (c=b=a=k=0; k<AX_BC_RMAX*AX_BC_RMESH; k++)
    {
        AX_bc[iw][k].radius = (k+0.5) / AX_BC_RMESH;
        AX_CircleWindowScan (x0, y0, AX_bc[iw][k].radius,
                             AX_FREEZE_WIDTH, AX_FREEZE_WIDTH,
                             AX_AOP_top(iw), AX_BOP_top(iw));
        AX_bc[iw][k].c = AX_bc[iw][0].c + c;
        AX_bc[iw][k].b = AX_bc[iw][0].b + b;
        AX_bc[iw][k].B = AX_bc[iw][0].B + b;
        AX_bc[iw][k].light_influence = AX_bc[iw][0].light_influence + b;
        AX_bc[iw][k].a = AX_bc[iw][0].a + a;
        AX_bc[iw][k].A = AX_bc[iw][0].A + a;
        /* enhance the weak */
        if (AX_BOP_TOP(iw,0).offset == 1)
        {
            AX_BOP_TOP(iw,0).offset = 2;
            AX_BOP_TOP(iw,1).p.i = i0;
            AX_BOP_TOP(iw,1).p.j = j0;
            AX_AOP_TOP(iw,0).b.offset = 1;
        }
        width = (AX_bc[iw][k].radius <= AX_BC_rmax[iw]) ?
            AX_size[iw].width : AX_FREEZE_WIDTH;
        for (m=1; m<AX_AOP_TOP(iw,0).b.offset; m++)
        {
            AX_bc[iw][k].A[m-1].i = AX_AOP_TOP(iw,m).b.p.i - i0;
            AX_bc[iw][k].A[m-1].j = AX_AOP_TOP(iw,m).b.p.j - j0;
            AX_bc[iw][k].a[m-1].c1 = AX_AOP_TOP(iw,m).c.area;
            AX_bc[iw][k].a[m-1].offset = AX_bc[iw][k].A[m-1].i +
                AX_bc[iw][k].A[m-1].j * width;
        }
        AX_bc[iw][k].na = AX_AOP_TOP(iw,0).b.offset - 1;
        a += AX_bc[iw][k].na;
        if (a > a0) pe("AX_plugin_BC_module: a0 = %d not enough\n", a0);
        for (m=1; m<AX_BOP_TOP(iw,0).offset; m++)
        {
            AX_bc[iw][k].B[m-1].i = AX_BOP_TOP(iw,m).p.i - i0;
            AX_bc[iw][k].B[m-1].j = AX_BOP_TOP(iw,m).p.j - j0;
            AX_bc[iw][k].light_influence[m-1] = AX_B_light_influence_and_z
                ( AX_bc[iw][k].B[m-1].i, AX_bc[iw][k].B[m-1].j,
                  AX_bc[iw][k].radius, &AX_bc[iw][k].b[m-1].c1 );
            AX_bc[iw][k].b[m-1].offset = AX_bc[iw][k].B[m-1].i +
                AX_bc[iw][k].B[m-1].j * width;
        }
        AX_bc[iw][k].nb = AX_BOP_TOP(iw,0).offset - 1;
        b += AX_bc[iw][k].nb;
        if (b > b0) pe("AX_plugin_BC_module: b0 = %d not enough\n", b0);
        if (AX_bc[iw][k].nb == 0)
        { /* give the weak a mask */
            jw = 0;
            AX_bc[iw][k].nc = 1;
            AX_bc[iw][k].c[0] = 0;
            AX_bc[iw][k].c[1] = 1;
        }
        else
        {
            jw = AX_BOP_TOP(iw,1).p.j - j0;
            AX_bc[iw][k].nc = AX_BOP_TOP(iw,AX_bc[iw][k].nb).p.j
                - AX_BOP_TOP(iw,1).p.j + 1;
            for (m=0; m<AX_bc[iw][k].nc; m++)
            {
                AX_bc[iw][k].c[2*m]   =  HUGE_INT;
                AX_bc[iw][k].c[2*m+1] = -HUGE_INT;
            }
            for (m=1; m<AX_BOP_TOP(iw,0).offset; m++)
            {
                if (AX_BOP_TOP(iw,m).p.i-i0 < AX_bc[iw][k].
                    c[2*(AX_BOP_TOP(iw,m).p.j-j0-jw)])
                    AX_bc[iw][k].
                        c[2*(AX_BOP_TOP(iw,m).p.j-j0-jw)]
                        = AX_BOP_TOP(iw,m).p.i-i0;
                if (AX_BOP_TOP(iw,m).p.i-i0 >= AX_bc[iw][k].
                    c[2*(AX_BOP_TOP(iw,m).p.j-j0-jw)+1])
                    AX_bc[iw][k].
                        c[2*(AX_BOP_TOP(iw,m).p.j-j0-jw)+1] =
                        AX_BOP_TOP(iw,m).p.i-i0 + 1;
            }
        }
        for (m=0; m<AX_bc[iw][k].nc; m++)
        {
            AX_bc[iw][k].c[2*m]   += (m+jw) * width;
            AX_bc[iw][k].c[2*m+1] += (m+jw) * width;
        }
        c += 2 * AX_bc[iw][k].nc;
        if (c > c0) pe("AX_plugin_BC_module: c0 = %d not enough\n", c0);
    }
  exit:
    /* printf ("better ratios c=%f b=%f a=%f -> size=%d bytes\n", */
    /* c0_ratio*c/c0, b0_ratio*b/b0, a0_ratio*a/a0, c0*sizeof(int) + */
    /* b0*(sizeof(AX_Float)+sizeof(AX_BC_Unit)+sizeof(AX_PAIR)) + */
    /* a0*(sizeof(AX_BC_Unit)+sizeof(AX_PAIR)) + */
    /* INT(AX_BC_BP_RATIO * AX_size[iw].width * */
    /* AX_size[iw].height * sizeof(AX_3D_pixel)) ); */
    MALLOC ( AX_plugin_BC_module, AX_BC_bp[iw],
             INT(AX_BC_BP_RATIO *AX_size[iw].width * AX_size[iw].height),
             AX_3D_pixel );
    AX_BC_bp_top [iw] = 0;
    AX_Module_Is_PluggedIn [iw] [MODULE_BC] = 1;
    return ( c0*sizeof(int) +
             b0*(sizeof(AX_Float)+sizeof(AX_BC_Unit)+sizeof(AX_PAIR)) +
             a0*(sizeof(AX_BC_Unit)+sizeof(AX_PAIR)) +
             INT(AX_BC_BP_RATIO * AX_size[iw].width *
                 AX_size[iw].height * sizeof(AX_3D_pixel) ) );
} /* end AX_plugin_BC_module() */


/* For the kth cached ball, given that the original offset is */
/* with respect to "from_width", reassign it to "to_width".   */
static void AX_BC_ReOffset (int iw, int k, int from_width, int to_width)
{
    int i,j,m,shift,rmax;
    if (from_width == to_width) return;
    if (2 * AX_bc[iw][k].radius + 1 > from_width)
        pe ("AX_BC_ReOffset: impossible from_width, sir\n");
    if (2 * AX_bc[iw][k].radius + 1 > to_width)
        pe ("AX_BC_ReOffset: impossible to_width, sir\n");
    rmax = AX_bc[iw][k].radius + 0.5;
    shift = rmax * from_width + rmax;
    for (m=0; m<AX_bc[iw][k].na; m++)
    {
        AX_bc[iw][k].a[m].offset += shift;
        i = AX_bc[iw][k].a[m].offset % from_width - rmax;
        j = AX_bc[iw][k].a[m].offset / from_width - rmax;
        AX_bc[iw][k].a[m].offset = i + j * to_width;
    }
    for (m=0; m<AX_bc[iw][k].nb; m++)
    {
        AX_bc[iw][k].b[m].offset += shift;
        i = AX_bc[iw][k].b[m].offset % from_width - rmax;
        j = AX_bc[iw][k].b[m].offset / from_width - rmax;
        AX_bc[iw][k].b[m].offset = i + j * to_width;
    }
    for (m=0; m<2*AX_bc[iw][k].nc; m++)
    {
        AX_bc[iw][k].c[m] += shift;
        i = AX_bc[iw][k].c[m] % from_width - rmax;
        j = AX_bc[iw][k].c[m] / from_width - rmax;
        AX_bc[iw][k].c[m] = i + j * to_width;
    }
    return;
} /* end AX_BC_ReOffset() */


/*
  Resize AX_bc[iw][]: AX_BC_offshift[iw] (used on balls
  with radius <= AX_BC_rmax[iw]) had always been
  (AX_size[iw].width + 1) * AX_BC_rmax[iw] whether
  AX_BC_rmax[iw] is AX_BC_RMAX or not, therefore we
  can deduce old AX_size[iw].width from AX_BC_offshift[iw]
  and AX_BC_rmax[iw]. Also, we guarantee that
  2 * AX_BC_rmax[iw] + 1 <= AX_size[iw].width. If
  AX_BC_RMAX satisfy that, fine; else let
  AX_BC_rmax[iw] = [ (AX_size[iw].width-1)/2 ].
  For ball cache with radius > AX_BC_rmax[iw],
  we will preserve its information using a constant
  width scale = 2*AX_BC_RMAX+1 but will not use it
  in this session. In next sessions we may thaw
  some ball cache and re-arm with new AX_size[iw].width.
*/

/* Update ball cache information according to new AX_size[iw] */
int AX_resize_BC_module (int iw)
{
    int k, from_width, new_rmax;
    AX_AssertWinID ("AX_resize_BC_module", iw);
    if ( ! AX_Module_Is_PluggedIn [iw] [MODULE_BC] )
        pe ("AX_resize_BC_module:\n"
            "you are crazy, this module is not plugged in\n");
    from_width = AX_BC_offshift[iw] / AX_BC_rmax[iw] - 1;
    if (from_width == AX_size[iw].width) return(1);
    new_rmax = (AX_size[iw].width - 1) / 2;
    if (new_rmax >= AX_BC_RMAX) new_rmax = AX_BC_RMAX;
    AX_BC_offshift[iw] = (AX_size[iw].width + 1) * new_rmax;
    for (k=0; k<AX_BC_RMAX*AX_BC_RMESH; k++)
    {
        if ( AX_bc[iw][k].radius <= AX_BC_rmax[iw] )
        {
            if ( AX_bc[iw][k].radius <= new_rmax )
                AX_BC_ReOffset (iw, k, from_width, AX_size[iw].width);
            else
                AX_BC_ReOffset (iw, k, from_width, AX_FREEZE_WIDTH);
        }
        else
        {
            if ( AX_bc[iw][k].radius <= new_rmax )
                AX_BC_ReOffset (iw, k, AX_FREEZE_WIDTH, AX_size[iw].width);
        }
    }
    AX_BC_rmax[iw] = new_rmax;
    REALLOC(AX_resize_BC_module, AX_BC_bp[iw],
            INT(AX_BC_BP_RATIO * AX_size[iw].width * AX_size[iw].height),
            AX_3D_pixel);
    AX_BC_bp_top [iw] = 0;
    return(1);
} /* end AX_resize_BC_module() */


/* Free ball cache allocations */
int AX_plugout_BC_module (int iw)
{
    int jw;
    AX_AssertWinID ("AX_plugout_BC_module", iw);
    if ( ! AX_Module_Is_PluggedIn [iw] [MODULE_BC] )
        pe ("AX_plugout_BC_module:\n"
            "you are crazy, this module is not plugged in\n");
    Free(AX_BC_bp[iw]);
    Free(AX_bc[iw][0].c);
    Free(AX_bc[iw][0].b);
    Free(AX_bc[iw][0].a);
    AX_Module_Is_PluggedIn [iw] [MODULE_BC] = 0;
    for ( jw = 0; jw < AX_MAXWIN; jw++ )
        if ( AX_Module_Is_PluggedIn [jw] [MODULE_BC] ) return(1);
    Free(AX_bc[iw][0].B);
    Free(AX_bc[iw][0].A);
    Free(AX_bc[iw][0].light_influence);
    return(1);
} /* end AX_plugout_BC_module() */


/****************/
/* Ball Support */
/****************/


/* Realloc n_balls */
void AX_3D_Balls_Realloc (AX_3D_Balls *B, int n_balls)
{
    B->n_balls = n_balls;
    REALLOC (AX_3D_Balls_Realloc, B->BALL, n_balls, AX_3D_Ball);
    REALLOC (AX_3D_Balls_Realloc, B->idx,  n_balls, int);
    REALLOC (AX_3D_Balls_Realloc, B->status,  n_balls, unsigned char);
    return;
} /* end AX_3D_Balls_Realloc() */


/* Realloc n_balls followed by x0,y0,z0,radius,r,g,b,... interface */
void AX_3D_Balls_Recreate (AX_3D_Balls *B, int n_balls, ...)
{
    register int i;
    va_list ap;
    AX_3D_Balls_Realloc (B, n_balls);
    va_start(ap, n_balls);
    for (i=0; i<n_balls; i++)
    {
        B->BALL[i].x[0] = va_arg(ap, double);
        B->BALL[i].x[1] = va_arg(ap, double);
        B->BALL[i].x[2] = va_arg(ap, double);
        B->BALL[i].radius = va_arg(ap, double);
        AX_3D_AssignRGB (B->BALL[i], va_arg(ap, double), va_arg(ap, double),
                         va_arg(ap, double));
    }
    va_end(ap);
    return;
} /* end AX_3D_Balls_Recreate() */


/* free allocated memory and set NULL */
void AX_3D_Balls_Free (AX_3D_Balls *B)
{
    Free (B->status);
    Free (B->idx);
    Free (B->BALL);
    return;
} /* end AX_3D_Balls_Free() */


#define QUICKSORT_M  7
static void AX_sort_3D_Balls (AX_3D_Balls *B)
{
    register AX_3D_Ball *ball;
    register int i, j, tmp, *idx;
    register AX_Float a;
    int idxt, ir=B->n_shown, k, l=1, jstack=0;
    int istack[128];
    idx = B->idx - 1;
    ball = B->BALL;
    for (;;)
    {
	if (ir-l < QUICKSORT_M)
	{
	    for (j=l+1; j<=ir; j++)
	    {
		idxt = idx[j];
		a = ball[idxt].vx[2];
		for (i=j-1; i>=l; i--)
		{
		    if (ball[idx[i]].vx[2] <= a) break;
		    idx[i+1]=idx[i];
		}
		idx[i+1]=idxt;
	    }
	    if (jstack == 0) break;
	    ir = istack[jstack--];
	    l = istack[jstack--];
	}
	else
	{
	    k = (l+ir) >> 1;
	    SWAP(idx[k],idx[l+1],tmp);
	    if (ball[idx[l]].vx[2] > ball[idx[ir]].vx[2])
                SWAP(idx[l],idx[ir],tmp);
	    if (ball[idx[l+1]].vx[2] > ball[idx[ir]].vx[2])
                SWAP(idx[l+1],idx[ir],tmp);
	    if (ball[idx[l]].vx[2] > ball[idx[l+1]].vx[2])
                SWAP(idx[l],idx[l+1],tmp);
	    i = l+1;
	    j = ir;
	    idxt = idx[l+1];
	    a = ball[idxt].vx[2];
	    for (;;)
	    {
		do i++; while (ball[idx[i]].vx[2] < a);
		do j--; while (ball[idx[j]].vx[2] > a);
		if (j < i) break;
		SWAP (idx[i],idx[j],tmp);
	    }
	    idx[l+1]=idx[j];
	    idx[j]=idxt;
	    jstack += 2;
	    if (ir-i+1 >= j-l)
	    {
		istack[jstack] = ir;
		istack[jstack-1] = i;
		ir = j-1;
	    }
	    else
	    {
		istack[jstack] = j-1;
		istack[jstack-1] = l;
		l = i;
	    }
	}
    }
    return;
} /* end AX_sort_3D_Balls() */
#undef QUICKSORT_M


/* Calculate coordinates in viewport frame and determine visibility */
void AX_3D_Balls_Zdet (int iw, AX_3D_Balls *B)
{
    register AX_3D_Ball *ball=NULL;
    register int i, m;
    register AX_Float tmp=0;
    int j, fp_activated, activated[AX_3D_MAX_FILTER_PLANE];
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
                        if ( AX_V3DOT(AX_3D[iw].fp[j].dx, ball->x) >
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
        if (i >=0 ) /* viewpoint ejection */
        {
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
                if (m >= 1000)
                    pe ("AX_3D_Balls_Zdet: viewpoint ejection failed\n");
            }
        }
        else break;
    }
    AX_sort_3D_Balls(B);
    return;
} /* end AX_3D_Balls_Zdet() */


/* register the balls in z-buffer */
void AX_3D_Balls_Zpaper (int iw, AX_3D_Balls *B)
{
    register AX_3D_Ball *ball;
    register int m,k,w,v;
    register AX_Float *base, *BASE, tmp;
    register AX_BC_Unit *b;
    register AX_PAIR *bb;
    register AX_3D_pixel *bp;
    int i, n_shown;
    AX_Float edge_influence, lr,lg,lb, light_influence, vradius;
    B->bp_occupied = AX_BC_bp_top[iw];
    n_shown = B->n_shown;
    /* the balls are presumably sorted in vx[2] */
    for ( B->n_shown=i=0; i<n_shown; i++ )
    {
        ball = B->BALL + B->idx[i];
        if ( B->status[B->idx[i]] & AX_3D_BALL_CACHED )  
        { /* use cached scans */
            base = AX_zmap[iw] + AX_OFF(iw, ball->i0, ball->j0);
            if ( B->status[B->idx[i]] & AX_3D_BALL_COMPLETE )
            {
                tmp = ball->vx[2] - ball->radius;
                for (m=ball->vs.bc->nc; m--;)
                    for (BASE=base+ball->vs.bc->c[2*m];
                         BASE<base+ball->vs.bc->c[2*m+1]; BASE++)
                        if (tmp < *BASE) goto next1;
                continue;
              next1:
                b = ball->vs.bc->b;
                for (m=ball->vs.bc->nb; m--;)
                {
                    tmp = ball->vx[2] - b[m].c1 * ball->radius;
                    if (tmp < base[b[m].offset]) base[b[m].offset] = tmp;
                }
            }
            else
            {  /* AX_3D_BALL_CLIPPED */
                k = 0;
                b = ball->vs.bc->b;
                bb = ball->vs.bc->B;
                for (m=ball->vs.bc->nb; m--;)
                {
                    w = ball->i0 + bb[m].i;
                    v = ball->j0 + bb[m].j;
                    if (AX_IJINWIN(iw,w,v))
                    { /* pixel-based clipping */
                        BASE = base + b[m].offset;
                        tmp = ball->vx[2] - b[m].c1 * ball->radius;
                        if ( tmp < *BASE )
                        {
                            *BASE = tmp;
                            k++;
                        }
                    }
                }
                if (k==0) continue;
            }
        }
        else  /* scan it on the fly */
        {
            vradius = AX_3D[iw].k / ball->vx[2] * ball->radius;
            AX_CircleWindowScan
                ( ball->vx[0], ball->vx[1], vradius,
                  AX_size[iw].width, AX_size[iw].height,
                  AX_AOP_top(iw), AX_BOP_top(iw) );
            ball->vs.bp = bp = AX_BC_bp [iw] + AX_BC_bp_top [iw];
            for (k=m=1; m<AX_AOP_TOP(iw,0).b.offset; m++)
            {
                bp[k].offset = AX_OFF
                    (iw, AX_AOP_TOP(iw,m).b.p.i, AX_AOP_TOP(iw,m).b.p.j);
                if ( ball->vx[2] < AX_zmap [iw] [bp[k].offset] )
                {
                    bp[k].c1 = AX_AOP_TOP(iw,m).c.area;
                    k++;
                }
            }
            bp[0].offset = k;
            AX_3D_LRGB(AX_BS_TOP*ball->r, AX_BS_TOP*ball->g,
                       AX_BS_TOP*ball->b, lr,lg,lb);
            tmp = ball->vx[2] - ball->radius;
            for (k++,m=1; m<AX_BOP_TOP(iw,0).offset; m++)
            {
                bp[k].offset = AX_OFF
                    (iw, AX_BOP_TOP(iw,m).p.i, AX_BOP_TOP(iw,m).p.j);
                if ( tmp > AX_zmap [iw] [bp[k].offset] ) continue;
                light_influence = AX_B_light_influence_and_z
                    ( AX_BOP_TOP(iw,m).p.i + 0.5 - ball->vx[0],
                      AX_BOP_TOP(iw,m).p.j + 0.5 - ball->vx[1],
                      vradius, &edge_influence );
                bp[k].c1 = ball->vx[2] - edge_influence * ball->radius;
                if ( bp[k].c1 < AX_zmap [iw] [bp[k].offset] )
                {
                    edge_influence *= AX_BS_TOPI;
                    edge_influence += AX_BS_EDGE;
                    AX_zmap [iw] [bp[k].offset] = bp[k].c1;
                    bp[k].carrier = AX_Colorcarrier
                        (edge_influence*ball->r + light_influence*lr,
                         edge_influence*ball->g + light_influence*lg,
                         edge_influence*ball->b + light_influence*lb);
                    k++;
                }
            }
            if (k > 2)
            {
                bp[bp[0].offset].offset = k - bp[0].offset;
                AX_BC_bp_top [iw] += k;
            }
            else continue;
        }
        B->idx[B->n_shown++] = B->idx[i];
    }
    B->bp_occupied = AX_BC_bp_top [iw] - B->bp_occupied;
    return;
} /* end AX_3D_Balls_Zpaper() */


/* Draw the balls to pixmap (block pixels only) according to z-buffer */
void AX_3D_Balls_ZBprint (int iw, AX_3D_Balls *B)
{
    register AX_3D_Ball *ball;
    register int m,w,v;
    register AX_Float *base, edge_influence,lr,lg,lb, light_influence;
    register AX_Shm p;
    register AX_BC_Unit *b;
    register AX_3D_pixel *bp;
    register AX_PAIR *bb;
    int i;

    AX_C(

        for (i=B->n_shown; i--;)
            if (B->status[B->idx[i]] & AX_3D_BALL_CACHED)
            {  /* use cached ball scans */
                ball = B->BALL + B->idx[i];
                AX_3D_LRGB(AX_BS_TOP*ball->r, AX_BS_TOP*ball->g,
                           AX_BS_TOP*ball->b, lr,lg,lb);
                m = AX_OFF(iw, ball->i0, ball->j0);
                base = AX_zmap[iw] + m;
                p.uc = AX_mem[iw].uc + m;
                b = ball->vs.bc->b;
                if (B->status[B->idx[i]] & AX_3D_BALL_COMPLETE)
                {
                    for (m=ball->vs.bc->nb; m--; )
                        if ( ball->vx[2] - b[m].c1 * ball->radius <
                             base[b[m].offset] * (1+AX_TINY) )
                        {
                            edge_influence = AX_BS_EDGE + AX_BS_TOPI * b[m].c1;
                            light_influence = ball->vs.bc->light_influence[m];
                            p.uc[b[m].offset] = AX_Colorcarrier
                                ( edge_influence * ball->r +
                                  light_influence * lr,
                                  edge_influence * ball->g +
                                  light_influence * lg,
                                  edge_influence * ball->b +
                                  light_influence * lb);
                        }
                }
                else  /* AX_3D_BALL_CLIPPED */
                {
                    bb = ball->vs.bc->B;
                    for (m=ball->vs.bc->nb; m--;)
                    {
                        w = ball->i0 + bb[m].i;
                        v = ball->j0 + bb[m].j;
                        if (AX_IJINWIN(iw,w,v) &&
                            ( ball->vx[2] - b[m].c1 * ball->radius <
                              base[b[m].offset] * (1+AX_TINY) ) )
                        {
                            edge_influence = AX_BS_EDGE + AX_BS_TOPI * b[m].c1;
                            light_influence = ball->vs.bc->light_influence[m];
                            p.uc[b[m].offset] = AX_Colorcarrier
                                ( edge_influence * ball->r +
                                  light_influence * lr,
                                  edge_influence * ball->g +
                                  light_influence * lg,
                                  edge_influence * ball->b +
                                  light_influence * lb);
                        }
                    }
                }
            }
            else  /* in ready-to-deploy pixels */
            {
                ball = B->BALL + B->idx[i];
                bp = ball->vs.bp + ball->vs.bp[0].offset + 1;
                for (m=bp[-1].offset-1; m--;)
                    if (bp[m].c1 <= AX_zmap [iw] [bp[m].offset])
                        AX_mem[iw].uc[bp[m].offset] = bp[m].carrier;
            },

        for (i=B->n_shown; i--;)
            if (B->status[B->idx[i]] & AX_3D_BALL_CACHED)
            {  /* use cached ball scans */
                ball = B->BALL + B->idx[i];
                AX_3D_LRGB(AX_BS_TOP*ball->r, AX_BS_TOP*ball->g,
                           AX_BS_TOP*ball->b, lr,lg,lb);
                m = AX_OFF(iw, ball->i0, ball->j0);
                base = AX_zmap[iw] + m;
                p.i2 = AX_mem[iw].i2 + m;
                b = ball->vs.bc->b;
                if (B->status[B->idx[i]] & AX_3D_BALL_COMPLETE)
                {
                    for (m=ball->vs.bc->nb; m--; )
                        if ( ball->vx[2] - b[m].c1 * ball->radius <
                             base[b[m].offset] * (1+AX_TINY) )
                        {
                            edge_influence = AX_BS_EDGE + AX_BS_TOPI * b[m].c1;
                            light_influence = ball->vs.bc->light_influence[m];
                            p.i2[b[m].offset] = AX_Colorcarrier
                                ( edge_influence * ball->r +
                                  light_influence * lr,
                                  edge_influence * ball->g +
                                  light_influence * lg,
                                  edge_influence * ball->b +
                                  light_influence * lb);
                        }
                }
                else  /* AX_3D_BALL_CLIPPED */
                {
                    bb = ball->vs.bc->B;
                    for (m=ball->vs.bc->nb; m--;)
                    {
                        w = ball->i0 + bb[m].i;
                        v = ball->j0 + bb[m].j;
                        if (AX_IJINWIN(iw,w,v) &&
                            ( ball->vx[2] - b[m].c1 * ball->radius <
                              base[b[m].offset] * (1+AX_TINY) ) )
                        {
                            edge_influence = AX_BS_EDGE + AX_BS_TOPI * b[m].c1;
                            light_influence = ball->vs.bc->light_influence[m];
                            p.i2[b[m].offset] = AX_Colorcarrier
                                ( edge_influence * ball->r +
                                  light_influence * lr,
                                  edge_influence * ball->g +
                                  light_influence * lg,
                                  edge_influence * ball->b +
                                  light_influence * lb);
                        }
                    }
                }
            }
            else  /* in ready-to-deploy pixels */
            {
                ball = B->BALL + B->idx[i];
                bp = ball->vs.bp + ball->vs.bp[0].offset + 1;
                for (m=bp[-1].offset-1; m--;)
                    if (bp[m].c1 <= AX_zmap [iw] [bp[m].offset])
                        AX_mem[iw].i2[bp[m].offset] = bp[m].carrier;
            },

        for (i=B->n_shown; i--;)
            if (B->status[B->idx[i]] & AX_3D_BALL_CACHED)
            {  /* use cached ball scans */
                ball = B->BALL + B->idx[i];
                AX_3D_LRGB(AX_BS_TOP*ball->r, AX_BS_TOP*ball->g,
                           AX_BS_TOP*ball->b, lr,lg,lb);
                m = AX_OFF(iw, ball->i0, ball->j0);
                base = AX_zmap[iw] + m;
                p.i4 = AX_mem[iw].i4 + m;
                b = ball->vs.bc->b;
                if (B->status[B->idx[i]] & AX_3D_BALL_COMPLETE)
                {
                    for (m=ball->vs.bc->nb; m--; )
                        if ( ball->vx[2] - b[m].c1 * ball->radius <
                             base[b[m].offset] * (1+AX_TINY) )
                        {
                            edge_influence = AX_BS_EDGE + AX_BS_TOPI * b[m].c1;
                            light_influence = ball->vs.bc->light_influence[m];
                            p.i4[b[m].offset] = AX_Colorcarrier
                                ( edge_influence * ball->r +
                                  light_influence * lr,
                                  edge_influence * ball->g +
                                  light_influence * lg,
                                  edge_influence * ball->b +
                                  light_influence * lb);
                        }
                }
                else  /* AX_3D_BALL_CLIPPED */
                {
                    bb = ball->vs.bc->B;
                    for (m=ball->vs.bc->nb; m--;)
                    {
                        w = ball->i0 + bb[m].i;
                        v = ball->j0 + bb[m].j;
                        if (AX_IJINWIN(iw,w,v) &&
                            ( ball->vx[2] - b[m].c1 * ball->radius <
                              base[b[m].offset] * (1+AX_TINY) ) )
                        {
                            edge_influence = AX_BS_EDGE + AX_BS_TOPI * b[m].c1;
                            light_influence = ball->vs.bc->light_influence[m];
                            p.i4[b[m].offset] = AX_Colorcarrier
                                ( edge_influence * ball->r +
                                  light_influence * lr,
                                  edge_influence * ball->g +
                                  light_influence * lg,
                                  edge_influence * ball->b +
                                  light_influence * lb);
                        }
                    }
                }
            }
            else  /* in ready-to-deploy pixels */
            {
                ball = B->BALL + B->idx[i];
                bp = ball->vs.bp + ball->vs.bp[0].offset + 1;
                for (m=bp[-1].offset-1; m--;)
                    if (bp[m].c1 <= AX_zmap [iw] [bp[m].offset])
                        AX_mem[iw].i4[bp[m].offset] = bp[m].carrier;
            }
                
        );
    return;
} /* end AX_3D_Balls_ZBprint() */


/* draw the balls to pixmap (alias pixels only) according to z-buffer */
void AX_3D_Balls_ZAprint (int iw, AX_3D_Balls *B)
{
    register AX_3D_Ball *ball;
    register int m,w,v;
    register AX_PAIR *aa;
    register AX_Shm p;
    register AX_Float *base, tmp;
    register AX_BC_Unit *a;
    register AX_3D_pixel *bp;
    register AX_Carrier edge_carrier;
    int i;
    AX_C(
        
        for (i=B->n_shown; i--;)
        {
            ball = B->BALL + B->idx[i];
            edge_carrier = AX_Colorcarrier
                (AX_BS_EDGE*ball->r, AX_BS_EDGE*ball->g, AX_BS_EDGE*ball->b);
            tmp = ball->vx[2];
            if (B->status[B->idx[i]] & AX_3D_BALL_CACHED)
            {  /* use cached ball scans */
                m = AX_OFF(iw, ball->i0, ball->j0);
                base = AX_zmap[iw] + m;
                p.uc = AX_mem[iw].uc + m;
                a = ball->vs.bc->a;
                if (B->status[B->idx[i]] & AX_3D_BALL_COMPLETE)
                {
                    for (m=ball->vs.bc->na; m--;)
                        if (tmp < base[a[m].offset])
                            p.uc[a[m].offset] = AX_mix
                                ( p.uc[a[m].offset], edge_carrier, a[m].c1 );
                }
                else
                { /* AX_3D_BALL_CLIPPED */
                    aa = ball->vs.bc->A;
                    for (m=ball->vs.bc->na; m--;)
                    {
                        w = ball->i0 + aa[m].i;
                        v = ball->j0 + aa[m].j;
                        if (AX_IJINWIN(iw,w,v) &&
                            (tmp < base[a[m].offset]))
                            p.uc[a[m].offset] = AX_mix
                                ( p.uc[a[m].offset], edge_carrier, a[m].c1 );
                    }
                }
            }
            else  /* in ready-to-deploy pixels */
            {
                bp = ball->vs.bp + 1;
                for (m=bp[-1].offset-1; m--;)
                    if (tmp <= AX_zmap[iw][bp[m].offset])
                        AX_mem[iw].uc[bp[m].offset] = AX_mix
                            (AX_mem[iw].uc[bp[m].offset], edge_carrier,
                             bp[m].c1);
            }
        },
                
        for (i=B->n_shown; i--;)
        {
            ball = B->BALL + B->idx[i];
            edge_carrier = AX_Colorcarrier
                (AX_BS_EDGE*ball->r, AX_BS_EDGE*ball->g, AX_BS_EDGE*ball->b);
            tmp = ball->vx[2];
            if (B->status[B->idx[i]] & AX_3D_BALL_CACHED)
            {  /* use cached ball scans */
                m = AX_OFF(iw, ball->i0, ball->j0);
                base = AX_zmap[iw] + m;
                p.i2 = AX_mem[iw].i2 + m;
                a = ball->vs.bc->a;
                if (B->status[B->idx[i]] & AX_3D_BALL_COMPLETE)
                {
                    for (m=ball->vs.bc->na; m--;)
                        if (tmp < base[a[m].offset])
                            p.i2[a[m].offset] = AX_mix
                                ( p.i2[a[m].offset], edge_carrier, a[m].c1 );
                }
                else
                { /* AX_3D_BALL_CLIPPED */
                    aa = ball->vs.bc->A;
                    for (m=ball->vs.bc->na; m--;)
                    {
                        w = ball->i0 + aa[m].i;
                        v = ball->j0 + aa[m].j;
                        if (AX_IJINWIN(iw,w,v) &&
                            (tmp < base[a[m].offset]))
                            p.i2[a[m].offset] = AX_mix
                                ( p.i2[a[m].offset], edge_carrier, a[m].c1 );
                    }
                }
            }
            else  /* in ready-to-deploy pixels */
            {
                bp = ball->vs.bp + 1;
                for (m=bp[-1].offset-1; m--;)
                    if (tmp <= AX_zmap[iw][bp[m].offset])
                        AX_mem[iw].i2[bp[m].offset] = AX_mix
                            (AX_mem[iw].i2[bp[m].offset], edge_carrier,
                             bp[m].c1);
            }
        },

        for (i=B->n_shown; i--;)
        {
            ball = B->BALL + B->idx[i];
            edge_carrier = AX_Colorcarrier
                (AX_BS_EDGE*ball->r, AX_BS_EDGE*ball->g, AX_BS_EDGE*ball->b);
            tmp = ball->vx[2];
            if (B->status[B->idx[i]] & AX_3D_BALL_CACHED)
            {  /* use cached ball scans */
                m = AX_OFF(iw, ball->i0, ball->j0);
                base = AX_zmap[iw] + m;
                p.i4 = AX_mem[iw].i4 + m;
                a = ball->vs.bc->a;
                if (B->status[B->idx[i]] & AX_3D_BALL_COMPLETE)
                {
                    for (m=ball->vs.bc->na; m--;)
                        if (tmp < base[a[m].offset])
                            p.i4[a[m].offset] = AX_mix
                                ( p.i4[a[m].offset], edge_carrier, a[m].c1 );
                }
                else
                { /* AX_3D_BALL_CLIPPED */
                    aa = ball->vs.bc->A;
                    for (m=ball->vs.bc->na; m--;)
                    {
                        w = ball->i0 + aa[m].i;
                        v = ball->j0 + aa[m].j;
                        if (AX_IJINWIN(iw,w,v) &&
                            (tmp < base[a[m].offset]))
                            p.i4[a[m].offset] = AX_mix
                                ( p.i4[a[m].offset], edge_carrier, a[m].c1 );
                    }
                }
            }
            else  /* in ready-to-deploy pixels */
            {
                bp = ball->vs.bp + 1;
                for (m=bp[-1].offset-1; m--;)
                    if (tmp <= AX_zmap[iw][bp[m].offset])
                        AX_mem[iw].i4[bp[m].offset] = AX_mix
                            (AX_mem[iw].i4[bp[m].offset], edge_carrier,
                             bp[m].c1);
            }
        }
        
        );
    AX_BC_bp_top[iw] -= B->bp_occupied;
    return;
} /* end AX_3D_Balls_ZAprint() */


/* Find out the index of the topmost ball at a given pixel */
/* after Zpaper/ZBprint(); return -1 if there is no match. */
int AX_3D_Balls_Zgrab (int iw, AX_3D_Balls *B, int gi, int gj)
{
    register AX_3D_Ball *ball;
    register int m,goffset,k;
    register AX_Float *base, vradius, x2, y2, r2;
    register AX_BC_Unit *b;
    int i;
    if (AX_IJOUWIN(iw,gi,gj)) return(-1);
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
                        return (B->idx[i]);
                b = ball->vs.bc->a;
                for (m=ball->vs.bc->na; m--;)
                    if ((b[m].offset == k) &&
                        (ball->vx[2] < base[b[m].offset]))
                        return (B->idx[i]);
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
                        return (B->idx[i]);
                }
                for (m=1; m<AX_AOP_TOP(iw,0).b.offset; m++)
                {
                    if ( AX_OFF(iw, AX_AOP_TOP(iw,m).b.p.i,
                                AX_AOP_TOP(iw,m).b.p.j) != goffset )
                        continue;
                    if ( ball->vx[2] <= AX_zmap[iw][goffset] )
                        return (B->idx[i]);
                }
            }
        }
    }
    return(-1);
} /* end AX_3D_Balls_Zgrab() */


#ifdef _3D_Balls_TEST
#define z0      10.
#define n_balls 4
#define PNG_FILENAME "/tmp/balls.png"
#define EPS_FILENAME "/tmp/balls.eps"
int main (int argc, char *argv[])
{
    AX_3D_Balls B[1] = {0};
    AX_3D_Balls_Recreate (b, n_balls,
                          -1.5, 0., z0, 0.73,  0.8, 0.2, 0.2,
                          1.5,  0., z0, 0.435, 0.8, 0.8, 0.8,
                          0.,  1.5, z0, 0.655, 0.3, 0.3, 0.3,
                          0., -1.5, z0, 0.75,  0.2, 0.2, 0.8);
    AX_openwindow (getpid(), NULL, AX_DEFWIDTH, AX_DEFHEIGHT);
    AX_PLUGIN_3D_module(0);
    AX_PLUGIN_BC_module(0);
    AX_namedmop(0,AX_WHITE);
    AX_3D_Balls_Zdet (0, b);
    AX_3D_Balls_Zpaper (0, b);
    AX_3D_Balls_ZBprint (0, b);
    AX_3D_Balls_ZAprint (0, b);
    AX_dump(0); AX_show(0);
    AXPressKey(0);
    AX_save_pixmap_as_PNG(0, PNG_FILENAME);
    printf ("Image saved on \"%s\".\n", PNG_FILENAME);
    AX_save_pixmap_as_EPS(0, EPS_FILENAME);
    printf ("Image saved on \"%s\".\n", EPS_FILENAME);
    AX_3D_Balls_Free(b);
    AX_closewindow(0);
    Press_return();
    return (0);
}
#endif /* _3D_Balls_TEST */


#ifdef _3D_Balls_Bench_TEST
#define z0        1.
#define n_frames  100
#define n_balls   1000
int main (int argc, char *argv[])
{
    register int i,j;
    double velocity [3*n_balls];
    AX_3D_Balls B[1] = {0};
    AX_3D_Balls_Realloc (B, n_balls);
    for (i=0; i<n_balls; i++)
    {
        AX_V3FRANDOM(B->BALL[i].x);
        B->BALL[i].x[2] += z0;
        B->BALL[i].radius = Frandom() / cbrt(n_balls*10.);
        B->BALL[i].r = Frandom();
        B->BALL[i].g = Frandom();
        B->BALL[i].b = Frandom();
        AX_V3FRANDOM(&velocity[3*i]);
    }
    AX_openwindow (getpid(), NULL, AX_DEFWIDTH, AX_DEFHEIGHT);
    AX_plugin_3D_module(0);
    AX_PLUGIN_BC_module(0);
    start_chronometer();
    for (j=0; j<n_frames; j++)
    {
        AX_namedmop(0, AX_RED);
        AX_clearzmap(0);
        for (i=0; i<n_balls; i++)
            AX_V3ADDDIV (B->BALL[i].x, &velocity[3*i], n_frames);
        AX_3D_Balls_Zdet (0, B);
        AX_3D_Balls_Zpaper (0, B);
        AX_3D_Balls_ZBprint (0, B);
        AX_3D_Balls_ZAprint (0, B);
        AX_dump(0); AX_show(0);
    }
    stop_chronometer();
    printf ("3D_Balls_Bench: %d frames of %d balls took %s.\n",
            n_frames, n_balls, earthtime(stopped_usertime()));
    AX_3D_Balls_Free(B);
    AX_closewindow(0);
    Press_return();
    return (0);
}
#endif /* _3D_Balls_Bench_TEST */

/* mightyjuli 16 bit, float,  12.0 s */
/* mightyjuli 16 bit, double, 14.9 s */
/* mightyjuli 32 bit, float,  12.7 s */
/* mightyjuli 32 bit, double, 15.3 s */


/*********************/
/* Ellipsoid Support */
/*********************/

/* Realloc n_ellipsoids */
void AX_3D_Ellipsoids_Realloc (AX_3D_Ellipsoids *E, int n_ellipsoids)
{
    int i;
    E->n_ellipsoids = E->n_shown = n_ellipsoids;
    REALLOC (AX_3D_Ellipsoids_Realloc, E->ELLIPSOID, n_ellipsoids,
             AX_3D_Ellipsoid);
    REALLOC (AX_3D_Ellipsoids_Realloc, E->idx,  n_ellipsoids, int);
    for (i=0; i<n_ellipsoids; i++) E->idx[i] = i;
    return;
} /* end AX_3D_Ellipsoids_Realloc() */


/* Realloc n_ellipsoids followed by x0,y0,z0,G,r,g,b,... interface */
void AX_3D_Ellipsoids_Recreate (AX_3D_Ellipsoids *E, int n_ellipsoids, ...)
{
    int i;
    AX_Float (*G)[3];
    va_list ap;
    AX_3D_Ellipsoids_Realloc (E, n_ellipsoids);
    va_start(ap, n_ellipsoids);
    for (i=0; i<n_ellipsoids; i++)
    {
        E->ELLIPSOID[i].x[0] = va_arg(ap, double);
        E->ELLIPSOID[i].x[1] = va_arg(ap, double);
        E->ELLIPSOID[i].x[2] = va_arg(ap, double);
        G = (AX_Float (*)[3]) va_arg(ap, AX_Float *);
        AX_M3EQV (G, E->ELLIPSOID[i].G);
        AX_3D_AssignRGB (E->ELLIPSOID[i], va_arg(ap, double),
                         va_arg(ap, double), va_arg(ap, double));
    }
    va_end(ap);
    return;
} /* end AX_3D_Ellipsoids_Recreate() */


/* free allocated memory and set NULL */
void AX_3D_Ellipsoids_Free (AX_3D_Ellipsoids *E)
{
    Free (E->ELLIPSOID);
    Free (E->idx);
    return;
} /* end AX_3D_Ellipsoids_Free() */


#define QUICKSORT_M  7
static void AX_sort_3D_Ellipsoids (AX_3D_Ellipsoids *E)
{
    register int i,j,tmp;
    int idxt, ir=E->n_shown, k, l=1, jstack=0;
    int istack[128];
    register AX_Float a;
    E->idx--;
    for (;;)
    {
	if (ir-l < QUICKSORT_M)
	{
            for (j=l+1; j<=ir; j++)
	    {
		idxt = E->idx[j];
		a = E->ELLIPSOID[idxt].vx[2];
		for (i=j-1; i>=l; i--)
		{
		    if (E->ELLIPSOID[E->idx[i]].vx[2] <= a) break;
		    E->idx[i+1]=E->idx[i];
		}
		E->idx[i+1]=idxt;
	    }
	    if (jstack == 0) break;
	    ir = istack[jstack--];
	    l = istack[jstack--];
	}
	else
	{
	    k = (l+ir) >> 1;
	    SWAP(E->idx[k],E->idx[l+1],tmp);
	    if (E->ELLIPSOID[E->idx[l]].vx[2] >
                E->ELLIPSOID[E->idx[ir]].vx[2])
                SWAP(E->idx[l],E->idx[ir],tmp);
	    if (E->ELLIPSOID[E->idx[l+1]].vx[2] >
                E->ELLIPSOID[E->idx[ir]].vx[2])
                SWAP(E->idx[l+1],E->idx[ir],tmp);
	    if (E->ELLIPSOID[E->idx[l]].vx[2] >
                E->ELLIPSOID[E->idx[l+1]].vx[2])
                SWAP(E->idx[l],E->idx[l+1],tmp);
	    i = l+1;
	    j = ir;
	    idxt = E->idx[l+1];
	    a = E->ELLIPSOID[idxt].vx[2];
	    for (;;)
	    {
		do i++; while (E->ELLIPSOID[E->idx[i]].vx[2] < a);
		do j--; while (E->ELLIPSOID[E->idx[j]].vx[2] > a);
		if (j < i) break;
		SWAP (E->idx[i],E->idx[j],tmp);
	    }
	    E->idx[l+1]=E->idx[j];
	    E->idx[j]=idxt;
	    jstack += 2;
	    if (ir-i+1 >= j-l)
	    {
		istack[jstack] = ir;
		istack[jstack-1] = i;
		ir = j-1;
	    }
	    else
	    {
		istack[jstack] = j-1;
		istack[jstack-1] = l;
		l = i;
	    }
	}
    }
    E->idx++;
    return;
} /* end AX_sort_3D_Ellipsoids() */
#undef QUICKSORT_M


#define AX_3D_Ellipsoid_Determinant(a,b,c,d,e,f,x) ( \
  (a)*(x)[0]*(x)[0]+2*(b)*(x)[0]*(x)[1]+2*(c)*(x)[0]*(x)[2]+ \
  (d)*(x)[1]*(x)[1]+2*(e)*(x)[1]*(x)[2]+(f)*(x)[2]*(x)[2] )

#define AX_SOLVE_quadratic(A,B,C,D,xm) { D=(B)*(B)-(A)*(C); \
  if (D>0) xm = -(B)/(A)-sqrt(D)/ABS(A); else xm = -(B)/(A); }

/* Calculate coordinates in viewport frame and determine visibility */
void AX_3D_Ellipsoids_Zdet (int iw, AX_3D_Ellipsoids *E)
{
    register AX_3D_Ellipsoid *ellipsoid;
    register int m;
    AX_Float dx[3], X0[2], UU[3][3], K;
    int i;
    for (m=0;;m++)
    {
        for (E->n_shown=i=0; i<E->n_ellipsoids; i++)
        {
            ellipsoid = E->ELLIPSOID + i;
            AX_V3sub (ellipsoid->x, AX_3D[iw].x, dx);
            AX_M3mulV3 (AX_3D[iw].V, dx, ellipsoid->vx);
            if ((ellipsoid->vx[2] <= 0) ||
                (ellipsoid->vx[2] > AX_3D[iw].zcut)) continue;
            AX_M3MUL (AX_3D[iw].V, ellipsoid->G, UU);
            ellipsoid->a = AX_V3DOT(UU[0], AX_3D[iw].V[0]);
            ellipsoid->B = AX_V3DOT(UU[0], AX_3D[iw].V[1]);
            ellipsoid->c = AX_V3DOT(UU[0], AX_3D[iw].V[2]);
            ellipsoid->d = AX_V3DOT(UU[1], AX_3D[iw].V[1]);
            ellipsoid->e = AX_V3DOT(UU[1], AX_3D[iw].V[2]);
            ellipsoid->f = AX_V3DOT(UU[2], AX_3D[iw].V[2]);
            ellipsoid->D = AX_3D_Ellipsoid_Determinant
                (ellipsoid->a,ellipsoid->B,ellipsoid->c,
                 ellipsoid->d,ellipsoid->e,ellipsoid->f,ellipsoid->vx);
            if (ellipsoid->D<1) break; /* ejection */
            K = AX_3D[iw].k / ellipsoid->vx[2];
            X0[0] = ellipsoid->vx[0] * K + AX_3D[iw].wx;
            X0[1] = ellipsoid->vx[1] * K + AX_3D[iw].wy;
            AX_Ellipseassign
                (X0[0], X0[1],
                 (ellipsoid->a-ellipsoid->c*ellipsoid->c/ellipsoid->f)/K/K,
                 (ellipsoid->B-ellipsoid->c*ellipsoid->e/ellipsoid->f)/K/K,
                 (ellipsoid->d-ellipsoid->e*ellipsoid->e/ellipsoid->f)/K/K,
                 ellipsoid->p);
            if (AX_EllipseIntersectWindow
                (ellipsoid->p, AX_size[iw].width, AX_size[iw].height))
                E->idx[E->n_shown++] = i;
        }
        if (i < E->n_ellipsoids) /* viewpoint ejection */
        {
            if (ellipsoid->D == 0)
            {
                K = (1+AX_TINY) / sqrt(ellipsoid->f);
                AX_V3ADDMUL (AX_3D[iw].x, K, AX_3D[iw].V[2]);
            }
            else
            {
                K = (1+AX_TINY) / sqrt(ellipsoid->D) - 1;
                AX_V3SUBMUL (AX_3D[iw].x, K, dx);
            }
            if (m >= 3)
            {
                AX_V3sub (AX_3D[iw].x, E->ELLIPSOID[i].x, dx);
                AX_3D[iw].x[0] += Frandom() * dx[0];
                AX_3D[iw].x[1] += Frandom() * dx[1];
                AX_3D[iw].x[2] += Frandom() * dx[2];
                if (m >= 1000)
                    pe ("AX_3D_Ellipsoids_Zdet: viewpoint ejection failed\n");
            }
        }
        else break;
    }
    AX_sort_3D_Ellipsoids(E);
    return;
} /* end AX_3D_Ellipsoids_Zdet() */


/* Generate ellipsoid scans and register them in the z-buffer */
void AX_3D_Ellipsoids_Zpaper (int iw, AX_3D_Ellipsoids *E)
{
    register AX_3D_Ellipsoid *ellipsoid;
    register int m,n;
    int i;
    AX_Float dx[3], U[3], W[3], B,C,D,K;
    AX_Float edge_influence, lr,lg,lb, light_influence;
    for (i=0; i<E->n_shown; i++)
    {
        ellipsoid = E->ELLIPSOID + E->idx[i];
        K = AX_3D[iw].k / ellipsoid->vx[2];
        AX_EllipseWindowScan
            (ellipsoid->p, AX_size[iw].width, AX_size[iw].height,
             AX_AOP_top(iw), AX_BOP_top(iw));
        /* allocate pixel storage: will free it in Zprint() */
        MALLOC (AX_3D_Ellipsoids_Zpaper, ellipsoid->img,
                AX_AOP_TOP(iw,0).b.offset + AX_BOP_TOP(iw,0).offset,
                AX_3D_Pixel);
        for (n=m=1; m<AX_AOP_TOP(iw,0).b.offset; m++)
        {
            ellipsoid->img[n].offset = AX_OFF
                (iw, AX_AOP_TOP(iw,m).b.p.i, AX_AOP_TOP(iw,m).b.p.j);
            dx[0] = (AX_AOP_TOP(iw,m).b.p.i + 0.5 - ellipsoid->p->x0) / K;
            dx[1] = (AX_AOP_TOP(iw,m).b.p.j + 0.5 - ellipsoid->p->y0) / K;
            B = ellipsoid->c*dx[0] + ellipsoid->e*dx[1];
            C = ellipsoid->a*dx[0]*dx[0] + 2*ellipsoid->B*dx[0]*dx[1] +
                ellipsoid->d*dx[1]*dx[1];
            AX_SOLVE_quadratic(ellipsoid->f,B,C,ellipsoid->D,dx[2]);
            ellipsoid->img[n].z = ellipsoid->vx[2] + dx[2];
            if (ellipsoid->img[n].z < AX_zmap[iw][ellipsoid->img[n].offset])
            {
                ellipsoid->img[n].c.area = AX_AOP_TOP(iw,m).c.area;
                n++;
            }
        }
        ellipsoid->img[0].offset = (n++);
        AX_3D_LRGB (AX_ES_TOP*ellipsoid->r, AX_ES_TOP*ellipsoid->g,
                    AX_ES_TOP*ellipsoid->b, lr, lg, lb);
        for (m=1; m<AX_BOP_TOP(iw,0).offset; m++)
        {
            ellipsoid->img[n].offset = AX_OFF
                (iw, AX_BOP_TOP(iw,m).p.i, AX_BOP_TOP(iw,m).p.j);
            dx[0] = (AX_BOP_TOP(iw,m).p.i + 0.5 - ellipsoid->p->x0) / K;
            dx[1] = (AX_BOP_TOP(iw,m).p.j + 0.5 - ellipsoid->p->y0) / K;
            B = ellipsoid->c*dx[0] + ellipsoid->e*dx[1];
            C = ellipsoid->a*dx[0]*dx[0] + 2*ellipsoid->B*dx[0]*dx[1] +
                ellipsoid->d*dx[1]*dx[1] - 1;
            AX_SOLVE_quadratic(ellipsoid->f,B,C,ellipsoid->D,dx[2]);
            ellipsoid->img[n].z = ellipsoid->vx[2] + dx[2];
            if (ellipsoid->img[n].z < AX_zmap[iw][ellipsoid->img[n].offset])
            {
                U[0] = ellipsoid->a*dx[0] + ellipsoid->B*dx[1] +
                    ellipsoid->c*dx[2];
                U[1] = ellipsoid->B*dx[0] + ellipsoid->d*dx[1] +
                    ellipsoid->e*dx[2];
                U[2] = ellipsoid->c*dx[0] + ellipsoid->e*dx[1] +
                    ellipsoid->f*dx[2];
                AX_V3NORMALIZE (U, D);
                AX_V3add(ellipsoid->vx, dx, W);
                edge_influence = -AX_V3DOT(U,W) / AX_V3LENGTH(W);
                light_influence = edge_influence * AX_3D_LDOT(U);
                if (light_influence < 0) light_influence = 0;
                edge_influence = AX_ES_EDGE + AX_ES_TOPI * edge_influence;
                ellipsoid->img[n].c.carrier = AX_Colorcarrier
                    ( edge_influence * ellipsoid->r + light_influence * lr,
                      edge_influence * ellipsoid->g + light_influence * lg,
                      edge_influence * ellipsoid->b + light_influence * lb );
                AX_zmap[iw][ellipsoid->img[n].offset] = ellipsoid->img[n].z;
                n++;
            }
        }
        ellipsoid->img[ellipsoid->img[0].offset].offset =
            n - ellipsoid->img[0].offset;
    }
    return;
} /* end AX_3D_Ellipsoids_Zpaper() */


/* draw the ellipsoids to pixmap (block pixels only) according to z-buffer */
void AX_3D_Ellipsoids_ZBprint (int iw, AX_3D_Ellipsoids *E)
{
    register AX_3D_Ellipsoid *ellipsoid;
    register int m;
    int i;
    AX_C(
        
        for (i=0; i<E->n_shown; i++)
        {
            ellipsoid = E->ELLIPSOID + E->idx[i];
            for (m=ellipsoid->img[0].offset+1;
                 m<ellipsoid->img[0].offset+ellipsoid->img
                     [ellipsoid->img[0].offset].offset; m++)
                if (ellipsoid->img[m].z <=
                    AX_zmap[iw][ellipsoid->img[m].offset])
                    AX_mem[iw].uc[ellipsoid->img[m].offset] =
                        ellipsoid->img[m].c.carrier;
        },

        for (i=0; i<E->n_shown; i++)
        {
            ellipsoid = E->ELLIPSOID + E->idx[i];
            for (m=ellipsoid->img[0].offset+1;
                 m<ellipsoid->img[0].offset+ellipsoid->img
                     [ellipsoid->img[0].offset].offset; m++)
                if (ellipsoid->img[m].z <=
                    AX_zmap[iw][ellipsoid->img[m].offset])
                    AX_mem[iw].i2[ellipsoid->img[m].offset] =
                        ellipsoid->img[m].c.carrier;
        },

        for (i=0; i<E->n_shown; i++)
        {
            ellipsoid = E->ELLIPSOID + E->idx[i];
            for (m=ellipsoid->img[0].offset+1;
                 m<ellipsoid->img[0].offset+ellipsoid->img
                     [ellipsoid->img[0].offset].offset; m++)
                if (ellipsoid->img[m].z <=
                    AX_zmap[iw][ellipsoid->img[m].offset])
                    AX_mem[iw].i4[ellipsoid->img[m].offset] =
                        ellipsoid->img[m].c.carrier;
        }
        
        );
    return;
} /* end AX_3D_Ellipsoids_ZBprint() */


/* draw the ellipsoids to pixmap (alias pixels only) according to z-buffer */
void AX_3D_Ellipsoids_ZAprint (int iw, AX_3D_Ellipsoids *E)
{
    register AX_3D_Ellipsoid *ellipsoid;
    register int m;
    register AX_Carrier edge_carrier;
    int i;
    AX_C(

        for (i=E->n_shown-1; i>=0; i--)
        {
            ellipsoid = E->ELLIPSOID + E->idx[i];
            edge_carrier = AX_Colorcarrier
                (AX_ES_EDGE*ellipsoid->r,AX_ES_EDGE*ellipsoid->g,
                 AX_ES_EDGE*ellipsoid->b);
            for (m=1; m<ellipsoid->img[0].offset; m++)
                if (ellipsoid->img[m].z <=
                    AX_zmap[iw][ellipsoid->img[m].offset])
                    AX_mem[iw].uc[ellipsoid->img[m].offset] = AX_mix
                        ( AX_mem[iw].uc[ellipsoid->img[m].offset],
                          edge_carrier, ellipsoid->img[m].c.area );
            free (ellipsoid->img);
        },

        for (i=E->n_shown-1; i>=0; i--)
        {
            ellipsoid = E->ELLIPSOID + E->idx[i];
            edge_carrier = AX_Colorcarrier
                (AX_ES_EDGE*ellipsoid->r, AX_ES_EDGE*ellipsoid->g,
                 AX_ES_EDGE*ellipsoid->b);
            for (m=1; m<ellipsoid->img[0].offset; m++)
                if (ellipsoid->img[m].z <=
                    AX_zmap[iw][ellipsoid->img[m].offset])
                    AX_mem[iw].i2[ellipsoid->img[m].offset] = AX_mix
                        ( AX_mem[iw].i2[ellipsoid->img[m].offset],
                          edge_carrier, ellipsoid->img[m].c.area );
            free (ellipsoid->img);
        },

        for (i=E->n_shown-1; i>=0; i--)
        {
            ellipsoid = E->ELLIPSOID + E->idx[i];
            edge_carrier = AX_Colorcarrier
                (AX_ES_EDGE*ellipsoid->r,AX_ES_EDGE*ellipsoid->g,
                 AX_ES_EDGE*ellipsoid->b);
            for (m=1; m<ellipsoid->img[0].offset; m++)
                if (ellipsoid->img[m].z <=
                    AX_zmap[iw][ellipsoid->img[m].offset])
                    AX_mem[iw].i4[ellipsoid->img[m].offset] = AX_mix
                        ( AX_mem[iw].i4[ellipsoid->img[m].offset],
                          edge_carrier, ellipsoid->img[m].c.area );
            free (ellipsoid->img);
        }
        
        );
    return;
} /* end AX_3D_Ellipsoids_ZAprint() */


#ifdef _3D_Ellipsoids_TEST
#define Z0        7
#define boxwidth  3.
#define boxheight 3.
#define boxthickness 3.
#define n_ellipsoids 2
#define PNG_FILENAME "/tmp/ellipsoids.png"
#define EPS_FILENAME "/tmp/ellipsoids.eps"
int main (int argc, char *argv[])
{
    AX_3D_Ellipsoids E[1]={0};
    AX_3D_Define_PW(L);
    AX_Float H[3][3],R[3][3],A[3][3],G0[3][3],G1[3][3];
    AX_M3RANDOMROTATION(R,H[0][0]);
    AX_M3TRANSPOSE(R,H);
    AX_M3MUL(R,H,A);
    AX_M3DIAGONAL(boxwidth,boxheight,boxthickness,A);
    AX_M3MUL(A,R,H);
    AX_M3IDENTITY(G0);
    G0[0][0] = 0.5;
    G0[1][1] = 0.7;
    AX_M3IDENTITY(G1);
    G1[1][0] = G1[0][1] = 0.3;
    AX_3D_Ellipsoids_Recreate (E, n_ellipsoids,
                               -1.5, 0., 6., G0, 1., 0., 0.,
                               1.5, 0., 6., G1, .5, .5, .5 );
    AX_3D_PW_Assign(-boxwidth/2, -boxheight/2, Z0, H, 0, 0, 0, L);
    AX_openwindow (getpid(), NULL, AX_DEFWIDTH, AX_DEFHEIGHT);
    AX_PLUGIN_3D_module(0);
    AX_namedmop(0,AX_WHITE);
    AX_3D_Ellipsoids_Zdet (0, E);
    AX_3D_Ellipsoids_Zpaper (0, E);
    AX_3D_Ellipsoids_ZBprint (0, E);
    AX_3D_Ellipsoids_ZAprint (0, E);
    AX_3D_Lines_Zprint (0, L);
    AX_dump(0); AX_show(0);
    AX_openwindow (getpid(), NULL, AX_DEFWIDTH, AX_DEFHEIGHT);
    AX_PLUGIN_3D_module(1);
    AX_namedmop(1,AX_WHITE);
    AX_3D_Ellipsoids_Zdet (1, E);
    AX_3D_Ellipsoids_Zpaper (1, E);
    AX_3D_Ellipsoids_ZBprint (1, E);
    /* AX_3D_Ellipsoids_ZAprint (1, E); */
    AX_3D_Lines_Zprint (1, L);
    AX_dump(1); AX_show(1);
    AX_save_pixmap_as_PNG(0, PNG_FILENAME);
    printf ("Image saved on \"%s\".\n", PNG_FILENAME);
    AX_save_pixmap_as_EPS(0, EPS_FILENAME);
    printf ("Image saved on \"%s\".\n", EPS_FILENAME);
    Press_return();
    AX_closewindow(0);
    Press_return();
    AX_closewindow(1);
    Press_return();
    return (0);
}
#endif /* _3D_Ellipsoids_TEST */


/********************/
/* Cylinder Support */
/********************/

/* Realloc n_cylinders */
void AX_3D_Cylinders_Realloc (AX_3D_Cylinders *C, int n_cylinders)
{
    C->n_cylinders = n_cylinders;
    REALLOC(AX_3D_Cylinders_Realloc, C->CYLINDER, n_cylinders, AX_3D_Cylinder);
    Free (C->stack);
    Free (C->cylinder);
    return;
} /* end AX_3D_Cylinders_Realloc() */


/* Realloc n_cylinders followed by x0,y0,z0,x1,y1,z1,radius,r,g,b,... */
void AX_3D_Cylinders_Recreate (AX_3D_Cylinders *C, int n_cylinders, ...)
{
    int i;
    va_list ap;
    double x1[3];
    AX_3D_Cylinders_Realloc (C, n_cylinders);
    va_start(ap, n_cylinders);
    for (i=0; i<n_cylinders; i++)
    {
        C->CYLINDER[i].x0[0] = va_arg(ap, double);
        C->CYLINDER[i].x0[1] = va_arg(ap, double);
        C->CYLINDER[i].x0[2] = va_arg(ap, double);
        x1[0] = va_arg(ap, double);
        x1[1] = va_arg(ap, double);
        x1[2] = va_arg(ap, double);
        AX_3D_CYLINDER_ARM_AXIS (C->CYLINDER[i], x1);
        C->CYLINDER[i].radius = va_arg(ap, double);
        AX_3D_AssignRGB (C->CYLINDER[i], va_arg(ap, double),
                         va_arg(ap, double), va_arg(ap, double));
    }
    va_end(ap);
    return;
} /* end AX_3D_Cylinders_Recreate() */


/* Realloc n_cylinders followed by x0[],x1[],radius,r,g,b,... */
void AX_3D_Cylinders_RECREATE (AX_3D_Cylinders *C, int n_cylinders, ...)
{
    int i;
    AX_Float *x;
    va_list ap;
    AX_3D_Cylinders_Realloc (C, n_cylinders);
    va_start(ap, n_cylinders);
    for (i=0; i<n_cylinders; i++)
    {
        x = va_arg(ap, AX_Float *);
        AX_V3EQV (x, C->CYLINDER[i].x0);
        x = va_arg(ap, AX_Float *);
        AX_3D_CYLINDER_ARM_AXIS (C->CYLINDER[i], x);
        C->CYLINDER[i].radius = va_arg(ap, double);
        AX_3D_AssignRGB (C->CYLINDER[i], va_arg(ap, double),
                         va_arg(ap, double), va_arg(ap, double));
    }
    va_end(ap);
    return;
} /* end AX_3D_Cylinders_RECREATE() */


/* free allocated memory and set NULL */
void AX_3D_Cylinders_Free (AX_3D_Cylinders *C)
{
    Free (C->stack);
    Free (C->cylinder);
    Free (C->CYLINDER);
    return;
} /* end AX_3D_Cylinders_Free() */


#define AX_3D_CYLINDER_REIL  0.40
/* Gets intercepted before passing the z-plane */
#define AX_3D_CYLINDER_VEIL  0.01
/* #define AX_3D_CYLINDER_VEIL  0.001 */

#define AX_3D_CYLINDERS_MINAREA2      1.0
#define AX_3D_CYLINDERS_STACK_RATIO   2.1

#define AX_3D_CYLINDER_MESH  12
const AX_Float table[2*AX_3D_CYLINDER_MESH] =
{0.866025,0.500000, 0.500000,0.866025, 0.,1.,
 -0.500000,0.866025, -0.866025,0.500000, -1.,0.,
 -0.866025,-0.500000, -0.500000,-0.866025, 0.,-1.,
 0.500000,-0.866025, 0.866025,-0.500000, 1.,0.};

#define AX_SOLVE_QUADRATIC(A,B,C,D,xm,xp) { (D)=(B)*(B)-(A)*(C); \
  if ((D)>0) { (D)=sqrt(D)/ABS(A); (xm)=-(D)-(B)/(A); (xp)=(D)-(B)/(A); } \
  else (xm)=(xp)=-(B)/(A); }

/* Calculate coordinates in viewport frame and determine visibility */
void AX_3D_Cylinders_Zdet (int iw, AX_3D_Cylinders *C)
{
    register AX_3D_Cylinder *cylinder;
    register AX_3D_Cylinder_Power *power;
    register double tmp=0;
    int i, m;
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
                        if ( ( AX_V3DOT(AX_3D[iw].fp[j].dx, cylinder->x0) >
                               AX_3D[iw].fp[j].d0-AX_TINY ) ||
                             ( AX_V3DOT(AX_3D[iw].fp[j].dx, x1) >
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
        if (i < C->n_cylinders)
        {  /* viewpoint ejection */
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
            if (m>=1000)
                pe ("AX_3D_Cylinders_Zdet: viewpoint ejection failed\n");
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


/* Generate cylinder scans and register them in the z-buffer */
void AX_3D_Cylinders_Zpaper(int iw, AX_3D_Cylinders *C)
{
    register AX_3D_Cylinder_Power *power;
    register AX_3D_Cylinder *cylinder;
    register AX_3D_Pixel *img;
    register int m,n;
    register double tmp;
    int i, stack_top, n_power, stack_max;
    double UU[3], xx[3], X0[2], A, B, edge_influence, lr,lg,lb,
        light_influence;
    AX_3D_Cylinder **cylinders;
    stack_max = INT ( AX_3D_CYLINDERS_STACK_RATIO *
                      AX_size[iw].width * AX_size[iw].height );
    REALLOC( AX_3D_Cylinders_Zpaper, C->stack, stack_max, AX_3D_Pixel );
    n_power = C->n_power;
    for (C->n_power=stack_top=i=0; i<n_power; i++)
    {
        power = C->power + C->idx[i];
        cylinder = C->cylinder[C->idx[i]];
        AX_PolygonScan (power->p, AX_AOP_top(iw), AX_BOP_top(iw));
        if (AX_BOP_TOP(iw,0).offset == 1) continue;
        if (stack_top+AX_AOP_TOP(iw,0).b.offset+AX_BOP_TOP(iw,0).offset
            >= stack_max) break;
        img = C->stack + stack_top;
        for (n=m=1; m<AX_AOP_TOP(iw,0).b.offset; m++)
        {
            img[n].offset = AX_OFF
                (iw, AX_AOP_TOP(iw,m).b.p.i, AX_AOP_TOP(iw,m).b.p.j);
            X0[0] = AX_AOP_TOP(iw,m).b.p.i + 0.5 - AX_3D[iw].wx;
            X0[1] = AX_AOP_TOP(iw,m).b.p.j + 0.5 - AX_3D[iw].wy;
            tmp = (power->U[0]*X0[0]+power->U[1]*X0[1]) / AX_3D[iw].k +
                power->U[2];
            A = (X0[0]*X0[0]+X0[1]*X0[1])/power->dx[2]+1-tmp*tmp;
            B = (-power->y0[0]*X0[0]-power->y0[1]*X0[1])/AX_3D[iw].k-
                power->y0[2]-power->dx[0]*tmp;
            AX_SOLVE_quadratic(A,B,power->dx[1],tmp,img[n].z);
            if (img[n].z < AX_zmap[iw][img[n].offset])
            {
                xx[0] = X0[0] / AX_3D[iw].k * img[n].z;
                xx[1] = X0[1] / AX_3D[iw].k * img[n].z;
                xx[2] = img[n].z;
                AX_V3sub (xx, power->y0, UU);
                tmp = AX_V3DOT (UU, power->U);
                if ((tmp<0) || (tmp>power->U[3])) continue;
                img[n].c.area = AX_AOP_TOP(iw,m).c.area;
                n++;
            }
        }
        img[0].offset = n++;
        AX_3D_LRGB (AX_CS_TOP*cylinder->r, AX_CS_TOP*cylinder->g,
                    AX_CS_TOP*cylinder->b, lr,lg,lb);
        for (m=1; m<AX_BOP_TOP(iw,0).offset; m++)
        {
            img[n].offset = AX_OFF
                (iw, AX_BOP_TOP(iw,m).p.i, AX_BOP_TOP(iw,m).p.j);
            X0[0] = AX_BOP_TOP(iw,m).p.i + 0.5 - AX_3D[iw].wx;
            X0[1] = AX_BOP_TOP(iw,m).p.j + 0.5 - AX_3D[iw].wy;
            tmp = (power->U[0]*X0[0]+power->U[1]*X0[1])
                /AX_3D[iw].k+power->U[2];
            A = (X0[0]*X0[0]+X0[1]*X0[1])/power->dx[2]+1-tmp*tmp;
            B = (-power->y0[0]*X0[0]-power->y0[1]*X0[1])/AX_3D[iw].k-
                power->y0[2]-power->dx[0]*tmp;
            AX_SOLVE_quadratic(A,B,power->dx[1],tmp,img[n].z);
            if (img[n].z < AX_zmap[iw][img[n].offset])
            {
                xx[0] = X0[0] / AX_3D[iw].k * img[n].z;
                xx[1] = X0[1] / AX_3D[iw].k * img[n].z;
                xx[2] = img[n].z;
                AX_V3sub (xx, power->y0, UU);
                tmp = AX_V3DOT(UU, power->U);
                if ((tmp<0) || (tmp>power->U[3])) continue;
                AX_V3SUBMUL(UU,tmp,power->U);
                edge_influence = -AX_V3DOT(UU,xx) / cylinder->radius
                    / AX_V3LENGTH(xx);
                light_influence = edge_influence * AX_3D_LDOT(UU)
                    / cylinder->radius;
                if (light_influence < 0) light_influence = 0;
                edge_influence = AX_CS_EDGE + AX_CS_TOPI * edge_influence;
                if ( edge_influence < 0 ) edge_influence = 0;
                img[n].c.carrier = AX_Colorcarrier
                    ( edge_influence * cylinder->r + light_influence * lr,
                      edge_influence * cylinder->g + light_influence * lg,
                      edge_influence * cylinder->b + light_influence * lb );
                AX_zmap[iw][img[n].offset] = img[n].z;
                n++;
            }
        }
        if (n > 2)
        {
            img[img[0].offset].offset = n - img[0].offset;
            cylinder->img = img;
            stack_top += n;
            C->idx[C->n_power++] = C->idx[i];
        }
    }
    Free(C->power);
    REALLOC( AX_3D_Cylinders_Zpaper, C->stack, stack_top, AX_3D_Pixel );
    MALLOC (AX_3D_Cylinders_Zpaper, cylinders, C->n_power, AX_3D_Cylinder *);
    for (i=0; i<C->n_power; i++) cylinders[i] = C->cylinder[C->idx[i]];
    Free(C->idx);
    Free(C->cylinder);
    C->cylinder = cylinders;
    return;
} /* end AX_3D_Cylinders_Zpaper() */


/* draw the cylinders to pixmap (block pixels only) according to z-buffer */
void AX_3D_Cylinders_ZBprint (int iw, AX_3D_Cylinders *C)
{
    register AX_3D_Pixel *img;
    register int m;
    int i;
    AX_C(

        for (i=0; i<C->n_power; i++)
        {
            img = C->cylinder[i]->img;
            for (m=img[0].offset+1;
                 m<img[0].offset+img[img[0].offset].offset; m++)
                if (img[m].z <= AX_zmap[iw][img[m].offset])
                    AX_mem[iw].uc[img[m].offset] = img[m].c.carrier;
        },

        for (i=0; i<C->n_power; i++)
        {
            img = C->cylinder[i]->img;
            for (m=img[0].offset+1;
                 m<img[0].offset+img[img[0].offset].offset; m++)
                if (img[m].z <= AX_zmap[iw][img[m].offset])
                    AX_mem[iw].i2[img[m].offset] = img[m].c.carrier;
        },

        for (i=0; i<C->n_power; i++)
        {
            img = C->cylinder[i]->img;
            for (m=img[0].offset+1;
                 m<img[0].offset+img[img[0].offset].offset; m++)
                if (img[m].z <= AX_zmap[iw][img[m].offset])
                    AX_mem[iw].i4[img[m].offset] = img[m].c.carrier;
        }
        
        );
    return;
} /* end AX_3D_Cylinders_ZBprint() */


/* draw the cylinders to pixmap (alias pixels only) according to z-buffer */
void AX_3D_Cylinders_ZAprint (int iw, AX_3D_Cylinders *C)
{
    register AX_3D_Pixel *img;
    register int m;
    register AX_Carrier edge_carrier;
    int i;
    AX_C(
        
        for (i=C->n_power-1; i>=0; i--)
        {
            img = C->cylinder[i]->img;
            edge_carrier = AX_Colorcarrier
                (AX_ES_EDGE*C->cylinder[i]->r,
                 AX_ES_EDGE*C->cylinder[i]->g,
                 AX_ES_EDGE*C->cylinder[i]->b);
            for (m=1; m<img[0].offset; m++)
                if (img[m].z <= AX_zmap[iw][img[m].offset])
                    AX_mem[iw].uc[img[m].offset] = AX_mix
                        ( AX_mem[iw].uc[img[m].offset],
                          edge_carrier, img[m].c.area );
        },

        for (i=C->n_power-1; i>=0; i--)
        {
            img = C->cylinder[i]->img;
            edge_carrier = AX_Colorcarrier
                (AX_ES_EDGE*C->cylinder[i]->r,
                 AX_ES_EDGE*C->cylinder[i]->g,
                 AX_ES_EDGE*C->cylinder[i]->b);
            for (m=1; m<img[0].offset; m++)
                if (img[m].z <= AX_zmap[iw][img[m].offset])
                    AX_mem[iw].i2[img[m].offset] = AX_mix
                        ( AX_mem[iw].i2[img[m].offset],
                          edge_carrier, img[m].c.area );
        },

        for (i=C->n_power-1; i>=0; i--)
        {
            img = C->cylinder[i]->img;
            edge_carrier = AX_Colorcarrier
                (AX_ES_EDGE*C->cylinder[i]->r,
                 AX_ES_EDGE*C->cylinder[i]->g,
                 AX_ES_EDGE*C->cylinder[i]->b);
            for (m=1; m<img[0].offset; m++)
                if (img[m].z <= AX_zmap[iw][img[m].offset])
                    AX_mem[iw].i4[img[m].offset] = AX_mix
                        ( AX_mem[iw].i4[img[m].offset],
                          edge_carrier, img[m].c.area );
        }
        
        );
    return;
} /* end AX_3D_Cylinders_ZAprint() */


/* Find out the index of the topmost cylinder at a given pixel */
/* after Zpaper/ZBprint(); return -1 if there is no match.     */
int AX_3D_Cylinders_Zgrab (int iw, AX_3D_Cylinders *C, int gi, int gj)
{
    register AX_3D_Pixel *img;
    register int m,goffset;
    int i;
    if (AX_IJOUWIN(iw,gi,gj)) return(-1);
    goffset = AX_OFF(iw,gi,gj);
    for (i=0; i<C->n_power; i++)
    {
        img = C->cylinder[i]->img;
        for (m=1; m<img[0].offset; m++)
            if (goffset == img[m].offset)
                if (img[m].z <= AX_zmap[iw][goffset])
                    return (C->cylinder[i]-C->CYLINDER);
        for (m=img[0].offset+1;
             m<img[0].offset+img[img[0].offset].offset; m++)
            if (goffset == img[m].offset)
                if (img[m].z <= AX_zmap[iw][goffset])
                    return (C->cylinder[i]-C->CYLINDER);
    }
    return(-1);
} /* end AX_3D_Cylinders_Zgrab() */


#ifdef _3D_Cylinders_TEST
#define n_balls        4
#define n_cylinders    2
#define PNG_FILENAME  "/tmp/cylinders.png"
#define EPS_FILENAME  "/tmp/cylinders.eps"
int main (int argc, char *argv[])
{
    AX_3D_Balls B[1] = {0};
    AX_3D_Cylinders C[1] = {0};
    AX_3D_Balls_Recreate (B, n_balls,
                          -3., 0., 10., 1.5, 1.,  0.,  0.,
                          3.,  0., 10.,  2., 0.5, 0.5, 0.5,
                          0., -3., 10.,  2., 1.,  0.,  0.,
                          0.,  3., 10.,  2., 0.5, 0.5, 0.5);
    AX_3D_Cylinders_RECREATE (C, n_cylinders,
                              B->BALL[0].x, B->BALL[1].x, 0.51, 0., 0.5, 0.5,
                              B->BALL[2].x, B->BALL[3].x, 0.51, 0., 0.5, 0.5);
    AX_openwindow (getpid(), NULL, AX_DEFWIDTH, AX_DEFHEIGHT);
    AX_PLUGIN_3D_module(0);
    /* AX_PLUGIN_BC_module(0); */
    AX_namedmop(0,AX_WHITE);
    /* AX_3D_Balls_Zdet (0, b); */
    AX_3D_Cylinders_Zdet (0, C);
    /* AX_3D_Balls_Zpaper (0, b); */
    AX_3D_Cylinders_Zpaper (0, C);
    /* AX_3D_Balls_ZBprint (0, b); */
    AX_3D_Cylinders_ZBprint (0, C);
    /* AX_3D_Balls_ZAprint (0, b); */
    AX_3D_Cylinders_ZAprint (0, C);
    AX_dump(0); AX_show(0);
    AXPressKey(0);
    AX_save_pixmap_as_PNG(0, PNG_FILENAME);
    printf ("Image saved on \"%s\".\n", PNG_FILENAME);
    AX_save_pixmap_as_EPS(0, EPS_FILENAME);
    printf ("Image saved on \"%s\".\n", EPS_FILENAME);
    AX_3D_Balls_Free(B);
    AX_3D_Cylinders_Free(C);
    AX_closewindow(0);
    Press_return();
    return(0);
}
#endif /* _3D_Cylinders_TEST */
