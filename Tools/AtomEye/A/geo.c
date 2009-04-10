/*******************************************/
/* Threaded Atomistic Configuration Viewer */
/*                                         */
/* Fri Mar 31 2000 Ju Li <liju99@mit.edu>  */
/*******************************************/

#include "A.h"

double *geo[MAX_GEO_MEASURES] = {0};
GeoList geolist[MAX_GEO_MEASURES] = {
    {"shear strain", &evaluate_shear_strain, 1},
    {"central symm", &evaluate_central_symm, 0},
};
double avg_M[3][3], avg_shear_strain, avg_coordination;
short *coordination=NULL, coordination_crystal;
int coordination_hist[ATOM_COORDINATION_MAX+1];
int shear_strain_subtract_mean=0, central_symm_neighbormax=-1;

    
static void addsort_distance2 (double *d2, int *p, double distance2)
{
    register int m,n;
    if ((*p) < coordination_crystal)
    {
        for (m=0; m<(*p); m++)
            if (distance2 < d2[m])
            {
                for (n=(*p)++; n>m; n--) d2[n] = d2[n-1];
                d2[m] = distance2;
                return;
            }
        d2[(*p)++] = distance2;
        return;
    }
    for (m=0; m<coordination_crystal; m++)
        if (distance2 < d2[m])
        {
            for (n=coordination_crystal-1; n>m; n--) d2[n] = d2[n-1];
            d2[m] = distance2;
            return;
        }
    return;
} /* end addsort_distance2() */


void evaluate_shear_strain (double **shear_strain)
{
    register int i,j;
    double tmp[DIMENSION],dx[DIMENSION+1],metric[6];
    int *participant;
    double *distance2;
    REALLOC (evaluate_shear_strain, *shear_strain, 6*np, double);
    bzero ((void *)(*shear_strain), (6*np)*sizeof(double));
    REALLOC (evaluate_shear_strain, coordination, np, short);
    bzero ((void *)(coordination), np*sizeof(short));
    for (i=0; i<np; i++)
        for (j=N->idx[i]; j<N->idx[i+1]; j++)
        {
            coordination[i]++;
            coordination[N->list[j]]++;
        }
    for (i=0; i<=ATOM_COORDINATION_MAX; i++) coordination_hist[i] = 0;
    for (i=0; i<np; i++)
        if (coordination[i] > ATOM_COORDINATION_MAX)
            pe ("evaluate_shear_strain:\n"
                "ATOM_COORDINATION_MAX = %d exceeded\n",ATOM_COORDINATION_MAX);
        else coordination_hist[coordination[i]]++;
    j = -1;
    for (i=0; i<=ATOM_COORDINATION_MAX; i++)
        if (coordination_hist[i] > j)
        { /* the most populous coordination number */
            j = coordination_hist[i];
            coordination_crystal = i;
        }
    CALLOC (evaluate_shear_strain, participant, np, int);
    MALLOC (evaluate_shear_strain, distance2, coordination_crystal*np,
            double);
    for (i=0; i<np; i++)
        for (j=N->idx[i]; j<N->idx[i+1]; j++)
        {
            V3SUB (&s[DIMENSION*N->list[j]], &s[DIMENSION*i], tmp);
            V3ImagE (tmp);
            V3M3LENGTH2 (tmp, H, dx);
            addsort_distance2 (&distance2[coordination_crystal*i],
                               &participant[i], dx[3]);
            addsort_distance2 (&distance2[coordination_crystal*N->list[j]],
                               &participant[N->list[j]], dx[3]);
        }
    dx[0] = dx[1] = 0;
    for (i=0; i<np; i++)
    {
        if (participant[i] == coordination_crystal)
        {
            for (j=0; j<participant[i]; j++)
                dx[1] += distance2[coordination_crystal*i+j];
            dx[0] += coordination_crystal;
        }
        if (participant[i] > 0)
            distance2[coordination_crystal*i] =
                distance2[coordination_crystal*i+participant[i]-1] + SMALL;
        else distance2[coordination_crystal*i] = SMALL;
        participant[i] = 0;
    }
    dx[3] = dx[1] / dx[0] / DIMENSION;
    for (i=0; i<np; i++)
        for (j=N->idx[i]; j<N->idx[i+1]; j++)
        {
            V3SUB (&s[DIMENSION*N->list[j]], &s[DIMENSION*i], tmp);
            V3ImagE (tmp);
            V3mM3 (tmp, H, dx);
            if ( (participant[i] < coordination_crystal) &&
                 (V3LENGTH2(dx) < distance2[coordination_crystal*i]) )
            {
                (*shear_strain)[6*i]   += dx[0]*dx[0] / dx[3];
                (*shear_strain)[6*i+1] += dx[0]*dx[1] / dx[3];
                (*shear_strain)[6*i+2] += dx[0]*dx[2] / dx[3];
                (*shear_strain)[6*i+3] += dx[1]*dx[1] / dx[3];
                (*shear_strain)[6*i+4] += dx[1]*dx[2] / dx[3];
                (*shear_strain)[6*i+5] += dx[2]*dx[2] / dx[3];
                participant[i]++;
            }
            if ( (participant[N->list[j]] < coordination_crystal) &&
                 (V3LENGTH2(dx) < distance2[coordination_crystal*N->list[j]]) )
            {
                (*shear_strain)[6*N->list[j]]   += dx[0]*dx[0] / dx[3];
                (*shear_strain)[6*N->list[j]+1] += dx[0]*dx[1] / dx[3];
                (*shear_strain)[6*N->list[j]+2] += dx[0]*dx[2] / dx[3];
                (*shear_strain)[6*N->list[j]+3] += dx[1]*dx[1] / dx[3];
                (*shear_strain)[6*N->list[j]+4] += dx[1]*dx[2] / dx[3];
                (*shear_strain)[6*N->list[j]+5] += dx[2]*dx[2] / dx[3];
                participant[N->list[j]]++;
            }
        }
    free (distance2);
    for (j=0; j<6; j++) metric[j] = 0;
    for (avg_coordination=i=0; i<np; i++)
    {
        avg_coordination += coordination[i];
        for (j=0; j<6; j++)
        {
            if (participant[i] > 0) (*shear_strain)[6*i+j] /= participant[i];
            metric[j] += (*shear_strain)[6*i+j];
        }
    }
    free(participant);
    avg_coordination /= np;
    for (j=0; j<6; j++) metric[j] /= np;
    for (i=0; i<np; i++)
    {
        if (shear_strain_subtract_mean)
        {
            avg_M[0][0]               = (*shear_strain)[6*i]   - metric[0];
            avg_M[0][1] = avg_M[1][0] = (*shear_strain)[6*i+1] - metric[1];
            avg_M[0][2] = avg_M[2][0] = (*shear_strain)[6*i+2] - metric[2];
            avg_M[1][1]               = (*shear_strain)[6*i+3] - metric[3];
            avg_M[1][2] = avg_M[2][1] = (*shear_strain)[6*i+4] - metric[4];
            avg_M[2][2]               = (*shear_strain)[6*i+5] - metric[5];
        }
        else
        {
            avg_M[0][0]               = (*shear_strain)[6*i];
            avg_M[0][1] = avg_M[1][0] = (*shear_strain)[6*i+1];
            avg_M[0][2] = avg_M[2][0] = (*shear_strain)[6*i+2];
            avg_M[1][1]               = (*shear_strain)[6*i+3];
            avg_M[1][2] = avg_M[2][1] = (*shear_strain)[6*i+4];
            avg_M[2][2]               = (*shear_strain)[6*i+5];
        }
        (*shear_strain)[i] = SymmetricM3MisesInvariant(avg_M) / 2;
    }
    REALLOC (main, *shear_strain, np, double);
    avg_M[0][0]               = metric[0];
    avg_M[0][1] = avg_M[1][0] = metric[1];
    avg_M[0][2] = avg_M[2][0] = metric[2];
    avg_M[1][1]               = metric[3];
    avg_M[1][2] = avg_M[2][1] = metric[4];
    avg_M[2][2]               = metric[5];
    avg_shear_strain = SymmetricM3MisesInvariant(avg_M) / 2;
    return;
} /* end evaluate_shear_strain() */


bool change_shear_strain_subtract_mean (int iw)
{
    shear_strain_subtract_mean =
        1 - shear_strain_subtract_mean;
    evaluate_shear_strain (geo+GEO_SHEAR_STRAIN);
    n[iw].color_mode = COLOR_MODE_AUXILIARY;
    n[iw].auxiliary_idx = CONFIG_MAX_AUXILIARY + GEO_SHEAR_STRAIN;
    return (color_encode_auxiliary(iw));
} /* end change_shear_strain_subtract_mean() */


static double compute_central_symm (int n, double *dxlist)
{
    register int i,j;
    int k=0,m,*idx;
    double c=0, *d2, *dx, d2sum=0, tmp[DIMENSION], d2min;
    if (n <= 0) return(0);
    if (n == 1) return(1);
    m = MIN( n, central_symm_neighbormax );
    m = m / 2 * 2;
    MALLOC(compute_central_symm, d2,  n, double);
    MALLOC(compute_central_symm, dx,  DIMENSION*n, double);
    MALLOC(compute_central_symm, idx, n, int);
    for (i=0; i<n; i++) d2[i] = dxlist[(DIMENSION+1)*i+DIMENSION];
    qsort_glibc (n, d2, idx, USE_NEW_IDX);
    for (i=0; i<m; i++)
    {
        d2sum += dxlist[(DIMENSION+1)*idx[i]+DIMENSION];
        V3EQV( dxlist+(DIMENSION+1)*idx[i], dx+DIMENSION*i);
        idx[i] = 1;
    }
    for (i=0; i<m; i++)
        if (idx[i])
        {
            d2min = SINGLE_PRECISION_INFINITY;
            for (j=i+1; j<m; j++)
                if (idx[j])
                {
                    V3ADD( dx+DIMENSION*i, dx+DIMENSION*j, tmp );
                    d2[j] = V3LENGTH2( tmp );
                    if (d2[j] < d2min)
                    {
                        k = j;
                        d2min = d2[j];
                    }
                }                
            c += d2min;
            idx[k] = 0;
        }
    Free(idx);
    Free(dx);
    Free(d2);
    return(c / d2sum / 2);
} /* end compute_central_symm() */


void evaluate_central_symm (double **central_symm)
{
    register int i,j;
    int n;
    double tmp[DIMENSION], *dx=NULL;
    Neighborlist M[1];
    if (central_symm_neighbormax <= 0)
        central_symm_neighbormax = coordination_crystal;
    /* it must be an even number */
    central_symm_neighbormax = central_symm_neighbormax / 2 * 2;
    Neighborlist_create_nonpairwise_compressed_image
        (Config_Aapp_to_Alib, ct, &tp, N, M);
    REALLOC_ZERO (evaluate_central_symm, *central_symm, np, double);
    if (central_symm_neighbormax <= 0) return;
    n = Neighborlist_max_records(Config_Aapp_to_Alib, M);
    REALLOC (evaluate_central_symm, dx, (DIMENSION+1)*n, double);
    for (i=0; i<np; i++)
    {
        n = 0;
        for (j=M->idx[i]; j<M->idx[i+1]; j++)
        {
            V3SUB (&s[DIMENSION*M->list[j]], &s[DIMENSION*i], tmp);
            V3ImagE (tmp);
            V3M3LENGTH2 (tmp, H, dx + (DIMENSION+1)*n);
            n++;
        }
        if (n != coordination[i])
            pe("evaluate_central_symm: checksum failed!\n");
        else (*central_symm)[i] = compute_central_symm(n,dx);
    }
    Free(dx);
    Neighborlist_free_nonpairwise_image (M);
    return;
} /* end evaluate_central_symm() */


bool change_central_symm_neighbormax (int iw)
{
    char danswer[MAX_FILENAME_SIZE], *answer;
    int new_central_symm_neighbormax;
    xterm_get_focus(iw); clear_stdin_buffer();
    sprintf (danswer, "%d", central_symm_neighbormax);
    answer = readline_gets("\ncentral_symm_neighbormax",danswer);
    xterm_release_focus(iw);
    sscanf(answer, "%d", &new_central_symm_neighbormax);
    if (new_central_symm_neighbormax <= 0)
        new_central_symm_neighbormax = coordination_crystal;
    /* it must be an even number */
    new_central_symm_neighbormax = new_central_symm_neighbormax / 2 * 2;
    if ( new_central_symm_neighbormax == central_symm_neighbormax )
        return (FALSE);
    central_symm_neighbormax = new_central_symm_neighbormax;
    evaluate_central_symm (geo+GEO_CENTRAL_SYMM);
    n[iw].color_mode = COLOR_MODE_AUXILIARY;
    n[iw].auxiliary_idx = CONFIG_MAX_AUXILIARY + GEO_CENTRAL_SYMM;
    return(color_encode_auxiliary(iw));
} /* end change_central_symm_neighbormax() */


void geo_clear_has_evaluated_flags()
{
    int i;
    for (i=0; i<MAX_GEO_MEASURES; i++)
        geolist[i].has_evaluated = 0;
    return;
} /* end geo_clear_has_evaluated_flags() */


void geo_set_should_evaluate_flag (i)
{
    geolist[i].should_evaluate = 1;
    return;
} /* end geo_set_should_evaluate_flag() */


void evaluate_geo_measures()
{
    int i;
    for (i=0; i<MAX_GEO_MEASURES; i++)
        if ( geolist[i].should_evaluate &&
             (!geolist[i].has_evaluated) )
        {
            geolist[i].fun(geo + i);
            geolist[i].has_evaluated = 1;
        }
    return;
} /* end evaluate_geo_measures() */
