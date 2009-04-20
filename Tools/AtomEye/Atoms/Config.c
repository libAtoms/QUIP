/*************************************************/
/* Atoms: -llapack -lblas -lm                    */
/*        -lVecMat3 -lVecMat -lScalar -lIO       */
/*                                               */
/* Physical constants, macros, and configuration */
/* operators for atomistic simulations.          */
/*                                               */
/* Dec. 12, 1999  Ju Li <liju99@mit.edu>         */
/*************************************************/

#include "Atoms.h"

#include <netcdf.h>
#include <sys/stat.h>

#define netcdf_check(s) if ((retval = (s))) pe("NetCDF Error: %s %d %s\n", __FILE__, __LINE__, nc_strerror(retval));

#define N_A 6.0221479e23
#define ELEM_CHARGE 1.60217653e-19
#define MASS_CONVERT 1.0e7 / (N_A * ELEM_CHARGE)

#include <xyz_netcdf.h>


/****************************************************/
/* Minimal Specification of atomistic configuration */
/****************************************************/

#define BUFFERSIZE MAX(TERMSIZE,COMMENT_SIZE)
/* fprintf to a stream but adding comment escape sequence at front */
void comment (FILE *out, char *fmt, ...)
{
    va_list ap;
    static char *buffer=NULL, *start, *end, *p, *q;
    if ( buffer )
    { /* been allocated */
      beginning:
        for (p=start,q=fmt; (p<end) && (*q!=EOS); p++,q++)
            *p = *q;
        if (*q == EOS)
        {
            if (p<end) *p=EOS; else goto reallocation;
            va_start (ap, fmt);
            vfprintf (out, buffer, ap);
            return;
        }
        else
        {
          reallocation:
            while ((*q)!=EOS) q++;
            free(buffer);
            buffer = (char *) malloc(COMMENT_SIZE+(q-fmt));
            end = buffer + (COMMENT_SIZE+(q-fmt));
        }
    }
    else
    {
        buffer = (char *) malloc(BUFFERSIZE);
        end = buffer + BUFFERSIZE;
    }
    start = buffer + COMMENT_CHAR;
    for (p=buffer,q=COMMENT; p<start; p++,q++) *p = *q;
    goto beginning;
} /* end comment() */
#undef BUFFERSIZE

#ifdef _comment_TEST
int main (int argc, char *argv[])
{
    comment(stdout, "this is %s\n", "bad");
    comment(stdout, "this is %s\n", "worse");
    return (0);
}
#endif /* _comment_TEST */


void Config_analyze_H (Alib_Declare_Config, FILE *out)
{
    double tmp[3][3];
    M3MULTIPLY (ulength_IN_A, H, tmp);
    if (out) Mfprintf(out, "H = %3M|| %.6lf %.6lf %.6lf | = "
                      "%3M|| %.6lf %.6lf %.6lf | A\n", H, tmp);
    return;
} /* end Config_analyze_H() */


/* Volume, number density and mass density  */
/* printouts; return mass density [g/cm^3]. */
double Config_analyze_density (Alib_Declare_Config, FILE *out)
{
    register int i;
    register double volume, totalmass;
    volume = M3DETERMINANT(H);
    if (volume < 0) volume = -volume;
    totalmass = 0;
    for (i=(*np); i--;) totalmass += (*mass)[i];
    if (out)
    {
        fprintf (out, "supercell volume = %.15g = %.15g A^3\n",
                 volume, volume * uvolume_IN_A3);
        fprintf (out, "avg. atomic volume = %.15g = %.15g A^3\n",
                 volume / (*np), volume * uvolume_IN_A3 / (*np));
        fprintf (out, "atomic number density = %.15g = %.15g A^-3\n",
                 (*np) / volume, (*np) / (volume * uvolume_IN_A3));
        fprintf (out, "avg. mass density = %.15g g/cm^3\n",
                 totalmass * umass_IN_G / volume / uvolume_IN_CM3);
    }
    return (totalmass * umass_IN_G / volume / uvolume_IN_CM3);
} /* end Config_analyze_density() */


/* Chemical species and isotope analysis printouts; returns statistics */
ConfigSpecies *Config_analyze_species (Alib_Declare_Config, FILE *out)
{
    static ConfigSpecies root = {0,0,NULL};
    ConfigSpecies *p;
    register int i,j;
    double totalweight;
    root.counter = 0;
    for (i=0; i<(*np); i++)
    {
        for (j=0,p=&root; j<root.counter; p=p->next,j++)
            if ( ( (*mass)[i] == (*mass)[p->next->first] ) &&
                 ( !strcmp(SYMBOL(i),SYMBOL(p->next->first)) ) )
            {
                p->next->counter++;
                break;
            }
        if (j == root.counter)
        {
            if (p->next == NULL)
                p->next = (ConfigSpecies *)
                    malloc(sizeof(ConfigSpecies));
            p->next->first = i;
            p->next->counter = 1;
            p->next->next = NULL;
            root.counter++;
        }
    }
    if ( out )
    {
        for (totalweight=0.,j=0,p=&root; j<root.counter; p=p->next,j++)
            totalweight += (*mass)[p->next->first] * p->next->counter;
        fprintf (out, "-------------- Species Summary --------------\n");
        fprintf (out, "Type  Mass[amu]   Count  Abundance  Wt. Pct.\n");
        for (j=0,p=&root; j<root.counter; p=p->next,j++)
            fprintf (out, " %2s    %7.3f%9d   %6.2f%%   %6.2f%%\n",
                     SYMBOL(p->next->first),
                     (*mass)[p->next->first] * umass_IN_AMU,
                     p->next->counter, p->next->counter*100./(*np),
                     100*(*mass)[p->next->first]*p->next->counter/totalweight);
        fprintf (out, "---------------------------------------------\n");
    }
    return (&root);
} /* end Config_analyze_species() */


/* find the nearest neighbor to atom w and return their bond length */
int Xtal_analyze_nearest_bond
(Alib_Declare_Config, int w, double *nearest_bond_IN_ulength, FILE *out)
{
    int i,j,k,n,nearest=-1;
    double ds[3],dx[3];
    *nearest_bond_IN_ulength = DOUBLE_PRECISION_INFINITY;
    for (i=-1; i<=1; i++)
        for (j=-1; j<=1; j++)
            for (k=-1; k<=1; k++)
                for (n=0; n<(*np); n++)
                {
                    if ((i==0)&&(j==0)&&(k==0)&&(n==w)) continue;
                    ds[0] = i + (*s)[DIMENSION*n]   - (*s)[DIMENSION*w];
                    ds[1] = j + (*s)[DIMENSION*n+1] - (*s)[DIMENSION*w+1];
                    ds[2] = k + (*s)[DIMENSION*n+2] - (*s)[DIMENSION*w+2];
                    V3mM3 (ds, H, dx);
                    ds[0] = V3LENGTH2(dx);
                    if (ds[0] < *nearest_bond_IN_ulength)
                    {
                        *nearest_bond_IN_ulength = ds[0];
                        nearest = n;
                    }
                }
    *nearest_bond_IN_ulength = sqrt(*nearest_bond_IN_ulength);
    if (out)
        fprintf (out, "For the %s atom, its nearest neighbor is a %2s atom\n"
                 "bonding distance = %.15g = %.15g A.\n", word_for_order(w),
                 SYMBOL(nearest), *nearest_bond_IN_ulength,
                 *nearest_bond_IN_ulength * ulength_IN_A);
    return (nearest);
} /* Xtal_analyze_nearest_bond() */


/* Save configuration in line-based CFG ASCII file */
void Config_save (Alib_Declare_Config, bool PBC, char *fname)
{
    register int i,j;
    V3 mys;
    /* safely open a file for (over)writing */
    FILE *out = WOpen(fname);
    /* this must be the first line */
    fprintf (out, "Number of particles = %d\n", (*np));
    comment (out, "(required) this must be the first line\n\n");
    fprintf (out, "A = 1.0 Angstrom (basic length-scale)\n");
    comment (out, "(optional) basic length-scale: default A = 1.0 "
             "[Angstrom]\n\n");
    for (i=0; i<DIMENSION; i++)
    {
        for (j=0; j<DIMENSION; j++)
            fprintf ( out, "H0(%1d,%1d) = %.15g A\n", i+1, j+1,
                      H[i][j] * ulength_IN_A );
        comment(out, "(required) this is the supercell's %s edge, in A\n\n",
                word_for_order(i+1));
    }

    for (i=0; i<DIMENSION; i++)
        for (j=0; j<DIMENSION; j++)
            fprintf (out, "Transform(%1d,%1d) = %.15g\n",
                     i+1, j+1, DOUBLE(KroneckerDelta(i,j)));
    comment (out, "(optional) apply additional transformation on H0:  "
             "H = H0 * Transform;\n");
    comment (out, "default = Identity matrix.\n\n");

    for (i=0; i<DIMENSION; i++)
        for (j=i; j<DIMENSION; j++)
            fprintf (out, "eta(%1d,%1d) = %.15g\n", i+1, j+1, 0.);
    comment (out, "(optional) apply additional Lagrangian strain on H0:\n");
    comment (out, "H = H0 * sqrt(Identity_matrix + 2 * eta);\n");
    comment (out, "default = zero matrix.\n\n");
    
    comment (out, "ENSUING ARE THE ATOMS, EACH ATOM DESCRIBED BY A ROW\n");
    comment (out, "1st entry is atomic mass in a.m.u.\n");
    comment (out, "2nd entry is the chemical symbol (max %d chars)\n\n",
             SYMBOL_CHAR);
    for (i=0; i<DIMENSION; i++)
        comment (out, "%s entry is reduced coordinate s%1d (dimensionless)\n",
                 word_for_order(3+i), i+1);
    comment (out, "real coordinates x = s * H,  "
             "x, s are 1x%1d row vectors\n\n", DIMENSION);
    for (i=0; i<DIMENSION; i++)
        comment (out, "%s entry is d(s%1d)/dt in basic rate-scale R\n",
                 word_for_order(3+DIMENSION+i), i+1);
    fprintf (out, "R = 1.0 [ns^-1]\n");
    comment (out, "(optional) basic rate-scale: default R = 1.0 [ns^-1]\n\n");
    for (i=0; i<(*np); i++)
    {
        fprintf (out, "%.15g", (*mass)[i] * umass_IN_AMU);
        fprintf (out, " %2s",  SYMBOL(i));
        if (PBC) V3TRIM ( &(*s)[DIMENSION*i], mys );
        else V3EQV ( &(*s)[DIMENSION*i], mys );
        for (j=0; j<DIMENSION; j++)
            zfprintf (out, " %.15g", mys[j] );
        for (j=0; j<DIMENSION; j++)
            zfprintf (out, " %.15g", (*s1)[DIMENSION*i+j] / utime_IN_NS);
        fcr(out);
    }
    fcr(out);
    reset_scratch();
    fprintf(fp_scratch, "Analysis of this configuration:\n\n");
    Config_analyze_density (Config_Alib_to_Alib, fp_scratch);
    fcr(fp_scratch);
    Config_analyze_motion (Config_Alib_to_Alib, fp_scratch);
    fcr(fp_scratch);
    Config_analyze_species (Config_Alib_to_Alib, fp_scratch);
    dump_scratched_to (out, COMMENT);
    kill_scratch();
    Zclose (out, fname);
    return;
} /* end Config_save() */


/************************************************************/
/* CONFIG is more flexible than Config:                     */
/* 1. Atomic mass and chemical symbol (1st and 2nd entry in */
/*    Config) can be omitted if equal to the previous atom. */
/* 2. s1,s2,s3 (3rd,4th,5th entry in Config) can be stored  */
/*    in arbitrary precisions.                              */
/* 3. d(s1)/dt,d(s2)/dt,d(s3)/dt (6,7,8th entry in Config)  */
/*    may be omitted or stored in arbitrary precisions.     */
/* 4. New properties of arbitrary precision can be added.   */
/************************************************************/
void vCONFIG_save
(Alib_Declare_Config, bool PBC, char *fname, char *s_formats,
 char *velocity_formats, int num_auxiliary, va_list ap)
{
    register int i,j;
    int entry_count = 0;
    double oldmass;
    char oldsymbol[SYMBOL_SIZE];
    V3 mys;
    struct auxiliary
    {
        char *description;
        char *format;
        char *pointer;
        int bytes_separation;
    } aux[CONFIG_MAX_AUXILIARY];
    SimpleStatistics ss;
    /* safely open a file for (over)writing */
    FILE *out = WOpen(fname);

    if (num_auxiliary > CONFIG_MAX_AUXILIARY)
        pe ("vCONFIG_save: num_auxiliary=%d > CONFIG_MAX_AUXILIARY=%d.\n",
            num_auxiliary, CONFIG_MAX_AUXILIARY);

    for (i=0; i<num_auxiliary; i++)
    {
        aux[i].description = va_arg(ap, char *); /* "s11 [GPa]" */
        aux[i].format = va_arg(ap, char *); /* "%.15g" */
        aux[i].pointer = va_arg(ap, char *); /* (char *)s */
        aux[i].bytes_separation = va_arg(ap, int); /* sizeof(double) * 6 */
    }
    va_end (ap);
     
    /* this must be the first line */
    fprintf (out, "Number of particles = %d\n", (*np));
    comment (out, "(required) this must be the first line\n\n");
    fprintf (out, "A = 1.0 Angstrom (basic length-scale)\n");
    comment (out, "(optional) basic length-scale: default A = 1.0 "
             "[Angstrom]\n\n");
    for (i=0; i<DIMENSION; i++)
    {
        for (j=0; j<DIMENSION; j++)
            fprintf ( out, "H0(%1d,%1d) = %.15g A\n", i+1, j+1,
                      H[i][j] * ulength_IN_A );
        comment(out, "(required) this is the supercell's %s edge, in A\n\n",
                word_for_order(i+1));
    }

    for (i=0; i<DIMENSION; i++)
        for (j=0; j<DIMENSION; j++)
            fprintf (out, "Transform(%1d,%1d) = %.15g\n",
                     i+1, j+1, DOUBLE(KroneckerDelta(i,j)));
    comment (out, "(optional) apply additional transformation on H0:  "
             "H = H0 * Transform;\n");
    comment (out, "default = Identity matrix.\n\n");

    for (i=0; i<DIMENSION; i++)
        for (j=i; j<DIMENSION; j++)
            fprintf (out, "eta(%1d,%1d) = %.15g\n", i+1, j+1, 0.);
    comment (out, "(optional) apply additional Lagrangian strain on H0:\n");
    comment (out, "H = H0 * sqrt(Identity_matrix + 2 * eta);\n");
    comment (out, "default = zero matrix.\n\n");
    
    comment (out, "Each atom is described by a row:\n\n");
    for (i=0; i<DIMENSION; i++)
        comment (out, "%s entry is reduced coordinate s%1d (dimensionless)\n",
                 word_for_order(entry_count+1+i), i+1);
    entry_count += DIMENSION;
    comment (out, "real coordinates x = s * H,  "
             "x, s are 1x%1d row vectors\n\n", DIMENSION);

    if (velocity_formats==NULL)
    { /* do not save velocity */
        fputs(".NO_VELOCITY.\n", out);
        comment (out, "Atom velocities are deemed irrelevant for\n");
        comment (out, "this configuration so they are not stored.\n\n");
    }
    else
    {
        for (i=0; i<DIMENSION; i++)
            comment (out, "%s entry is d(s%1d)/dt in basic rate-scale R\n",
                     word_for_order(entry_count+1+i), i+1);
        entry_count += DIMENSION;
        fprintf (out, "R = 1.0 [ns^-1]\n");
        comment (out, "(optional) basic rate-scale: "
                 "default R = 1.0 [ns^-1]\n\n");
    }

    entry_count += num_auxiliary;
    fprintf (out, "entry_count = %d\n\n", entry_count);

    for (i=0; i<num_auxiliary; i++)
    {
        CalculateSimpleStatistics
            (*np, aux[i].pointer, aux[i].bytes_separation, IOVAL_DOUBLE, &ss);
        fprintf (out, "auxiliary[%d] = %s\n", i, aux[i].description);
        comment (out, "%s entry: [%g (atom %d), %g (atom %d)],\n",
                 word_for_order(entry_count-num_auxiliary+1+i),
                 ss.min, ss.idx_min, ss.max, ss.idx_max);
        comment (out, "average = %g, std.dev. = %g\n\n",
                 ss.average, ss.standard_deviation);
    }

    comment (out, "These properties are piece-wise uniform:\n");
    fprintf (out, "%.15g\n", (*mass)[0] * umass_IN_AMU);
    oldmass = (*mass)[0];
    comment (out, "(required) atomic mass in a.m.u.\n");
    fprintf (out, "%2s\n",  SYMBOL(0));
    strcpy (oldsymbol, SYMBOL(0));
    comment (out, "(required) chemical symbol (max %d chars)\n",
             SYMBOL_CHAR);

    for (i=0; i<(*np); i++)
    { /* both CONFIG and Config respect index */
        if ( (*mass)[i] != oldmass )
        {
            oldmass = (*mass)[i];
            fprintf (out, "%.15g\n", (*mass)[i] * umass_IN_AMU);
        }
        if ( DIFFERENT_SYMBOL(oldsymbol,SYMBOL(i)) )
        {
            safe_symbol(SYMBOL(i),oldsymbol);
            fprintf (out, "%2s\n", SYMBOL(i));
        }
        if (PBC) V3TRIM ( &(*s)[DIMENSION*i], mys );
        else V3EQV ( &(*s)[DIMENSION*i], mys );
        zfprintf (out, s_formats, V3E(mys));
        V3DIV ( &(*s1)[DIMENSION*i], utime_IN_NS, mys );
        if (velocity_formats!=NULL)
            zfprintf (out, velocity_formats, V3E(mys));
        for (j=0; j<num_auxiliary; j++)
            zfprintf (out, aux[j].format,
                      IOVAL(aux[j].pointer+i*aux[j].bytes_separation,
                            IOVAL_DOUBLE));
        fcr(out);
    }
    fcr(out);
    reset_scratch();
    fprintf(fp_scratch, "Analysis of this configuration:\n\n");
    Config_analyze_density (Config_Alib_to_Alib, fp_scratch);
    fcr(fp_scratch);
    Config_analyze_motion (Config_Alib_to_Alib, fp_scratch);
    fcr(fp_scratch);
    Config_analyze_species (Config_Alib_to_Alib, fp_scratch);
    dump_scratched_to (out, COMMENT);
    kill_scratch();
    Zclose (out, fname);
    return;
} /* end vCONFIG_save() */


void CONFIG_save
(Alib_Declare_Config, bool PBC, char *fname, char *s_formats,
 char *velocity_formats, int num_auxiliary, ...)
{
    va_list ap;
    va_start (ap, num_auxiliary);
    vCONFIG_save (Config_Alib_to_Alib, PBC, fname, s_formats,
                  velocity_formats, num_auxiliary, ap);
    return;
} /* end CONFIG_save() */


/* PBC enforced */
void CONFIG_SAVE
(Alib_Declare_Config, char *fname, char *s_formats,
 char *velocity_formats, int num_auxiliary, ...)
{
    va_list ap;
    va_start (ap, num_auxiliary);
    vCONFIG_SAVE (Config_Alib_to_Alib, fname, s_formats,
                  velocity_formats, num_auxiliary, ap);
    return;
} /* end CONFIG_SAVE() */


/* no velocity */
void CONFIG_NV_save
(Alib_Declare_Config, bool PBC, char *fname, char *s_formats,
 int num_auxiliary, ...)
{
    va_list ap;
    va_start (ap, num_auxiliary);
    vCONFIG_NV_save(Config_Alib_to_Alib, PBC, fname,
                    s_formats, num_auxiliary, ap);
    return;
} /* end CONFIG_NV_save() */


/* no velocity & PBC enforced */
void CONFIG_NV_SAVE
(Alib_Declare_Config, char *fname, char *s_formats, int num_auxiliary, ...)
{
    va_list ap;
    va_start (ap, num_auxiliary);
    vCONFIG_NV_SAVE(Config_Alib_to_Alib, fname, s_formats, num_auxiliary, ap);
    return;
} /* end CONFIG_NV_SAVE() */


/* no auxiliary */
void CONFig_save(Alib_Declare_Config, bool PBC, char *fname,
                 char *s_formats, char *velocity_formats)
{
    va_list ap;
    vCONFIG_save (Config_Alib_to_Alib, PBC, fname, s_formats,
                  velocity_formats, 0, ap);
    return;
} /* end CONFig_save() */


/***********************************************************************/
/* Save configuration in Protein Data Bank format so it can be viewed  */
/* in external viewers: you can choose different origin for the saved  */
/* PDB configuration by nonzero new_origin_s0,s1,s2, and PBC specifies */
/* whether s[] should then be fold into [0,1)x[0,1)x[0,1). More info:  */
/* http://rcsb.org/pdb/docs/format/pdbguide2.2/guide2.2_frame.html.    */
/***********************************************************************/
void Config_save_as_pdb
(Alib_Declare_Config, double new_origin_s0, double new_origin_s1,
 double new_origin_s2, bool PBC, char *fname)
{
    register int i;
    double mys[3], myx[3], HH[3][3];
    Crystallographic X;
    FILE *out = wOpen(fname);
    fprintf (out, "HEADER    libAtoms:Config_save_as_pdb; %d atoms\n", (*np));
    /* Save H[][] in standard orientation: will be reloaded as such then */
    H_to_Crystallographic(H,X);
    Crystallographic_to_H(X,HH);
    fprintf (out, "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %11s%4d\n",
             X.a*ulength_IN_A, X.b*ulength_IN_A, X.c*ulength_IN_A,
             X.alpha, X.beta, X.gamma, "P 1        ", 1);
    for (i=0; i<(*np); i++)
    {
        mys[0] = (*s)[DIMENSION*i]   - new_origin_s0;
        mys[1] = (*s)[DIMENSION*i+1] - new_origin_s1;
        mys[2] = (*s)[DIMENSION*i+2] - new_origin_s2;
        if (PBC) V3TRIM (mys, mys);
        V3mM3 ( mys, HH, myx );
        fprintf (out, "ATOM  %5d %2s%16s%8.3f%8.3f%8.3f\n",
                 i+1, SYMBOL(i), "           1    ", myx[0]*ulength_IN_A,
                 myx[1]*ulength_IN_A, myx[2]*ulength_IN_A);
    }
    fclose(out);
    return;
} /* end Config_save_as_pdb() */


/* JL -> .pdb converter */
#ifdef _c2p
int main (int argc, char *argv[])
{
    Aapp_Define_Config;
    char read[TERMSIZE], write[2*TERMSIZE];
    char *r, *p, *q;
    if (argc == 1)
    { /* no input argument */
        printf ("Convert config file of JL's format: ");
        for (r=read; r<read+TERMSIZE-1; r++)
            if ((*r=(char)getc(stdin))=='\n') break;
        *r = (char)0;
    }
    else snprintf (read, TERMSIZE, "%s", argv[1]);
    for (r=read; (*r==' ')||(*r=='\t'); r++);
    for (q=r; *q!=(char)0; q++)
        if ((*q==' ')||(*q=='\t'))
        {
            *q=(char)0;
            break;
        }
    Config_load (r, NULL, Config_Aapp_to_Alib);
    sprintf (write, "%s.pdb", r);
    if (argc <= 2)
    { /* interactive mode */
        printf("Write to (default=\"%s\"): ", write);
        q = write + strlen(write);
        for (p=q; p<write+2*TERMSIZE-5; p++)
            if ((*p=(char)getc(stdin))=='\n') break;
        if (p==q) q = write;
        *p = (char)0;
        /* add file suffix if without */
        if (strcasecmp(p-4,".pdb")) sprintf(p,".pdb");
    }
    else q = argv[2];
    for (; (*q==' ')||(*q=='\t'); q++);
    Config_save_as_PDB (Config_Aapp_to_Alib, q);
    printf ("\"%s\" converted -> \"%s\".\n", r, q);
    return (0);
}
#endif /* _c2p */


/* .pdb -> JL converter */
#ifdef _p2c
int main (int argc, char *argv[])
{
    Aapp_Define_Config;
    char read[TERMSIZE], write[2*TERMSIZE];
    char *r, *p, *q;
    if (argc == 1)
    { /* no input argument */
        printf ("Convert config file of PDB format: ");
        for (r=read; r<read+TERMSIZE-1; r++)
            if ((*r=(char)getc(stdin))=='\n') break;
        *r = (char)0;
    }
    else snprintf (read, TERMSIZE, "%s", argv[1]);
    for (r=read; (*r==' ')||(*r=='\t'); r++);
    for (q=r; *q!=(char)0; q++)
        if ((*q==' ')||(*q=='\t'))
        {
            *q=(char)0;
            break;
        }
    Config_Load (r, NULL, Config_Aapp_to_Alib);
    sprintf (write, "%s.cfg", r);
    if (argc <= 2)
    { /* interactive mode */
        printf("Write to (default=\"%s\"): ", write);
        q = write + strlen(write);
        for (p=q; p<write+2*TERMSIZE-5; p++)
            if ((*p=(char)getc(stdin))=='\n') break;
        if (p==q) q = write;
        *p = (char)0;
        /* add file suffix if without */
        if (strcasecmp(p-4,".cfg")) sprintf(p,".cfg");
    }
    else q = argv[2];
    for (; (*q==' ')||(*q=='\t'); q++);
    Config_save (Config_Aapp_to_Alib, TRUE, q);
    printf ("\"%s\" converted -> \"%s\".\n", r, q);
    return (0);
}
#endif /* _p2c */


/* JL -> .xyz converter */
#ifdef _cfg2xyz
int main (int argc, char *argv[])
{
    Aapp_Define_Config;
    char read[TERMSIZE], write[2*TERMSIZE];
    char *r, *p, *q;
    int i;
    FILE *fp;
    V3 myx;
    if (argc == 1)
    { /* no input argument */
        printf ("Convert config file of JL's format: ");
        for (r=read; r<read+TERMSIZE-1; r++)
            if ((*r=(char)getc(stdin))=='\n') break;
        *r = (char)0;
    }
    else snprintf (read, TERMSIZE, "%s", argv[1]);
    for (r=read; (*r==' ')||(*r=='\t'); r++);
    for (q=r; *q!=(char)0; q++)
        if ((*q==' ')||(*q=='\t'))
        {
            *q=(char)0;
            break;
        }
    Config_load (r, NULL, Config_Aapp_to_Alib);
    sprintf (write, "%s.xyz", r);
    if (argc <= 2)
    { /* interactive mode */
        printf("Write to (default=\"%s\"): ", write);
        q = write + strlen(write);
        for (p=q; p<write+2*TERMSIZE-5; p++)
            if ((*p=(char)getc(stdin))=='\n') break;
        if (p==q) q = write;
        *p = (char)0;
        /* add file suffix if without */
        if (strcasecmp(p-4,".xyz")) sprintf(p,".xyz");
    }
    else q = argv[2];
    for (; (*q==' ')||(*q=='\t'); q++);

    fp = WOpen(q);
    fprintf (fp, "%d\n", np);
    fprintf (fp, "Supercell 1st edge=[%.15g %.15g %.15g], "
             "2nd edge=[%.15g %.15g %.15g], 3rd edge=[%.15g %.15g %.15g] A\n",
             M3e(H));
    for (i=0; i<np; i++)
    {
        /* V3TriM( &s[DIMENSION*i] ); */
        /* XMakemol prefers center of mass at [0 0 0] */
        V3ImagE( &s[DIMENSION*i] );
        V3mM3( &s[DIMENSION*i], H, myx );
        fprintf (fp, "%2s %.15g %.15g %.15g\n", SYM(i),V3E(myx));
    }
    Zclose(fp, q);
    printf ("\"%s\" converted -> \"%s\".\n", r, q);
    return (0);
}
#endif /* _cfg2xyz */


/* Save the current config, including allocation state, to a stack member */
void Config_push (Alib_Declare_Config, ConfigStack *cs)
{
    M3EQV (H, cs->h);
    cs->NP     = (*np);
    cs->symBOL = (*symbol);
    cs->MASS   = (*mass);
    cs->S      = (*s);
    cs->S1     = (*s1);
    cs->ULENGTH_in_A = (ulength_IN_A);
    cs->UMASS_in_AMU = (umass_IN_AMU);
    cs->UTIME_in_NS  = (utime_IN_NS);
    return;
} /* end Config_push() */


/* Restore the config, including allocation state, from a stack member */
void Config_pop (ConfigStack *cs, Alib_Declare_Config)
{
    M3EQV (cs->h, H);
    (*np)      = cs->NP;
    (*symbol)  = cs->symBOL;
    (*mass)    = cs->MASS;
    (*s)       = cs->S;
    (*s1)      = cs->S1;
    return;
} /* end Config_pop() */


/* Allocate new memory, then clone the configuration stored in "cs" */
void Config_clone (ConfigStack *cs, Alib_Declare_Config)
{
    M3EQV (cs->h, H);
    (*np) = cs->NP;
    Config_alloc (Config_Alib_to_Alib);
    BEQV ( SYMBOL_SIZE*(*np), cs->symBOL,  (*symbol) );
    VEQV ( (*np),             cs->MASS,    (*mass)   );
    VEQV ( DIMENSION*(*np),   cs->S,       (*s)      );
    VEQV ( DIMENSION*(*np),   cs->S1,      (*s1)     );
    return;
} /* end Config_clone() */


/* Save the current config, including allocation state, to a stack member */
void CONFIG_push (Alib_Declare_CONFIG, ConfigStack *cs)
{
    M3EQV(H, cs->h);
    cs->NP     = (*np);
    cs->symBOL = (*symbol);
    cs->MASS   = (*mass);
    cs->S      = (*s);
    cs->S1     = (*s1);
    cs->ULENGTH_in_A = (*ulength_IN_A);
    cs->UMASS_in_AMU = (*umass_IN_AMU);
    cs->UTIME_in_NS  = (*utime_IN_NS);
    return;
} /* end CONFIG_push() */


/* Restore the config, including allocation state, from a stack member */
void CONFIG_pop (ConfigStack *cs, Alib_Declare_CONFIG)
{
    M3EQV(cs->h, H);
    (*np)      = cs->NP;
    (*symbol)  = cs->symBOL;
    (*mass)    = cs->MASS;
    (*s)       = cs->S;
    (*s1)      = cs->S1;
    (*ulength_IN_A) = cs->ULENGTH_in_A;
    (*umass_IN_AMU) = cs->UMASS_in_AMU;
    (*utime_IN_NS)  = cs->UTIME_in_NS;
    return;
} /* end CONFIG_pop() */


/* Free memory pointers stored in a stack member */
void CONFIG_erase (ConfigStack *cs)
{
    Free (cs->symBOL);
    Free (cs->MASS);
    Free (cs->S);
    Free (cs->S1);
    return;
} /* end CONFIG_erase() */


/* Allocate new memory pointers symbol, mass, s, s1 for current np */
void Config_alloc (Alib_Declare_Config)
{
    MALLOC(Config_alloc, *symbol, SYMBOL_SIZE*(*np), char);
    MALLOC(Config_alloc, *mass,   *np,               double);
    MALLOC(Config_alloc, *s,      DIMENSION*(*np),   double);
    MALLOC(Config_alloc, *s1,     DIMENSION*(*np),   double);
    return;
} /* end Config_alloc() */


/* NULL-safe release of memory pointers symbol, mass, s, s1 */
void Config_free (Alib_Declare_Config)
{
    Free (*symbol);
    Free (*mass);
    Free (*s);
    Free (*s1);
    return;
} /* end Config_free() */


/* NULL-safe (re)allocation of symbol, mass, s, s1 for current np */
void Config_realloc (Alib_Declare_Config)
{
    REALLOC(Config_realloc, *symbol, SYMBOL_SIZE*(*np), char);
    REALLOC(Config_realloc, *mass,   *np,               double);
    REALLOC(Config_realloc, *s,      DIMENSION*(*np),   double);
    REALLOC(Config_realloc, *s1,     DIMENSION*(*np),   double);
    return;
} /* end Config_realloc() */


/** configuration transformations **/

/* particle-conservative transformations */

/* s0:=s0+s0_trans, s1:=s1+s1_trans, s2:=s2+s2_trans for all atoms */
void Config_translate
(double s0_trans, double s1_trans, double s2_trans, Alib_Declare_Config)
{
    register int i;
    for (i=(*np); i--;)
    {
        (*s)[DIMENSION*i]   += s0_trans;
        (*s)[DIMENSION*i+1] += s1_trans;
        (*s)[DIMENSION*i+2] += s2_trans;
    }
    return;
} /* end Config_translate() */


#ifdef _trans
int main (int argc, char *argv[])
{
    Aapp_Define_Config;
    double s0_trans, s1_trans, s2_trans;
    if (argc != 6)
    {
        printf ("Purpose: translate configuration by [ds0 ds1 ds2].\n");
        printf ("Usage: %s in_file ds0 ds1 ds2 out_file\n", argv[0]);
        return (1);
    }
    printf ("Loading \"%s\"...\n\n", argv[1]);
    Config_Load (argv[1], stdout, Config_Aapp_to_Alib);
    s0_trans = atof(argv[2]);
    s1_trans = atof(argv[3]);
    s2_trans = atof(argv[4]);
    Config_translate (s0_trans, s1_trans, s2_trans, Config_Aapp_to_Alib);
    Config_save (Config_Aapp_to_Alib, TRUE, argv[5]);
    printf ("\"%s\" with [%g %g %g] boost -> \"%s\".\n",
            argv[1], s0_trans, s1_trans, s2_trans, argv[5]);
    return (0);
}
#endif /* _trans */


/* Set the center of mass to S0,S1,S2; return the amount of shift in s[] */
double *Config_set_cm (double S0, double S1, double S2, Alib_Declare_Config)
{
    register int i;
    register double totalmass = 0;
    static double s_shift[3];
    V3ZERO (s_shift);
    for (i=(*np); i--;)
    {
        totalmass += (*mass)[i];
        V3ADDmuL( (*mass)[i], (*s)+DIMENSION*i, s_shift);
    }
    V3MuL (-1./totalmass, s_shift);
    s_shift[0] += S0;
    s_shift[1] += S1;
    s_shift[2] += S2;
    Config_Translate (s_shift, Config_Alib_to_Alib);
    return (s_shift);
} /* end Config_set_cm() */


/* achieve x' = x*R by H' = H*R  */
void Config_rotate_via_H (Alib_Declare_Config, double R[3][3])
{
    double TMP[3][3];
    M3MUL (H, R, TMP);
    M3EQV (TMP, H);
    return;
} /* end Config_rotate_via_H() */


/* Achieve x':=x*R by s':=s*H*R*H^-1 (you have to know what you are doing) */
void Config_rotate_via_s
(Alib_Declare_Config, double S0, double S1, double S2, double R[3][3])
{
    register int i;
    V3 s0;
    double HI[3][3], TMP[3][3], RR[3][3];
    V3ASSIGN (S0,S1,S2, s0);
    M3INV (H, HI, TMP[0][0]);
    M3MUL (H, R, TMP);
    M3MUL (TMP, HI, RR);
    for (i=(*np); i--;)
    {
        V3SUB( &((*s)[DIMENSION*i]), s0, TMP[1] );
        V3mM3( TMP[1], RR, TMP[0] );
        V3ADD( s0, TMP[0], &((*s)[DIMENSION*i]) );
        V3MM3( &((*s1)[DIMENSION*i]), RR, TMP[0] );
    }
    return;
} /* end Config_rotate_via_s() */


/* s,H -> snew,Hnew: (snew-s0) * Hnew := (s-s0) * H */
void Config_transplant
(Alib_Declare_Config, double S0, double S1, double S2, double Hnew[3][3])
{
    register int i;
    V3 s0;
    M3 HnewI, TMP, RR;
    V3ASSIGN (S0,S1,S2, s0);
    M3INV (Hnew, HnewI, TMP[0][0]);
    M3MUL (H, HnewI, RR);
    for (i=(*np); i--;)
    {
        V3SUB( &((*s)[DIMENSION*i]), s0, TMP[1] );
        V3mM3( TMP[1], RR, TMP[0] );
        V3ADD( s0, TMP[0], &((*s)[DIMENSION*i]) );
        V3MM3( &((*s1)[DIMENSION*i]), RR, TMP[0] );
    }
    M3EQV (Hnew, H);
    return;
} /* end Config_transplant() */


#define DIRECTMAX 3.
/* Fold all atoms' s0 into [S0,S0+1), s1 into [S1,S1+1), s2 into [S2,S2+1). */
void Config_fold_into_pbc (Alib_Declare_Config, double S0,double S1,double S2)
{
    register int i;
    register double T0, T1, T2;
    double DS0, DS1, DS2;
    if ( (fabs(S0) > DIRECTMAX) ||
         (fabs(S1) > DIRECTMAX) ||
         (fabs(S2) > DIRECTMAX) )
    {
        DS0 = floor(S0);
        DS1 = floor(S1);
        DS2 = floor(S2);
        S0 -= DS0;
        S1 -= DS1;
        S2 -= DS2;
    }
    else DS0 = DS1 = DS2 = 0;
    T0 = S0 + 1;
    T1 = S1 + 1;
    T2 = S2 + 1;
    for (i=(*np); i--;)
    {
        while ((*s)[DIMENSION*i]   <  S0) (*s)[DIMENSION*i]++;
        while ((*s)[DIMENSION*i]   >= T0) (*s)[DIMENSION*i]--;
        while ((*s)[DIMENSION*i+1] <  S1) (*s)[DIMENSION*i+1]++;
        while ((*s)[DIMENSION*i+1] >= T1) (*s)[DIMENSION*i+1]--;
        while ((*s)[DIMENSION*i+2] <  S2) (*s)[DIMENSION*i+2]++;
        while ((*s)[DIMENSION*i+2] >= T2) (*s)[DIMENSION*i+2]--;
    }
    if ( (DS0 != 0) || (DS1 != 0) || (DS2 != 0) )
        Config_translate (DS0, DS1, DS2, Config_Alib_to_Alib);
    return;
} /* end Config_fold_into_pbc() */
#undef DIRECTMAX


/* Swap atoms i and j in index */
void Config_swapatom (Alib_Declare_Config, int i, int j)
{
    int k;
    char c;
    register double tmp;
    for (k=0; k<SYMBOL_SIZE; k++)
        SWAP ( SYMBOL(i)[k], SYMBOL(j)[k], c );
    SWAP ( (*mass)[i], (*mass)[j], tmp );
    V3SWAP ( &(*s)[DIMENSION*i],  &(*s)[DIMENSION*j],  tmp );
    V3SWAP ( &(*s1)[DIMENSION*i], &(*s1)[DIMENSION*j], tmp );
    return;
} /* end Config_swapatom() */


/* Open up a rift in "direction" by expanding that cell */
/* edge by "ratio" while shrinking the s[] by as much.  */
void Config_openrift (Alib_Declare_Config, int direction, double ratio)
{
    register int i;
    V3MuL ( ratio, H[direction] );
    for (i=(*np); i--;) (*s)[DIMENSION*i+direction] /= ratio;
    return;
} /* end Config_openrift() */


/* Same as Config_openrift() but with s[]'s staying at the center */
void Config_Openrift (Alib_Declare_Config, int direction, double ratio)
{
    register int i;
    V3MuL ( ratio, H[direction] );
    for (i=(*np); i--;)
        (*s)[DIMENSION*i+direction] =
            0.5 + ( (*s)[DIMENSION*i+direction] - 0.5 ) / ratio;
    return;
} /* end Config_Openrift() */


/* Perturb each atom in each H(i,:) direction by uniformly */
/* distributed random (-s_amplitude[i]/2, s_amplitude[i]/2). */
void Config_perturb_s (Alib_Declare_Config, double s_amplitude[3])
{
    register int i;
    V3 ds;
    for (i=(*np); i--;)
    {
        V3FRANDOM(ds);
        V3ProD(ds, s_amplitude);
        V3AdD (ds, &((*s)[DIMENSION*i]));
    }
    return;
} /* end Config_perturb_s() */


/* Same as Config_perturb_s() except center of mass does not move. */
void Config_Perturb_s (Alib_Declare_Config, double s_amplitude[3])
{
    register int i;
    register double cc, totalmass;
    V3 ds, dcm;
    M3 HI;
    M3INV (H, HI, cc);
    V3ZERO(dcm);
    totalmass = 0;
    for (i=(*np); i--;)
    {
        V3FRANDOM(ds);
        V3ProD(ds, s_amplitude);
        V3AdD (ds, &((*s)[DIMENSION*i]));
        totalmass += ((*mass)[i]);
        V3ADDmuL (((*mass)[i]), ds, dcm); 
    }
    if (totalmass > 0)
    {
        V3DiV (dcm, totalmass);
        for (i=(*np); i--;)
            V3SuB (&((*s)[DIMENSION*i]), dcm);
    }
    return;
} /* end Config_Perturb_s() */


#ifdef _perturb_s
int main (int argc, char *argv[])
{
    Aapp_Define_Config;
    char *out_fname;
    double s_amplitude[3]={0};
    TimeRandomize();
    if (argc == 6)
    {
        out_fname = argv[5];
        s_amplitude[2] = atof(argv[4]);
        s_amplitude[1] = atof(argv[3]);
        s_amplitude[0] = atof(argv[2]);
    }
    else if (argc == 5)
    {
        out_fname = argv[4];
        s_amplitude[1] = atof(argv[3]);
        s_amplitude[0] = atof(argv[2]);
    }
    else if (argc == 4)
    {
        out_fname = argv[3];
        s_amplitude[0] = atof(argv[2]);
    }
    else
    {
        printf ("\nPurpose: perturb each atom in H(i,:) direction\n"
                "(i=1..3) by uniformly distributed random amplitude(s)\n"
                "(-s_amplitude[i]/2, s_amplitude[i]/2). The center\n"
                "of mass is kept fixed.\n\n");
        printf ("Usage: %s in_fname 0.01 out_fname (shake in 1)\n"
                "       %s in_fname 0.01 0.02 out_fname (shake in 1,2)\n"
                "       %s in_fname 0 0.03 0.005 out_fname (shake in 2,3)\n\n",
                argv[0], argv[0], argv[0]);
        return (1);
    }
    printf ("Loading \"%s\"...\n\n", argv[1]);
    CONFIG_LOAD (argv[1], Config_Aapp_to_Alib);
    Config_Perturb_s (Config_Aapp_to_Alib, s_amplitude);
    Config_save (Config_Aapp_to_Alib, TRUE, out_fname);
    printf ("Shake in h1[] by %g A, h2[] by %g A, h3[] by %g A:\n",
            V3LENGTH(H[0]) * fabs(s_amplitude[0])/2,
            V3LENGTH(H[1]) * fabs(s_amplitude[1])/2,
            V3LENGTH(H[2]) * fabs(s_amplitude[2])/2 );
    printf ("\"%s\" -> \"%s\".\n", argv[1], out_fname);
    return (0);
}
#endif /* _perturb_s */


/* Perturb every atom in each Cartesian direction by uniformly */
/* distributed random (-x_amplitude[i]/2, x_amplitude[i]/2).   */
void Config_perturb_x (Alib_Declare_Config, double x_amplitude[3])
{
    register int i;
    double cc;
    V3 dx, ds;
    M3 HI;
    M3INV (H, HI, cc);
    for (i=(*np); i--;)
    {
        dx[0] = FRANDOM() * x_amplitude[0];
        dx[1] = FRANDOM() * x_amplitude[1];
        dx[2] = FRANDOM() * x_amplitude[2];
        V3mM3 (dx, HI, ds);
        V3AdD (ds, &((*s)[DIMENSION*i]));
    }
    return;
} /* end Config_perturb_x() */


/* Same as Config_perturb_x() except center of mass does not move */
void Config_Perturb_x (Alib_Declare_Config, double x_amplitude[3])
{
    register int i;
    double cc, totalmass;
    V3 dx, ds, dcm;
    M3 HI;
    M3INV (H, HI, cc);
    V3ZERO(dcm);
    totalmass = 0;
    for (i=(*np); i--;)
    {
        dx[0] = FRANDOM() * x_amplitude[0];
        dx[1] = FRANDOM() * x_amplitude[1];
        dx[2] = FRANDOM() * x_amplitude[2];
        V3mM3 (dx, HI, ds);
        V3AdD (ds, &((*s)[DIMENSION*i]));
        totalmass += ((*mass)[i]);
        V3ADDmuL (((*mass)[i]), ds, dcm); 
    }
    V3DiV (dcm, totalmass);
    for (i=(*np); i--;)
        V3SuB (&((*s)[DIMENSION*i]), dcm);
    return;
} /* end Config_Perturb_x() */


#ifdef _perturb_x
int main (int argc, char *argv[])
{
    Aapp_Define_Config;
    char *out_fname;
    double x_amplitude[3]={0};
    TimeRandomize();
    if (argc == 6)
    {
        out_fname = argv[5];
        x_amplitude[2] = atof(argv[4]);
        x_amplitude[1] = atof(argv[3]);
        x_amplitude[0] = atof(argv[2]);
    }
    else if (argc == 5)
    {
        out_fname = argv[4];
        x_amplitude[1] = atof(argv[3]);
        x_amplitude[0] = atof(argv[2]);
    }
    else if (argc == 4)
    {
        out_fname = argv[3];
        x_amplitude[0] = atof(argv[2]);
    }
    else
    {
        printf ("\nPurpose: perturb each atom in Cartesian directions [A] "
                "by uniformly\ndistributed random amplitude "
                "(-x_amplitude[i]/2, x_amplitude[i]/2).\n"
                "The center of mass is kept fixed.\n\n");
        printf ("Usage: %s in_fname 0.1 out_fname (shake in x)\n"
                "       %s in_fname 0.1 0.06 out_fname (shake in x,y)\n"
                "       %s in_fname 0 0.06 0.1 out_fname (shake in y,z)\n\n",
                argv[0], argv[0], argv[0]);
        return (1);
    }
    printf ("Loading \"%s\"...\n\n", argv[1]);
    CONFIG_LOAD (argv[1], Config_Aapp_to_Alib);
    Config_Perturb_x (Config_Aapp_to_Alib, x_amplitude);
    Config_save (Config_Aapp_to_Alib, TRUE, out_fname);
    printf ("Shake each atom in x by <|%g| A, y by <|%g| A, z by <|%g| A:\n",
            fabs(x_amplitude[0])/2,
            fabs(x_amplitude[1])/2,
            fabs(x_amplitude[2])/2 );
    printf ("\"%s\" -> \"%s\".\n", argv[1], out_fname);
    return (0);
}
#endif /* _perturb_x */


#define valid_atom(i) ( ((i)>=0) && ((i)<(*np)) )

/* Convert tagged atoms to the kind specified   */
/* behind; return the number of converted atoms */
int Config_convert_taglist_to (char taglist[], Alib_Declare_Config,
                               char *atom_symbol, double atom_mass_IN_AMU)
{
    int i, atoms_converted;
    for (atoms_converted=i=0; i<(*np); i++)
        if (taglist[i])
        {
            safe_symbol ( atom_symbol, SYMBOL(i) );
            (*mass)[i]  = atom_mass_IN_AMU / umass_IN_AMU;
            atoms_converted++;
        }
    return (atoms_converted);
} /* end Config_convert_taglist_to() */


/* particle-reducing transformations */

/* Get rid of atoms with taglist[i]!=0; return the number of atoms lost */
int Config_rid_of_taglist (char taglist[], Alib_Declare_Config)
{
    int i,j,k,loss;
    ConfigStack cs[1];
    CONFIG_PUSH(cs);
    for (loss=i=0; i<(*np); i++) if (taglist[i]) loss++;
    (*np) -= loss;
    Config_alloc(Config_Alib_to_Alib);
    for (i=j=0; i<(*np)+loss; i++)
        if (!taglist[i])
            Config_RETRIEVE(cs,i,Config_Alib_to_Alib,j,k);
    CONFIG_erase (cs);
    return (loss);
} /* end Config_rid_of_taglist() */


/* Delete atoms outside of s0[] * H + [0..1, 0..1, 0..1) * Hparallelepiped */
/* Return the number of atoms retained.                                    */
int Config_keep_atoms_in_parallelepiped
(V3 s0, M3 Hparallelepiped, Alib_Declare_Config)
{
    register int i;
    V3 mys, myx;
    M3 HIparallelepiped;
    char *taglist = IOALLOC(*np);
    M3inv (Hparallelepiped, HIparallelepiped);
    for (i=(*np); i--;)
    {
        V3TriM( &(*s)[DIMENSION*i] );
        V3SUB( &(*s)[DIMENSION*i], s0, mys );
        V3mM3( mys, H, myx );
        V3mM3( myx, HIparallelepiped, mys);
        taglist[i] =
            (mys[0]<-TINY) || (mys[0]>1-TINY) ||
            (mys[1]<-TINY) || (mys[1]>1-TINY) ||
            (mys[2]<-TINY) || (mys[2]>1-TINY);
    }
    Config_rid_of_TAGLIST(taglist,Config_Alib_to_Alib,i);
    return(*np);
} /* end Config_keep_atoms_in_parallelepiped() */


/* Same as Config_keep_atoms_in_parallelepiped except */
/* H[][] is also changed to Hparallelepiped[][].      */
int Config_Keep_atoms_in_parallelepiped
(V3 s0, M3 Hparallelepiped, Alib_Declare_Config)
{
    register int i;
    V3 mys, myx;
    M3 HIparallelepiped;
    char *taglist = IOALLOC(*np);
    M3inv (Hparallelepiped, HIparallelepiped);
    for (i=(*np); i--;)
    {
        V3TriM( &(*s)[DIMENSION*i] );
        V3SUB( &(*s)[DIMENSION*i], s0, mys );
        V3mM3( mys, H, myx );
        V3mM3( myx, HIparallelepiped, mys);
        taglist[i] =
            (mys[0]<-TINY) || (mys[0]>1-TINY) ||
            (mys[1]<-TINY) || (mys[1]>1-TINY) ||
            (mys[2]<-TINY) || (mys[2]>1-TINY);
    }
    Config_rid_of_TAGLIST(taglist,Config_Alib_to_Alib,i);
    for (i=(*np); i--;)
    {
        V3SUB( &(*s)[DIMENSION*i], s0, mys );
        V3mM3( mys, H, myx );
        V3mM3( myx, HIparallelepiped, &(*s)[DIMENSION*i]);
    }
    M3EQV( Hparallelepiped, H );
    return(*np);
} /* end Config_Keep_atoms_in_parallelepiped() */


#ifdef _period_reduce
int main (int argc, char *argv[])
{
    int c;
    double sa, sb, s0[3]={0}, Hparallelepiped[3][3];
    char *out_fname;
    Aapp_Define_Config;
    if (argc != 6)
    {
        printf ("\nPurpose: reduce a periodically repeated configuration\n"
                "         (inverse of mul).\n\n");
        printf ("Usage: %s in_file c sa sb out_file\n\n", argv[0]);
        printf ("       c=[0..2]: edge vector to be reduced.\n");
        printf ("       sa,sb=[0..1]: range in that edge vector we keep.\n\n");
        return (1);
    }
    printf ("Loading \"%s\"...\n\n", argv[1]);
    Config_Load (argv[1], stdout, Config_Aapp_to_Alib);
    c  = atoi(argv[2]);
    sa = atof(argv[3]);
    sb = atof(argv[4]);
    out_fname = argv[5];
    s0[c] = sa;
    M3EQV (H, Hparallelepiped);
    V3MuL (sb-sa, Hparallelepiped[c]);
    Config_Keep_atoms_in_parallelepiped
        (s0, Hparallelepiped, Config_Aapp_to_Alib);
    Config_save (Config_Aapp_to_Alib, TRUE, out_fname);
    printf ("Range [%g %g] of axis %d -> \"%s\".\n", sa, sb, c, out_fname);
    return (0);
}
#endif /* _period_reduce */


/* Get rid of those atoms on the blacklist; return the number of atoms lost */
int Config_rid_of_blacklist
(int how_many, int their_index[], Alib_Declare_Config)
{
    int i;
    char *taglist;
    taglist = IOALLOC(*np);
    for (i=0; i<how_many; i++)
        if (valid_atom(their_index[i])) taglist[their_index[i]] = 1;
    Config_rid_of_TAGLIST(taglist,Config_Alib_to_Alib,i);
    return (i);
} /* end Config_rid_of_blacklist() */


/* Same as Config_rid_of_blacklist with a different input interface */
int Config_decimate (Alib_Declare_Config, int how_many, ...)
{
    int i, *blacklist;
    va_list ap;
    va_start(ap,how_many);
    blacklist = Ialloc(how_many);
    for (i=0; i<how_many; i++) blacklist[i] = va_arg(ap, int);
    Config_rid_of_BLACKLIST(how_many,blacklist,Config_Alib_to_Alib,how_many);
    va_end (ap);
    return(how_many);
} /* end Config_decimate() */


/* Randomly get rid of prescribed amount of atoms; return "how_many" */
int Config_randomly_decimate (int how_many, Alib_Declare_Config)
{
    register int i,j,tmp;
    int *blacklist;
    blacklist = Ialloc(*np);
    Sequentially_index((*np), blacklist);
    /* This way may be slower than directly working with taglist but is */
    /* guaranteed to finish in predictable time when how_many is large. */
    for (i=0; i<how_many; i++)
    {
        j = Fran(i,(*np)-1);
        SWAP(blacklist[i], blacklist[j], tmp);
    }
    Config_rid_of_BLACKLIST(how_many,blacklist,Config_Alib_to_Alib,how_many);
    return(how_many);
} /* end Config_randomly_decimate() */

#ifdef _Config_decimate_TEST
#define CONFIG_FILE "/tmp/a"
#define PDB_FILE    "/tmp/a.pdb"
int main (int argc, char *argv[])
{
    Aapp_Define_Config;
    Config_rebuild_Xtal (Config_Aapp_to_Alib, 1,1,1, Xtal_abstracts+XTAL_FCC,
                         "Ar", -1., -1.);
    printf ("lose %d atoms\n", Config_DECIMATE(0,Config_Aapp_to_Alib));
    Config_SAVE (Config_Aapp_to_Alib, CONFIG_FILE);
    printf ("config saved on \"%s\".\n", CONFIG_FILE);
    Config_save_as_PDB (Config_Aapp_to_Alib, PDB_FILE);
    printf ("config saved on \"%s\".\n", PDB_FILE);
    return (0);
}
#endif /* _Config_decimate_TEST */

bool Config_punch_use_pbc = TRUE;
/* bool Config_punch_use_pbc = FALSE; */

/* punch a real space void |dx| <= r[A] centered at S0,S1,S2 */
int Config_punch_real_space_void
(double S0, double S1, double S2, double r_IN_A, Alib_Declare_Config)
{
    register int i;
    double mys[3], myx[3];
    char *taglist = IOALLOC(*np);
    /* convert r2 to ULENGTH squared */
    r_IN_A *= r_IN_A * A2_IN_uarea;
    for (i=(*np); i--;)
    {
        mys[0] = (*s)[DIMENSION*i]   - S0;
        mys[1] = (*s)[DIMENSION*i+1] - S1;
        mys[2] = (*s)[DIMENSION*i+2] - S2;
        if (Config_punch_use_pbc) V3ImagE( mys );
        V3mM3( mys, H, myx );
        taglist[i] = (V3LENGTH2(myx) <= r_IN_A);
    }
    return(Config_rid_of_TAGLIST(taglist,Config_Alib_to_Alib,i));
} /* end Config_punch_real_space_void() */


/* punch a s-space void |ds| <= sr centered at S0,S1,S2 */
int Config_punch_s_space_void
(double S0, double S1, double S2, double sr, Alib_Declare_Config)
{
    register int i;
    double mys[3];
    char *taglist = IOALLOC(*np);
    sr *= sr;
    for (i=(*np); i--;)
    {
        mys[0] = (*s)[DIMENSION*i]   - S0;
        mys[1] = (*s)[DIMENSION*i+1] - S1;
        mys[2] = (*s)[DIMENSION*i+2] - S2;
        if (Config_punch_use_pbc) V3ImagE( mys );
        taglist[i] = (V3LENGTH2(mys) <= sr);
    }
    return(Config_rid_of_TAGLIST(taglist,Config_Alib_to_Alib,i));
} /* end Config_punch_s_space_void() */


/* punch a real space ellipsoid dx * E[A^-2] * dx' <= 1 centered at S0,S1,S2 */
int Config_punch_real_space_ellipsoid
(double S0, double S1, double S2, double E_IN__A2[3][3], Alib_Declare_Config)
{
    register int i;
    double E[3][3], mys[3], myx[3];
    char *taglist = IOALLOC(*np);
    M3MULTIPLY(uarea_IN_A2, E_IN__A2, E);
    for (i=(*np); i--;)
    {
        mys[0] = (*s)[DIMENSION*i]   - S0;
        mys[1] = (*s)[DIMENSION*i+1] - S1;
        mys[2] = (*s)[DIMENSION*i+2] - S2;
        if (Config_punch_use_pbc) V3ImagE( mys );
        V3mM3( mys, H, myx );
        taglist[i] = (V3ADOT(myx,E,myx,mys) <= 1);
    }
    return(Config_rid_of_TAGLIST(taglist,Config_Alib_to_Alib,i));
} /* end Config_punch_real_space_ellipsoid() */


#ifdef _punchcylinder
int main (int argc, char *argv[])
{
    int c;
    double s0[3], E[3][3], ra_in_A, rb_in_A;
    char *out_fname;
    Aapp_Define_Config;
    if ( (argc != 7) && (argc != 8) )
    {
        printf ("Purpose: punch a cylindrical void through configuration.\n");
        printf ("Usage: %s in_file c sa sb ra_in_A <rb_in_A> out_file\n",
                argv[0]);
        printf ("c [0..2]: cylinder axis index.\n");
        printf ("sa,sb [0..1]: two other reduced coordinates of\n"
                "              the cylinder center.\n");
        printf ("ra_in_A,rb_in_A: radii in sa,sb direction.\n");
        return (1);
    }
    printf ("Loading \"%s\"...\n\n", argv[1]);
    Config_Load (argv[1], stdout, Config_Aapp_to_Alib);
    c = atoi(argv[2]);
    ra_in_A = atof(argv[5]);
    if (argc == 7)
    {
        rb_in_A = ra_in_A;
        out_fname = argv[6];
    }
    else  /* argc == 8 */
    {
        rb_in_A = atof(argv[6]);
        out_fname = argv[7];
    }
    if (c == 0)
    {
        V3ASSIGN (0, atof(argv[3]), atof(argv[4]), s0);
        M3diagonal (0, 1/SQUARE(ra_in_A), 1/SQUARE(rb_in_A), E);
    }
    else if (c == 1)
    {
        V3ASSIGN (atof(argv[3]), 0, atof(argv[4]), s0);
        M3diagonal (1/SQUARE(ra_in_A), 0, 1/SQUARE(rb_in_A), E);
    }
    else
    {
        V3ASSIGN (atof(argv[3]), atof(argv[4]), 0, s0);
        M3diagonal (1/SQUARE(ra_in_A), 1/SQUARE(rb_in_A), 0, E);
    }
    Config_punch_real_space_ellipsoid
        (s0[0], s0[1], s0[2], E, Config_Aapp_to_Alib);
    Config_save (Config_Aapp_to_Alib, TRUE, out_fname);
    printf ("cylindrical void of ra = %g, rb = %g [A] created in axis-%d\n"
            "centered at (%g %g %g) -> \"%s\".\n", ra_in_A, rb_in_A,
            c, s0[0], s0[1], s0[2], out_fname);
    return (0);
}
#endif /* _punchcylinder */


/*************************************************************************/
/* Same as Config_punch_real_space_ellipsoid() except E_IN__A2[3][3] is  */
/* deduced from one of the following inputs: 1,ax,ay,az,a_IN_A; 2,ax,ay, */
/* az,a_IN_A,bx,by,bz,b_in_A; 3,ax,ay,az,a_IN_A,bx,by,bz,b_in_A,c_in_A.  */
/* In which 1 means a fissure of halfwidth a_IN_A; 2 means an elliptical */
/* crack of radii a_IN_A, b_in_A; 3 means an true ellipsoid. If any of   */
/* the _IN_A's is <= 0, it is taken to be +infinity. And, ax,ay,az don't */
/* have to be unit vector, and bx,by,bz don't have to be perpendicular   */
/* to a - the normalization and orthogonalization are taken care of by   */
/* the subroutine. Notice that in order to produce a through fissure or  */
/* crack under PBC, one should use H[i][] as ax,ay,az, bx,by,bz, etc.    */
/*************************************************************************/
int Config_punch_real_space_Ellipsoid
(Alib_Declare_Config, double S0, double S1, double S2, int n_axis, ...)
{
    double A[3][3], AT[3][3], L[3][3], E_IN__A2[3][3];
    va_list ap;
    va_start(ap, n_axis);
    M3ZERO(A);
    M3ZERO(L);
    if (n_axis >= 1)
    {
        A[0][0] = va_arg(ap, double);
        A[0][1] = va_arg(ap, double);
        A[0][2] = va_arg(ap, double);
        V3NORMALIZE(A[0], L[0][0]);
        L[0][0] = va_arg(ap, double);
        if (L[0][0] <= 0) L[0][0] = 0;
        else L[0][0] = 1. / SQUARE(L[0][0]);
    }
    if (n_axis >= 2)
    {
        A[1][0] = va_arg(ap, double);
        A[1][1] = va_arg(ap, double);
        A[1][2] = va_arg(ap, double);
        L[1][1] = V3DOT(A[1], A[0]);
        V3SUBmuL (A[1], L[1][1], A[0]);
        V3NORMALIZE (A[1], L[1][1]);
        L[1][1] = va_arg(ap, double);
        if (L[1][1] <= 0) L[1][1] = 0;
        else L[1][1] = 1. / SQUARE(L[1][1]);
    }
    if (n_axis >= 3)
    {
        V3CROSS(A[0], A[1], A[2]);
        L[2][2] = va_arg(ap, double);
        if (L[2][2] <= 0) L[2][2] = 0;
        else L[2][2] = 1. / SQUARE(L[2][2]);
    }
    M3TRANSPOSE(A, AT);
    M3MUL(AT, L, E_IN__A2);
    M3MUL(E_IN__A2, A, L);
    M3EQV(L, E_IN__A2);
    va_end (ap);
    return( Config_punch_real_space_ellipsoid
            (S0,S1,S2,E_IN__A2, Config_Alib_to_Alib) );
} /* end Config_punch_real_space_Ellipsoid() */


/* Delete atoms more than halfway in real-space distance */
/* from point "in" to point "in+ds". No PBC is assumed.  */
int Config_vcut (double in_S0, double in_S1, double in_S2,
                 double ds0, double ds1, double ds2, Alib_Declare_Config)
{
    register int i;
    double s0[3], ds[3], dx[4], sn[3];
    char *taglist;
    V3ASSIGN (in_S0, in_S1, in_S2, s0);
    V3ASSIGN (  ds0,   ds1,   ds2, ds);
    V3mM3( ds, H, dx );
    dx[3] = V3LENGTH2( dx );
    M3mV3 (H, dx, sn);
    V3MuL (2./dx[3], sn);
    MALLOC ( Config_vcut, taglist, *np, char );
    for (i=(*np); i--;)
    {
        V3SUB ( &((*s)[DIMENSION*i]), s0, ds );
        taglist[i] = (V3DOT(ds,sn) > 1);
    }
    return(Config_rid_of_TAGLIST(taglist,Config_Alib_to_Alib,i));
} /* end Config_vcut() */


#ifdef _vcut
int main (int argc, char *argv[])
{
    int under_pbc, np_old;
    double s0[3], ds[3];
    Aapp_Define_Config;
    if (argc != 10)
    {
        printf ("Purpose: cut configuration beyond (and including)\n"
                "(S0,S1,S2) in the direction dx[] = ds[] * H[][].\n");
        printf ("Usage: %s in_file S0 S1 S2 ds0 ds1 ds2 out_file under_pbc\n",
                argv[0]);
        return (1);
    }
    s0[0] = atof(argv[2]);
    s0[1] = atof(argv[3]);
    s0[2] = atof(argv[4]);
    ds[0] = atof(argv[5]);
    ds[1] = atof(argv[6]);
    ds[2] = atof(argv[7]);
    under_pbc = atoi(argv[9]);
    printf ("Loading \"%s\" (assuming %sPBC)...\n\n", argv[1],
            under_pbc? "" : "no ");
    Config_Load (argv[1], NULL, Config_Aapp_to_Alib);
    np_old = np;
    if (under_pbc) Config_fold_into_PBC (Config_Aapp_to_Alib);
    Config_VCUT( s0[0],s0[1],s0[2], ds[0],ds[1],ds[2], Config_Aapp_to_Alib);
    Config_save(Config_Aapp_to_Alib, FALSE, argv[8]);
    printf ("s0 = (%g %g %g), ds = (%g %g %g) cut -> \"%s\".\n",
            s0[0],s0[1],s0[2], ds[0],ds[1],ds[2], argv[8]);
    printf ("np_old = %d, np = %d\n", np_old, np);
    return (0);
}
#endif /* _vcut */


/* Delete atoms from S[] in the ds[] direction. No PBC is assumed. */
int Config_SCUT (double S0, double S1, double S2,
                 double ds0, double ds1, double ds2, Alib_Declare_Config)
{
    M3 HOLD;
    M3EQV (H, HOLD);
    M3IDENTITY (H);
    Config_VCUT (S0,S1,S2,ds0,ds1,ds2,Config_Alib_to_Alib);
    M3EQV (HOLD, H);
    return (0);
} /* end Config_vcut() */


#ifdef _scut
int main (int argc, char *argv[])
{
    int under_pbc, np_old;
    double s0[3], ds[3];
    Aapp_Define_Config;
    if (argc != 10)
    {
        printf ("Purpose: scut configuration beyond (and including)\n"
                "(S0,S1,S2) in the direction ds[].\n");
        printf ("Usage: %s in_file S0 S1 S2 ds0 ds1 ds2 out_file under_pbc\n",
                argv[0]);
        return (1);
    }
    s0[0] = atof(argv[2]);
    s0[1] = atof(argv[3]);
    s0[2] = atof(argv[4]);
    ds[0] = atof(argv[5]);
    ds[1] = atof(argv[6]);
    ds[2] = atof(argv[7]);
    under_pbc = atoi(argv[9]);
    printf ("Loading \"%s\" (assuming %sPBC)...\n\n", argv[1],
            under_pbc? "" : "no ");
    Config_Load (argv[1], NULL, Config_Aapp_to_Alib);
    np_old = np;
    if (under_pbc) Config_fold_into_PBC (Config_Aapp_to_Alib);
    Config_SCUT( s0[0],s0[1],s0[2], ds[0],ds[1],ds[2], Config_Aapp_to_Alib);
    Config_save(Config_Aapp_to_Alib, FALSE, argv[8]);
    printf ("s0 = (%g %g %g), ds = (%g %g %g) cut -> \"%s\".\n",
            s0[0],s0[1],s0[2], ds[0],ds[1],ds[2], argv[8]);
    printf ("np_old = %d, np_new = %d\n", np_old, np);
    return (0);
}
#endif /* _scut */


#ifdef _div
int main (int argc, char *argv[])
{
    int i;
    Aapp_Define_Config;
    double nc[3];
    if (argc != 6)
    {
        printf ("Purpose: save [0,1/nc0) x [0,1/nc1) x [0,1/nc2) portion of a "
                "configuration.\n");
        printf ("Usage: %s in_file nc0 nc1 nc2 out_file\n", argv[0]);
        return (1);
    }
    printf ("Loading \"%s\"...\n\n", argv[1]);
    Config_Load (argv[1], stdout, Config_Aapp_to_Alib);
    Config_fold_into_PBC (Config_Aapp_to_Alib);
    nc[0] = atof(argv[2]);
    nc[1] = atof(argv[3]);
    nc[2] = atof(argv[4]);
    Config_SCUT (1./nc[0],0,0, 1.,0,0, Config_Aapp_to_Alib);
    Config_SCUT (0,1./nc[1],0, 0,1.,0, Config_Aapp_to_Alib);
    Config_SCUT (0,0,1./nc[2], 0,0,1., Config_Aapp_to_Alib);
    for (i=DIMENSION; i--;)
    {
        if (nc[i]<1) pe ("nc[%d] (now=%g) must >= 1", i, nc[i]);
        V3DiV (H[i], nc[i]);
    }
    for (i=np; i--;) V3ProD (s+3*i, nc);
    Config_save (Config_Aapp_to_Alib, TRUE, argv[5]);
    printf ("1/%g x 1/%g x 1/%g of \"%s\" -> \"%s\".\n", nc[0], nc[1], nc[2],
            argv[1], argv[5]);
    return (0);
}
#endif /* _div */


/* Particle-increasing transformations */

/* Add an atom to config */
int Config_add_atom
(int index_after_insertion, char *atom_symbol, double atom_mass_IN_AMU,
 double atom_s0, double atom_s1, double atom_s2, double atom_s10_IN__NS,
 double atom_s11_IN__NS, double atom_s12_IN__NS, Alib_Declare_Config)
{
    int i,k;
    ConfigStack cs[1];
    CONFIG_PUSH(cs);
    (*np)++;
    if (OUW(index_after_insertion,*np))
        pe ("Config_add_atom: index_after_insertion = %d not between [0,%d]\n",
            index_after_insertion, (*np)-1);
    Config_alloc(Config_Alib_to_Alib);
    for (i=0; i<index_after_insertion;)
        Config_RETRIEVE(cs,i,Config_Alib_to_Alib,i,k);
    safe_symbol ( atom_symbol, SYMBOL(index_after_insertion) );
    (*mass)[index_after_insertion]  = atom_mass_IN_AMU / umass_IN_AMU;
    (*s)[DIMENSION*index_after_insertion]   = atom_s0;
    (*s)[DIMENSION*index_after_insertion+1] = atom_s1;
    (*s)[DIMENSION*index_after_insertion+2] = atom_s2;
    (*s1)[DIMENSION*index_after_insertion]   = atom_s10_IN__NS * utime_IN_NS;
    (*s1)[DIMENSION*index_after_insertion+1] = atom_s11_IN__NS * utime_IN_NS;
    (*s1)[DIMENSION*index_after_insertion+2] = atom_s12_IN__NS * utime_IN_NS;
    for (i=index_after_insertion+1; i<(*np);)
        Config_RETRIEVE(cs,i-1,Config_Alib_to_Alib,i,k);
    CONFIG_erase (cs);
    return(1);
} /* end Config_add_atom() */

#ifdef _Config_add_atom_TEST
#define NC           4
#define CONFIG_FILE "/tmp/a"
int main (int argc, char *argv[])
{
    Aapp_Define_Config;
    Config_rebuild_Xtal (Config_Aapp_to_Alib, NC,NC,NC,
                         Xtal_abstracts+XTAL_FCC, "Ar", -1., -1.);
    printf ("adding %d atoms\n",
            Config_add_atom (np, "Ar", ATOM_MASS_IN_AMU(Z_Ar),
                             s[0]+.25/NC, s[1]+.25/NC, s[2]+.25/NC,
                             0.,0.,0., Config_Aapp_to_Alib));
    Config_SAVE (Config_Aapp_to_Alib, CONFIG_FILE);
    printf ("config saved on \"%s\".\n", CONFIG_FILE);
    return (0);
}
#endif /* _Config_add_atom_TEST */


/* Concatenate "cs" configuration to the end of current configuration, */
/* keeping H[][] of the current configuration. "cs" will be erased.    */
int Config_cat (ConfigStack *cs, Alib_Declare_Config)
{
    int i,j,k;
    (*np) += cs->NP;
    Config_realloc (Config_Alib_to_Alib);
    j = (*np) - cs->NP;
    for (i=0; i<cs->NP; i++)
        Config_RETRIEVE (cs, i, Config_Alib_to_Alib, j, k);
    CONFIG_erase (cs);
    return (1);
} /* end Config_cat() */


/* Make n0 x n1 x n2 stacked copies of the current configuration */
void Config_multiply (int n0, int n1, int n2, Alib_Declare_Config)
{
    int oldnp, i, j, k, m, n, w;
    ConfigStack cs[1];
    CONFIG_PUSH (cs);
    oldnp = (*np);
    (*np) *= n0 * n1 * n2;
    Config_alloc (Config_Alib_to_Alib);
    n = 0;
    for (i=0; i<n0; i++)
        for (j=0; j<n1; j++)
            for (k=0; k<n2; k++)
                for (m=0; m<oldnp; m++)
                {
                    Config_RETRIEVE(cs,m,Config_Alib_to_Alib,n,w);
                    (*s)[DIMENSION*n-3] = ((*s)[DIMENSION*n-3] + i) / n0;
                    (*s)[DIMENSION*n-2] = ((*s)[DIMENSION*n-2] + j) / n1;
                    (*s)[DIMENSION*n-1] = ((*s)[DIMENSION*n-1] + k) / n2;
                    (*s1)[DIMENSION*n-3] /= n0;
                    (*s1)[DIMENSION*n-2] /= n1;
                    (*s1)[DIMENSION*n-1] /= n2;
                }
    for (k=0; k<CONFIG_num_auxiliary; k++)
    {
        REALLOC(Config_multiply, CONFIG_auxiliary[k], (*np), double);
        for (i=1; i<n0*n1*n2; i++)
            VEQV(oldnp, CONFIG_auxiliary[k], CONFIG_auxiliary[k]+i*oldnp);
    }
    V3MuL(n0,H[0]);
    V3MuL(n1,H[1]);
    V3MuL(n2,H[2]);
    CONFIG_erase (cs);
    return;
} /* end Config_multiply() */


#ifdef _mul
int main (int argc, char *argv[])
{
    Aapp_Define_Config;
    int nc[3];
    if (argc != 6)
    {
        printf ("Purpose: make nc0 x nc1 x nc2 stacked copies of a "
                "configuration.\n");
        printf ("Usage: %s in_file nc0 nc1 nc2 out_file\n", argv[0]);
        return (1);
    }
    printf ("Loading \"%s\"...\n\n", argv[1]);
    Config_Load (argv[1], stdout, Config_Aapp_to_Alib);
    nc[0] = atoi(argv[2]);
    nc[1] = atoi(argv[3]);
    nc[2] = atoi(argv[4]);
    Config_multiply (nc[0], nc[1], nc[2], Config_Aapp_to_Alib);
    Config_save (Config_Aapp_to_Alib, TRUE, argv[5]);
    printf ("%d x %d x %d x \"%s\" -> \"%s\".\n", nc[0], nc[1], nc[2],
            argv[1], argv[5]);
    return (0);
}
#endif /* _mul */


/** Configuration Builders **/

/* The following crystal abstracts are provided as a public service */
Xtal_Abstract Xtal_abstracts [XTAL_ABSTRACTS_MAX] = {
    /* XTAL_fcc: */
    {1, 1, 90,90,90, 1,1, {{0,.5,.5},{.5,0,.5},{.5,.5,0}}, {0,0,0}, {0}},
    /* XTAL_FCC: */
    {4, 1, 90,90,90, 1,1, {{1,0,0},{0,1,0},{0,0,1}},
     {0,0,0,0,.5,.5,.5,0,.5,.5,.5,0}, {0,0,0,0}},
    /* XTAL_Fcc: x=[1 1 -2]/2, y=[1 1 1], z=[1 -1 0]/2 */
    {6, 1, 90,90,90, SQRT2,1/SQRT3, {{1,0,0},{0,1,0},{0,0,1}},
     {0,0,0, .5,0,.5, THIRD,THIRD,0, 5*SIXTH,THIRD,.5,
      SIXTH,2*THIRD,.5, 2*THIRD,2*THIRD,0,}, {0,0,0,0,0,0}},
    /* XTAL_bcc: */
    {1, 1, 90,90,90, 1,1, {{-.5,.5,.5},{.5,-.5,.5},{.5,.5,-.5}}, {0,0,0}, {0}},
    /* XTAL_BCC: */
    {2, 1, 90,90,90, 1,1, {{1,0,0},{0,1,0},{0,0,1}}, {0,0,0,.5,.5,.5}, {0,0}},
    /* XTAL_Bcc: x=[1 1 1]/2, y=[0 -1 1], z=[2 -1 -1] */
    {6, 1, 90,90,90, SQRT8/SQRT3,SQRT8, {{1,0,0},{0,1,0},{0,0,1}},
     {0,0,0, 2*THIRD,0,THIRD, THIRD,0,2*THIRD,
      THIRD,HALF,SIXTH, 0,HALF,HALF,  2*THIRD,HALF,5*SIXTH}, {0,0,0,0,0,0}},
    /* XTAL_bCC: x=[1 1 -2], y=[-1 1 0], z=[1 1 1]/2 */
    {6, 1, 90,90,90, 1./SQRT3,1./SQRT8, {{1,0,0},{0,1,0},{0,0,1}},
     {0,0,0, THIRD,0,THIRD, 2*THIRD,0,2*THIRD,
      HALF,HALF,0, 5*SIXTH,HALF,THIRD, SIXTH,HALF,2*THIRD}, {0,0,0,0,0,0}},
    /* XTAL_dia: */
    {2, 1, 90,90,90, 1,1, {{0,.5,.5},{.5,0,.5},{.5,.5,0}},
     {0,0,0,.25,.25,.25}, {0,0}},
    /* XTAL_Dia: x=[1 1 -2], y=[1 1 1], z=[1 -1 0] */
    {12, 1, 90,90,90, SQRT2,1/SQRT3, {{1,0,0},{0,1,0},{0,0,1}},
     {0,0,0, .5,0,.5,
      0,FOURTH,0, .5,FOURTH,.5,
      THIRD,THIRD,0, 5*SIXTH,THIRD,.5,
      THIRD,THIRD+FOURTH,0, 5*SIXTH,THIRD+FOURTH,.5,
      SIXTH,2*THIRD,.5, 2*THIRD,2*THIRD,0,
      SIXTH,2*THIRD+FOURTH,.5, 2*THIRD,2*THIRD+FOURTH,0},
     {0,0,0,0,0,0,0,0,0,0,0,0}},
    /* XTAL_DiA: */
    {12, 1, 90,90,90, SQRT6,SQRT3, {{1,0,0},{0,1,0},{0,0,1}},
     {0,0,0, .5,0,.5,   0,3*FOURTH,0, .5,3*FOURTH,.5,
      .5,THIRD,SIXTH,   0,THIRD,2*THIRD,
      .5,THIRD-FOURTH,SIXTH, 0,THIRD-FOURTH,2*THIRD,
      0,2*THIRD,THIRD, .5,2*THIRD,5*SIXTH,
      0,2*THIRD-FOURTH,THIRD, .5,2*THIRD-FOURTH,5*SIXTH},
     {0,0,0,0,0,0,0,0,0,0,0,0}},
    /* XTAL_DIA: */
    {8, 1, 90,90,90, 1,1, {{1,0,0},{0,1,0},{0,0,1}},
     {0,0,0,.25,.25,.25, 0,.5,.5,.25,.75,.75, .5,0,.5,.75,.25,.75,
      .5,.5,0,.75,.75,.25}, {0,0,0,0,0,0,0,0}},
    /* XTAL_sc, XTAL_SC */
    {1, 1, 90,90,90, 1,1, {{1,0,0},{0,1,0},{0,0,1}}, {0,0,0}, {0}},
    /* XTAL_rho */
    {0},
    /* XTAL_RHO */
    {0},
    /* XTAL_hcp */
    {2, 1, 90,90,120, 1,-SQRT8/SQRT3, {{1,0,0},{0,1,0},{0,0,1}},
     {0,0,0, 2./3,1./3,.5}, {0,0}},
    /* XTAL_Hcp: x=[1 -1 0 0], y=[0 0 0 1], z=[-1 -1 2 0]/3: XTAL_Fcc analog */
    {4, 1, 90,90,90, -SQRT8/3,1/SQRT3, {{1,0,0},{0,1,0},{0,0,1}},
     {0,0,0, .5,0,.5, THIRD,.5,0, 5*SIXTH,.5,.5}, {0,0,0,0,0,0}},
    /* XTAL_HCP */
    {4, 1, 90,90,90, SQRT3,-SQRT8/SQRT3, {{1,0,0},{0,1,0},{0,0,1}},
     {0,0,0, .5,.5,0, .5,1./6,.5, 0,2./3,.5}, {0,0,0,0}},
    /* XTAL_gra2H */
    {4, 1, 90,90,120, 1,-2.72, {{1,0,0},{0,1,0},{0,0,1}},
     {0,0,0, 2./3,1./3,0, 2./3,1./3,.5, 1./3,2./3,.5}, {0,0,0,0}},
    /* XTAL_GRA2H */
    {0},
    /* XTAL_gra3R */
    {0},
    /* XTAL_GRA3R */
    {0},
    /* XTAL_tet, XTAL_TET */
    {1, 1, 90,90,90, 1,-1, {{1,0,0},{0,1,0},{0,0,1}}, {0,0,0}, {0}},
    /* XTAL_ort, XTAL_ORT */
    {1, 1, 90,90,90, -1,-1, {{1,0,0},{0,1,0},{0,0,1}}, {0,0,0}, {0}},
    /* XTAL_zns */
    {2, 2, 90,90,90, 1,1, {{0,.5,.5},{.5,0,.5},{.5,.5,0}},
     {0,0,0, .25,.25,.25}, {0,1}},
    /* XTAL_ZNS */
    {8, 2, 90,90,90, 1,1, {{1,0,0},{0,1,0},{0,0,1}},
     {0,0,0, 0,.5,.5, .5,0,.5, .5,.5,0,
      .25,.25,.25, .25,.75,.75, .75,.25,.75, .75,.75,.25}, {0,0,0,0,1,1,1,1}},
    /* XTAL_nacl */
    {2, 2, 90,90,90, 1,1, {{0,.5,.5},{.5,0,.5},{.5,.5,0}},
     {0,0,0, .5,.5,.5}, {0,1}},
    /* XTAL_NACL */
    {8, 2, 90,90,90, 1,1, {{1,0,0},{0,1,0},{0,0,1}},
     {0,0,0, 0,.5,.5, .5,0,.5, .5,.5,0, .5,0,0, .5,.5,.5, 0,0,.5, 0,.5,0},
     {0,0,0,0,1,1,1,1}},
    /* XTAL_cscl, XTAL_CSCL */
    {2, 2, 90,90,90, 1,1, {{1,0,0},{0,1,0},{0,0,1}}, {0,0,0,.5,.5,.5}, {0,1}},
    /* XTAL_wc */
    {2, 2, 90,90,120, 1,-5.361/5.492, {{1,0,0},{0,1,0},{0,0,1}},
     {0,0,0, 2./3,1./3,.5}, {0,1}},
    /* XTAL_WC */
    {4, 2, 90,90,90, SQRT3,-5.361/5.492, {{1,0,0},{0,1,0},{0,0,1}},
     {0,0,0, .5,.5,0, .5,1./6,.5, 0,2./3,.5}, {0,0,1,1}},
    /* XTAL_ni3al */
    {4, 2, 90,90,90, 1,1, {{1,0,0},{0,1,0},{0,0,1}},
     {0,0,0, 0,.5,.5, .5,0,.5, .5,.5,0}, {1,0,0,0}},
    /* XTAL_Ni3Al: x=[1 1 -2], y=[1 1 1], z=[1 -1 0] */
    {24, 2, 90,90,90, 1/SQRT2,1/SQRT3, {{1,0,0},{0,1,0},{0,0,1}},
     {0,0,0, .5,0,.5,
      0,0,.5, .25,0,.25, .25,0,.75, .5,0,1, .75,0,.25, .75,0,.75,
      1./3,1./3,0, 5./6,1./3,.5,
      1./3,1./3,.5, 7./12,1./3,.25, 7./12,1./3,.75, 5./6,1./3,1,
      1./12,1./3,.25, 1./12,1./3,.75,
      1./6,2./3,.5, 2./3,2./3,0,
      1./6,2./3,1., 5./12,2./3,.75, 5./12,2./3,.25, 2./3,2./3,.5,
      11./12,2./3,.75, 11./12,2./3,.25},
     {1,1,0,0,0,0,0,0,
      1,1,0,0,0,0,0,0,
      1,1,0,0,0,0,0,0}},
};


#define ATOM_RADIUS_IN_A(Z) ATOM_EMPIRICAL_RADIUS_IN_A(Z)
/******************************************************************/
/* Easy way to build a Xtal of n0 x n1 x n2 unit cells; each Xtal */
/* Abstract has its own argument list in the following fashion:   */
/* [(symbol_0,mass_0_IN_AMU),...], a_IN_A,(b_IN_A, c_IN_A); in    */
/* which only the symbols need to be verbatim. mass_i_IN_AMU, if  */
/* given a negative value, will be filled in from the periodic    */
/* table. a_IN_A, if given a negative value, will be deduced from */
/* the atomic radii. For symmetries with a=b=c, do not give b,c;  */
/* otherwise specify the free b or c or both in that order. If    */
/* given b_IN_A(c_IN_A) < 0, we will fill in the standard ratios. */
/******************************************************************/
void Config_rebuild_Xtal
(Alib_Declare_Config, int n0, int n1, int n2, Xtal_Abstract *st, ...)
{
    register int i,j;
    char *symbol_name, *q;
    double *mass_in_amu, T[3][3], nearest_bond_IN_ulength;
    Crystallographic X;
    va_list ap;
    va_start(ap,st);
    symbol_name = IOALLOC(st->nkd*SYMBOL_SIZE);
    mass_in_amu = Valloc(st->nkd);
    for (i=0; i<st->nkd; i++)
    { /* build the lookup table */
        q = va_arg(ap, char *);
        safe_symbol(q,Symbol(symbol_name,i));
        mass_in_amu[i] = va_arg(ap, double);
        if ( mass_in_amu[i] < 0 )
        {
            for (j=1; j<=MENDELEYEV_MAX; j++)
                if (!strcmp(Symbol(symbol_name,i),ATOM_SYMBOL(j))) break;
            if (j > MENDELEYEV_MAX)
                pe("Config_rebuild_Xtal: mass of \"%2s\" not found in table;\n"
                   "please specify it explicitly in amu.\n",
                   Symbol(symbol_name,i));
            else mass_in_amu[i] = ATOM_MASS_IN_AMU(j);
        }
    }
    (*np) = st->npa;
    Config_realloc (Config_Alib_to_Alib);
    M3INV(st->h, T, nearest_bond_IN_ulength);
    for (i=0; i<(*np); i++)
    { /* from lookup table to individual atoms */
        if ((st->kind[i]<0) || (st->kind[i]>=st->nkd))
            pe("Config_rebuild_Xtal: abstract st->kind[%d]=%d illegal "
               "(st->nkd=%d).\n", i, st->kind[i], st->nkd);
        strcpy(SYMBOL(i),Symbol(symbol_name,st->kind[i]));
        (*mass)[i] = mass_in_amu[st->kind[i]] / umass_IN_AMU;
        V3mM3(&(st->s[DIMENSION*i]), T, &((*s)[DIMENSION*i]));
        V3ZERO(&((*s1)[DIMENSION*i]));
    }
    X.a = va_arg(ap, double);
    if (X.a < 0)
    { /* deduce the nearest neighbor distance */
        CrystallographicAssign( 1, fabs(st->ba_ratio), fabs(st->ca_ratio),
                                st->alpha, st->beta, st->gamma, X );
        Crystallographic_to_H(X,T);
        M3MUL(st->h,T,H);
        i = Xtal_analyze_nearest_bond (Config_Alib_to_Alib, 0,
                                       &nearest_bond_IN_ulength, NULL);
        for (j=1; j<=MENDELEYEV_MAX; j++)
            if (!strcmp(SYMBOL(i), ATOM_SYMBOL(j))) break;
        if ( (j > MENDELEYEV_MAX) || (ATOM_RADIUS_IN_A(j) == NVL) )
            pe ("Config_rebuild_Xtal: charge radius of \"%2s\" not found "
                "in table\nplease specify a_in_A explicitly.\n", SYMBOL(i));
        else X.a = ATOM_RADIUS_IN_A(j);
        i = 0;
        for (j=1; j<=MENDELEYEV_MAX; j++)
            if (!strcmp(SYMBOL(i), ATOM_SYMBOL(j))) break;
        if ( (j > MENDELEYEV_MAX) || (ATOM_RADIUS_IN_A(j) == NVL) )
            pe ("Config_rebuild_Xtal: charge radius of \"%2s\" not found "
                "in table\nplease specify a_in_A explicitly.\n", SYMBOL(i));
        X.a += ATOM_RADIUS_IN_A(j);
        X.a /= nearest_bond_IN_ulength;
    }
    X.a /= ulength_IN_A;
    if (st->ba_ratio > 0) X.b = X.a * st->ba_ratio;
    else
    {
        X.b = va_arg(ap, double);
        if ( X.b < 0 ) X.b = - X.a * st->ba_ratio;
    }
    if (st->ca_ratio > 0) X.c = X.a * st->ca_ratio;
    else
    {
        X.c = va_arg(ap, double);
        if ( X.c < 0 ) X.c = - X.a * st->ca_ratio;
    }
    CrystallographicAssign (X.a, X.b, X.c, st->alpha, st->beta, st->gamma, X);
    Crystallographic_to_H(X,T);
    M3MUL(st->h,T,H);
    va_end (ap);
    free(symbol_name);
    free(mass_in_amu);
    Config_SET_CM(Config_Alib_to_Alib);
    Config_multiply(n0, n1, n2, Config_Alib_to_Alib);
    return;
} /* end Config_rebuild_Xtal() */


#ifdef _Config_rebuild_Xtal_TEST
#define ULENGTH_IN_A 2.
#define CONFIG_FILE "/tmp/a"
#define PDB_FILE    "/tmp/a.pdb"
int main (int argc, char *argv[])
{
    Aapp_Define_Config;
    /* Config_rebuild_Xtal (Config_Aapp_to_Alib, 2,2,2, */
    /* Xtal_abstracts+XTAL_DIA, "Si", -1., -1.); */
    /* Config_rebuild_Xtal (Config_Aapp_to_Alib, 2,2,2, */
    /* Xtal_abstracts+XTAL_hcp, "Zr", -1., -1.); */
    /* Config_rebuild_Xtal (Config_Aapp_to_Alib, 2,2,2, */
    /* Xtal_abstracts+XTAL_NACL, */
    /* "Zr",-1., "C",-1., -1.); */
    Config_rebuild_Xtal (Config_Aapp_to_Alib, 2,2,2,
                       Xtal_abstracts+XTAL_ZNS,
                       "Si",-1., "C",-1., -1.);
    /* Config_rebuild_Xtal (Config_Aapp_to_Alib, 8,8,8, */
    /* Xtal_abstracts+XTAL_DIA, "Si", -1., -1.); */
    /* printf ("lose %d atoms\n", Config_punch_real_space_Ellipsoid */
    /* (Config_Aapp_to_Alib, 0.5, 0.5, 0.5,  2, */
    /* H[0][0], H[0][1], H[0][2], 10., */
    /* H[1][0], H[1][1], H[1][2], 5.) ); */
    /* Config_multiply(2, 1, 2, Config_Aapp_to_Alib); */
    Config_SAVE (Config_Aapp_to_Alib, CONFIG_FILE);
    printf ("config saved on \"%s\".\n", CONFIG_FILE);
    Config_save_as_PDB (Config_Aapp_to_Alib, PDB_FILE);
    printf ("config saved on \"%s\".\n", PDB_FILE);
    return (0);
}
#endif /* _Config_rebuild_Xtal_TEST */


/***************************************************************/
/* A general config like those loaded from PDB files has H[][] */
/* and s[], but s[] may not be in [0,1). Sometime that will    */
/* cause problem, like in constructing neighborlist. This      */
/* function redefines H[][],s[] to make s[] in [0,1).          */
/***************************************************************/
void Config_to_bounding_box_config (Alib_Declare_Config, FILE *info)
{
    register int i,j;
    double my[DIMENSION], bound[4][DIMENSION], HOLD[DIMENSION][DIMENSION],
        HI[DIMENSION][DIMENSION];
    for (j=0; j<DIMENSION; j++)
    { /* lower and upper bounds */
        bound[0][j] =  DOUBLE_PRECISION_INFINITY;
        bound[1][j] = -DOUBLE_PRECISION_INFINITY;
    }
    for (i=0; i<*np; i++)
    {
        V3mM3( &(*s)[DIMENSION*i], H, my );
        for (j=0; j<DIMENSION; j++)
        {
            if (my[j] < bound[0][j]) bound[0][j] = my[j];
            if (my[j] > bound[1][j]) bound[1][j] = my[j];
        }
    }
    Fprintf(info, "The orthonormal bounding box (tight) is\n"
            "[%.4f, %.4f] x [%.4f, %.4f] x [%.4f, %.4f],\n",
            bound[0][0],bound[1][0], bound[0][1],bound[1][1],
            bound[0][2],bound[1][2]);
    for (j=0; j<DIMENSION; j++)
    {
        bound[2][j] = (bound[1][j]+bound[0][j])/2;
        bound[3][j] =  bound[1][j]-bound[0][j];
    }
    /* a little margin */
    my[0] = cbrt(bound[3][0]*bound[3][1]*bound[3][2]/(*np));
    if (my[0]==0) my[0] = 1;
    Fprintf(info, "and we will extend each side by %.4f,\n", my[0]);
    V3aDd (bound[3], my[0]);
    M3EQV (H, HOLD);
    M3DIAGONAL (bound[3], H);
    V3MULSUB (0.5,bound[3], bound[2], bound[2]);
    Fprintf(info, "then shift atoms by (%.4f %.4f %.4f) in real space.\n\n",
            bound[2][0], bound[2][1], bound[2][2]);
    M3inv(H,HI);
    for (i=0; i<*np; i++)
    {
        V3mM3( &(*s)[DIMENSION*i], HOLD, my );
        V3AdD( bound[2], my );
        V3mM3( my, HI, &(*s)[DIMENSION*i] );
    }
    Fprintf(info, "WARNING: H[][] and s[] are modified to fit [0,1),\n");
    Fprintf(info, "WARNING: this might not be the right thing to do!!\n\n");
    Config_analyze_H (Config_Alib_to_Alib, info); Fcr(info);
    Config_analyze_density (Config_Alib_to_Alib, info);
    return;
} /* end Config_to_bounding_box_config() */


/* PDB Tags Popularity Index: */
#define PDB_HEADER_ATOM    0
#define PDB_HEADER_HETATM  1
#define PDB_HEADER_CRYST1  2
#define PDB_HEADER_HEADER  3
#define PDB_HEADER_MAX     4
/* J-P. Chang needs atom number > 99999. "ATOM  "->"ATOM" */
char *PDB_HEADERS[PDB_HEADER_MAX]={"ATOM","HETATM","CRYST1","HEADER"};
/* char *PDB_HEADERS[PDB_HEADER_MAX]={"ATOM  ","HETATM","CRYST1","HEADER"}; */
#define PDB_JPCHANG_ATOM_EXTENSION 2
/* J-P. Chang needs atom number > 99999. "ATOM  "->"ATOM" */
/* can never catch them all, so forget about tag checking */
/* http://rcsb.org/pdb/docs/format/pdbguide2.2/guide2.2_frame.html */
#define PDB_HEADER_FTNOTE
#define PDB_HEADER_REMARK
#define PDB_HEADER_SOURCE
#define PDB_HEADER_AUTHOR
#define PDB_HEADER_REVDAT
#define PDB_LINESIZE 256
/* actual linesize should be 83, 256 is to accommodate HTML headers */
/*******************************************************************/
/* Load atomistic configuration from Protein Data Bank format. If  */
/* "info" is not NULL, we will report loading progress to "info".  */
/* If there is no CRYST1 tag in the file, the program will try to  */
/* assign H[][] as an upright bounding box to the atoms.           */
/*******************************************************************/
void Config_load_from_pdb (char *fname, FILE *info, Alib_Declare_Config)
{
    int i,j,k,rebuild=0,must_pbc=0,ext;
    char linebuffer[PDB_LINESIZE],buf[PDB_LINESIZE];
    double my[DIMENSION],HI[DIMENSION][DIMENSION],bound[4][DIMENSION];
    Crystallographic X;
    FILE *in;
    Fprintf(info, "Loading configuration from file \"%s\":\n", fname);
    in = ROpen(fname);
    /* first find out the number of atoms and their bounding box */
    for (i=0; i<DIMENSION; i++)
    { /* lower and upper bounds */
        bound[0][i] =  DOUBLE_PRECISION_INFINITY;
        bound[1][i] = -DOUBLE_PRECISION_INFINITY;
    }
    *np = k = 0;
    do
    {
        i = freadline_matchheader (in, PDB_HEADER_MAX, PDB_HEADERS,
                                   PDB_LINESIZE, linebuffer);
        switch (i)
        {
            case PDB_HEADER_HEADER:
            case -PDB_HEADER_HEADER-5:
                sscanf(linebuffer, "%d %s", &j, buf);
                if (!strcasecmp(buf,"MUST_PBC")) must_pbc=1;
                break;
            case  PDB_HEADER_ATOM:
            case -PDB_HEADER_ATOM-5:
                ext = PDB_JPCHANG_ATOM_EXTENSION;
                goto interpret1;
            case  PDB_HEADER_HETATM:
            case -PDB_HEADER_HETATM-5:
                ext = 0;
          interpret1:
                sscanf(linebuffer, "%d", &j);
                sscanf(linebuffer+ext+24, "%lf %lf %lf", my, my+1, my+2);
                if (j > *np) *np = j;
                k++;
                for (j=0; j<DIMENSION; j++)
                {
                    if (my[j] < bound[0][j]) bound[0][j] = my[j];
                    if (my[j] > bound[1][j]) bound[1][j] = my[j];
                }
        }
    } while (i>=-2);
    if (k != *np)
    {
        Fprintf(info, "Warning: non-compact atom indices: "
                "rebuild indices\n\n");
        *np = k;
        rebuild = 1;
    }
    if ((*np) == 0)
        pe("Config_load_from_pdb: file \"%s\" contains\n"
           "no \"%s\" tags.\n", fname, PDB_HEADERS[PDB_HEADER_ATOM]);
    Fprintf(info, "%d atoms found.. allocating memory.\n\n", *np);
    /* (Re)allocate memory */
    Config_realloc (Config_Alib_to_Alib);
    rewind (in);
    do
    { /* look for CRYST1 tag */
        i = freadline_matchheader (in, PDB_HEADER_MAX, PDB_HEADERS,
                                   PDB_LINESIZE, linebuffer);
        switch (i)
        {
            case  PDB_HEADER_CRYST1:
            case -PDB_HEADER_CRYST1-5:
                sscanf (linebuffer, "%lf %lf %lf %lf %lf %lf",
                        &X.a, &X.b, &X.c, &X.alpha, &X.beta, &X.gamma);
                i = -PDB_HEADER_MAX-5;
        }
    } while (i>=-2);
    V3ZERO (bound[2]);
    if (must_pbc || (i == -PDB_HEADER_MAX-5))
        Crystallographic_to_H(X,H);
    else
    { /* deduce the bounding box */
        for (k=0; k<DIMENSION; k++)
        {
            bound[2][k] = (bound[1][k]+bound[0][k])/2;
            bound[3][k] =  bound[1][k]-bound[0][k];
        }
        /* a little margin */
        my[0] = cbrt(bound[3][0]*bound[3][1]*bound[3][2]/(*np));
        if (my[0]==0) my[0] = 1;
        V3aDd (bound[3], my[0]);
        V3MULSUB (0.5,bound[3], bound[2], bound[2]);
        M3DIAGONAL (bound[3], H);
        Fprintf (info, "Warning: no CRYST1 tag exists -> use atom bounding\n"
                 "box as H[][], which might lead to wrong PBC. To bound them,"
                 "\nall atoms are shifted by (%g %g %g) A.\n\n", bound[2][0],
                 bound[2][1], bound[2][2]);
    }
    M3inv(H,HI);
    rewind (in);
    j = -1;
    do
    { /* now let's get the config */
        i = freadline_matchheader (in, PDB_HEADER_MAX, PDB_HEADERS,
                                   PDB_LINESIZE, linebuffer);
        switch (i)
        {
            case  PDB_HEADER_ATOM:
            case -PDB_HEADER_ATOM-5:
                ext = PDB_JPCHANG_ATOM_EXTENSION;
                goto interpret2;
            case  PDB_HEADER_HETATM:
            case -PDB_HEADER_HETATM-5:
                ext = 0;
          interpret2:
                sscanf(linebuffer, "%d", &k);
                if (rebuild) j++; else j=k-1;
                safe_symbol ( &linebuffer[ext+6], SYMBOL(j) );
                /* Antonino Romano's file: 15-16 */
                if ((*(SYMBOL(j))==' ') && (*(SYMBOL(j)+1)==' '))
                    safe_symbol ( &linebuffer[ext+8], SYMBOL(j) );
                /* nametags like "AC8", "NC2*" */
                if ( (SYMBOL(j)[1]=='C') &&
                     ((SYMBOL(j)[0]=='A') || (SYMBOL(j)[0]=='N')) )
                    SYMBOL(j)[0] = ' ';
                sscanf(linebuffer+ext+24, "%lf %lf %lf", my, my+1, my+2);
                V3AdD (bound[2], my);
                V3mM3 (my, HI, &((*s)[DIMENSION*j]));
                if (IN((*s)[DIMENSION*j],   -EPS,0)) (*s)[DIMENSION*j]  =0;
                if (IN((*s)[DIMENSION*j+1], -EPS,0)) (*s)[DIMENSION*j+1]=0;
                if (IN((*s)[DIMENSION*j+2], -EPS,0)) (*s)[DIMENSION*j+2]=0;
                if (IN((*s)[DIMENSION*j],  1,1+EPS)) (*s)[DIMENSION*j]  =1-EPS;
                if (IN((*s)[DIMENSION*j+1],1,1+EPS)) (*s)[DIMENSION*j+1]=1-EPS;
                if (IN((*s)[DIMENSION*j+2],1,1+EPS)) (*s)[DIMENSION*j+2]=1-EPS;
                if (must_pbc)
                {
                    V3TriM(&((*s)[DIMENSION*j]));
                    /* if ((*s)[DIMENSION*j]  ==1) (*s)[DIMENSION*j]  =0; */
                    /* if ((*s)[DIMENSION*j+1]==1) (*s)[DIMENSION*j+1]=0; */
                    /* if ((*s)[DIMENSION*j+2]==1) (*s)[DIMENSION*j+2]=0; */
                }
                /* " A" nametag means unidentified Atom in PDB file */
                for (k=0; k<=MENDELEYEV_MAX; k++)
                    if (!strcasecmp(SYMBOL(j),ATOM_SYMBOL(k))) break;
                if (k > MENDELEYEV_MAX)
                    pe("Config_load_from_pdb: file \"%s\"\n"
                       "contains symbol \"%s\" in line\n\"%s\"\n"
                       "that is not in our periodic table.\n",
                       fname, SYMBOL(j), linebuffer);
                else
                {
                    (*mass)[j] = ATOM_MASS_IN_AMU(k) / umass_IN_AMU;
                    /* PDB guys like to use upper case */
                    SYMBOL(j)[1] = ATOM_SYMBOL(k)[1];
                }
        }
    } while(i>=-2);
    M3DividE (H, ulength_IN_A);
    bzero((void *)(*s1), DIMENSION*(*np)*sizeof(double) );
    if (info != NULL)
    {
        Config_analyze_H       (Config_Alib_to_Alib, info); fcr(info);
        Config_analyze_density (Config_Alib_to_Alib, info); fcr(info);
        Config_analyze_species (Config_Alib_to_Alib, info); fcr(info);
    }
    fclose(in);
    return;
} /* end Config_load_from_pdb() */

#ifdef _Config_load_from_pdb_TEST
#define ULENGTH_IN_A  2.
/* #define LOAD_FILE "/Home/Archive/NetApp/Mime/DNA.pdb.z" */
/* from http://chemdept.uwsp.edu/pdbs/ */
#define LOAD_FILE "/Home/Archive/NetApp/Mime/1D66.pdb.bz2"
/* from http://alpha2.bmc.uu.se/~tom/pdb_browse.html   */
/* and http://beta.rcsb.org/pdb/pe/explorer/pe_tut.htm */
#define PDB_FILE    "/tmp/b.pdb"
#define CONFIG_FILE "/tmp/b"
int main (int argc, char *argv[])
{
    Aapp_Define_Config;
    Config_load_from_PDB (LOAD_FILE, Config_Aapp_to_Alib);
    Config_SAVE (Config_Aapp_to_Alib, CONFIG_FILE);
    printf ("config saved on \"%s\".\n", CONFIG_FILE);
    Config_save_as_PDB (Config_Aapp_to_Alib, PDB_FILE);
    printf ("config saved on \"%s\".\n", PDB_FILE);
    return (0);
}
#endif /* _Config_load_from_pdb_TEST */


/* possible auxiliary properties stored in CONFIG file: */
int CONFIG_num_auxiliary = 0;
double *CONFIG_auxiliary[CONFIG_MAX_AUXILIARY] = {0};
char CONFIG_auxiliary_name[CONFIG_MAX_AUXILIARY][TERMSIZE] = {{0}};
char CONFIG_auxiliary_unit[CONFIG_MAX_AUXILIARY][TERMSIZE] = {{0}};

/* Free all auxiliary properties */
void Config_free_auxiliary()
{
    int k;
    for (k=0; k<CONFIG_MAX_AUXILIARY; k++)
    {
        Free(CONFIG_auxiliary[k]);
        CONFIG_auxiliary_name[k][0] = EOS;
        CONFIG_auxiliary_unit[k][0] = EOS;
    }
    CONFIG_num_auxiliary = 0;
    return;
} /* end Config_free_auxiliary() */


/* Popularity Index: */
#define CFG_HEADER_COMMENT   0
#define CFG_HEADER_NUMBER    1
#define CFG_HEADER_A         2
#define CFG_HEADER_H0        3
#define CFG_HEADER_TRANSFORM 4
#define CFG_HEADER_NV        5
#define CFG_HEADER_R         6
#define CFG_HEADER_ENTRY     7
#define CFG_HEADER_AUXILIARY 8
#define CFG_HEADER_ETA       9
#define CFG_HEADER_MAX       10
char *CFG_HEADERS[CFG_HEADER_MAX] =
{"#", "Number of particles =", "A =", "H0", "Transform",
 ".NO_VELOCITY.", "R =", "entry_count =", "auxiliary", "eta"};
#define CFG_LINESIZE 1024

/* Load atomistic configuration from Ju Li's CFG ASCII file */
void Config_load (char *fname, FILE *info, Alib_Declare_Config)
{
    int i,j,k,n;
    char linebuffer[CFG_LINESIZE], as[SYMBOL_SIZE+1]={' '};
    char format[4*(CONFIG_MAX_AUXILIARY+DIMENSION+DIMENSION)];
    char oldsymbol[SYMBOL_SIZE]={0}, *ptr;
    double H0[3][3]=ZeroM3, Transform[3][3]=IdentityM3, eta[3][3]=ZeroM3;
    double A=1., R=1., oldmass=DOUBLE_PRECISION_INFINITY;
    bool no_velocity=FALSE, entry_count=0;
    FILE *in;

    Config_free_auxiliary();
    Fprintf(info, "Loading configuration from \"%s\":\n", fname);
    in = ROpen(fname);
  readline:
    i = freadline_matchheader (in, CFG_HEADER_MAX, CFG_HEADERS,
                               CFG_LINESIZE, linebuffer);
    switch(i)
    {
        case CFG_HEADER_COMMENT:
            goto readline;
        case  CFG_HEADER_NUMBER:
        case -CFG_HEADER_NUMBER-5:
            sscanf (linebuffer, "%d", np);
            if ((*np)<=0) pe("Config_load: np=%d <= 0.\n", (*np));
            Fprintf(info, "%d atoms found.. reallocating memory\n\n", *np);
            Config_realloc (Config_Alib_to_Alib);
            break;
        default:
            pe("Config_load: file \"%s\"\n"
               "does not look like Ju Li's CFG format because it does\n"
               "not have \"%s\" as its first line.\n", fname,
               CFG_HEADERS[CFG_HEADER_NUMBER]);
    }
    n = 0;
    do
    {
        i = freadline_matchheader (in, CFG_HEADER_MAX, CFG_HEADERS,
                                   CFG_LINESIZE, linebuffer);
        switch(i)
        {
            case CFG_HEADER_A:
                sscanf (linebuffer, "%lf", &A);
                break;
            case CFG_HEADER_H0:
                sscanf(linebuffer, "(%d,%d) =", &j, &k);
                sscanf(linebuffer+7, "%lf", &(H0[j-1][k-1]));
                break;
            case CFG_HEADER_TRANSFORM:
                sscanf(linebuffer, "(%d,%d) =", &j, &k);
                sscanf(linebuffer+7, "%lf", &(Transform[j-1][k-1]));
                break;
            case CFG_HEADER_ETA:
                sscanf(linebuffer, "(%d,%d) =", &j, &k);
                sscanf(linebuffer+7, "%lf", &(eta[j-1][k-1]));
                eta[k-1][j-1] = eta[j-1][k-1];
                break;
            case CFG_HEADER_NV:
                no_velocity = TRUE;
                break;
            case CFG_HEADER_R:
                sscanf (linebuffer, "%lf", &R);
                break;
            case CFG_HEADER_ENTRY: /* this is the signature of CONFIG format */
                sscanf (linebuffer, "%d", &entry_count);
                j = entry_count - DIMENSION - (no_velocity?0:DIMENSION);
                if (j > CONFIG_MAX_AUXILIARY)
                    pe ("Config_load: CONFIG_MAX_AUXILIARY = %d exceeded.\n",
                        CONFIG_MAX_AUXILIARY);
                for (k=0; k<j; k++)
                    REALLOC (Config_load, CONFIG_auxiliary[k], *np, double);
                for (k=j; k<CONFIG_num_auxiliary; k++)
                {
                    Free (CONFIG_auxiliary[k]);
                    CONFIG_auxiliary_name[k][0] = EOS;
                    CONFIG_auxiliary_unit[k][0] = EOS;
                }
                CONFIG_num_auxiliary = j;
                for (i=0; i<entry_count; i++)
                {
                    format[4*i]   = '%';
                    format[4*i+1] = 'l';
                    format[4*i+2] = 'f';
                    format[4*i+3] = (i < entry_count-1) ? ' ' : EOS;
                }
                break;
            case CFG_HEADER_AUXILIARY:
                sscanf (linebuffer, "[%d]", &k);
                if (OUW(k,CONFIG_num_auxiliary))
                    pe ("Config_load: CONFIG_num_auxiliary = %d exceeded.\n",
                        CONFIG_num_auxiliary);
                ptr = strchr(linebuffer,'=') + 2;
                strncpy(CONFIG_auxiliary_name[k], ptr, TERMSIZE);
                CONFIG_auxiliary_name[k][TERMSIZE-1] = EOS;
                ptr = strchr(CONFIG_auxiliary_name[k], '[');
                if (ptr != NULL)
                {
                    *(ptr-1) = EOS;
                    strncpy(CONFIG_auxiliary_unit[k], ptr+1, TERMSIZE);
                    ptr = strchr(CONFIG_auxiliary_unit[k], ']');
                    if (ptr != NULL) *ptr = EOS;
                }
                else CONFIG_auxiliary_unit[k][0] = EOS;
                break;
            case -2:
            case -3: /* atom data block */
                if ( (*oldsymbol==EOS) &&
                     (oldmass==DOUBLE_PRECISION_INFINITY) )
                { /* the first time encountering a naked sentence */
                    M3Multiply (A/ulength_IN_A, H0);
                    M3MUl (H0, Transform, H);
                    pure_deform (H0, eta, H);
                    R /= cbrt(fabs(M3DETERMINANT(H)/M3DETERMINANT(H0)));
                }
                if (entry_count == 0)
                { /* Config format */
                    fseek(in, -strlen(linebuffer)-1, SEEK_CUR);
                    for (n=0; n<(*np); n++)
                    {
                        fscanf (in, "%lf %s %lf %lf %lf %lf %lf %lf\n",
                                &((*mass)[n]), as+1,
                                &((*s)[DIMENSION*n]),
                                &((*s)[DIMENSION*n+1]),
                                &((*s)[DIMENSION*n+2]),
                                &((*s1)[DIMENSION*n]),
                                &((*s1)[DIMENSION*n+1]),
                                &((*s1)[DIMENSION*n+2]));
                        for (j=2; as[j]!=EOS; j++);
                        for (k=j-SYMBOL_CHAR; k<=j; k++)
                            SYMBOL(n)[k+SYMBOL_CHAR-j] = as[k];
                        (*mass)[n] /= umass_IN_AMU;
                        (*s1)[DIMENSION*n]   *= utime_IN_NS * R;
                        (*s1)[DIMENSION*n+1] *= utime_IN_NS * R;
                        (*s1)[DIMENSION*n+2] *= utime_IN_NS * R;
                    }                    
                }
                else
                { /* CONFIG format */
                  check_singleton:
                    ptr = blank_advance(linebuffer);
                    if (*blank_advance(nonblank_advance(ptr))==EOS)
                    { /* it is a singleton */
                        if (ISALPHA(*ptr))
                        { /* likely to be symbol input */
                            safe_symbol(ptr,SYMBOL(n));
                            if (*oldsymbol)
                            {
                                safe_symbol(ptr,oldsymbol);
                            }
                            else
                            {
                                safe_symbol(ptr,oldsymbol);
                                break;
                            }
                        }
                        else
                        { /* likely to be mass input */
                            sscanf(ptr, "%lf", &((*mass)[n]));
                            (*mass)[n] /= umass_IN_AMU;
                            if (oldmass != DOUBLE_PRECISION_INFINITY)
                                oldmass = (*mass)[n];
                            else
                            {
                                oldmass = (*mass)[n];
                                break;
                            }
                        }
                        goto getnextline;
                    }
                    else
                    { /* it isn't a singleton */
                        safe_symbol(oldsymbol,SYMBOL(n));
                        (*mass)[n] = oldmass;
                    }
                    if (no_velocity)
                    {
                        if ( sscanf(ptr, format,
                                    &((*s)[DIMENSION*n]),
                                    &((*s)[DIMENSION*n+1]),
                                    &((*s)[DIMENSION*n+2]),
                                    CONFIG_auxiliary[0]+n,
                                    CONFIG_auxiliary[1]+n,
                                    CONFIG_auxiliary[2]+n,
                                    CONFIG_auxiliary[3]+n,
                                    CONFIG_auxiliary[4]+n,
                                    CONFIG_auxiliary[5]+n,
                                    CONFIG_auxiliary[6]+n,
                                    CONFIG_auxiliary[7]+n,
                                    CONFIG_auxiliary[8]+n,
                                    CONFIG_auxiliary[9]+n,
                                    CONFIG_auxiliary[10]+n,
                                    CONFIG_auxiliary[11]+n,
                                    CONFIG_auxiliary[12]+n,
                                    CONFIG_auxiliary[13]+n,
                                    CONFIG_auxiliary[14]+n,
                                    CONFIG_auxiliary[15]+n,
                                    CONFIG_auxiliary[16+0]+n,
                                    CONFIG_auxiliary[16+1]+n,
                                    CONFIG_auxiliary[16+2]+n,
                                    CONFIG_auxiliary[16+3]+n,
                                    CONFIG_auxiliary[16+4]+n,
                                    CONFIG_auxiliary[16+5]+n,
                                    CONFIG_auxiliary[16+6]+n,
                                    CONFIG_auxiliary[16+7]+n,
                                    CONFIG_auxiliary[16+8]+n,
                                    CONFIG_auxiliary[16+9]+n,
                                    CONFIG_auxiliary[16+10]+n,
                                    CONFIG_auxiliary[16+11]+n,
                                    CONFIG_auxiliary[16+12]+n,
                                    CONFIG_auxiliary[16+13]+n,
                                    CONFIG_auxiliary[16+14]+n,
                                    CONFIG_auxiliary[16+15]+n,
                                    CONFIG_auxiliary[32+0]+n,
                                    CONFIG_auxiliary[32+1]+n,
                                    CONFIG_auxiliary[32+2]+n,
                                    CONFIG_auxiliary[32+3]+n,
                                    CONFIG_auxiliary[32+4]+n,
                                    CONFIG_auxiliary[32+5]+n,
                                    CONFIG_auxiliary[32+6]+n,
                                    CONFIG_auxiliary[32+7]+n,
                                    CONFIG_auxiliary[32+8]+n,
                                    CONFIG_auxiliary[32+9]+n,
                                    CONFIG_auxiliary[32+10]+n,
                                    CONFIG_auxiliary[32+11]+n,
                                    CONFIG_auxiliary[32+12]+n,
                                    CONFIG_auxiliary[32+13]+n,
                                    CONFIG_auxiliary[32+14]+n,
                                    CONFIG_auxiliary[32+15]+n) != entry_count )
                            pe("Config_load: incomplete row at atom %d.\n",n);
                        V3ZERO (&((*s1)[DIMENSION*n]));
                    }
                    else
                    {
                        if ( sscanf(ptr, format,
                                    &((*s)[DIMENSION*n]),
                                    &((*s)[DIMENSION*n+1]),
                                    &((*s)[DIMENSION*n+2]),
                                    &((*s1)[DIMENSION*n]),
                                    &((*s1)[DIMENSION*n+1]),
                                    &((*s1)[DIMENSION*n+2]),
                                    CONFIG_auxiliary[0]+n,
                                    CONFIG_auxiliary[1]+n,
                                    CONFIG_auxiliary[2]+n,
                                    CONFIG_auxiliary[3]+n,
                                    CONFIG_auxiliary[4]+n,
                                    CONFIG_auxiliary[5]+n,
                                    CONFIG_auxiliary[6]+n,
                                    CONFIG_auxiliary[7]+n,
                                    CONFIG_auxiliary[8]+n,
                                    CONFIG_auxiliary[9]+n,
                                    CONFIG_auxiliary[10]+n,
                                    CONFIG_auxiliary[11]+n,
                                    CONFIG_auxiliary[12]+n,
                                    CONFIG_auxiliary[13]+n,
                                    CONFIG_auxiliary[14]+n,
                                    CONFIG_auxiliary[15]+n,
                                    CONFIG_auxiliary[16+0]+n,
                                    CONFIG_auxiliary[16+1]+n,
                                    CONFIG_auxiliary[16+2]+n,
                                    CONFIG_auxiliary[16+3]+n,
                                    CONFIG_auxiliary[16+4]+n,
                                    CONFIG_auxiliary[16+5]+n,
                                    CONFIG_auxiliary[16+6]+n,
                                    CONFIG_auxiliary[16+7]+n,
                                    CONFIG_auxiliary[16+8]+n,
                                    CONFIG_auxiliary[16+9]+n,
                                    CONFIG_auxiliary[16+10]+n,
                                    CONFIG_auxiliary[16+11]+n,
                                    CONFIG_auxiliary[16+12]+n,
                                    CONFIG_auxiliary[16+13]+n,
                                    CONFIG_auxiliary[16+14]+n,
                                    CONFIG_auxiliary[16+15]+n,
                                    CONFIG_auxiliary[32+0]+n,
                                    CONFIG_auxiliary[32+1]+n,
                                    CONFIG_auxiliary[32+2]+n,
                                    CONFIG_auxiliary[32+3]+n,
                                    CONFIG_auxiliary[32+4]+n,
                                    CONFIG_auxiliary[32+5]+n,
                                    CONFIG_auxiliary[32+6]+n,
                                    CONFIG_auxiliary[32+7]+n,
                                    CONFIG_auxiliary[32+8]+n,
                                    CONFIG_auxiliary[32+9]+n,
                                    CONFIG_auxiliary[32+10]+n,
                                    CONFIG_auxiliary[32+11]+n,
                                    CONFIG_auxiliary[32+12]+n,
                                    CONFIG_auxiliary[32+13]+n,
                                    CONFIG_auxiliary[32+14]+n,
                                    CONFIG_auxiliary[32+15]+n) != entry_count )
                            pe("Config_load: incomplete row at atom %d.\n",n);
                        (*s1)[DIMENSION*n]   *= utime_IN_NS * R;
                        (*s1)[DIMENSION*n+1] *= utime_IN_NS * R;
                        (*s1)[DIMENSION*n+2] *= utime_IN_NS * R;
                    }
                    if ((++n) == (*np)) break;
                  getnextline:
                    if (fgets(linebuffer,CFG_LINESIZE,in))
                    {
                        for (ptr=linebuffer; *ptr!=EOS; ptr++)
                            if ( (*ptr=='\n') || (*ptr=='\r') )
                            {
                                *ptr = EOS;
                                break;
                            }
                        goto check_singleton;
                    }
                    pe ("Config_load: premature ending at atom %d.\n", n);
                }
        }
    } while (i>=-2);
    if (info != NULL)
    {
        Config_analyze_H       (Config_Alib_to_Alib, info); fcr(info);
        Config_analyze_density (Config_Alib_to_Alib, info); fcr(info);
        Config_analyze_species (Config_Alib_to_Alib, info); fcr(info);
    }
    fclose(in);
    return;
} /* end Config_load() */

/* Load atomistic configuration from NetCDF trajectory */
void Config_load_netcdf (char *fname, FILE *info, Alib_Declare_Config)
{
  static int n_aux = 0;
  static int ncid = -1;
  static int frame_dim_id, spatial_dim_id, atom_dim_id, cell_spatial_dim_id, 
    cell_angular_dim_id, label_dim_id;
  static int cell_lengths_var_id = -1, cell_angles_var_id = -1, coordinates_var_id= -1, species_var_id = -1, 
    mass_var_id = -1, Z_var_id = -1;
  static int species_var_dims = -1, mass_var_dims = -1, Z_var_dims = -1;
  static int aux_prop_var_id[CONFIG_MAX_AUXILIARY], aux_prop_dim[CONFIG_MAX_AUXILIARY], 
    aux_prop_type[CONFIG_MAX_AUXILIARY];
  static double aux_scale_factors[CONFIG_MAX_AUXILIARY];
  static double coordinates_scale_factor = 1.0, cell_lengths_scale_factor = 1.0;
  static char oldfname[512] = "";
  static int frame = 0;
  static size_t n_frame, n_atom, n_label, nfilter;

  int i,j,k,n;
  int retval, nvars, inc;
  size_t start1[1], count1[1], start2[2], count2[2],start3[3],count3[3];
  double cell_lengths[3], cell_angles[3], x[3], sp[3], scale_factor;
  double *coords;
  char *species, *nfname, *framestr;
  char buf1[512], buf2[512], buf3[512], namebuf[512], varname[NC_MAX_NAME+1], linebuffer[CFG_LINESIZE];
  nc_type vartype;
  int ndims, dimids[NC_MAX_VAR_DIMS], natts;
  M3 g;
  FILE *filter; 
  int *tmpint, *tmpint3, *atomfilter, gotfilter; 
  double *tmpreal, *tmpreal3;

  Fprintf(info, "Loading NetCDF configuration \"%s\"\n", fname);
  
  strcpy(namebuf, fname);
  framestr = namebuf;
  nfname = strsep(&framestr, ":");
  if (framestr == NULL) {
    strcpy(buf2, "first");
    framestr = buf2;
  }

  strcpy(buf3,nfname);
  strcat(buf3,".filter");
  gotfilter = (filter = fopen(buf3,"r")) != NULL;
  
  //  fprintf(info, "After split: framestr=\"%s\", nfname=\"%s\", fname=\"%s\"\n", framestr, nfname, fname);

  // Is this the same file we used last time?
  if (ncid != -1 && (strcmp(nfname, oldfname) != 0 || strcmp(framestr, "reload") == 0)) {
    Fprintf(info, "Closing NetCDF file \"%s\"\n", oldfname);
    netcdf_check(nc_close(ncid));
    ncid = -1;
  }
     
  // Open file and get dimension and variable information
  if (ncid == -1 ) {
    Config_free_auxiliary();
    Fprintf(info, "Opening NetCDF file \"%s\"\n", nfname);
#ifdef NETCDF4
    netcdf_check(nc_open(nfname, NC_NOWRITE, &ncid));
#else
    netcdf_check(nc_open(nfname, NC_64BIT_OFFSET | NC_NOWRITE, &ncid));
#endif
    strcpy(oldfname, nfname);

    // Inquire dimensions
    netcdf_check(nc_inq_dimid(ncid, "frame", &frame_dim_id));
    netcdf_check(nc_inq_dimid(ncid, "spatial", &spatial_dim_id));
    netcdf_check(nc_inq_dimid(ncid, "atom", &atom_dim_id));
    netcdf_check(nc_inq_dimid(ncid, "cell_spatial", &cell_spatial_dim_id));
    netcdf_check(nc_inq_dimid(ncid, "cell_angular", &cell_angular_dim_id));
    netcdf_check(nc_inq_dimid(ncid, "label", &label_dim_id));

    // Sizes of dimensions
    netcdf_check(nc_inq_dimlen(ncid, frame_dim_id, &n_frame));
    netcdf_check(nc_inq_dimlen(ncid, atom_dim_id, &n_atom));
    netcdf_check(nc_inq_dimlen(ncid, label_dim_id, &n_label));

    n_aux = 0;

    // Loop over all variables in file
    netcdf_check(nc_inq_nvars(ncid, &nvars));
    for (i=0; i<nvars; i++) {
      netcdf_check(nc_inq_var(ncid, i, varname, &vartype, &ndims, dimids, &natts));

      if (nc_get_att_double(ncid, i, "scale_factor", &scale_factor) != NC_NOERR) {
	scale_factor = 1.0;
      }

      if (strcmp(varname,"cell_lengths") == 0) {
	cell_lengths_var_id = i;
	cell_lengths_scale_factor = scale_factor;
	continue;
      }

      if (strcmp(varname,"cell_angles") == 0) {
	cell_angles_var_id = i;
	continue;
      }

      if (strcmp(varname,"coordinates") == 0) {
	coordinates_var_id = i;
	coordinates_scale_factor = scale_factor;
	continue;
      }

      if (strcmp(varname,"species") == 0) {
	species_var_id = i;
	species_var_dims = ndims;
	continue;
      }

      if (strcmp(varname,"mass") == 0) {
	mass_var_id = i;
	mass_var_dims = ndims;
	continue;
      }

/*       if (strcmp(varname,"velocities") == 0) { */
/* 	velocity_var_id = i; */
/* 	continue; */
/*       } */

      if (strcmp(varname,"Z") == 0) {
	Z_var_id = i;
	Z_var_dims = ndims;
	continue;
      }
								
      // Now we're onto to aux props
      // Skip unless it's a global or per-frame and per-atom property
      if (ndims == 1) {
	if (dimids[0] != atom_dim_id)  continue;
      } else if (ndims == 2) {
	if (dimids[0] != frame_dim_id || dimids[1] != atom_dim_id)  continue;
      } else if (ndims == 3) {
	if (dimids[0] != frame_dim_id || dimids[1] != atom_dim_id || dimids[2] != spatial_dim_id)  continue;
      } else
	continue;

      // Only look at real and int properties (for now)
      if (!(vartype == NC_DOUBLE || vartype == NC_FLOAT || vartype == NC_INT)) continue;

      if (nc_get_att_text(ncid, i, "units", buf1) == NC_NOERR) {
	strcpy(CONFIG_auxiliary_unit[n_aux], buf1);
	CONFIG_auxiliary_unit[n_aux][TERMSIZE-1] = EOS;
      }

      if (ndims == 1) {
	fprintf(info, "Config_load_netcdf: single column global aux prop %s varid %d\n", varname, i);
	strcpy(CONFIG_auxiliary_name[n_aux], varname);
	aux_prop_var_id[n_aux] = i;
	aux_prop_dim[n_aux] = -2;
	aux_prop_type[n_aux] = vartype;
	aux_scale_factors[n_aux] = scale_factor;
	n_aux++;
	if (n_aux == CONFIG_MAX_AUXILIARY) {
	  fprintf(info,"Config_load: Warning CONFIG_MAX_AUXILIARY = %d exceeded, skipping remaining properties\n",
		  CONFIG_MAX_AUXILIARY);
	  n_aux = CONFIG_MAX_AUXILIARY;
	  goto AUX_PROPS_DONE;
	}
      } else if (ndims == 2) {
	fprintf(info, "Config_load_netcdf: single column aux prop %s with varid %d\n", varname, i);
	strcpy(CONFIG_auxiliary_name[n_aux], varname);
	aux_prop_var_id[n_aux] = i;
	aux_prop_dim[n_aux] = -1;
	aux_prop_type[n_aux] = vartype;
	aux_scale_factors[n_aux] = scale_factor;
	n_aux++;
	if (n_aux == CONFIG_MAX_AUXILIARY) {
	  fprintf(info,"Config_load: Warning CONFIG_MAX_AUXILIARY = %d exceeded, skipping remaining properties\n",
		  CONFIG_MAX_AUXILIARY);
	  n_aux = CONFIG_MAX_AUXILIARY;
	  goto AUX_PROPS_DONE;
	}
      } else {
	fprintf(info, "Config_load_netcdf: 3-column aux prop %s with varid %d\n", varname, i);
	for (j=0; j<3; j++) {
	  sprintf(CONFIG_auxiliary_name[n_aux], "%s%d", varname, j);
	  aux_prop_var_id[n_aux] = i;
	  aux_prop_dim[n_aux] = j;
	  aux_prop_type[n_aux] = vartype;
	  aux_scale_factors[n_aux] = scale_factor;
	  n_aux++;
	  if (n_aux == CONFIG_MAX_AUXILIARY) {
	    fprintf(info,"Config_load: Warning CONFIG_MAX_AUXILIARY = %d exceeded, skipping remaining properties\n", 
		    CONFIG_MAX_AUXILIARY);
	    n_aux = CONFIG_MAX_AUXILIARY;
	    goto AUX_PROPS_DONE;
	  }
	}
      }
    }
      
  AUX_PROPS_DONE:
    // Check no essential variables are missing
    if (coordinates_var_id == -1) pe("Config_load_netcdf: no coordinates variable in file %s\n", nfname);
    if (species_var_id == -1 && Z_var_id == -1) pe("Config_load_netcdf: neither species nor Z variables in file %s\n", nfname);
    if (cell_lengths_var_id == -1) pe("Config_load_netcdf: no cell_lengths variable in file %s\n", nfname);
    if (cell_angles_var_id == -1) pe("Config_load_netcdf: no cell_angles variable in file %s\n", nfname);

    // Optional variables
    if (mass_var_id == -1) fprintf(info, "Config_load_netcdf: warning: no mass variable in file %s\n", nfname);
    /*     if (velocity_var_id == -1) fprintf(info, "Config_load_netcdf: warning: no velocities variable in file %s\n", nfname); */
    
    fprintf(info, "Config_load_netcdf: got %d frames, %d atoms and %d aux props\n", n_frame, n_atom, n_aux);
  }

  if (sscanf(framestr, "%d", &frame) == 1) {
    // Check frame is in range
    if ((frame < 0) || (frame >= n_frame)) {
      pe("Config_load_netcdf: frame should be in range 0 <= frame <= %d\n", n_frame-1);
    }
  }
  else if (strstr(framestr, "first")) {
    frame = 0;
  } 
  else if (strstr(framestr, "last")) {
    frame = n_frame-1;
  } 
  else if (strstr(framestr, "forward")) {
    if (sscanf(framestr, "forward:%d", &inc) != 1) inc = 1;
    if ((frame += inc) >= n_frame) frame = 0;
  } else if (strstr(framestr, "backward")) {
    if (sscanf(framestr, "backward:%d", &inc) != 1) inc = 1;
    if ((frame -= inc) <= -1) frame = n_frame-1;
  }
  else if (strstr(framestr, "reload")) {
    fprintf(info, "Config_load_netcdf: reloaded file, remaining at frame %d/%d\n", frame, n_frame);
  }
  else {
    pe("Config_load_netcdf: unknown framestr %s\n", framestr);
  }

  fprintf(info, "Config_load_netcdf: loading frame %d\n", frame);
  strcpy(fname, nfname);
  sprintf(buf2, ":%%0%dd", (int)ceil(log10(n_frame)));
  sprintf(buf1, buf2, frame);
  strcat(fname, buf1);

  // Check for .filter file
  if (gotfilter) {
    atomfilter = (int *)malloc(n_atom*sizeof(int));
    for (i=0; i<n_atom; i++) atomfilter[i] = 0;
    nfilter = 0;
    while (fgets(linebuffer,CFG_LINESIZE,filter)) {
      if (sscanf(linebuffer, "%d", &i) != 1) 
	pe("Error reading atom list: %s", linebuffer);
      if(i < 0 || i >= n_atom) pe("Atom filter atom %d out of range 0 < i < %d\n", i, n_atom);
      atomfilter[i] = 1;
      nfilter++;
    }
    fclose(filter);
    *np = nfilter;
    fprintf(info, "Config_load_netcdf: filter %s applied. Selected %d/%d atoms.\n", buf3, nfilter, n_atom);
    strcpy(CONFIG_auxiliary_name[n_aux],"original_index");
    strcpy(CONFIG_auxiliary_unit[n_aux],"");
  }
  else {
    // Include all atoms
    *np = n_atom;
  }

  // Reallocate to correct size
  Config_realloc (Config_Alib_to_Alib);
  
  // Allocate correct number of aux props
  for (k=0; k<n_aux+gotfilter; k++)
    REALLOC (Config_load, CONFIG_auxiliary[k], n_atom, double);
  for (k=n_aux+gotfilter; k<CONFIG_num_auxiliary; k++)
    {
      Free (CONFIG_auxiliary[k]);
      CONFIG_auxiliary_name[k][0] = EOS;
      CONFIG_auxiliary_unit[k][0] = EOS;
    }
  CONFIG_num_auxiliary = n_aux+gotfilter;

  // Load lattice
  start2[0] = frame;
  start2[1] = 0;
  count2[0] = 1;
  count2[1] = 3;
  netcdf_check(nc_get_vara_double(ncid, cell_lengths_var_id, start2, count2, cell_lengths));
  netcdf_check(nc_get_vara_double(ncid, cell_angles_var_id, start2, count2, cell_angles));

  cell_lengths[0] *= cell_lengths_scale_factor;
  cell_lengths[1] *= cell_lengths_scale_factor;
  cell_lengths[2] *= cell_lengths_scale_factor;
    
  fprintf(info, "Config_load_netcdf: cell_lengths=[%f %f %f]\n", cell_lengths[0],
	  cell_lengths[1], cell_lengths[2]);
  fprintf(info, "Config_load_netcdf: cell_angles=[%f %f %f]\n", cell_angles[0],
	  cell_angles[1], cell_angles[2]);

  lattice_abc_to_xyz(cell_lengths, cell_angles, H);
  M3inv (H, g);

/*   fprintf(info, "Config_load_netcdf: H = [[%10.2f %10.2f %10.2f]\n", H[0][0], */
/* 	  H[0][1], H[0][2]);  */
/*   fprintf(info, "                         [%10.2f %10.2f %10.2f]\n", H[1][0], */
/* 	  H[1][1], H[1][2]);  */
/*   fprintf(info, "                         [%10.2f %10.2f %10.2f]]\n", H[2][0], */
/* 	  H[2][1], H[2][2]);  */
/*   fprintf(info, "Config_load_netcdf: g = [[%10.2f %10.2f %10.2f]\n", g[0][0], */
/* 	  g[0][1], g[0][2]); */
/*   fprintf(info, "                         [%10.2f %10.2f %10.2f]\n", g[1][0], */
/* 	  g[1][1], g[1][2]); */
/*   fprintf(info, "                         [%10.2f %10.2f %10.2f]]\n", g[2][0], */
/* 	  g[2][1], g[2][2]); */
    

  // Temporary arrays for netcdf reading
  tmpint = (int *)malloc(sizeof(int)*n_atom);
  tmpreal = (double *)malloc(sizeof(double)*n_atom);
  tmpreal3 = (double *)malloc(3*sizeof(double)*n_atom);
  tmpint3 = (int *)malloc(3*sizeof(int)*n_atom);

    
  // Species
  if (species_var_id != -1) {
    start3[0] = frame;
    start3[1] = 0;
    start3[2] = 0;
    count3[0] = 1;
    count3[1] = n_atom;
    count3[2] = n_label;
    species = (char *)malloc(n_atom*n_label*sizeof(char));
    netcdf_check(nc_get_vara_text(ncid, species_var_id, start3, count3, species));
    n = 0;
    for (i=0; i<n_atom; i++) {
      if (gotfilter && !atomfilter[i]) continue;
      //fprintf(info, "Config_load_netcdf: species[%d] = %c%c\n", i, species[i*n_label], species[i*n_label+1]);
      safe_symbol((species+i*n_label),SYMBOL(n));
      n++;
    }
    free(species);
  } else {
    // get species from Z
    if (Z_var_dims == 1) {
      // Z is defined globally, i.e. once in the nc file
      start1[0] = 0;
      count1[0] = n_atom;
      netcdf_check(nc_get_vara_int(ncid, Z_var_id, start1, count1, tmpint));
    } else if (Z_var_dims == 2) {
      // Z is defined per frame
      start2[0] = frame;
      start2[1] = 0;
      count2[0] = 1;
      count2[1] = n_atom;
      netcdf_check(nc_get_vara_int(ncid, Z_var_id, start2, count2, tmpint));
    } else {
      pe("Config_load_netcdf: 'Z' dims = %i (!= 1,2).", Z_var_dims);
    }

    n = 0;
    for (i=0; i<n_atom; i++) {
      if (gotfilter && !atomfilter[i]) continue;
      safe_symbol(ATOM_SYMBOL(tmpint[i]),SYMBOL(n));
      n++;
    }
  }

  // Mass
  if (mass_var_id != -1) {
    start2[0]=frame;
    start2[1]=0;
    count2[0]=1;
    count2[1]=n_atom;

    if (gotfilter) {
      netcdf_check(nc_get_vara_double(ncid, mass_var_id, start2, count2, tmpreal));
      n = 0;
      for (i=0; i<n_atom; i++) {
	if (!atomfilter[i]) continue;
	(*mass)[n] = tmpreal[i];
	n++;
      }
    }
    else
      netcdf_check(nc_get_vara_double(ncid, mass_var_id, start2, count2, *mass));

    // convert from libAtoms units to AtomEye units
    for (n=0; n<*np; n++) {
      //fprintf(info, "Config_load_netcdf: mass[%d] = %f\n", n, (*mass)[n]);
      (*mass)[n] /= (MASS_CONVERT*umass_IN_AMU);
    }
  }else {
    // Fill in masses from SYMBOL
    for (n=0; n<*np; n++)
      (*mass)[n] = ATOM_MASS_IN_AMU(Search_atom_by_symbol(SYMBOL(n)));
  }

  // Coordinates
  start3[0] = frame;
  start3[1] = 0;
  start3[2] = 0;
  count3[0] = 1;
  count3[1] = n_atom;
  count3[2] = 3;
  coords = (double *)malloc(n_atom*3*sizeof(double));
  netcdf_check(nc_get_vara_double(ncid, coordinates_var_id, start3, count3, coords));
  n = 0;
  for (i=0; i<n_atom; i++){
    if (gotfilter && !atomfilter[i]) continue;
    /*       fprintf(info, "Config_load_netcdf: x[%d] = [%10.2f %10.2f %10.2f]\n", */
    /* 	      i, coords[i*3+0],coords[i*3+1],coords[i*3+2]); */
    x[0] = coords[i*3]*coordinates_scale_factor;
    x[1] = coords[i*3+1]*coordinates_scale_factor;
    x[2] = coords[i*3+2]*coordinates_scale_factor;
    V3mulM3(x,g,sp);
    (*s)[3*n] = sp[0]; (*s)[3*n+1] = sp[1]; (*s)[3*n+2] = sp[2];
    // Fill in original_index property
    if (gotfilter) *(CONFIG_auxiliary[n_aux]+n) = (double)i;
    /*       fprintf(info, "Config_load_netcdf: s[%d] = [%10.2f %10.2f %10.2f]\n", */
    /* 	      i, (*s)[3*i],(*s)[3*i+1], (*s)[3*i+2]); */
    n++;
  }

  // Velocities
/*   if (velocity_var_id != -1) { */
/*     start3[0] = frame; */
/*     start3[1] = 0; */
/*     start3[2] = 0; */
/*     count3[0] = 1; */
/*     count3[1] = n_atom; */
/*     count3[2] = 3; */
/*     netcdf_check(nc_get_vara_double(ncid, velocity_var_id, start3, count3, coords)); */
/*     for (i=0; i<n_atom; i++) { */
/*       x[0] = coords[i*3]; x[1] = coords[i*3+1]; x[2] = coords[i*3+2]; */
/*       V3mulM3(x,g,sp); */
/*       (*s1)[3*i] = sp[0]; (*s1)[3*i+1] = sp[1]; (*s1)[3*i+2] = sp[2]; */
/*     } */
/*   } */

  free(coords);

  // Auxiliary properties 
  if (n_aux != 0) {
    fprintf(info, "Config_load_netcdf: loading %d aux props\n", n_aux);
    //    for (i=0; i<n_aux; i++) {
    i = 0;
    while (i < n_aux) {
      if (aux_prop_type[i] == NC_DOUBLE || aux_prop_type[i] == NC_FLOAT) {
	//fprintf(info, "%d %d\n", i, aux_prop_dim[i]);

	if (aux_prop_dim[i] == -1) {
	  //fprintf(info, "per frame definition");
	  start2[0] = frame;
	  start2[1] = 0;
	  count2[0] = 1;
	  count2[1] = n_atom;
	  netcdf_check(nc_get_vara_double(ncid, aux_prop_var_id[i], start2, count2, tmpreal));
	} else if (aux_prop_dim[i] == -2) {
	  //fprintf(info, "per file definition");
	  start1[0] = 0;
	  count1[0] = n_atom;
	  netcdf_check(nc_get_vara_double(ncid, aux_prop_var_id[i], start1, count1, tmpreal));
	} else {
	  //fprintf(info, "per frame definition (3-dim)");
	  start3[0] = frame;
	  start3[1] = 0;
	  start3[2] = 0;
	  count3[0] = 1;
	  count3[1] = n_atom;
	  count3[2] = 3;
	  netcdf_check(nc_get_vara_double(ncid, aux_prop_var_id[i], start3, count3, tmpreal3));
	}

	
	if (aux_prop_dim[i] < 0) {
	  n = 0;
	  for (j=0; j<n_atom; j++) {
	    if (!gotfilter || atomfilter[j]) {
	      *(CONFIG_auxiliary[i]+n) = tmpreal[j]*aux_scale_factors[i];
	      n++;
	    }
	  }
	  i++;
	} else {
	  n = 0;
	  for (j=0; j<n_atom; j++) {
	    if (!gotfilter || atomfilter[j]) {
	      *(CONFIG_auxiliary[i]+n)   = tmpreal3[3*j+0]*aux_scale_factors[i];
	      *(CONFIG_auxiliary[i+1]+n) = tmpreal3[3*j+1]*aux_scale_factors[i];
	      *(CONFIG_auxiliary[i+2]+n) = tmpreal3[3*j+2]*aux_scale_factors[i];
	      n++;
	    }
	  }
	  i += 3;
	}

      } else {
	//fprintf(info, "integer: %d %d\n", i, aux_prop_dim[i]);

	if (aux_prop_dim[i] == -1) {
	  start2[0] = frame;
	  start2[1] = 0;
	  count2[0] = 1;
	  count2[1] = n_atom;
	  netcdf_check(nc_get_vara_int(ncid, aux_prop_var_id[i], start2, count2, tmpint));
	} else if (aux_prop_dim[i] == -2) {
	  start1[0] = 0;
	  count1[0] = n_atom;
	  netcdf_check(nc_get_vara_int(ncid, aux_prop_var_id[i], start1, count1, tmpint));
	} else {
	  start3[0] = frame;
	  start3[1] = 0;
	  start3[2] = 0;
	  count3[0] = 1;
	  count3[1] = n_atom;
	  count3[2] = 3;
	  netcdf_check(nc_get_vara_int(ncid, aux_prop_var_id[i], start3, count3, tmpint3));
	}

	if (aux_prop_dim[i] < 0) {
	  n = 0;
	  for (j=0; j<n_atom; j++) {
	    if (!gotfilter || atomfilter[j]) {
	      *(CONFIG_auxiliary[i]+n) = (double)tmpint[j];
	      n++;
	    }
	  }
	  i++;
	} else {
	  n = 0;
	  for (j=0; j<n_atom; j++) {
	    if (!gotfilter || atomfilter[j]) {
	      *(CONFIG_auxiliary[i]+n)   = (double)tmpint3[3*j+0]*aux_scale_factors[i];
	      *(CONFIG_auxiliary[i+1]+n) = (double)tmpint3[3*j+1]*aux_scale_factors[i];
	      *(CONFIG_auxiliary[i+2]+n) = (double)tmpint3[3*j+2]*aux_scale_factors[i];
	      n++;
	    }
	  }
	  i += 3;	  
	}
      }
    }
  }

  if (gotfilter) free(atomfilter);
  free(tmpint);  free(tmpint3);
  free(tmpreal); free(tmpreal3);
  
  if (info != NULL) {
    Config_analyze_H       (Config_Alib_to_Alib, info); fcr(info);
    Config_analyze_density (Config_Alib_to_Alib, info); fcr(info);
    Config_analyze_species (Config_Alib_to_Alib, info); fcr(info);
  }
  
  return;
} /* end Config_load_netcdf() */


/* Load atomistic configuration from extended XYZ file */
void Config_load_xyz (char *fname, FILE *info, Alib_Declare_Config)
{
  int i,n, entry_count,j,k,ncols;
  char sym[3];
  char linebuffer[CFG_LINESIZE], buf1[CFG_LINESIZE], buf2[CFG_LINESIZE], buf3[CFG_LINESIZE];
  char format[3+4*(CONFIG_MAX_AUXILIARY+DIMENSION+DIMENSION)];
  char fields[3*CONFIG_MAX_AUXILIARY][TERMSIZE];
  char *p, *p1, *framestr, *nfname;
  long frames[XYZ_MAX_FRAMES];
  int  atoms[XYZ_MAX_FRAMES];
  FILE *in, *filter;
  M3 g;
  double x[3], sp[3];
  static int frame;
  int tmpi, inc, n_frame, mass_aux_index = -1, gotfilter;
  int *atomfilter, nfilter, nxyz;

  Fprintf(info, "Loading XYZ configuration \"%s\"\n", fname);

  strcpy(buf1, fname);
  framestr = buf1;
  nfname = strsep(&framestr, ":");
  if (framestr == NULL) {
    strcpy(buf2, "first");
    framestr = buf2;
  }
  n_frame = xyz_find_frames(nfname, frames, atoms);

  strcpy(buf3,nfname);
  strcat(buf3,".filter");
  gotfilter = (filter = fopen(buf3,"r")) != NULL;

  Config_free_auxiliary();
  
  if (sscanf(framestr, "%d", &frame) == 1) {
    // Check frame is in range
    if ((frame < 0) || (frame >= n_frame)) {
      pe("Config_load_xyz: frame should be in range 0 <= frame <= %d\n", n_frame-1);
    }
  }
  else if (strstr(framestr, "first")) {
    frame = 0;
  }
  else if (strstr(framestr, "last")) {
    frame = n_frame-1;
  }
  else if (strstr(framestr, "forward")) {
    if (sscanf(framestr, "forward:%d", &inc) != 1) inc = 1;
    if ((frame += inc) >= n_frame) frame = 0;
  } else if (strstr(framestr, "backward")) {
    if (sscanf(framestr, "backward:%d", &inc) != 1) inc = 1;
    if ((frame -= inc) <= -1) frame = n_frame-1;
  }
  else if (strstr(framestr, "reload")) {
    fprintf(info, "Config_load_xyz: reloaded file, remaining at frame %d/%d\n", frame, n_frame);
  }
  else {
    pe("Config_load_xyz: unknown framestr %s\n", framestr);
  }

  in = ROpen(nfname);

  fprintf(info, "Config_load_xyz: loading frame %d\n", frame);
  strcpy(fname, nfname);
  sprintf(buf2, ":%%0%dd", (int)ceil(log10(n_frame)));
  sprintf(buf1, buf2, frame);
  strcat(fname, buf1);

  fseek(in, frames[frame], SEEK_SET);

  if (!fgets(linebuffer,CFG_LINESIZE,in)) pe("Config_load_xyz: premature file ending - expecting number of atoms\n");
  if (sscanf(linebuffer, "%d", &nxyz) != 1)  pe("Config_load_xyz: first line must be number of atoms\n");

  // Read comment line, which contains 'Lattice=' and 'Properties=' keys
  if (!fgets(linebuffer,CFG_LINESIZE,in)) pe("Config_load_xyz: premature file ending - expecting comment line\n");
  fprintf(info, "Config_load_xyz: comment line %s\n", linebuffer);

  // Read lattice, fortran ordering
  p = strstr(linebuffer, "Lattice=\"")+9;
  if (p == NULL) pe("Config_load_xyz: missing lattice\n");
  strcpy(buf1, p);
  p = strstr(buf1, "\"");
  if (p == NULL) pe("Config_load_xyz: missing endquote in lattice\n");
  *p = '\0';

  if (sscanf(buf1, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", 
	     &H[0][0], &H[0][1], &H[0][2], 
	     &H[1][0], &H[1][1], &H[1][2],
	     &H[2][0], &H[2][1], &H[2][2]) != 9)
	     pe("Config_load_xyz: error reading lattice from string %s", buf1);
  M3inv (H,g);

  // Scan properties string
  p = strstr(linebuffer, "Properties=")+11;
  if (p == NULL) pe("Config_load_xyz: missing Properties key\n");
  if (p[0] == '\"') {
    p++;
    p1 = strstr(p, "\"");
    if (p1 != NULL) pe("Config_load_xyz: missing endquote in properties\n");
    *p1 = '\0';
  }
  else {
    p1 = strstr(p, " ");
    if (p1 != NULL) *p1 = '\0';
  }

  k = 0;
  while ((p1 = strsep(&p, ":")) != NULL) {
    strcpy(fields[k++], p1);
  }

  // Check first column is species and next three are position
  if (strcmp(fields[0], "species") != 0 || strcmp(fields[1], "S") != 0 || sscanf(fields[2], "%d", &tmpi) != 1 || tmpi != 1)
    pe("Config_load_xyz: first property must be species\n");
  if (strcmp(fields[3], "pos") != 0 || strcmp(fields[4], "R") != 0 || sscanf(fields[5], "%d", &tmpi) != 1 || tmpi != 3)
    pe("Config_load_xyz: second property must be pos %s %s %s\n",fields[3],fields[4],fields[5]);

  // Process rest of aux props
  entry_count = 0;
  for (i=2; i<k/3; i++) {
/*     fprintf(info, "Config_load_xyz: got property %s:%s:%s\n", fields[3*i], fields[3*i+1], fields[3*i+2]); */
    
    // Abort if we've got any non-numeric fields
    if (strcmp(fields[3*i+1], "S") == 0) pe("Config_load_xyz: string properties not supported\n");
    if (strcmp(fields[3*i+1], "L") == 0) pe("Config_load_xyz: logical properties not supported\n");

    if (sscanf(fields[3*i+2], "%d", &ncols) != 1) pe("Config_load_xyz: Bad column count %s", fields[3*i+2]);
    if (strcmp(fields[3*i], "mass") == 0 && strcmp(fields[3*i+1],"R") == 0 && ncols == 1)
      mass_aux_index = entry_count;

    if (ncols == 1) {
      if (entry_count == CONFIG_MAX_AUXILIARY) {
	fprintf(info,"Config_load: Warning CONFIG_MAX_AUXILIARY = %d exceeded, skipping remaining properties\n", 
		CONFIG_MAX_AUXILIARY);
	entry_count = CONFIG_MAX_AUXILIARY;
	goto XYZ_AUX_PROPS_DONE;
      }
      strcpy(CONFIG_auxiliary_name[entry_count], fields[3*i]);
      entry_count++;
    }
    else
      for (j=0; j<ncols; j++) {
	if (entry_count == CONFIG_MAX_AUXILIARY) {
	  fprintf(info,"Config_load: Warning CONFIG_MAX_AUXILIARY = %d exceeded, skipping remaining properties\n", 
		  CONFIG_MAX_AUXILIARY);
	  entry_count = CONFIG_MAX_AUXILIARY;
	  goto XYZ_AUX_PROPS_DONE;
	}
	sprintf(CONFIG_auxiliary_name[entry_count], "%s%d", fields[3*i], j);
	entry_count++;
      }
  }

 XYZ_AUX_PROPS_DONE:
  entry_count += 4;  // For symbol and positions
  fprintf(info, "Config_load_xyz: entry_count = %d\n", entry_count);

  // Build format string - everything apart from species is read as a double
  format[0] = '%'; 
  format[1]='s'; 
  format[2]=' ';
  for (i=1; i<entry_count; i++)
    {
      format[3+4*(i-1)]   = '%';
      format[3+4*(i-1)+1] = 'l';
      format[3+4*(i-1)+2] = 'f';
      format[3+4*(i-1)+3] = (i < entry_count-1) ? ' ' : EOS;
    }

  // Check for .filter file
  if (gotfilter) {
    atomfilter = (int *)malloc(nxyz*sizeof(int));
    for (i=0; i<nxyz; i++) atomfilter[i] = 0;
    nfilter = 0;
    while (fgets(linebuffer,CFG_LINESIZE,filter)) {
      if (sscanf(linebuffer, "%d", &i) != 1) 
	pe("Error reading atom list: %s", linebuffer);
      if(i < 0 || i >= nxyz) pe("Atom filter atom %d out of range 0 < i < %d\n", i, nxyz);
      atomfilter[i] = 1;
      nfilter++;
    }
    strcpy(CONFIG_auxiliary_name[entry_count-4],"original_index");
    strcpy(CONFIG_auxiliary_unit[entry_count-4],"");
    fclose(filter);
    *np = nfilter;
    fprintf(info, "Config_load_xyz: filter %s applied. Selected %d/%d atoms.\n", buf3, nfilter, nxyz);
  }
  else {
    // Include all atoms
    *np = nxyz;
  }

  // Reallocate to correct size
  fprintf(info, "Config_load_xyz: got %d atoms\n", *np);
  Config_realloc (Config_Alib_to_Alib);

  // Allocate correct number of aux props
  for (k=0; k<entry_count-4+gotfilter; k++)
    REALLOC (Config_load, CONFIG_auxiliary[k], *np, double);
  for (k=entry_count-4+gotfilter; k<CONFIG_num_auxiliary; k++)
    {
      Free (CONFIG_auxiliary[k]);
      CONFIG_auxiliary_name[k][0] = EOS;
      CONFIG_auxiliary_unit[k][0] = EOS;
    }
  CONFIG_num_auxiliary = entry_count-4+gotfilter;

  // Now it's just one line per atom
  n = 0;
  for (i=0; i < nxyz; i++) {
    if (!fgets(linebuffer,CFG_LINESIZE,in)) pe("Config_load_xyz: premature file ending at atom %d\n",n);
    if (gotfilter && !atomfilter[i]) continue;
    if (sscanf(linebuffer, format, 
	       sym, &x[0], &x[1], &x[2], 
	       CONFIG_auxiliary[0]+n,
	       CONFIG_auxiliary[1]+n,
	       CONFIG_auxiliary[2]+n,
	       CONFIG_auxiliary[3]+n,
	       CONFIG_auxiliary[4]+n,
	       CONFIG_auxiliary[5]+n,
	       CONFIG_auxiliary[6]+n,
	       CONFIG_auxiliary[7]+n,
	       CONFIG_auxiliary[8]+n,
	       CONFIG_auxiliary[9]+n,
	       CONFIG_auxiliary[10]+n,
	       CONFIG_auxiliary[11]+n,
	       CONFIG_auxiliary[12]+n,
	       CONFIG_auxiliary[13]+n,
	       CONFIG_auxiliary[14]+n,
	       CONFIG_auxiliary[15]+n,
	       CONFIG_auxiliary[16+0]+n,
	       CONFIG_auxiliary[16+1]+n,
	       CONFIG_auxiliary[16+2]+n,
	       CONFIG_auxiliary[16+3]+n,
	       CONFIG_auxiliary[16+4]+n,
	       CONFIG_auxiliary[16+5]+n,
	       CONFIG_auxiliary[16+6]+n,
	       CONFIG_auxiliary[16+7]+n,
	       CONFIG_auxiliary[16+8]+n,
	       CONFIG_auxiliary[16+9]+n,
	       CONFIG_auxiliary[16+10]+n,
	       CONFIG_auxiliary[16+11]+n,
	       CONFIG_auxiliary[16+12]+n,
	       CONFIG_auxiliary[16+13]+n,
	       CONFIG_auxiliary[16+14]+n,
	       CONFIG_auxiliary[16+15]+n,
	       CONFIG_auxiliary[32+0]+n,
	       CONFIG_auxiliary[32+1]+n,
	       CONFIG_auxiliary[32+2]+n,
	       CONFIG_auxiliary[32+3]+n,
	       CONFIG_auxiliary[32+4]+n,
	       CONFIG_auxiliary[32+5]+n,
	       CONFIG_auxiliary[32+6]+n,
	       CONFIG_auxiliary[32+7]+n,
	       CONFIG_auxiliary[32+8]+n,
	       CONFIG_auxiliary[32+9]+n,
	       CONFIG_auxiliary[32+10]+n,
	       CONFIG_auxiliary[32+11]+n,
	       CONFIG_auxiliary[32+12]+n,
	       CONFIG_auxiliary[32+13]+n,
	       CONFIG_auxiliary[32+14]+n,
	       CONFIG_auxiliary[32+15]+n) != entry_count)
      pe("Error reading line %s", linebuffer);
    V3ZERO (&((*s1)[DIMENSION*n])); // No velocities
    V3mulM3(x,g,sp);                // Convert to fractional coordinates
    (*s)[3*n] = sp[0]; (*s)[3*n+1] = sp[1]; (*s)[3*n+2] = sp[2];
    safe_symbol(sym,SYMBOL(n));
    if (gotfilter) *(CONFIG_auxiliary[CONFIG_num_auxiliary-1]+n) = (double)i; // original_index
    n++;
  }

  // Set the mass, either from symbol or from mass aux prop if there is one
  if (mass_aux_index == -1) {
    for (i=0; i<*np; i++)
      (*mass)[i] = ATOM_MASS_IN_AMU(Search_atom_by_symbol(SYMBOL(i)));
  }
  else {
    // convert from libAtoms units to AtomEye units
    for (i=0; i<*np; i++)
      (*mass)[i] = *(CONFIG_auxiliary[mass_aux_index]+i)/(MASS_CONVERT*umass_IN_AMU);
  }

  if (info != NULL)
    {
      Config_analyze_H       (Config_Alib_to_Alib, info); fcr(info);
      Config_analyze_density (Config_Alib_to_Alib, info); fcr(info);
      Config_analyze_species (Config_Alib_to_Alib, info); fcr(info);
    }

  fclose(in);
  if (gotfilter) free(atomfilter);
  return;
} /* end Config_load_xyz() */



#ifdef _Config_load_TEST
#define ULENGTH_IN_A  2.
/* #define LOAD_FILE   "/tmp/b" */
/* bzip2 -c /tmp/b > /tmp/b.bz2 */
/* #define LOAD_FILE   "/tmp/b.bz2" */
#define LOAD_FILE "/asm/home/liju99/Ar/Data/config.amorphous"
#define PDB_FILE    "/tmp/c.pdb"
#define CONFIG_FILE "/tmp/c"
int main (int argc, char *argv[])
{
    Aapp_Define_Config;
    Config_LOAD(LOAD_FILE, Config_Aapp_to_Alib);
    Config_SAVE (Config_Aapp_to_Alib, CONFIG_FILE);
    printf ("config saved on \"%s\".\n", CONFIG_FILE);
    Config_save_as_PDB (Config_Aapp_to_Alib, PDB_FILE);
    printf ("config saved on \"%s\".\n", PDB_FILE);
    return (0);
}
#endif /* _Config_load_TEST */


/* Try to guess file format from filename suffix and   */
/* then load configuration using appropriate function. */
/* Return the file format loaded.                      */
int Config_Load (char *fname, FILE *info, Alib_Declare_Config)
{
    char *p, *q;
    Fprintf(info, "Guessing file format of \"%s\":\n",
            fname);
    p = strrchr(fname, '/');
    if (p == NULL) p = fname;
    q = strchr(p, '.');
    if (q != NULL)
    {
        p = strrstr(q, "cfg");
        if (p != NULL)
        {
            if (strstr(p, "pdb") || strstr(p, "PDB"))
	      goto PDB;
	    if (strstr(p, "nc") || strstr(p, "NC"))
	      goto NETCDF;
	    if (strstr(p, "xyz") || strstr(p, "XYZ"))
	      goto XYZ;
        }
        else {
	  if (strstr(q, "pdb") || strstr(q, "PDB"))
            goto PDB;
	  if (strstr(q, "nc") || strstr(q, "NC"))
	    goto NETCDF;
	  if (strstr(q, "xyz") || strstr(q, "XYZ"))
	    goto XYZ;
	}
    }
    Fprintf(info, "should be Ju Li's CFG format.\n\n");
    Config_load(fname, info, Config_Alib_to_Alib);
    return(CONFIG_CFG_LOADED);
  PDB:
    Fprintf(info, "should be Protein Data Bank format.\n\n");
    Config_load_from_pdb(fname, info, Config_Alib_to_Alib);
    return(CONFIG_PDB_LOADED);
/*  NETCDF: */
/*     Fprintf(info, "should be NetCDF trajectory.\n\n"); */
/*     Config_load_netcdf(fname, info, Config_Alib_to_Alib); */
/*     return(CONFIG_CFG_LOADED); // We want folding into primitive cell, not bounding box */
/*  XYZ: */
/*     Fprintf(info, "should be XYZ format\n\n"); */
/*     Config_load_xyz(fname, info, Config_Alib_to_Alib); */
/*     return(CONFIG_CFG_LOADED); // We want folding into primitive cell, not bounding box */
 NETCDF:
 XYZ:
    Fprintf(info, "calling Config_load_libatoms_filename\n\n");
    Config_load_libatoms_filename(fname, info, Config_Alib_to_Alib);
    return(CONFIG_CFG_LOADED); // We want folding into primitive cell, not bounding box
    
    
} /* end Config_Load() */

Atoms *config_libatoms;

/* Load atomistic configuration from libAtoms-compatible Atoms C structure */
void Config_load_libatoms (Atoms *atoms, FILE *info, Alib_Declare_Config)
{
  int i,n, entry_count,j,k, species_idx, pos_idx;
  M3 g;
  double x[3], sp[3];
  int naux, mass_aux_index = -1;

  Config_free_auxiliary();
  
  *np = atoms->n_atom;

  for (i=0; i<3; i++)
    for (j=0; j<3; j++)
      H[i][j] = atoms->lattice[j][i];
  M3inv (H,g);

  species_idx = atoms_find_property(atoms, "species");
  if (species_idx == -1) pe("No species in atoms");
  pos_idx = atoms_find_property(atoms, "pos");
  if (pos_idx == -1) pe("No pos in atoms");

  // Process rest of aux props
  entry_count = 0;
  for (i=0; i < atoms->n_property; i++) {
    if (i == species_idx || i == pos_idx) continue;

    if (atoms->property_type[i] == PROPERTY_STR) {
      Fprintf(info, "warning: ignoring string property %s\n", atoms->property_name[i]);
      continue;
    }

    if (strcmp(atoms->property_name[i], "mass") == 0 && 
	atoms->property_type[i] == PROPERTY_REAL && 
	atoms->property_ncols[i] == 1)
      mass_aux_index = entry_count;

    if (atoms->property_ncols[i] == 1) {
      strcpy(CONFIG_auxiliary_name[entry_count], atoms->property_name[i]);
      entry_count++;
    }
    else
      for (j=0; j<atoms->property_ncols[i]; j++) {
	sprintf(CONFIG_auxiliary_name[entry_count], "%s%d", atoms->property_name[i], j);
	entry_count++;
      }
  }

  entry_count += 4;  // For symbol and positions
  Fprintf(info, "Config_load_libatoms: entry_count = %d\n", entry_count);

  // Reallocate to correct size
  Fprintf(info, "Config_load_libatoms: got %d atoms\n", *np);
  Config_realloc (Config_Alib_to_Alib);

  // Allocate correct number of aux props
  for (k=0; k<entry_count-4; k++)
    REALLOC (Config_load, CONFIG_auxiliary[k], *np, double);
  for (k=entry_count-4; k<CONFIG_num_auxiliary; k++)
    {
      Free (CONFIG_auxiliary[k]);
      CONFIG_auxiliary_name[k][0] = EOS;
      CONFIG_auxiliary_unit[k][0] = EOS;
    }
  CONFIG_num_auxiliary = entry_count-4;

  for (n=0; n < *np; n++) {

    V3ZERO (&((*s1)[DIMENSION*n])); // No velocities
    for (i=0; i<3; i++) {
      x[i] = property_real(atoms, pos_idx, i, n);
    }
    V3mulM3(x,g,sp);                // Convert to fractional coordinates
    (*s)[3*n] = sp[0]; (*s)[3*n+1] = sp[1]; (*s)[3*n+2] = sp[2];
    safe_symbol(&property_str(atoms, species_idx, 0, n),SYMBOL(n));

    // Auxiliary properties
    naux = 0;
    for (i=0; i<atoms->n_property; i++) {
      if (i == species_idx || i == pos_idx) continue;

      for (j=0; j < atoms->property_ncols[i]; j++) {
	//printf("%s[%d,%d] type %d\n", atoms->property_name[i], n, j, atoms->property_type[i]);
	switch (atoms->property_type[i]) {
	case(PROPERTY_INT):
	  *(CONFIG_auxiliary[naux]+n) = (double)property_int(atoms, i, j, n);
	  break;
	
	case(PROPERTY_REAL):
	  *(CONFIG_auxiliary[naux]+n) = property_real(atoms, i, j, n);
	  break;

	case(PROPERTY_LOGICAL):
	  *(CONFIG_auxiliary[naux]+n) = (double)property_logical(atoms, i, j, n);
	  break;

	case(PROPERTY_STR):
	  break;
	}
      }
      naux += atoms->property_ncols[i];
    }
  }

  // Set the mass, either from symbol or from mass aux prop if there is one
  if (mass_aux_index == -1) {
    for (i=0; i<*np; i++)
      (*mass)[i] = ATOM_MASS_IN_AMU(Search_atom_by_symbol(SYMBOL(i)));
  }
  else {
    // convert from libAtoms units to AtomEye units
    for (i=0; i<*np; i++)
      (*mass)[i] = *(CONFIG_auxiliary[mass_aux_index]+i)/(MASS_CONVERT*umass_IN_AMU);
  }

  return;
} /* end Config_load_libatoms() */

void Config_load_libatoms_filename(char *fname, FILE *info, Alib_Declare_Config)
{
  char *p, *q, *framestr, *nfname, buf1[CFG_LINESIZE], buf2[CFG_LINESIZE], 
    buf3[CFG_LINESIZE], linebuffer[CFG_LINESIZE];
  int netcdf, xyz, gotfilter;
  int n_frame, i;
  static int frame;
  Atoms at;
  int *atomlist, natomlist;  
  int offset, nc_in;
  FILE *fxyz, *atomfile;
  int retval, inc;

  atoms_init(&at);

  xyz = 0;
  netcdf = 0;
  offset = 0;
    
  p = strrchr(fname, '/');
  if (p == NULL) p = fname;
  q = strchr(p, '.');
  if (q != NULL) {
    if (strstr(q, "nc") || strstr(q, "NC"))
      netcdf = 1;
    if (strstr(q, "xyz") || strstr(q, "XYZ"))
      xyz = 1;
  }
  if (!xyz && !netcdf) pe("Don't know how to read from file %s", fname);

  strcpy(buf1, fname);
  framestr = buf1;
  nfname = strsep(&framestr, ":");
  if (framestr == NULL) {
    strcpy(buf2, "first");
    framestr = buf2;
  }

  atomlist = NULL;
  natomlist = 0;

  strcpy(buf3,nfname);
  strcat(buf3,".filter");
  gotfilter = ((atomfile = fopen(buf3,"r")) != NULL);

  if (gotfilter) {
    natomlist = 0;  // Count lines in atomlist file
    while (fgets(linebuffer,LINESIZE,atomfile)) natomlist++;
      
    atomlist = malloc(natomlist*sizeof(int));
    if (atomlist == NULL) pe("Error allocating atomlist");
    fseek(atomfile, 0, SEEK_SET);
    for (i=0; i<natomlist; i++) {
      if (!fgets(linebuffer,LINESIZE,atomfile)) pe("Premature end of atom file");
      if (sscanf(linebuffer, "%d", &atomlist[i]) != 1) pe("Error reading line %d of atom list: %s", i+1, linebuffer);
      atomlist[i] -= offset;
    }
    fclose(atomfile);
  }

  if (xyz) {
    n_frame = xyz_find_frames(nfname, at.frames, at.atoms);
    if (n_frame == 0) pe("Error building frame index");
    Fprintf(stderr, "got xyz file with %d frames\n", n_frame);
    at.got_index = 1;
  }
  else {
    Fprintf(info, "Opening netcdf file %s\n", nfname);
#ifdef NETCDF4
    netcdf_check(nc_open(nfname, NC_NOWRITE, &nc_in));
#else
    netcdf_check(nc_open(nfname, NC_64BIT_OFFSET | NC_NOWRITE, &nc_in));
#endif
    n_frame = read_netcdf(nc_in, &at, 0, atomlist, natomlist, 1, 0, 1, 1, 0, 0.0);
    Fprintf(stderr, "got netcdf file with %d frames\n", n_frame);
    if (n_frame == 0) pe("Error reading netcdf");
    at.got_index = 1;
  }

  if (sscanf(framestr, "%d", &frame) == 1) {
    // Check frame is in range
    if ((frame < 0) || (frame >= n_frame)) {
      pe("Config_load_libatoms_filename: frame should be in range 0 <= frame <= %d\n", n_frame-1);
    }
  }
  else if (strstr(framestr, "first")) {
    frame = 0;
  }
  else if (strstr(framestr, "last")) {
    frame = n_frame-1;
  }
  else if (strstr(framestr, "forward")) {
    if (sscanf(framestr, "forward:%d", &inc) != 1) inc = 1;
    if ((frame += inc) >= n_frame) frame = 0;
  } else if (strstr(framestr, "backward")) {
    if (sscanf(framestr, "backward:%d", &inc) != 1) inc = 1;
    if ((frame -= inc) <= -1) frame = n_frame-1;
  }
  else if (strstr(framestr, "reload")) {
    Fprintf(info, "Config_load_libatoms_filename: reloaded file, remaining at frame %d/%d\n", frame, n_frame);
  }
  else {
    pe("Config_load_libatoms_filename: unknown framestr %s\n", framestr);
  }

  if (xyz) {
    if ((fxyz = fopen(nfname, "r")) == NULL) pe("Cannot open xyz file %s.\n", nfname);
    if (!read_xyz(fxyz, &at, atomlist, natomlist, frame, 0, 0, 1, 0)) 
      pe("Error reading frame %d from xyz file %s.\n", frame, nfname);
    fclose(fxyz);
  }
  else {
    if (!read_netcdf(nc_in, &at, frame, atomlist, natomlist, 0, 0, 1, 1, 0, 0.0))
      pe("Error reading frame %d from netcdf file %s.\n", frame, nfname);
    nc_close(nc_in);
  }

  Config_load_libatoms(&at, info, Config_Alib_to_Alib);
  atoms_free(&at);

  strcpy(fname, nfname);
  sprintf(buf2, ":%%0%dd", (int)ceil(log10(n_frame)));
  sprintf(buf1, buf2, frame);
  strcat(fname, buf1);
}


#ifdef _CONFIG_TEST
#define CONFIG_FILE "/tmp/aa"
#define CONFIG_FILE_BZ2 "/tmp/b.bz2"
#define CONFIG_FILE_GZ "/tmp/c.gz"
int main (int argc, char *argv[])
{
    int i;
    char *si = "Si";
    Aapp_Define_Config;
    Chemtab ct[1]={{0}};
    Tp *tp=NULL;
    Config_rebuild_Xtal (Config_Aapp_to_Alib, 1,1,1, Xtal_abstracts+XTAL_FCC,
                         " C", -1., -1.);
    for (i=0; i<np; i++)
    {
        s[DIMENSION*i] += 0.0321321434535176746;
        s[DIMENSION*i+1] += 0.0265365656576576324;
    }
    Config_motion_initialize (100., Config_Aapp_to_Alib);
    for (i=0; i<100; i++)
    {
        printf ("%d\n", i);
        CONFIG_NV_SAVE
            (Config_Aapp_to_Alib, CONFIG_FILE, "%.7g %.7g %.7g", 2,
             "s1 []", " %.7g", (char *)s,     sizeof(double)*DIMENSION,
             "s2 []", " %.7g", (char *)(s+1), sizeof(double)*DIMENSION);
    }
    printf ("** CONFIG saved on \"%s\" **\n\n", CONFIG_FILE);
    /* cat /tmp/a */
    Config_Load (CONFIG_FILE, stdout, Config_Aapp_to_Alib);
    Config_analyze_motion (Config_Aapp_to_Alib, stdout); cr();
    Config_motion_initialize (50., Config_Aapp_to_Alib);
    for (i=0; i<100; i++)
    {
        printf ("%d\n", i);
        CONFIG_SAVE
            (Config_Aapp_to_Alib, CONFIG_FILE_BZ2, "%.7g %.7g %.7g",
             " %.7g %.7g %.7g", 2,
             "s1 []", " %.7g", (char *)s,     sizeof(double)*DIMENSION,
             "s2 []", " %.7g", (char *)(s+1), sizeof(double)*DIMENSION);
    }
    printf ("** CONFIG saved on \"%s\" **\n\n", CONFIG_FILE_BZ2);
    /* bzcat /tmp/b.bz2 */
    Config_Load (CONFIG_FILE_BZ2, stdout, Config_Aapp_to_Alib);
    Config_analyze_motion (Config_Aapp_to_Alib, stdout); cr();
    /* for (i=0; i<np; i++) */
    /* printf ("%f %f\n", CONFIG_auxiliary[0][i], CONFIG_auxiliary[1][i]); */
    Config_free_auxiliary();
    safe_symbol(si,SYM(0)); mass[0] = 54. / UMASS_IN_AMU;
    safe_symbol(si,SYM(3)); mass[3] = 53. / UMASS_IN_AMU;
#ifndef ATOMEYE_LIB
    rebind_CT (Config_Aapp_to_Alib, "Si C", ct, &tp); cr();
#else
    rebind_ct (Config_Aapp_to_Alib, "Si C", ct, &tp, NULL); cr();
#endif
    Config_compact_index (ct, &tp, Config_Aapp_to_Alib);
    CONFIG_SAVE
        (Config_Aapp_to_Alib, CONFIG_FILE_GZ, "%.7g %.7g %.7g",
         " %.7g %.7g %.7g", 2,
         "s1 []", " %.7g", (char *)s,     sizeof(double)*DIMENSION,
         "s2 []", " %.7g", (char *)(s+1), sizeof(double)*DIMENSION);
    printf ("** CONFIG saved on \"%s\" **\n\n", CONFIG_FILE_GZ);
    /* gzip -c -d /tmp/c.gz */
    Config_Load (CONFIG_FILE_GZ, stdout, Config_Aapp_to_Alib);
    Config_analyze_motion (Config_Aapp_to_Alib, stdout); cr();
    return (0);
}
#endif /* _CONFIG_TEST */


/* Determine crystal structure expanded in  */
/* parallelepiped HN[][] with respect to x0 */
#ifdef _xtal
int main (int argc, char *argv[])
{
    int i,j,k,nc[3],*idx;
    double cc,x0[3],s0[3],*value;
    char *taglist;
    M3 HI,HN;
    Aapp_Define_Config;
    Crystallographic X;
    TermString buf;

    fscanf (stdin, "Old lattice:\n");
    for (i=0; i<3; i++)
        fscanf (stdin, "%lf %lf %lf\n", &H[i][0], &H[i][1], &H[i][2]);
    fscanf (stdin, "Number of atoms in old unit cell:\n");
    fscanf (stdin, "%d\n", &np);
    Config_alloc (Config_Aapp_to_Alib);
    fscanf (stdin, "symbol          s0              s1              s2\n");
    for (i=0; i<np; i++)
    {
        *(SYM(i))   = getc(stdin);
        *(SYM(i)+1) = getc(stdin);
        SAFE_SYMBOL(SYM(i));
        mass[i] = ATOM_MASS_IN_AMU(Search_atom_by_symbol(SYM(i)));
        fscanf (stdin, "%lf %lf %lf\n", &s[3*i], &s[3*i+1], &s[3*i+2]);
        V3ZERO (s1+3*i);
    }
    
    fscanf (stdin, "\nNew origin:\n");
    fscanf (stdin, "%lf %lf %lf\n", x0, x0+1, x0+2);
    M3INV (H,HI,cc);
    V3mM3 (x0,HI,s0);
    V3ImagE (s0);
    fscanf (stdin, "New lattice:\n");
    for (i=0; i<3; i++)
        fscanf (stdin, "%lf %lf %lf\n", &HN[i][0], &HN[i][1], &HN[i][2]);
    M3rows_to_cover_sphere (H, M3maxrowradius(HN), nc);
    V3aDd (nc,2);
    Config_multiply (2*nc[0]+1,2*nc[1]+1,2*nc[2]+1,Config_Aapp_to_Alib);

    s0[0] = (s0[0]+nc[0])/(2*nc[0]+1);
    s0[1] = (s0[1]+nc[1])/(2*nc[1]+1);
    s0[2] = (s0[2]+nc[2])/(2*nc[2]+1);
    Config_translate (-s0[0],-s0[1],-s0[2],Config_Aapp_to_Alib);
    Config_transplant (Config_Aapp_to_Alib,0.,0.,0.,HN);

    MALLOC(xtal, taglist, np, char);
    MALLOC(xtal, value, np, double);
    MALLOC(xtal, idx, np, int);
    for (j=i=0; i<np; i++)
        if ( (s[3*i]  <-TINY) || (s[3*i]  >1-TINY) ||
             (s[3*i+1]<-TINY) || (s[3*i+1]>1-TINY) ||
             (s[3*i+2]<-TINY) || (s[3*i+2]>1-TINY) )
            taglist[i] = 1;
        else
        {
            taglist[i] = 0;
            idx[j] = i;
            /* y has 1st priority because it is used for layering */
            value[i] = s[3*i+1] * 1e6 + s[3*i+2] * 1e3 + s[3*i];
            j++;
        }
    qsort_numerical_recipes (j,value,idx,USE_OLD_IDX);
    printf ("Number of atoms in the new unit cell = %d\n", j);
    for (i=0; i<j; i++)
    {
        k = idx[i];
        printf ("%2s %16.12f %16.12f %16.12f\n",
                SYM(k), s[3*k], s[3*k+1], s[3*k+2]);
    }
    Config_rid_of_taglist (taglist, Config_Aapp_to_Alib);
    free(taglist);
    free(value);
    free(idx);

    fscanf (stdin, "\nView with multiplicities and lattice constant [A]:\n");
    fscanf (stdin, "%d %d %d %lf\n", nc, nc+1, nc+2, &cc);
    M3MultiplY (cc,H);
    Config_multiply (nc[0],nc[1],nc[2],Config_Aapp_to_Alib);
    fscanf (stdin, "Save configuration on file (default=/tmp/config):\n%s\n",
            buf);
    if (!strcasecmp(buf,"default")) strcpy (buf, "/tmp/config");
    H_to_Crystallographic(H,X);
    Crystallographic_to_H(X,H);
    Config_SET_CM(Config_Aapp_to_Alib);
    Config_SAVE (Config_Aapp_to_Alib,buf);
    printf ("saved on \"%s\".\n", buf);
    Config_free(Config_Aapp_to_Alib);
    return (0);
}
#endif /* _xtal */


/* change s[] to its own image in [0,1)^3 */
void Config_TRIM (Alib_Declare_Config)
{
    VTRIM( DIMENSION*(*np), *s );
    return;
} /* end Config_TRIM() */


/* Create complete chemical disorder but maintaining the same kinetic energy */
void Config_chemical_randomize (Alib_Declare_Config)
{
    register int i, j;
    double *invariant = NULL;
    MALLOC (Config_chemical_randomize, invariant, DIMENSION*(*np), double);
    for (i=(*np); i--;)
        V3MUL ( sqrt((*mass)[i]), *s1+DIMENSION*i, invariant+DIMENSION*i );
    for (i=0; i<*np-1; i++)
    { /* randomly select a member in i..max */
        j = Fran(i, *np-1);
        swapd ( *mass+i, *mass+j, 1 );
        swapb ( *symbol+i*SYMBOL_SIZE, *symbol+j*SYMBOL_SIZE, SYMBOL_SIZE );
    }
    for (i=(*np); i--;)
        V3DIV ( invariant+DIMENSION*i, sqrt((*mass)[i]), *s1+DIMENSION*i );
    free (invariant);
    return;
} /* end Config_chemical_randomize() */


#ifdef _chemical_randomize
int main (int argc, char *argv[])
{
    Aapp_Define_Config;
    char *out_fname;
    TimeRandomize();
    if (argc == 3) out_fname = argv[2];
    else if (argc == 2) out_fname = argv[1];
    else
    {
        printf ("\nPurpose: Create complete chemical disorder \n"
                "         while maintaining the same kinetic energy.\n\n");
        printf ("Usage: %s in_fname (out_fname=in_fname)\n"
                "       %s in_fname out_fname\n\n", argv[0], argv[0]);
        return (1);
    }

    printf ("Loading \"%s\"...\n\n", argv[1]);
    CONFIG_LOAD (argv[1], Config_Aapp_to_Alib);
    Config_chemical_randomize (Config_Aapp_to_Alib);
    Config_save (Config_Aapp_to_Alib, TRUE, out_fname);
    printf ("\"%s\" -> \"%s\".\n", argv[1], out_fname);
    return (0);
}
#endif /* _chemical_randomize */


#ifdef _affine_deform
int main (int argc, char *argv[])
{
    Aapp_Define_Config;
    char *out_fname;
    M3 J, tmp;
    
    if (argc == 9+3) out_fname = argv[9+2];
    else if (argc == 9+2) out_fname = argv[1];
    else
    {
        printf ("\nPurpose: Hnew[][] = Hold[][] * J[][]\n\n");
        printf ("Usage: %s in_fname J11 J12 J13 J21 J22 J23 J31 J32 J33\n"
                "                                  (out_fname=in_fname)\n"
                "       %s in_fname J11 J12 J13 J21 J22 J23 J31 J32 J33 "
                "out_fname\n\n", argv[0], argv[0]);
        return (1);
    }

    printf ("Loading \"%s\"...\n\n", argv[1]);
    CONFIG_LOAD (argv[1], Config_Aapp_to_Alib);

    J[0][0] = atof(argv[2]);
    J[0][1] = atof(argv[3]);
    J[0][2] = atof(argv[4]);

    J[1][0] = atof(argv[5]);
    J[1][1] = atof(argv[6]);
    J[1][2] = atof(argv[7]);

    J[2][0] = atof(argv[8]);
    J[2][1] = atof(argv[9]);
    J[2][2] = atof(argv[10]);

    M3MUl (H, J, tmp);

    Config_save (Config_Aapp_to_Alib, TRUE, out_fname);
    printf ("\"%s\" -> \"%s\".\n", argv[1], out_fname);
    return (0);
}
#endif /* _affine_deform */


#ifdef _infectH
int main (int argc, char *argv[])
{
    M3 H0;
    Aapp_Define_Config;
    TermString in_file, out_file;
    if (argc < 3)
    {
        printf ("Purpose: make out_file have the same H[][] as in_file but\n"
                "         keeping its internal coordinates\n");
        printf ("Usage: %s in_file out_file\n", argv[0]);
        return (1);
    }
    strcpy (in_file, argv[1]);
    strcpy (out_file, argv[2]);
    ft = stdout;
    Fprintf (ft, "Loading \"%s\"...\n\n", in_file);
    Config_Load (in_file, ft, Config_Aapp_to_Alib);
    M3EQV (H, H0);
    Fprintf (ft, "Loading \"%s\"...\n\n", out_file);
    Config_Load (out_file, ft, Config_Aapp_to_Alib);
    M3EQV (H0, H);
    CONFIG_NV_SAVE (Config_Aapp_to_Alib, out_file, "%.15g %.15g %.15g", 0);
    S3fPRmul(ft, "H = %M A\n", H0, ULENGTH_IN_A);
    printf ("-> \"%s\".\n", out_file);
    return (0);
}
#endif /* _infectH */


#ifdef _concat
int main (int argc, char *argv[])
{
    int i, dir, np2;
    double ratio;
    Aapp_Define_Config;
    TermString in_file_1, in_file_2, out_file;
    ConfigStack cs[1];
    M3 H1, H2;
    
    if (argc < 5)
    {
        printf ("Purpose: concatenate two configurations in "
                "one direction (select from 0,1,2)\n");
        printf ("Usage: %s cfg_file_1 cfg_file_2 direction out_file\n",
                argv[0]);
        return (1);
    }
    strcpy (in_file_1, argv[1]);
    strcpy (in_file_2, argv[2]);
    dir = atoi(argv[3]);
    if ( OUW(dir,DIMENSION) )
        pe ("Must select concatenation direction (now=%d) among 0,1,2\n", dir);
    strcpy (out_file, argv[4]);

    ft = stdout;
    Fprintf (ft, "Loading \"%s\"...\n\n", in_file_2);
    Config_Load (in_file_2, ft, Config_Aapp_to_Alib);
    M3EQV (H, H2);
    np2 = np;
    Config_push (Config_Aapp_to_Alib, cs);

    Config_alloc (Config_Aapp_to_Alib);
    Fprintf (ft, "Loading \"%s\"...\n\n", in_file_1);
    Config_Load (in_file_1, ft, Config_Aapp_to_Alib);
    M3EQV (H, H1);

    Config_cat (cs, Config_Aapp_to_Alib);

    ratio = DOUBLE(np2) / np;
    for (i=DIMENSION; i--;)
    {
        if ( i==dir )
            V3ADD (H1[i], H2[i], H[i]);
        else
            V3ADDMULMUL (1-ratio, H1[i], ratio, H2[i], H[i]);
    }

    ratio = V3DOT( H1[dir], H[dir] ) / V3DOT( H[dir], H[dir] );
    for (i=np; i--;)
        if ( i < np-np2 ) s[3*i+dir] *= ratio;
        else s[3*i+dir] = ratio + s[3*i+dir]*(1-ratio);
    
    CONFIG_NV_SAVE (Config_Aapp_to_Alib, out_file, "%.15g %.15g %.15g", 0);
    S3fPRmul(ft, "new H = %M A\n", H, ULENGTH_IN_A);
    printf ("%d + %d = %d atoms -> \"%s\".\n", np-np2, np2, np, out_file);

    return (0);
}
#endif /* _concat */
