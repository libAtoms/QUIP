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


/* Vienna Ab-initio Simulation Package configuration toolkit */

/* Load atomistic configuration from POTCAR and POSCAR files */
void Config_load_from_VASP
(char *POTCAR_fname, char *POSCAR_fname, FILE *info, Alib_Declare_Config)
{
    int i, j, k, num_species, species_numatoms[MENDELEYEV_MAX];
    char *p=NULL, *q=NULL, buf[FGETS_MAXSIZE], name[SYMBOL_SIZE];
    bool has_dynamics, is_cartesian;
    double a0_in_A, tmp[DIMENSION], HI[DIMENSION][DIMENSION];
    FILE *in;
    Fprintf(info, "Loading configuration from file \"%s\":\n", POSCAR_fname);
    in = ROpen( POSCAR_fname );
    FGets(buf, in);
    /* p = blank_advance(buf); */
    /* q = nonblank_advance(p); */
    /* *q = EOS; */
    Fprintf(info, "configuration name = \"%s\"\n", buf);
    FGets(buf, in);
    sscanf(buf, "%lf\n", &a0_in_A);
    FGets(buf, in);
    sscanf(buf, "%lf %lf %lf\n", &H[0][0], &H[0][1], &H[0][2]);
    FGets(buf, in);
    sscanf(buf, "%lf %lf %lf\n", &H[1][0], &H[1][1], &H[1][2]);
    FGets(buf, in);
    sscanf(buf, "%lf %lf %lf\n", &H[2][0], &H[2][1], &H[2][2]);
    M3MultiplY (a0_in_A / ulength_IN_A, H);
    M3INV (H, HI, tmp[0]);
    FGets(buf, in);
    *np = 0;
    num_species = 0;
    p = blank_advance(buf);
    while (1)
    {
        if ( (*p == EOS) || (*p == '!') ) break;
        q = nonblank_advance(p);
        sscanf(p, "%d", &species_numatoms[num_species]);
        *np += species_numatoms[num_species];
        p=blank_advance(q);
        num_species++;
    }
    Fprintf(info, "num_species = %d, np = %d;\n", num_species, *np);
    Config_realloc (Config_Alib_to_Alib);
    FGets(buf, in);
    has_dynamics = (strcasestr(buf, "dynamics")!=NULL);
    if (has_dynamics) FGets(buf, in);
    is_cartesian = (strcasestr(buf, "Cartesian")!=NULL);
    for (i=0; i<*np; i++)
    {
        FGets(buf, in);
        sscanf(buf, "%lf %lf %lf\n", &(*s)[DIMENSION*i],
               &(*s)[DIMENSION*i+1], &(*s)[DIMENSION*i+2]);
        if ( is_cartesian )
            V3MM3MUL ( &(*s)[DIMENSION*i], HI, a0_in_A/ulength_IN_A, tmp );
    }
    if (has_dynamics)
    {
        is_cartesian = (strcasestr(buf, "Cartesian")!=NULL);
        for (i=0; i<*np; i++)
        {
            FGets(buf, in);
            sscanf(buf, "%lf %lf %lf\n", &(*s1)[DIMENSION*i],
                   &(*s1)[DIMENSION*i+1], &(*s1)[DIMENSION*i+2]);
            if (is_cartesian )
            {
                V3MuL ( utime_IN_FS / ulength_IN_A, &(*s1)[DIMENSION*i] );
                V3MM3 ( &(*s1)[DIMENSION*i], HI, tmp );
            }
        }
    }
    fclose(in);
    Fprintf(info, "Loading species designation from file \"%s\":\n",
            POTCAR_fname);
    in = ROpen( POTCAR_fname );
    for (i=k=0; i<num_species; i++)
    {
        FGets(buf, in);
        safe_symbol(blank_advance(nonblank_advance(blank_advance(buf))),
                    name);
        do
        {
            FGets(buf, in);
        } while ((q=strstr(buf, "POMASS"))==NULL);
        sscanf( q, "POMASS =%lf", &tmp[0] );
        tmp[0] /= umass_IN_AMU;
        for (j=0; j<species_numatoms[i]; j++,k++)
        {
            safe_symbol(name, SYMBOL(k));
            (*mass)[k] = tmp[0];
        }
        Fprintf(info, "%d %c%c atoms assigned;\n", species_numatoms[i],
                name[0], name[1]);
        do
        {
            FGets(buf, in);
        } while ((q=strstr(buf, "End of Dataset"))==NULL);
    }
    fclose(in);
    return;
} /* end Config_load_from_VASP() */


#ifdef _vasp2cfg
int main (int argc, char *argv[])
{
    Aapp_Define_Config;
    char *POSCAR_fname="POSCAR", *POTCAR_fname="POTCAR";
    if (argc == 3) POSCAR_fname = argv[1];
    else if (argc == 4)
    {
        POTCAR_fname = argv[1];
        POSCAR_fname = argv[2];
    }
    else if (argc == 1)
    {
        printf ("\nPurpose: convert VASP POSCAR + "
                "POTCAR to Ju Li's configuration format.\n");
        printf ("see "JLCPL_URL"Graphics/A/ for details.\n\n");
        printf ("Usage: %s cfg_fname\n\t(POSCAR, POTCAR in ./)\n\n", argv[0]);
        printf ("   or, %s POSCAR_fname cfg_fname\n\t(POTCAR in ./)\n\n",
                argv[0]);
        printf ("   or, %s POTCAR_fname POSCAR_fname cfg_fname\n\n", argv[0]);
        return (1);
    }
    Config_load_from_VASP
        (POTCAR_fname, POSCAR_fname, stdout, Config_Aapp_to_Alib);
    Config_save (Config_Aapp_to_Alib, TRUE, argv[argc-1]);
    printf ("\"%s\" + \"%s\" -> \"%s\".\n",
            POTCAR_fname, POSCAR_fname, argv[argc-1]);
    return (0);
}
#endif /* _vasp2cfg */


#ifdef _cfg2vasp
int main (int argc, char *argv[])
{
    Aapp_Define_Config;
    Chemtab ct[1]={{0}};
    Tp *tp=NULL;
    FILE *in, *out;
    int i,j,k;
    char *potdir = "/usr/local/VASP/potpaw_GGA/elements", *PREC, *cmd,
        buf[FGETS_MAXSIZE];
    double kspacing_IN__A, HI[DIMENSION][DIMENSION], lengths[DIMENSION];
    if (argc == 4)
    {
        PREC = argv[2];
        kspacing_IN__A = atof(argv[3]);
    }
    else if (argc == 5)
    {
        PREC = argv[2];
        kspacing_IN__A = atof(argv[3]);
        potdir = argv[4];
    }
    else
    {
        printf("\nPurpose: create POSCAR\n"
               "                POTCAR\n"
               "                KPOINTS (if not in ./)\n"
               "                INCAR   (if not in ./)\n"
               "files for Vienna Ab-initio Simulation Package "
               "calculation.\n\n");
        printf ("Usage: %s cfg_fname PREC[Low|Med|High] kspacingxA\n"
                "\t(potdir=%s)\n\n", argv[0], potdir);
        printf ("   or, %s cfg_fname PREC[Low|Med|High] kspacingxA potdir\n"
                "\t(e.g. /usr/local/VASP/potpaw/elements)\n\n",
                argv[0]);
        return (1);
    }
    printf ("Loading \"%s\"...\n\n", argv[1]);
    CONFIG_LOAD (argv[1], Config_Aapp_to_Alib);
    rebind_CT (Config_Aapp_to_Alib, "", ct, &tp); cr();
    out = wOpen("POSCAR");
    fprintf (out, "%s\n", argv[1]);
    fprintf (out, "1.0\n");
    for (i=0; i<DIMENSION; i++)
    {
        for (j=0; j<DIMENSION; j++)
            fprintf ( out, "%.15g ", H[i][j] * ULENGTH_IN_A );
        fcr( out );
    }
    for (j=0; j<ct->t; j++)
        fprintf ( out, "%d ", ct->count[j] );
    fcr( out );
    fprintf (out, "Selective dynamics\n");
    fprintf (out, "Direct\n");
    for (j=0; j<ct->t; j++)
    {
        for (i=0; i<np; i++)
            if (tp[i] == j)
                fprintf ( out, "%.15g %.15g %.15g T T T\n",
                          s[DIMENSION*i], s[DIMENSION*i+1],
                          s[DIMENSION*i+2] );
    }
    fclose(out);
    printf ("\"%s\" -> \"POSCAR\";\n", argv[1]);
    out = wOpen("POTCAR");
    for (j=0; j<ct->t; j++)
    {
        cmd = str4(potdir,"/",blank_advance(SYM(ct->first[j])),"/POTCAR");
        in = ROpen( cmd );
        do
        {
            FGets(buf, in);
            fprintf (out, "%s\n", buf);
        } while (strstr(buf, "End of Dataset")==NULL);
        fclose(in);
        printf ("%d %s atoms:\n", ct->count[j],
                blank_advance(SYM(ct->first[j])));
        printf ("\"%s\" -> \"POTCAR\";\n", cmd);
    }
    fclose(out);
    if (!Fexists("KPOINTS"))
    {
        M3INV (H, HI, lengths[0]);
        M3columnlengths (HI, lengths);
        V3MuL (2*PI/ULENGTH_IN_A, lengths);
        out = wOpen("KPOINTS");
        i = ceil( lengths[0] / kspacing_IN__A);
        j = ceil( lengths[1] / kspacing_IN__A);
        k = ceil( lengths[2] / kspacing_IN__A);
        fprintf (out,
                 "Automatic mesh\n"
                 "0\n"
                 "Monkhorst-Pack\n"
                 "%d %d %d\n"
                 "0. 0. 0.\n",
                 i, j, k);
        fclose(out);
        printf ("k-spacing lengths < %g/A -> %d x %d x %d -> \"KPOINTS\";\n",
                kspacing_IN__A, i,j,k);
    }
    if (!Fexists("INCAR"))
    {
        out = wOpen("INCAR");
        fprintf (out,
                 "SYSTEM = %s\n"
                 "\n"
                 "Start Parameter for This Run:\n"
                 "NWRITE = 2\n"
                 "\n"
                 "Electronic Relaxation 1:\n"
                 "PREC   = %s\n"
                 "\n"
                 "Ionic Relaxation:\n"
                 "NSW    = 1000000\n"
                 "NBLOCK = 1\n"
                 "KBLOCK = 1\n"
                 "IBRION = 2\n"
                 "POTIM  = 0.5\n"
                 "ISIF   = 2\n"
                 "TEBEG  = 300.\n"
                 "SMASS  = 0\n"
                 "\n"
                 "DOS Related Values:\n"
                 "ISMEAR = 1\n"
                 "SIGMA  = 0.2\n"
                 "EMIN   = 10.\n"
                 "EMAX   = 0.\n"
                 "LELF   = T\n"
                 "LVTOT  = T\n"
                 "\n"
                 "Electronic Relaxation 2:\n"
                 "IALGO  = %d\n"
                 "LREAL  = %s\n"
                 "\n"
                 "http://cms.mpi.univie.ac.at/vasp/vasp/node95.html\n",
                 argv[1],
                 PREC,
                 (np>20) ?  48    :  38,
                 (np>20) ? "Auto" : ".FALSE."
                 );
        fclose(out);
        printf ("PREC = %s -> \"INCAR\";\n", PREC);
    }
    return (0);
}
#endif /* _cfg2vasp */


#ifdef _vasp2out
#define OUTCAR_HEADER_ENERGY   0
#define OUTCAR_HEADER_STRESS   1
#define OUTCAR_HEADER_CELL     2
#define OUTCAR_HEADER_ATOM     3
#define OUTCAR_HEADER_ION_LOOP 4   
#define OUTCAR_HEADER_COMPLETE 5
#define OUTCAR_HEADER_MAX      6
char *OUTCAR_HEADERS[OUTCAR_HEADER_MAX] =
{ "  energy  without entropy=",
  "  in kB",
  "      direct lattice vectors",
  " POSITION",
  "     LOOP+:",
  "                 Voluntary"
};
#define OUTCAR_LINESIZE 1024
#define draw(name,value) \
  fprintf( fp_draw, name" %e %.15g\n", drawtime, (double)(value) )
#define DRAW(NAME,value) draw( STR(NAME), value )
#define M3DRAW(NAME,A) (       draw(STR(NAME)"11",A[0][0]), \
  draw(STR(NAME)"12",A[0][1]), draw(STR(NAME)"13",A[0][2]), \
  draw(STR(NAME)"21",A[1][0]), draw(STR(NAME)"22",A[1][1]), \
  draw(STR(NAME)"23",A[1][2]), draw(STR(NAME)"31",A[2][0]), \
  draw(STR(NAME)"32",A[2][1]), draw(STR(NAME)"33",A[2][2]) )
#define M3DRAWMUL(NAME,A,factor) (      draw(STR(NAME)"11",A[0][0]*(factor)),\
  draw(STR(NAME)"12",A[0][1]*(factor)), draw(STR(NAME)"13",A[0][2]*(factor)),\
  draw(STR(NAME)"21",A[1][0]*(factor)), draw(STR(NAME)"22",A[1][1]*(factor)),\
  draw(STR(NAME)"23",A[1][2]*(factor)), draw(STR(NAME)"31",A[2][0]*(factor)),\
  draw(STR(NAME)"32",A[2][1]*(factor)), draw(STR(NAME)"33",A[2][2]*(factor)) )
#define SYMMAT_DRAW(NAME,A) ( \
  draw(STR(NAME)"11",A[0][0]), draw(STR(NAME)"12",A[0][1]), \
  draw(STR(NAME)"13",A[0][2]), draw(STR(NAME)"22",A[1][1]), \
  draw(STR(NAME)"23",A[1][2]), draw(STR(NAME)"33",A[2][2]) )
#define SYMMAT_DRAWMUL(NAME,A,factor) ( \
  draw(STR(NAME)"11",A[0][0]*(factor)), draw(STR(NAME)"12",A[0][1]*(factor)),\
  draw(STR(NAME)"13",A[0][2]*(factor)), draw(STR(NAME)"22",A[1][1]*(factor)),\
  draw(STR(NAME)"23",A[1][2]*(factor)), draw(STR(NAME)"33",A[2][2]*(factor)) )
int main (int argc, char *argv[])
{
    Aapp_Define_Config;
    Chemtab ct[1]={{0}};
    Tp *tp=NULL;
    char linebuffer[OUTCAR_LINESIZE];
    char *dirname;
    FILE *outcar, *fp_draw, *fp_cfg, *fp_scr_anim, *fp_screen;
    int i, j, k, step = 0, is_complete = 0;
    double drawtime=0, pote=0, tmp, x[3], *force=NULL,
        totalforcenorm, maxforcenorm;
    M3 stress, H0, HI, eta;
    if (argc != 2) 
    {
        /* printf ("\nPurpose: parse VASP OUTCAR (energy, cell geom., stress)\n" */
        /* "and save intermediate configurations w/ force-on-atoms\n" */
        /* "in Ju Li's Extended CFG format. The original INCAR, POSCAR,\n" */
        /* "POTCAR, KPOINTS inputs and DOSCAR, OUTCAR outputs will also\n" */
        /* "be copied for preservation.\n\n"); */
        /* printf ("Usage: %s output_directory\n\n", argv[0]); */
        /* return (1); */
        dirname = "out";
    }
    else dirname = argv[1];
    fp_screen = wOpen(str2(dirname, "/results.txt"));
    ftie (stdout, fp_screen);
    Fcr(ft);
    Config_load_from_VASP ("POTCAR", "POSCAR", ft, Config_Aapp_to_Alib);
    Fcr(ft); rebind_ct (Config_Aapp_to_Alib, "", ct, &tp, ft); Fcr(ft);
    dopen (dirname);
    outcar = rOpen("OUTCAR");
    fp_draw = wOpen(str2(dirname, "/draw.out"));
    fp_scr_anim = wOpen(str2(dirname, "/scr_anim"));
    fprintf (fp_scr_anim, "%d\n", 90);
    do
    {
        i = freadline_matchheader (outcar, OUTCAR_HEADER_MAX, OUTCAR_HEADERS,
                                   OUTCAR_LINESIZE, linebuffer);
        switch(i)
        {
            case OUTCAR_HEADER_STRESS:
                sscanf (linebuffer, "%lf %lf %lf %lf %lf %lf\n",
                        &stress[0][0], &stress[1][1], &stress[2][2],
                        &stress[0][1], &stress[1][2], &stress[0][2] );
                M3MultiplY (-1000*BAR_IN_GPA, stress);
                SYMMAT_DRAWMUL ( t, stress, 1 );
                fprintf (ft, "stress = %g %g %g %g %g %g GPa\n",
                         stress[0][0], stress[1][1], stress[2][2],
                         stress[1][2], stress[0][2], stress[0][1] );
                break;
            case OUTCAR_HEADER_ENERGY:
                sscanf (linebuffer, "%lf energy(sigma->0) = %lf\n",
                        &tmp, &pote);
                DRAW ( pote, pote );
                fprintf (ft, "step = %d, pote = %.15g eV\n", step, pote);
                break;
            case OUTCAR_HEADER_CELL:
                FGETS(linebuffer, OUTCAR_LINESIZE, outcar);
                sscanf(linebuffer,"%lf %lf %lf\n",&H[0][0],&H[0][1],&H[0][2]);
                FGETS(linebuffer, OUTCAR_LINESIZE, outcar);
                sscanf(linebuffer,"%lf %lf %lf\n",&H[1][0],&H[1][1],&H[1][2]);
                FGETS(linebuffer, OUTCAR_LINESIZE, outcar);
                sscanf(linebuffer,"%lf %lf %lf\n",&H[2][0],&H[2][1],&H[2][2]);
                M3DRAWMUL ( H, H, 1 );
                M3INV (H, HI, tmp);
                if (step == 0)
                {
                    M3EQV(H,H0);
                    if (pote!=0) S3fPR (ft, "H0 = %M [A];\n", H);
                }
                else if ( M3NE(H,H0) )
                {
                    S3fPR (ft, "H = %M [A];\n", H);
                    Lagrangian_strain(H0, H, eta);
                    S3fPR (ft, "eta = %M;\n", eta);
                }
                break;
            case OUTCAR_HEADER_ATOM:
                FGETS(linebuffer, OUTCAR_LINESIZE, outcar);
                REALLOC(main, force, DIMENSION*np, double);
                totalforcenorm = maxforcenorm = 0;
                for (i=0; i<np; i++)
                {
                    FGETS(linebuffer, OUTCAR_LINESIZE, outcar);
                    sscanf(linebuffer,"%lf %lf %lf %lf %lf %lf\n",
                           x, x+1, x+2,
                           force+DIMENSION*i,
                           force+DIMENSION*i+1,
                           force+DIMENSION*i+2);
                    V3mM3 ( x, HI, s+DIMENSION*i );
                    tmp = V3LENGTH2( force+DIMENSION*i );
                    totalforcenorm += tmp;
                    tmp = sqrt(tmp);
                    if (tmp > maxforcenorm) maxforcenorm = tmp;
                }
                totalforcenorm = sqrt( totalforcenorm );
                fprintf (ft, "total force norm = %.1e, "
                         "max force on one atom = %.1e eV/A\n",
                         totalforcenorm, maxforcenorm);
                sprintf(linebuffer, "%s/v%05d.cfg", dirname, step);
                fp_cfg = wopen(linebuffer);
                fprintf (fp_cfg, "Number of particles = %d\n", np);
                fprintf (fp_cfg, "H0(1,1) = %.15g A\n", H[0][0]);
                fprintf (fp_cfg, "H0(1,2) = %.15g A\n", H[0][1]);
                fprintf (fp_cfg, "H0(1,3) = %.15g A\n", H[0][2]);
                fprintf (fp_cfg, "H0(2,1) = %.15g A\n", H[1][0]);
                fprintf (fp_cfg, "H0(2,2) = %.15g A\n", H[1][1]);
                fprintf (fp_cfg, "H0(2,3) = %.15g A\n", H[1][2]);
                fprintf (fp_cfg, "H0(3,1) = %.15g A\n", H[2][0]);
                fprintf (fp_cfg, "H0(3,2) = %.15g A\n", H[2][1]);
                fprintf (fp_cfg, "H0(3,3) = %.15g A\n", H[2][2]);
                fprintf (fp_cfg, ".NO_VELOCITY.\n");
                fprintf (fp_cfg, "entry_count = %d\n", 3+3);
                fprintf (fp_cfg, "auxiliary[0] = fx [eV/A]\n");
                fprintf (fp_cfg, "auxiliary[1] = fy [eV/A]\n");
                fprintf (fp_cfg, "auxiliary[2] = fz [eV/A]\n");
                for (i=j=0; j<ct->t; j++)
                {
                    fprintf (fp_cfg, "%f\n%s\n", mass[i]*UMASS_IN_AMU, SYM(i));
                    for (k=0; k<ct->count[j]; k++,i++)
                        fprintf ( fp_cfg,
                                  "%.15g %.15g %.15g %.15g %.15g %.15g\n",
                                  s[DIMENSION*i],
                                  s[DIMENSION*i+1],
                                  s[DIMENSION*i+2],
                                  force[DIMENSION*i],
                                  force[DIMENSION*i+1],
                                  force[DIMENSION*i+2]);
                }
                fclose (fp_cfg);
                fprintf (ft, "OUTCAR -> \"%s\"\n", linebuffer);
                fprintf (fp_scr_anim, "./v%05d.cfg ./Jpg/v%05d.jpg\n",
                         step, step);
                break;
            case OUTCAR_HEADER_ION_LOOP:
                step ++;
                drawtime += 1;
                Fcr(ft);
                break;
            case OUTCAR_HEADER_COMPLETE:
                is_complete = 1;
                break;
        }
    } while (i>=-2);
    Free(force);
    fclose(fp_draw);
    fclose(fp_scr_anim);
    fclose(outcar);
    fprintf (ft, "%s\n\n", is_complete?
             "VASP run completed successfully.":
             "VASP run is NOT completed yet.");
    system(str3("cp INCAR ",   dirname, "/INCAR"));
    system(str3("cp KPOINTS ", dirname, "/KPOINTS"));
    system(str3("cp POSCAR ",  dirname, "/POSCAR"));
    system(str3("cp POTCAR ",  dirname, "/POTCAR"));
    system(str3("cp OUTCAR ",  dirname, "/OUTCAR"));
    system(str3("cp DOSCAR ",  dirname, "/DOSCAR"));
    fbrk();
    fclose(fp_screen);
    return (0);
}
#endif /* _vasp2out */


#ifdef _chg2xsf
#define MAXLINE      256
int main (int argc, char *argv[])
{
    Aapp_Define_Config;
    int i, j, lineno;
    char linebuf[MAXLINE], *CHGCAR_fname="CHGCAR", *POTCAR_fname="POTCAR";
    FILE *in, *out;
    double *s_shift, x_shift[3], thickness[3];

    if (argc == 3) CHGCAR_fname = argv[1];
    else if (argc == 4)
    {
        POTCAR_fname = argv[1];
        CHGCAR_fname = argv[2];
    }
    else if (argc == 1)
    {
        printf ("\nPurpose: convert VASP CHGCAR to XCrySDen XSF format:\n");
        printf ("see http://www.xcrysden.org/doc/XSF.html for details.\n\n");
        printf ("Usage: %s xsf_fname\n\t(CHGCAR, POTCAR in ./)\n\n", argv[0]);
        printf ("   or, %s CHGCAR_fname xsf_fname\n\t(POTCAR in ./)\n\n",
                argv[0]);
        printf ("   or, %s POTCAR_fname CHGCAR_fname xsf_fname\n\n", argv[0]);
        return (1);
    }
    Config_load_from_VASP
        (POTCAR_fname, CHGCAR_fname, stdout, Config_Aapp_to_Alib);

    M3rowthicknesses (H, thickness);
    printf ("Thickness[0] = %.15g [A]\n", thickness[0]);
    printf ("Thickness[1] = %.15g [A]\n", thickness[1]);
    printf ("Thickness[2] = %.15g [A]\n", thickness[2]);
    
    in = ROpen( CHGCAR_fname );
    out = fopen(argv[argc-1], "w");

    fprintf (out, "CRYSTAL\n");
    fprintf (out, "PRIMVEC\n");
    for (i=0; i<3; i++)
        fprintf(out,"%.15g %.15g %.15g\n", H[i][0], H[i][1], H[i][2]);
    fprintf(out, "PRIMCOORD\n");
    fprintf(out, "%d 1\n", np);

    s_shift = Config_SET_CM (Config_Aapp_to_Alib);
    V3mM3 (s_shift, H, x_shift);

    for (i=0; i<np; i++)
    {
        V3Trim (s+3*i);
        V3mM3 (s+3*i, H, s1+3*i);
        fprintf (out, "%2s: %.15g %.15g %.15g\n", SYM(i), V3E(s1+3*i));
    }
    fprintf(out, "\n");
    fprintf(out, "BEGIN_BLOCK_DATAGRID_3D\n");
    fprintf(out, "Electron Density (after x %.5g)\n", 1./M3VOLUME(H));
    fprintf(out, "DATAGRID_3D_#1\n");

    for (lineno=0; fgets(linebuf,MAXLINE,in)!=NULL; ++lineno)
    {
        if (lineno == (8+np))
        {
            fprintf(out, "%s", linebuf);
            fprintf(out, "%.15g %.15g %.15g\n", V3E(x_shift));
            for(i=0;i<3;i++)
                fprintf(out, "%.15g %.15g %.15g\n", H[i][0],H[i][1],H[i][2]);
        }
        if (lineno > (8+np))
        {            
            if (linebuf[0] == 'a') break;
            fprintf(out, "%s", linebuf);
        }
    }
    fprintf(out, "END_DATAGRID_3D\n");
    fprintf(out, "END_BLOCK_DATAGRID_3D\n");          
    fclose(out);
    fclose(in);

    printf ("\"%s\" + \"%s\" -> \"%s\".\n",
            POTCAR_fname, CHGCAR_fname, argv[argc-1]);
    return (0);
}
#endif /* _chg2xsf */
