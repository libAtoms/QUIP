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

/* Dacapo total energy program            */
/* http://www.fysik.dtu.dk/campos/Dacapo/ */

#ifdef _cfg2dacapo
int main (int argc, char *argv[])
{
    Aapp_Define_Config;
    Chemtab ct[1]={{0}};
    ConfigSpecies root[1]={{0}};
    char *fname_py;
    Tp *tp=NULL;
    FILE *in, *out;
    int i,j,k,len;
    char *potdir = "/usr/local/dacapo/DACAPOPATH";
    double kspacing_IN__A, HI[DIMENSION][DIMENSION], lengths[DIMENSION];

    if (argc == 3)
    {
        kspacing_IN__A = atof(argv[2]);
    }
    else if (argc == 4)
    {
        kspacing_IN__A = atof(argv[2]);
        potdir = argv[3];
    }
    else
    {
        printf("\nPurpose: create template python script for CAMPOS/Dacapo "
               "DFT calculation.\n\n");
        printf ("Usage: %s cfg_fname kspacingxA\n"
                "\t(potdir=%s)\n\n", argv[0], potdir);
        printf ("   or, %s cfg_fname kspacingxA potdir\n"
                "\t(e.g. ~/dacapo/DACAPOPATH)\n\n",
                argv[0]);
        return (1);
    }
    printf ("\nLoading \"%s\"...\n\n", argv[1]);
    CONFIG_LOAD (argv[1], Config_Aapp_to_Alib);
    rebind_CT (Config_Aapp_to_Alib, "", ct, &tp); cr();
    
    fname_py = str2(file_basename1(argv[1]), ".py");
    out = wOpen(fname_py);

    fprintf(out, "#!/usr/bin/python\n");
    fprintf (out, "# http://www.google.com/search?q=dacapoNETCDF&btnI=f\n");
    fprintf (out, "# http://www.google.com/search?q=dacapo_tutorial&btnI=f\n");
    fprintf (out, "# http://oldwww.fysik.dtu.dk/campos/FAQ/FAQ.html\n\n");

    fprintf (out, "# -------- packages -------- #\n\n");
    fprintf(out, "from Simulations.Dacapo import *\n");
    fprintf(out, "from Simulations.Dacapo.EigenState import EigenState\n");
    fprintf(out, "from Simulations.Dacapo.ListOfEigenStates "
            "import ListOfEigenStates\n\n");
    fprintf (out, "# -------- potentials -------- #\n\n");

    fprintf(out, "sim = Simulation()\n\n");
    for (j=0; j<ct->t; j++)
    {
        fprintf (out, "%s.SetProperty\(\"PseudopotentialPath\",\n",
                 blank_advance(SYM(ct->first[j])));
        fprintf (out, "\t\t\"%s/%2s_us_gga_7.3.4.pseudo\")\n",
                 potdir, blank_advance(SYM(ct->first[j])));
    }
    fcr (out); fcr (out);
    fprintf (out, "# -------- configurations -------- #\n\n");
    fprintf ( out, "# ucell:  unit is A\n");
    fprintf ( out, "ucell = BravaisLattice([\n");
    for (i=0; i<DIMENSION; i++)
    {
        fprintf ( out, "\t\t\t[ ");
        for (j=0; j<DIMENSION; j++)
        {
            if (j == DIMENSION-1)
                fprintf ( out, "%.15g ", H[i][j] * ULENGTH_IN_A );
            else
                fprintf ( out, "%.15g, ", H[i][j] * ULENGTH_IN_A );
        }
        if (i == DIMENSION-1)
            fprintf ( out, " ] ])\n");
        else
            fprintf ( out, " ],\n");
    }
    
    fprintf( out, "\nsim.atoms = ListOfAtoms([\n");
    for (i=0; i<np; i++)
    {
        fprintf (out, "\tAtom( type=%2s, position=Vector\n"
                 "([ %.15g, %.15g, %.15g ], ucell) )",
                 SYM(i), s[3*i], s[3*i+1], s[3*i+2]);
        if (i != np-1)
            fprintf (out, ",\n");
        else
            fprintf (out, "],\n");
    }
    fprintf (out, "\tunitcell=ucell)\n\n");
    fprintf (out, "# -------- input -------- #\n\n");
    fprintf (out, "sim.encut = PlaneWaveCutoff( 250.00)    # unit is eV\n");
    fprintf (out, "sim.eband = ElectronicBands(20)"
             "\t# change: electronsperatom*numberofatoms/2+10\n");
    fprintf (out, "sim.vxc = NetCDF.Entry(name=\"ExcFunctional\", "
             "value=\"PZ\")");
    fprintf (out, "\n# options: PZ VWN PW91 PBE revPBE RPBE");
    fcr (out);
    fprintf (out, "sim.sym = NetCDF.Entry(name=\"UseSymmetry\", "
             "value=\"Off\")");
    fprintf (out, "\n# options: Off Maximum");
    fcr (out);
    fcr (out);
    fprintf (out, "# -------- output -------- #\n\n");
    fprintf (out, "sim.debug = NetCDF.Entry(name= \"PrintDebugInfo\", "
             "value=\"MediumLevel\")\n");
    fprintf (out, "# options: Off MediumLevel HighLevel\n");
    fprintf (out,"sim.netcdf = NetCDF.Entry(name=\"NetCDFOutputControl\")\n");
    fprintf (out, "sim.netcdf.PrintEffPotential = \"Yes\"\n");
    fprintf (out, "sim.netcdf.PrintElsPotential = \"Yes\"\n");
    fprintf (out, "sim.netcdf.PrintChargeDensity = \"Yes\"\n");
    fprintf (out,"sim.pdos = NetCDF.Entry(name=\"PrintAtomProjectedDOS\")\n");
    fprintf (out, "sim.pdos.EnergyWindow = (-23., 1.)   # unit is eV\n");
    fprintf (out, "sim.pdos.EnergyWidth = 0.30          # unit is eV\n");
    fprintf (out, "sim.pdos.NumberEnergyPoints = 100\n");
    fprintf (out, "sim.pdos.CutoffRadius = 1            # unit is A\n");
    fcr (out);
    fprintf (out, "# -------- some default parameters -------- #\n\n");
    M3INV (H, HI, lengths[0]);
    M3columnlengths (HI, lengths);
    V3MuL (2*PI/ULENGTH_IN_A, lengths);
    i = ceil( lengths[0] / kspacing_IN__A);
    j = ceil( lengths[1] / kspacing_IN__A);
    k = ceil( lengths[2] / kspacing_IN__A);
    fprintf (out, "# sim.kpt = NetCDF.Entry(name=\"KpointSetup\", "
             "value=[%d,%d,%d])\n",i, j, k);
    fprintf (out, "# sim.kpt.gridtype = \"MonkhorstPack\"\n");
    fprintf (out, "# options: MonkhorstPack ChadiCohen\n");
    fprintf (out, "# sim.eminimization = "
             "NetCDF.Entry(name=\"ElectronicMinimization\")\n");
    fprintf (out, "# sim.eminimization.Method = \"eigsolve\"\n");
    fprintf (out, "# options: resmin eigsolve rmm-diis");
    fcr (out);
    fprintf (out, "# sim.eminimization.DiagonalizationsPerBand = 2");
    fcr (out);
    fprintf (out, "# sim.extracharge = NetCDF.Entry(name=\"ExtraCharge\", "
             "value = 0.)\n");
    fprintf (out, "# unit is electron/supercell\n");
    fprintf (out, "# sim.dyn = NetCDF.Entry(\"Dynamics\")\n");
    fprintf (out, "# sim.dyn.Type = \"Relaxation\"\n");
    fprintf (out, "# options: Static Relaxation ElasticBand "
             "ExternalIonMotion\n");
    fprintf (out, "# sim.dyn.Method = \"Quickmin_Atomwise\"\n");
    fprintf (out, "# options: Verlet Quickmin_Atomwise "
             "Quickmin_Sum BFGS ConjugateGradient\n");
    fprintf (out, "# sim.dyn.Step = 2.0");
    fprintf (out, "\t# unit is femto-second\n");
    fcr (out);
    fprintf (out, "# -------- execute -------- #\n\n");
    fprintf (out, "sim.Execute(outfile=\"%s.nc\", ",file_basename1(argv[1]));
    fprintf (out, "ascii=\"%s.txt\")\n", file_basename1(argv[1]));
    fprintf (out, "sim.UpdateFromNetCDFFile(\"%s.nc\")\n",
             file_basename1(argv[1]));
    fprintf (out, "sim.all = ListOfEigenStates()\n");
    fcr (out);
    fprintf (out, "# -------- terminal information -------- #\n\n");
    fprintf (out, "print \"Total potential energy = \", "
             "sim.atoms.GetTotalPotentialEnergy()\n");
    fprintf (out, "print \"Forces on atoms        = \", "
             "sim.atoms.GetCartesianForces()\n");

    fchmod (fileno(out), (mode_t)0755);
    fclose (out);
    
    printf ("\"%s\" -> \"%s.py\"\n\n", argv[1], file_basename1(argv[1]));
    
    return (0);
}
#endif /* _cfg2dacapo */
