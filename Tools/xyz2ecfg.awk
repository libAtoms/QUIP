#!/bin/awk -f
# H0 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# H0 X
# H0 X   libAtoms+QUIP: atomistic simulation library
# H0 X
# H0 X   Portions of this code were written by
# H0 X     Albert Bartok-Partay, Silvia Cereda, Gabor Csanyi, James Kermode,
# H0 X     Ivan Solt, Wojciech Szlachta, Csilla Varnai, Steven Winfield.
# H0 X
# H0 X   Copyright 2006-2010.
# H0 X
# H0 X   These portions of the source code are released under the GNU General
# H0 X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
# H0 X
# H0 X   If you would like to license the source code under different terms,
# H0 X   please contact Gabor Csanyi, gabor@csanyi.net
# H0 X
# H0 X   Portions of this code were written by Noam Bernstein as part of
# H0 X   his employment for the U.S. Government, and are not subject
# H0 X   to copyright in the USA.
# H0 X
# H0 X
# H0 X   When using this software, please cite the following reference:
# H0 X
# H0 X   http://www.libatoms.org
# H0 X
# H0 X  Additional contributions by
# H0 X    Alessio Comisso, Chiara Gattinoni, and Gianpietro Moras
# H0 X
# H0 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

function ceiling(x)
{
  return (x == int(x)) ? x : int(x)+1
}

BEGIN {
  command="basename " ARGV[1] " .xyz";
  command | getline cfgbase;
  frame=0;
  print "Counting frames..."
  command="grep -c Lattice " ARGV[1];
  command | getline nframes;
  print "Got " nframes " frames.";
  format=sprintf("%%0%dd",ceiling(log(nframes)/log(10.0)))

  mass["Ac"] = 227.000000;
  mass["Ag"] = 107.868200;
  mass["Al"] = 26.981540;
  mass["Am"] = 243.000000;
  mass["Ar"] = 39.948000;
  mass["As"] = 74.921600;
  mass["At"] = 210.000000;
  mass["Au"] = 196.966550;
  mass["B"] = 10.811000;
  mass["Ba"] = 137.327000;
  mass["Be"] = 9.012187;
  mass["Bh"] = 264.000000;
  mass["Bi"] = 208.980380;
  mass["Bk"] = 247.000000;
  mass["Br"] = 79.904000;
  mass["C"] = 12.010700;
  mass["Ca"] = 40.078000;
  mass["Cd"] = 112.411000;
  mass["Ce"] = 140.116000;
  mass["Cf"] = 251.000000;
  mass["Cl"] = 35.452700;
  mass["Cm"] = 247.000000;
  mass["Co"] = 58.933200;
  mass["Cr"] = 51.996100;
  mass["Cs"] = 132.905450;
  mass["Cu"] = 63.546000;
  mass["Db"] = 262.000000;
  mass["Ds"] = 271.000000;
  mass["Dy"] = 162.500000;
  mass["Er"] = 167.260000;
  mass["Es"] = 252.000000;
  mass["Eu"] = 151.964000;
  mass["F"] = 18.998400;
  mass["Fe"] = 55.845000;
  mass["Fm"] = 257.000000;
  mass["Fr"] = 223.000000;
  mass["Ga"] = 69.723000;
  mass["Gd"] = 157.250000;
  mass["Ge"] = 72.610000;
  mass["H"] = 1.007940;
  mass["He"] = 4.002600;
  mass["Hf"] = 178.490000;
  mass["Hg"] = 200.590000;
  mass["Ho"] = 164.930320;
  mass["Hs"] = 265.000000;
  mass["I"] = 126.904470;
  mass["In"] = 114.818000;
  mass["Ir"] = 192.217000;
  mass["K"] = 39.098300;
  mass["Kr"] = 83.800000;
  mass["La"] = 138.905500;
  mass["Li"] = 6.941000;
  mass["Lr"] = 262.000000;
  mass["Lu"] = 174.967000;
  mass["Md"] = 258.000000;
  mass["Mg"] = 24.305000;
  mass["Mn"] = 54.938050;
  mass["Mo"] = 95.940000;
  mass["Mt"] = 268.000000;
  mass["N"] = 14.006740;
  mass["Na"] = 22.989770;
  mass["Nb"] = 92.906380;
  mass["Nd"] = 144.240000;
  mass["Ne"] = 20.179700;
  mass["Ni"] = 58.693400;
  mass["No"] = 259.000000;
  mass["Np"] = 237.000000;
  mass["O"] = 15.999400;
  mass["Os"] = 190.230000;
  mass["P"] = 30.973760;
  mass["Pa"] = 231.035880;
  mass["Pb"] = 207.200000;
  mass["Pd"] = 106.420000;
  mass["Pm"] = 145.000000;
  mass["Po"] = 209.000000;
  mass["Pr"] = 140.907650;
  mass["Pt"] = 195.078000;
  mass["Pu"] = 244.000000;
  mass["Ra"] = 226.000000;
  mass["Rb"] = 85.467800;
  mass["Re"] = 186.207000;
  mass["Rf"] = 261.000000;
  mass["Rg"] = 272.000000;
  mass["Rh"] = 102.905500;
  mass["Rn"] = 222.000000;
  mass["Ru"] = 101.070000;
  mass["S"] = 32.066000;
  mass["Sb"] = 121.760000;
  mass["Sc"] = 44.955910;
  mass["Se"] = 78.960000;
  mass["Sg"] = 263.000000;
  mass["Si"] = 28.085500;
  mass["Sm"] = 150.360000;
  mass["Sn"] = 118.710000;
  mass["Sr"] = 87.620000;
  mass["Ta"] = 180.947900;
  mass["Tb"] = 158.925340;
  mass["Tc"] = 98.000000;
  mass["Te"] = 127.600000;
  mass["Th"] = 232.038100;
  mass["Ti"] = 47.867000;
  mass["Tl"] = 204.383300;
  mass["Tm"] = 168.934210;
  mass["U"] = 238.028900;
  mass["Uub"] = 285.000000;
  mass["Uuh"] = 292.000000;
  mass["Uup"] = 288.000000;
  mass["Uuq"] = 289.000000;
  mass["Uut"] = 284.000000;
  mass["V"] = 50.941500;
  mass["W"] = 183.840000;
  mass["Xe"] = 131.290000;
  mass["Y"] = 88.905850;
  mass["Yb"] = 173.040000;
  mass["Zn"] = 65.390000;
  mass["Zr"] = 91.224000;

}

{
  if (nframes == 1) {
    cfgfile=cfgbase ".cfg";
  } else {
    cfgfile=cfgbase sprintf(format,frame) ".cfg";
  }

  natoms=$1;
  printf "Frame "frame": " natoms" atoms              \r";

  print "Number of particles = "natoms > cfgfile;

  getline comment;
  print "# "comment >> cfgfile;

  match(comment,/Lattice="([^"]*)/,a);
  lat=a[1];
  split(lat,lattice);
   

  for (x=0; x<3; x++)
    for (y=0; y<3; y++) {
      r[x,y]=lattice[x*3+y+1];
      print "H0("x+1","y+1") = "r[x,y] >> cfgfile;
    }

  det = r[0,0]*(r[1,1]*r[2,2] - r[2,1]*r[1,2]) \
    - r[1,0]*(r[0,1]*r[2,2] - r[2,1]*r[0,2]) \
  + r[2,0]*(r[0,1]*r[1,2] - r[1,1]*r[0,2]);

  g[0,0] = ((r[1,1] * r[2,2]) -\
	    (r[1,2] * r[2,1]))/det;
  g[1,0] = ((r[0,2] * r[2,1]) -\
	    (r[0,1] * r[2,2]))/det;
  g[2,0] = ((r[0,1] * r[1,2]) -\
	    (r[0,2] * r[1,1]))/det;
  g[0,1] = ((r[1,2] * r[2,0]) -\
	    (r[1,0] * r[2,2]))/det;
  g[1,1] = ((r[0,0] * r[2,2]) -\
	     (r[0,2] * r[2,0]))/det;
  g[2,1] = ((r[0,2] * r[1,0]) -\
	     (r[0,0] * r[1,2]))/det;
  g[0,2] = ((r[1,0] * r[2,1]) -\
	     (r[1,1] * r[2,0]))/det;
  g[1,2] = ((r[0,1] * r[2,0]) -\
	    (r[0,0] * r[2,1]))/det;
  g[2,2] = ((r[0,0] * r[1,1]) -\
	     (r[0,1] * r[1,0]))/det;

  print ".NO_VELOCITY." >> cfgfile;

  match(comment,/Properties=([a-zA-Z0-9._:]*)/,a);
  props=a[1];
  nfields = split(props,a,":");

  delete lines;
  p=0;
  for (i = 7; i <= nfields; i++) {
    if (i % 3 == 1) name=a[i];

    if (i % 3 == 0) {
      for (j=0; j<a[i]; j++) {
        lines[p]="auxiliary["p"] = "name;
	p=p+1;
      }
    }
  }

  print "entry_count = "p+3 >> cfgfile;

  for (i=0; i<p; i++)
    print lines[i] >> cfgfile

  for (i=0; i<natoms; i++) {
    getline;
    
    print $1 >> cfgfile;
    print mass[$1] >> cfgfile;

    s[0]=$2;
    s[1]=$3;
    s[2]=$4;

    t[0] = g[0,0]*s[0]+g[0,1]*s[1]+g[0,2]*s[2];
    t[1] = g[1,0]*s[0]+g[1,1]*s[1]+g[1,2]*s[2];
    t[2] = g[2,0]*s[0]+g[2,1]*s[1]+g[2,2]*s[2];

    printf "%16.8f%16.8f%16.8f",t[0],t[1],t[2] >> cfgfile;

    for (j=5;j<=NF;j++) printf "%16.8f",$j >> cfgfile;
    printf "\n" >> cfgfile;
  }     

  close(cfgfile);
  frame=frame+1;
}

END {
  print "Done "frame" frames.                  "
}
