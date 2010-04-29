#!/bin/bash
# HJ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HJ X
# HJ X   libAtoms+QUIP: atomistic simulation library
# HJ X
# HJ X   Portions of this code were written by
# HJ X     Albert Bartok-Partay, Silvia Cereda, Gabor Csanyi, James Kermode,
# HJ X     Ivan Solt, Wojciech Szlachta, Csilla Varnai, Steven Winfield.
# HJ X
# HJ X   Copyright 2006-2010.
# HJ X
# HJ X   These portions of the source code are released under the GNU General
# HJ X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
# HJ X
# HJ X   If you would like to license the source code under different terms,
# HJ X   please contact Gabor Csanyi, gabor@csanyi.net
# HJ X
# HJ X   Portions of this code were written by Noam Bernstein as part of
# HJ X   his employment for the U.S. Government, and are not subject
# HJ X   to copyright in the USA.
# HJ X
# HJ X
# HJ X   When using this software, please cite the following reference:
# HJ X
# HJ X   http://www.libatoms.org
# HJ X
# HJ X  Additional contributions by
# HJ X    Alessio Comisso, Chiara Gattinoni, and Gianpietro Moras
# HJ X
# HJ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

# Driver script for DL_Poly module
# Invoke DL_Poly to get forces and energy of classical region

# Input: XYZ format file with additional column for number
# of neighbours of "other" species, i.e. number of O neighbours
# for Si atoms and number of O neighbours for Si atoms
# e.g.:
# Si x y z 2.0  has 2 oxygen neighbours
# Distance units are Angstrom

# Output: first line is total energy. Subsequent lines are in 
# format "Element Name Fx Fy Fz Q" where (Fx,Fy,Fz) is the force
# on atom and Q the charge of atom. Force units are eV/A and
# charge units multiples of e.

# Conversion contant from DL_POLY forces to eV/A
# Energy unit is 10 J/mol, so conv. factor is 10/(N_A*e)
FCONV=0.00010364269508285074

dl_poly=./DLPOLY_sio.Y
in_file=$1
out_file=$2
#two_body=./two_body
#three_body=./three_body
template=./FIELD.template
species_table=./species_table
charges=./dlpoly_charges

# Extract a property from an extended XYZ file and print it
# $1 : file name
# $2 : property name
# $3 : if set to 1, print species labels, otherwise don't
function print_property {
    awk 'NR == 2 {
  match($0,/Properties="?([^" ]*)/,a);
  nf = split(a[1],b,/:/);
  
  sum=0;
  for (i = 1; i <= nf; i+=3) {
    if(b[i] != "'$2'")
      { sum = sum + b[i+2]}
    else
      { begin=sum+1; end=sum+b[i+2]; break; }
  };
  n = 1;
  for (i = begin+1; i <= end +1; i++) {
     fields[n]=i;
     n++;
  }
  n_fields=n-1;
 }

NR > 2 {
  if ('$3'==1) printf "%s ",$1
  for (i = 1; i <= n_fields; i++) printf "%16.8f",$fields[i];
  printf "\n";
}' $1
}


rm -f $out_file $charges CONFIG FIELD REVCON OUTPUT

N=`awk 'NR==1 { print $1 }' $in_file` # Number of atoms
params=`head -2 $in_file | tail -1` # Parameter line

# Get the lattice from param line
lattice=`echo $params | awk '{match($0,/Lattice="?([^"]*)/,a); print a[1]}'`


# List of atomic species present in input file
species=`tail -$N $in_file | awk '{print $1}' | sort | uniq`

echo -n "dl_poly_driver: got $N atoms, species: "
# Report species counts
for s in $species; do
    echo -n $s `grep -c "^$s" $in_file` ' '
done
echo

# Check we've got species table entry for all species
for s in $species; do
    if ! grep -Eq "^$s" $species_table; then
	echo dl_poly_driver: No species table entry for element $s
	exit 1
    fi
done

# CONFIG file header
echo "dl_poly_driver" > CONFIG
echo "0      2" >> CONFIG
echo $lattice | awk '{printf "%16.8f%16.8f%16.8f\n",$1,$2,$3}' >> CONFIG
echo $lattice | awk '{printf "%16.8f%16.8f%16.8f\n",$4,$5,$6}' >> CONFIG
echo $lattice | awk '{printf "%16.8f%16.8f%16.8f\n",$7,$8,$9}' >> CONFIG

# FIELD file header
cat > FIELD <<EOF
LOTF DL_Poly driver
UNITS eV
molecular types  1
Molecule  1
nummols  1
EOF

# Atom lines for CONFIG and FIELD files
# Lookup charges and masses of each element using species_table
echo atoms $N >> FIELD
awk '\
BEGIN {
  while(getline < "'$species_table'") {
    mass[$1]    = $2
    charge[$1]  = $3
  }
}
NR > 2 {
  printf "%s\n%16.9f%20.9f%20.9f\n",$1,$2,$3,$4 >>"CONFIG"
  printf "%8s%12.4f%12.4f%10d%10d%10d\n",$1,mass[$1],charge[$1]*$5,1,0,NR-2 >>"FIELD"
  printf "%8s%12.4f\n", $1, charge[$1]*$5 >>"'$charges'"
}' $in_file
echo finish >> FIELD

# Append force fields from template
cat $template >> FIELD

# Add only the lines we need from two_body: those which
# define potentials between species that we have in input file
#two_body_name=`head -1 $two_body`
#awk '\
#NR==1 {
#  n = 0
#  split(species,s1)
#  for (s in s1) sp[s1[s]]=1
#}
#NR > 1 && $1 in sp && $2 in sp { 
#  lines[n++] = $0
#}
#END {
#  print name, n
#  for (i=0; i<n;i++) print lines[i]
#}' species="$species" name="$two_body_name" $two_body >> FIELD

# Now do the same for three body potential
#three_body_name=`head -1 $three_body`
#awk '\
#NR==1 {
#  n = 0
#  split(species,s1)
#  for (s in s1) sp[s1[s]]=1
#}
#NR > 1 && $1 in sp && $2 in sp && $3 in sp { 
#  lines[n++] = $0
#}
#END {
#  print name, n
#  for (i=0; i<n;i++) print lines[i]
#}' species="$species" name="$three_body_name" $three_body >> FIELD

# Finish FIELD file
echo close >> FIELD

# Run DL_Poly
energy=`$dl_poly | awk '{print $1}'`

if grep -q error OUTPUT; then
    echo dl_poly_driver: error running $dl_poly
    exit 1
fi

if [[ ! ( -f REVCON && -f OUTPUT) ]]; then
    echo dl_poly_driver: can\'t find REVCON and/or OUTPUT files
    exit 1
fi

runtime=`grep 'time elapsed' OUTPUT | tail -1 | awk '{print $6}'`

# First line of output is number of atoms
echo $N > $out_file

echo Lattice\=\"$lattice\" Properties=pos:R:3:force:R:3 energy=$energy >> $out_file

# Get forces from REVCON file, mapping pseudo species
# back to original species names. Also print charges
paste <(print_property $in_file pos 1) <(awk '\
BEGIN {
  i=0
  while (getline < "'$charges'") charge[i++]=$2
  i=0
}

NR > 5 && NR % 4 == 1 {
  printf "%16.8f%16.8f%16.8f\n",$1*FCONV,$2*FCONV,$3*FCONV
}' FCONV=$FCONV REVCON) >> $out_file

# Check we've got N forces and 1 energy
if [[ `wc -l $out_file | awk '{print $1}'` != $(($N+2)) ]]
then
    echo dl_poly_driver: Error parsing DL_POLY output files
    exit 1
fi

#nl $charges >> charge_log

echo dl_poly_driver: done, E\=`printf %.1f $energy` eV, time $runtime s
exit 0
