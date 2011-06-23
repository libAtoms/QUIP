#!/bin/bash
# HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HND X
# HND X   libAtoms+QUIP: atomistic simulation library
# HND X
# HND X   Portions of this code were written by
# HND X     Albert Bartok-Partay, Silvia Cereda, Gabor Csanyi, James Kermode,
# HND X     Ivan Solt, Wojciech Szlachta, Csilla Varnai, Steven Winfield.
# HND X
# HND X   Copyright 2006-2010.
# HND X
# HND X   Not for distribution
# HND X
# HND X   Portions of this code were written by Noam Bernstein as part of
# HND X   his employment for the U.S. Government, and are not subject
# HND X   to copyright in the USA.
# HND X
# HND X   When using this software, please cite the following reference:
# HND X
# HND X   http://www.libatoms.org
# HND X
# HND X  Additional contributions by
# HND X    Alessio Comisso, Chiara Gattinoni, and Gianpietro Moras
# HND X
# HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

# FilePot script that invokes MolPro to get forces and energy of a
# cluster. Molpro reads XYZ format, so only the name of the input xyz 
# file needs to be added to the MolPro template

# Input file is $1, output goes to $2. Both file are in extended
# XYZ format.

olddir=`pwd`
[[ -z "$molpro" ]] && molpro="molpro"                                # Path to MolPro executable
[[ -z "$molpro_command_file" ]] && molpro_command_file=$olddir/molpro_command_file # Template command input file for MolPro
[[ -z "$molpro_dir" ]] && molpro_dir=$olddir                         # Directory in which to run MolPro

test_mode=0                                                          # Set to 1 to test script without actually running MolPro


inpfile=$1
outfile=$2
stem=`basename $inpfile .xyz`

cd $molpro_dir

# Extract a property from an extended XYZ file and print it
# $1 : file name
# $2 : property name
# $3 : if set to 1, print species labels, otherwise don't
# $4 : if nonempty, multiply property by this
function print_property {
awk 'NR == 2 {
  match($0,/Properties="?([^" ]*)/,a);
  nf = split(a[1],b,/:/);
  
  sum=0;
  for (i = 1; i <= nf; i+=3) {
    if (b[i] != "'$2'")
      { sum = sum + b[i+2]; }
    else
      { begin=sum+1; end=sum+b[i+2]; break; }
  }
  n = 1;
  for (i = begin; i <= end; i++) {
     fields[n]=i;
     n++;
  }
  n_fields=n-1;
}

NR > 2 {
  if ('$3'==1) printf "%s ",$1;
  for (i = 1; i <= n_fields; i++) printf "%16.8f",$fields[i];
  printf "\n";
}' $1
}

# Make a working directory if necessary
if [[ ! -d $stem ]]; then
    mkdir $stem
fi
cd $stem

cp $olddir/$inpfile .

if [[ -f $olddir/*.usp ]]; then
    cp $olddir/*.usp .
fi

N=`head -1 $inpfile`                # Number of atoms
params=`head -2 $inpfile | tail -1` # Parameter line

# Get the lattice from param line
lattice=`echo $params | awk '{match($0,/Lattice="?([^"]*)/,a); print a[1]}'`

echo molpro_driver ${stem}: got $N atoms

# Amend Molpro command file
cp ${molpro_command_file} ${stem}_molpro
echo "geometry=${inpfile}" >> ${stem}_molpro


# Invoke MolPro
if [[ $test_mode == 1 ]]; then
    echo molpro_driver: test mode
else
    rm -f ${stem}_molpro.out ${stem}_molpro.xml ${olddir}/${stem}.out
    ${molpro} ${stem}_molpro
fi


# Extract information from output file
energy_Ha=grep -A 1 '<property name="Energy"' ${stem}_molpro.xml | tr -d '\n' | sed 's/.*value="\(.*\)".*/\1/'
energy=echo "${energy_Ha} * 27.2113961" | bc -ql

# First line of output is number of atoms
$N > ${stem}.out

# Second line of output is parameter line
echo Lattice\=\"$lattice\" Properties=\"species:S:1:pos:R:3:force:R:3\" energy\=$energy\" >> ${stem}.out

# Extract gradient
awk '/<gradient/,/<[/]gradient/ {if(!/>/){print -$1*27.2113961/0.529177249, -$2*27.2113961/0.529177249, -$3*27.2113961/0.529177249}}' ${stem}_molpro.xml > ${stem}_tmpforce

# Merge pos and gradient
print_property $inpfile pos 1 > ${stem}_tmppos
paste ${stem}_tmppos ${stem}_tmpforce >> ${stem}.out


# Save all  output
cat ${stem}_molpro.out >> ${stem}_molpro_output


# Copy output file
cp ${stem}.out $olddir/$outfile

echo molpro_driver ${stem}: done, E\=$energy eV
exit 0
