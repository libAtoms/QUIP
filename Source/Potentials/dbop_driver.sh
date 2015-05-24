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

# FilePot script that invokes DBOP to get forces and energy of a
# system. This script converts XYZ input file to the DBOP input file,

# Input file is $1, output goes to $2. Both file are in extended
# XYZ format.

olddir=`pwd`
[[ -z "$bopdir" ]] && bopdir=$olddir     # Directory in which to find BOP param files
[[ -z "$bopexe" ]] && bopexe=$olddir/dbop  # Path to executable

test_mode=0                                # Set to 1 to test script without actually running the force model

# name of seq. ; on BSD systems, it happens to be called jot
if which seq >& /dev/null ; then 
    SEQ=seq
elif which jot >& /dev/null ; then
    SEQ=jot
else
    echo "Cannot find \`seq\' or equivalent"
    exit
fi

xyzfile=$1
outfile=$2
stem=`basename $xyzfile .xyz`

cd "$olddir"

# Extract a property from an extended XYZ file and print it
# $1 : file name
# $2 : property name
# $3 : if set to 1, print species labels, otherwise don't
# $4 : if nonempty, multiply property by this
function print_property {
gawk 'NR == 2 {
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

cp "$olddir/$xyzfile" .


N=`head -1 $xyzfile`                # Number of atoms
params=`head -2 $xyzfile | tail -1` # Parameter line

# Get the lattice from param line
lattice=`echo $params | gawk '{match($0,/Lattice="?([^"]*)/,a); print a[1]}'`

echo dbop_driver ${stem}: got $N atoms
#echo dbop_driver ${stem}: got unit cell, $lattice

# we now need to invert the lattice matrix because BOP needs the atomic positions in fractional coords
invdet=`echo $lattice | awk '{printf "%24.16e", 1/($1*($5*$9-$6*$8)-$2*($4*$9-$6*$7)+$3*($4*$8-$5*$7));}'`
i1=`echo $lattice | awk -v id="$invdet" '{printf "%24.16e",  ($5*$9-$6*$8)*id;}'`
i4=`echo $lattice | awk -v id="$invdet" '{printf "%24.16e", -($4*$9-$6*$7)*id;}'`
i7=`echo $lattice | awk -v id="$invdet" '{printf "%24.16e",  ($4*$8-$5*$7)*id;}'`
i2=`echo $lattice | awk -v id="$invdet" '{printf "%24.16e", -($2*$9-$3*$8)*id;}'`
i5=`echo $lattice | awk -v id="$invdet" '{printf "%24.16e",  ($1*$9-$3*$7)*id;}'`
i8=`echo $lattice | awk -v id="$invdet" '{printf "%24.16e", -($1*$8-$2*$7)*id;}'`
i3=`echo $lattice | awk -v id="$invdet" '{printf "%24.16e",  ($2*$6-$3*$5)*id;}'`
i6=`echo $lattice | awk -v id="$invdet" '{printf "%24.16e", -($1*$6-$3*$4)*id;}'`
i9=`echo $lattice | awk -v id="$invdet" '{printf "%24.16e",  ($1*$5-$2*$4)*id;}'`
#echo dbop_driver ${stem}: inverse unit cell: ${i1} ${i2} ${i3} ${i4} ${i5} ${i6} ${i7} ${i8} ${i9}  

# Make BOP input file from input xyz file
echo "A" > input
echo $lattice | awk '{print $1,$2,$3}' >> input
echo $lattice | awk '{print $4,$5,$6}' >> input
echo $lattice | awk '{print $7,$8,$9}' >> input
echo "LEN" >> input
echo "1.0 1.0 1.0" >> input
echo "LATPAR" >> input
echo "3.1652" >> input
echo "ND" >> input
echo $N >> input
echo "D" >> input
print_property $xyzfile pos 1 | awk "{printf \"%24.16e%24.16e%24.16e %s %s\n\", \$2*${i1}+\$3*${i4}+\$4*${i7}, \$2*${i2}+\$3*${i5}+\$4*${i8}, \$2*${i3}+\$3*${i6}+\$4*${i9}, \$1, \"0.0\"}" >> input
echo "NINERT" >> input
echo "0" >> input
echo "DINERT" >> input
echo >> input


# Invoke BOP
if [[ $test_mode == 1 ]]; then
    echo dbop_driver: test mode
else
    rm -f quip.out out phonon.eng block.unwrapped block.out "${olddir}/${stem}.out"
    ln -sf "${bopdir}/BNDSCL.WW" .
    ln -sf "${bopdir}/BOP.atomdat" .
    ln -sf "${bopdir}/PARAM.BOP" .
    ln -sf "${bopdir}/SPLINE.WW" 
    ln -sf "${bopdir}/ENV.WW"  .
    
    ln -sf PARAM.BOP fort.8
    ln -sf SPLINE.WW fort.12
    ln -sf ENV.WW fort.14
    ln -sf input fort.31
    ln -sf out fort.9
    "${bopexe}" 2>&1 1>outerr
    rm fort.*
fi

# Extract information from quip.out file
energy=`grep -A 1 ENERGY quip.out | tail -1|tr -d ' '`

# No stress tensor available from BOP yet


# First line of output is number of atoms
echo $N > ${stem}.out

# Second line of output is parameter line
echo Lattice\=\"$lattice\" Properties=\"species:S:1:pos:R:3:force:R:3\" energy\=$energy >> ${stem}.out

# Extract forces
grep -A $N FORCES quip.out | tail -$N  > ${stem}_tmpforce

# Combine atomic positions and forces
print_property $xyzfile pos 1 > ${stem}_tmppos 
paste ${stem}_tmppos ${stem}_tmpforce >> ${stem}.out
    
# Save BOP raw output

#cat out >> ${stem}.out_log
#cat outerr >> ${stem}.outerr_log


if [[ `wc -l ${stem}.out | awk '{print $1}'` != $(($N+2)) ]]; then
    echo dbop_driver ${stem}: Error parsing BOP output file 
    exit 1
fi

# Check there are no crazily large forces
max_force=`print_property ${stem}.out force 0 | awk '{ print sqrt($1*$1 + $2*$2 + $3*$3) }' | sort -g -r | head -1`
echo dbop_driver ${stem}: max force is $max_force 

# Copy output file
cp ${stem}.out "$olddir/$outfile"

echo dbop_driver ${stem}: done, E\=$energy eV, time $ctime s
exit 0
