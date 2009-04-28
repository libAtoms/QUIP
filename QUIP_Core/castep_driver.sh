#!/bin/bash

# FilePot script that invokes CASTEP to get forces and energy of a
# cluster This script converts XYZ input file to castep cell file,
# then uses template .cell and .param files for the rest of the CASTEP
# parameters.  The atom number mangling necessary to convert between
# simple atom number 1..N and the (species, number within species)
# numbering that CASTEP uses is dealt with by this script.

# Input file is $1, output goes to $2. Both file are in extended
# XYZ format.

# name of seq. ; on BSD systems, it happens to be called jot
if which seq >& /dev/null ; then 
    SEQ=seq
elif which jot >& /dev/null ; then
    SEQ=jot
else
    echo "Cannot find \`seq\' or equivalent"
    exit
fi

olddir=`pwd`
[[ -z "$castep" ]] && castep=$olddir/castep  # Path to CASTEP executable
template=$olddir/castep_driver         # Template .cell and .param files
castep_dir=.                           # Directory in which to run CASTEP
use_check_files=1                      # Should we try to restart from .check files?
max_force_tol=99999.0                      # Max force that is considered reasonable:
                                       # if there are any larger forces we rerun CASTEP
test_mode=0                            # Set to 1 to test script without actually running castep

xyzfile=$1
outfile=$2
stem=`basename $xyzfile .xyz`

cd $castep_dir

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
    if(b[i] != "'$2'")
      { sum = sum + b[i+2]}
    else
      { begin=sum+1; end=sum+b[i+2]; break; }
  };
  n = 1;
  for (i = begin; i <= end; i++) {
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

# Make a working directory if necessary
if [[ ! -d $stem ]]
    then
    mkdir $stem
fi
cd $stem

cp $olddir/$xyzfile .

cp $olddir/*.usp .

N=`head -1 $xyzfile` # Number of atoms
params=`head -2 $xyzfile | tail -1` # Parameter line

# Get the lattice from param line
lattice=`echo $params | awk '{match($0,/Lattice="?([^"]*)/,a); print a[1]}'`

echo castep_driver ${stem}: got $N atoms

# Make CASTEP cell file from input xyz file
echo "%BLOCK LATTICE_CART" > ${stem}.cell
echo $lattice | awk '{print $1,$2,$3}' >> ${stem}.cell
echo $lattice | awk '{print $4,$5,$6}' >> ${stem}.cell
echo $lattice | awk '{print $7,$8,$9}' >> ${stem}.cell
echo "%ENDBLOCK LATTICE_CART" >> ${stem}.cell
echo >> ${stem}.cell
echo "%BLOCK POSITIONS_ABS" >> ${stem}.cell

# Extract positions from XYZ file
print_property $xyzfile pos 1 >> ${stem}.cell

echo "%ENDBLOCK POSITIONS_ABS" >> ${stem}.cell
echo >> ${stem}.cell

# Write castep input files from templates
cat ${template}.cell >> ${stem}.cell
cat ${template}.param > ${stem}.param

# If not the first time, then start from old wavefunctions
if [[ -r ${stem}.check && $use_check_files == 1 ]]
then
    echo "reuse: default" >> ${stem}.param
fi

for loop in 1 2 ; do # Loop at most twice: once reusing check file and once without

    # Invoke castep
    if [[ $test_mode == 1 ]]; then
	echo castep_driver: test mode
    else
	rm -f ${stem}.castep ${stem}.out ${stem}.0001.err ${olddir}/${stem}.out
	${castep} ${stem}
    fi

    n_failed=0
    while ! grep -q "Total time" ${stem}.castep; do
      echo castep_driver ${stem}: Error encountered while running CASTEP
      ((n_failed=$n_failed+1))
      if ((n_failed > 2)); then
	  echo castep_driver ${stem}: Too many failed restarts, aborting.
	  exit 1
      fi
    
      # If we were trying to reuse checkfile, try without this time
      cp ${stem}.param ${stem}.param.old
      awk '! /^reuse/' ${stem}.param.old > ${stem}.param
      rm ${stem}.param.old
    
      if [[ $test_mode == 1 ]]; then
	  echo castep_driver: test mode
      else
	  ${castep} ${stem}
      fi
    done

    # Extract information from .castep file
    ctime=`grep 'Total time' ${stem}.castep | awk '{print $4}'`

    if grep -i "task" ${stem}.param | grep -q -i 'geometry' ; then
	energy=`grep "Final Enthalpy" ${stem}.castep | awk '{print $5}'`
    elif ! grep "fix_occupancy" ${stem}.param | grep -q 'true'; then
	energy=`grep "Final energy" ${stem}.castep | awk '{print $4}'`
    elif ! grep "finite_basis_corr" ${stem}.param | grep -q 'none'; then
	energy=`grep "Total energy corrected for finite basis set" ${stem}.castep | 
                awk '{print $9}'`
    else
	energy=`grep "Final energy" ${stem}.castep | awk '{print $4}'`
    fi

    # First line of output is number of atoms
    echo $N > ${stem}.out

    # Second line of output is parameter line
    echo Lattice\=\"$lattice\" Properties=\"species:S:1:pos:R:3:force:R:3\" energy\=$energy >> ${stem}.out

    # Extract forces
    grep -A $(($N+5)) "***** Forces *****" ${stem}.castep | tail -$N | \
	awk '{print $2,$3,$4,$5,$6}' > ${stem}_tmpforce

    # Combine atomic positions and forces, converting from castep 
    # (species,species number) ordering to atom number ordering as we go
    
    if [[ $SEQ == "seq" ]] ; then
	numlist=`seq 1 $N`
    elif [[ $SEQ == "jot" ]] ; then
	numlist=`jot $N 1` 
    else
	echo "cannot interpret SEQ variable"
	exit 1
    fi

    if grep -i "task" ${stem}.param | grep -q -i 'geometry' ; then
	tail -$(($N*2+1)) ${stem}.geom | head -$N | awk -v factor=1 '{print $1, $3*factor,$4*factor,$5*factor}' > ${stem}_tmppos
    else
	print_property $xyzfile pos 1 > ${stem}_tmppos
    fi
    paste ${stem}_tmppos <(for i in $numlist
    do
       # Find species and species count of ith atom in XYZ file
       at=`awk 'NR > 2 { nz[$1] += 1 } NR-2 == '$i' {print $1, nz[$1]}' $xyzfile`
       # Look this up in CASTEP force list
       grep "^$at " ${stem}_tmpforce | awk '{printf "%16.8f%16.8f%16.8f\n",$3,$4,$5}'
    done) >> ${stem}.out

    

    # Save all castep output
    cat ${stem}.castep >> ${stem}_castep_output
    #rm ${stem}_tmpforce ${stem}_tmppos

    if [[ `wc -l ${stem}.out | awk '{print $1}'` != $(($N+2)) ]]
	then
	echo castep_driver ${stem}: Error parsing CASTEP output file 
	exit 1
    fi

    # Check there are no crazily large forces
    max_force=`print_property ${stem}.out force 0 | awk '{ print sqrt($1*$1 + $2*$2 + $3*$3) }' | sort -g -r | head -1`
    echo castep_driver ${stem}: max force is $max_force after iteration $loop
    if (( `echo "$max_force > $max_force_tol" | bc -l` == 1)); then
        # Repeat without using checkfile
	echo "castep_driver ${stem}: max force too large: repeating without check file"
	cp ${stem}.param ${stem}.param.old
	awk '! /^reuse/' ${stem}.param.old > ${stem}.param
	rm ${stem}.param.old
    else
	break
    fi
    
done



# Copy output file
cp ${stem}.out $olddir/$outfile

echo castep_driver ${stem}: done, E\=`printf %.1f $energy` eV, time $ctime s
exit 0
