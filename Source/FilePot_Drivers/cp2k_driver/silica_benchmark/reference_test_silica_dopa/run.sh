#!/bin/sh

#try to run CP2K
cp2k_program=$1
$1 cp2k.inp > cp2k.out

#succesful run?
if test $? -ne 0
then
   echo "CP2K aborted. Test failed."
   exit
fi

#do the forces and the energy match?
difference=`diff nvt-frc-1.xyz nvt-frc-1.xyz_ref | wc -l`
if test $difference -eq 0
then
   echo "Forces and energy are the same as in the reference file.  Test passed."
else
   echo "Force files are different.  Test failed."
   echo "See diff nvt-frc-1.xyz nvt-frc-1.xyz_ref"
fi
