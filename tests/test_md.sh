#!/bin/bash

test -e

if [ -z $QUIP_ROOT ]; then
   echo "$0: Need QUIP_ROOT defined"
   exit 1
fi
if [ -z $QUIP_ARCH ]; then
   echo "$0: Need QUIP_ARCH defined"
   exit 1
fi

mydir=`dirname $0`
bindir=$mydir/../build/$QUIP_ARCH

if [ ! -x $bindir/md ]; then
   (cd $QUIP_ROOT && make Programs) || exit 2
fi

TEST=test_md.sh

cat<<EOF > ${TEST}.in.xyz
8
Lattice="5.428835          0.000000          0.000000          0.000000          5.428835          0.000000          0.000000          0.000000          5.428835" Properties=species:S:1:pos:R:3
  Si      0.1000000      0.0000000      0.0000000
  Si      2.7144176      2.6144176      0.0000000
  Si      2.7144176      0.0000000      2.7144176
  Si      0.0000000      2.7144176      2.7144176
  Si      1.3572088      1.3572088      1.3572088
  Si      4.0716264      4.0716264      1.3572088
  Si      4.0716264      1.3572088      4.0716264
  Si      1.3572088      4.0716264      4.0716264
EOF

cat<<EOF2 > ${TEST}.traj_ref.xyz
Si 0.0950477455 0.0015295503 -0.0189412241
Si 2.7048749528 2.6311284359 0.0199746978
Si 2.7248335507 -0.0074628639 2.7298935304
Si 0.0306006379 2.7089328527 2.6801615812
Si 1.3545184816 1.3545294944 1.3545529375
Si 4.0990580111 4.0589759917 1.3301844425
Si 4.0959055785 1.3752201852 4.0672756875
Si 1.3571698576 4.0675026220 4.0499128352
EOF2

cat<<EOF3 > ${TEST}.out_ref
STAT        0.00      0.0000      0.0000    0.00E+00      0.00000000E+00      0.00000000E+00      0.00000000E+00     -0.34502725E+02     -0.34502725E+02
STAT        1.00      3.8651      0.0385    0.28E+01      0.00000000E+00      0.00000000E+00      0.39968204E-02     -0.34500324E+02     -0.34500324E+02
STAT        2.00     23.1251      0.2682    0.16E+02      0.00000000E+00      0.00000000E+00      0.23913155E-01     -0.34475903E+02     -0.34475903E+02
STAT        3.00     60.9994      0.8725    0.25E+02      0.00000000E+00      0.00000000E+00      0.63078286E-01     -0.34436119E+02     -0.34436119E+02
STAT        4.00    111.7440      1.9757    0.25E+02      0.00000000E+00      0.00000000E+00      0.11555223E+00     -0.34385532E+02     -0.34385532E+02
STAT        5.00    161.2551      3.5605    0.31E+02      0.00000000E+00      0.00000000E+00      0.16675065E+00     -0.34339609E+02     -0.34339609E+02
STAT        6.00    187.6323      5.3921    0.39E+02      0.00000000E+00      0.00000000E+00      0.19402676E+00     -0.34318264E+02     -0.34318264E+02
STAT        7.00    221.2641      7.5400    0.48E+02      0.00000000E+00      0.00000000E+00      0.22880466E+00     -0.34288324E+02     -0.34288324E+02
STAT        8.00    259.0816     10.0429    0.56E+02      0.00000000E+00      0.00000000E+00      0.26791102E+00     -0.34251151E+02     -0.34251151E+02
STAT        9.00    291.0613     12.8391    0.52E+02      0.00000000E+00      0.00000000E+00      0.30098061E+00     -0.34217136E+02     -0.34217136E+02
STAT       10.00    320.1778     15.8971    0.48E+02      0.00000000E+00      0.00000000E+00      0.33108938E+00     -0.34190243E+02     -0.34190243E+02
STAT       10.00    320.1778     15.8971    0.48E+02      0.00000000E+00      0.00000000E+00      0.33108938E+00     -0.34190243E+02     -0.34190243E+02
EOF3
# originally the test used rng_seed=1
# after the change in commit 8facdc861fa3e26f87ccb38bad4e64267bd9a8a3 this needed changing
# rng_seed=2065775975 gives 1 after 100 iterations
${MPIRUN} $bindir/md pot_init_args='{IP SW}' params_in_file=$QUIP_ROOT/share/Parameters/ip.parms.SW.xml trajectory_out_file=${TEST}.traj.xyz T=1500.0 dt=1.0 N_steps=10 rng_seed=2065775975 atoms_in_file=${TEST}.in.xyz > ${TEST}.out.raw
egrep '^ST' ${TEST}.out.raw > ${TEST}.out
tail -8 ${TEST}.traj.xyz | awk '{print $1" "$2" "$3" "$4}' > ${TEST}.traj_pos.xyz

diff=diff
which ndiff > /dev/null && diff=ndiff

error=0
echo -n "$0: "
nr=`wc -l < ${TEST}.traj_ref.xyz`
nt=`wc -l < ${TEST}.traj_pos.xyz`
if [[ $nr != $nt ]]; then
   echo -n "final pos line number mismatch $nr $nt "
else
   if $diff ${TEST}.traj_ref.xyz ${TEST}.traj_pos.xyz  > /dev/null 2>&1; then
      echo -n "final pos OK "
   else
      echo -n "final pos differs "
      error=1
   fi
fi

nr=`wc -l < ${TEST}.out`
nt=`wc -l < ${TEST}.out_ref`
if [[ $nr != $nt ]]; then
   echo -n "output line number mismatch $nr $nt"
else
   if $diff ${TEST}.out_ref ${TEST}.out > /dev/null 2>&1; then
      echo -n "E/T OK"
   else
      echo -n "E/T differs"
      error=1
   fi
fi

echo ""

[[ $error == 0 ]] && rm -f ${TEST}.*
exit $error
