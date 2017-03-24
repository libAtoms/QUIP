#!/bin/bash

set -e

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

TEST=test_md_velo.sh

cat<<EOF > ${TEST}.in.xyz
8
Lattice="5.428835          0.000000          0.000000          0.000000          5.428835          0.000000          0.000000          0.000000          5.428835" Properties=species:S:1:pos:R:3:velo:R:3
  Si      0.1000000      0.0000000      0.0000000       -0.00273261   0.00248437  -0.00193483
  Si      2.7144176      2.6144176      0.0000000       -0.00114229   0.00248983  -0.00046045
  Si      2.7144176      0.0000000      2.7144176       -0.00136687  -0.00313459  -0.00169086
  Si      0.0000000      2.7144176      2.7144176       -0.00267027   0.00378309  -0.00187122
  Si      1.3572088      1.3572088      1.3572088       -0.00012507  -0.00569576  -0.00238343
  Si      4.0716264      4.0716264      1.3572088        0.00623729  -0.00319620   0.00642180
  Si      4.0716264      1.3572088      4.0716264       -0.00087699   0.00292093   0.00186215
  Si      1.3572088      4.0716264      4.0716264        0.00267681   0.00034833   0.00005684
EOF

cat<<EOF2 > ${TEST}.traj_ref.xyz
Si 0.0706208553 0.0210087824 -0.0371004557
Si 2.6954892778 2.6530641578 0.0193577421
Si 2.7155704093 -0.0348569067 2.7159313280
Si 0.0127799060 2.7373868572 2.6666626085
Si 1.3522654873 1.3066579739 1.3334441315
Si 4.1515555815 4.0334823152 1.3832276592
Si 4.0857129422 1.4020548812 4.0810633934
Si 1.3780139134 4.0715581847 4.0504285120
EOF2

cat<<EOF3 > ${TEST}.out_ref
STAT        0.00    300.0002    300.0002    0.27E-14      0.00000000E+00      0.00000000E+00      0.31022607E+00     -0.34192499E+02     -0.34192499E+02
STAT        1.00    288.3695    299.8845    0.28E+01      0.00000000E+00      0.00000000E+00      0.29819894E+00     -0.34209869E+02     -0.34209869E+02
STAT        2.00    295.1148    299.8370    0.16E+02      0.00000000E+00      0.00000000E+00      0.30517411E+00     -0.34197779E+02     -0.34197779E+02
STAT        3.00    353.6651    300.3726    0.25E+02      0.00000000E+00      0.00000000E+00      0.36572018E+00     -0.34131658E+02     -0.34131658E+02
STAT        4.00    433.0721    301.6930    0.25E+02      0.00000000E+00      0.00000000E+00      0.44783385E+00     -0.34041472E+02     -0.34041472E+02
STAT        5.00    453.1315    303.1998    0.31E+02      0.00000000E+00      0.00000000E+00      0.46857703E+00     -0.34009793E+02     -0.34009793E+02
STAT        6.00    427.0636    304.4323    0.39E+02      0.00000000E+00      0.00000000E+00      0.44162058E+00     -0.34025902E+02     -0.34025902E+02
STAT        7.00    437.4332    305.7557    0.48E+02      0.00000000E+00      0.00000000E+00      0.45234366E+00     -0.34002837E+02     -0.34002837E+02
STAT        8.00    443.6212    307.1275    0.56E+02      0.00000000E+00      0.00000000E+00      0.45874252E+00     -0.33980628E+02     -0.33980628E+02
STAT        9.00    462.4982    308.6734    0.52E+02      0.00000000E+00      0.00000000E+00      0.47826304E+00     -0.33949816E+02     -0.33949816E+02
STAT       10.00    501.9522    310.5966    0.48E+02      0.00000000E+00      0.00000000E+00      0.51906185E+00     -0.33901334E+02     -0.33901334E+02
STAT       10.00    501.9522    310.5966    0.48E+02      0.00000000E+00      0.00000000E+00      0.51906185E+00     -0.33901334E+02     -0.33901334E+02
EOF3

# originally the test used rng_seed=1
# after the change in commit 8facdc861fa3e26f87ccb38bad4e64267bd9a8a3 this needed changing
# rng_seed=2065775975 gives 1 after 100 iterations
${MPIRUN} $bindir/md verbosity=ANAL pot_init_args='{IP SW}' param_filename=${QUIP_ROOT}/share/Parameters/ip.parms.SW.xml trajectory_filename=${TEST}.traj.xyz T=1500.0 dt=1.0 N_steps=10 rng_seed=2065775975 < ${TEST}.in.xyz > ${TEST}.out.raw
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
