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
Si 0.0706205374 0.0210087101 -0.0371005597
Si 2.6954892794 2.6530645599 0.0193579010
Si 2.7155705582 -0.0348569133 2.7159314325
Si 0.0127801891 2.7373867303 2.6666624999
Si 1.3522654917 1.3066580087 1.3334442519
Si 4.1515556249 4.0334823084 1.3832272757
Si 4.0857131657 1.4020548310 4.0810633005
Si 1.3780139693 4.0715580333 4.0504283863
EOF2

cat<<EOF3 > ${TEST}.out_ref
STAT        0.00    299.9967    299.9967    0.26E-14      0.00000000E+00      0.00000000E+00      0.31022046E+00     -0.34192505E+02     -0.34192505E+02
STAT        1.00    288.3662    299.8810    0.28E+01      0.00000000E+00      0.00000000E+00      0.29819357E+00     -0.34209874E+02     -0.34209874E+02
STAT        2.00    295.1115    299.8335    0.16E+02      0.00000000E+00      0.00000000E+00      0.30516876E+00     -0.34197785E+02     -0.34197785E+02
STAT        3.00    353.6617    300.3691    0.25E+02      0.00000000E+00      0.00000000E+00      0.36571438E+00     -0.34131664E+02     -0.34131664E+02
STAT        4.00    433.0686    301.6895    0.25E+02      0.00000000E+00      0.00000000E+00      0.44782743E+00     -0.34041478E+02     -0.34041478E+02
STAT        5.00    453.1283    303.1963    0.31E+02      0.00000000E+00      0.00000000E+00      0.46857074E+00     -0.34009800E+02     -0.34009800E+02
STAT        6.00    427.0608    304.4288    0.39E+02      0.00000000E+00      0.00000000E+00      0.44161488E+00     -0.34025908E+02     -0.34025908E+02
STAT        7.00    437.4305    305.7522    0.48E+02      0.00000000E+00      0.00000000E+00      0.45233800E+00     -0.34002844E+02     -0.34002844E+02
STAT        8.00    443.6188    307.1240    0.56E+02      0.00000000E+00      0.00000000E+00      0.45873713E+00     -0.33980635E+02     -0.33980635E+02
STAT        9.00    462.4964    308.6700    0.52E+02      0.00000000E+00      0.00000000E+00      0.47825809E+00     -0.33949822E+02     -0.33949822E+02
STAT       10.00    501.9507    310.5931    0.48E+02      0.00000000E+00      0.00000000E+00      0.51905701E+00     -0.33901341E+02     -0.33901341E+02
STAT       10.00    501.9507    310.5931    0.48E+02      0.00000000E+00      0.00000000E+00      0.51905701E+00     -0.33901341E+02     -0.33901341E+02
EOF3

# originally the test used rng_seed=1
# after the change in commit 8facdc861fa3e26f87ccb38bad4e64267bd9a8a3 this needed changing
# rng_seed=2065775975 gives 1 after 100 iterations
${MPIRUN} $bindir/md pot_init_args='{IP SW}' params_in_file=$QUIP_ROOT/share/Parameters/ip.parms.SW.xml trajectory_out_file=${TEST}.traj.xyz T=1500.0 dt=1.0 N_steps=10 rng_seed=2065775975 < ${TEST}.in.xyz > ${TEST}.out.raw
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
