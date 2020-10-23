#!/bin/bash

# set -e

if [ -z $QUIP_ROOT ]; then
   echo "$0: Need QUIP_ROOT defined"
   exit 1
fi
if [ -z $QUIP_ARCH ]; then
   echo "$0: Need QUIP_ARCH defined"
   exit 1
fi

TEST=`basename $0`

mydir=`dirname $0`
bindir=$mydir/../build/$QUIP_ARCH

if [ ! -x $bindir/convert ]; then
   (cd $QUIP_ROOT && make Programs) || exit 2
fi

cat<<EOF > ${TEST}.in_ok.xyz
1
Lattice="1 0 0 0 1 0 0 0 1" Properties=species:S:1:pos:R:3 scalar=1.0 vector3="1 1 1" matrix33="1 2 3 2 3 4 3 4 5" stress="1 1 1 0.1 0.1 0.1"
Si      0.1000000      0.0000000      0.0000000
EOF

error=0

echo -n "$0: "
$bindir/convert ${TEST}.in_ok.xyz ${TEST}.out.xyz
error=$(( $error + $? ))
echo "CUMULATIVE ERROR $error"

cat<<EOF > ${TEST}.in_bad.xyz
1
Lattice="1 0 0 0 1 0 0 0 1" Properties=species:S:1:pos:R:3 scalar=1.0 vector3="1 1 1" matrix33="1 2 3 2 3 4 3 4 5" stress="1 1 1 0.1 0.1 0.1" vector2="1.0 2.0"
Si      0.1000000      0.0000000      0.0000000
EOF

echo -n "$0: "
$bindir/convert ${TEST}.in_bad.xyz ${TEST}.out.xyz |& fgrep -q 'SYSTEM ABORT'
error=$(( $error + $? ))
echo "CUMULATIVE ERROR $error"

rm -f ${TEST}.*
exit $error
