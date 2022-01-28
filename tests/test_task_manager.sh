#!/bin/bash

set -e

if [ -z $QUIP_ARCH ]; then
   echo "$0: Need QUIP_ARCH defined"
   exit 1
fi

mydir=`dirname $0`
bindir=$mydir/../build/$QUIP_ARCH
QUIP_ROOT=${mydir}/..

if [ ! -x $bindir/test_task_manager ]; then
    (cd $QUIP_ROOT && make Programs) || exit 2
fi

$bindir/test_task_manager &> test_task_manager.log
grep "test_task_manager OK" test_task_manager.log &> /dev/null
error=$?

rm -f test_task_manager.log
exit $error
