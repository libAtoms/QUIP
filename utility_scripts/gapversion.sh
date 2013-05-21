#!/bin/bash

function git_date
{
   if [ -e $1 ]
   then
      DATE=$(git log -n1 --format="%at" $1)
      [ $? -eq 0 ] && echo $DATE || echo 0
   else
      echo 0
   fi
}

QUIP_ROOT=$(dirname $0)/..

GAP_FILES="QUIP_Core/IPModel_GAP.f95 GAP/gp_predict.f95 GAP-filler/clustering.f95 GAP-filler/gp_teach.f95 GAP-filler/teach_sparse_module.f95 GAP-filler/teach_sparse.f95"

for I in $GAP_FILES
do
   if [ -d ${QUIP_ROOT}/$(dirname $I) ]
   then
      cd ${QUIP_ROOT}/$(dirname $I)
      git_date $(basename $I)
      cd - >/dev/null
   fi
done | sort -n | tail -1
