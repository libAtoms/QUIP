#!/bin/bash

function git_date
{
   if [ -e $1 ]
   then
      DATE=$(git log -n1 --format="%at" $1)
      [[ $? -eq 0 && ! -z $DATE ]] && echo $DATE || echo 0
   else
      echo 0
   fi
}

QUIP_ROOT=$(dirname $0)/..
GAP_FILES="src/Potentials/IPModel_GAP.f95 \
   src/GAP/descriptors.f95 src/GAP/descriptors_wrapper.f95 \
   src/GAP/make_permutations_v2.f95 src/GAP/gp_predict.f95 \
   src/GAP-filler/clustering.f95 src/GAP-filler/gp_teach.f95 \
   src/GAP-filler/teach_sparse_module.f95 src/GAP-filler/teach_sparse.f95"

if [ -s "${QUIP_ROOT}/src/GAP/GAP_VERSION" ]; then
   echo -ne $(cat ${QUIP_ROOT}/src/GAP/GAP_VERSION)
   exit 0
elif [ -d ${QUIP_ROOT}/.git ]; then
   # first output latest file to stderr
   for I in $GAP_FILES
   do
      if [ -d ${QUIP_ROOT}/$(dirname $I) ]
      then
         cd ${QUIP_ROOT}/$(dirname $I)
         echo -n "$I "
         git_date $(basename $I)
         cd - >/dev/null
      fi
   done | sort -k2n | tail -1 | awk '{print $1}' 1>&2
   # now output time of latest file to stdout
   for I in $GAP_FILES
   do
      if [ -d ${QUIP_ROOT}/$(dirname $I) ]
      then
         cd ${QUIP_ROOT}/$(dirname $I)
         git_date $(basename $I)
         cd - >/dev/null
      fi
   done | sort -n | tail -1
else
   echo "0"
   exit 0
fi
