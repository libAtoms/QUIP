#!/bin/bash
#
# Determine the gap version from the git reported last changed date
# of file (or the file GAP_VERSION if it exists). The gapversion
# is the UNIX timestap of the most recently changed file.

GAP_FILES="src/Potentials/IPModel_GAP.f95 \
  src/GAP/descriptors.f95 src/GAP/descriptors_wrapper.f95 \
  src/GAP/make_permutations_v2.f95 src/GAP/gp_predict.f95 \
  src/GAP-filler/clustering.f95 src/GAP-filler/gp_teach.f95 \
  src/GAP-filler/teach_sparse_module.f95 src/GAP-filler/teach_sparse.f95"

QUIP_ROOT=$(dirname "$0")/..

# By default only show the date
VERBOSE=false

# Options for interactive use
while getopts 'v' opt; do
  case $opt in
    v)
      VERBOSE=true
      ;;
    \?)
      echo "Usage: $0"
      echo "'-v' to output name of newest file to stderr"
      exit 1
      ;;
  esac
done

function git_date
{
  if [ -e "$1" ]
  then
    DATE=$(git log -n1 --format="%at" "$1")
    [[ $? -eq 0 && ! -z $DATE ]] && echo "$DATE" || echo 0
  else
    echo 0
  fi
}

if [ -s "${QUIP_ROOT}/src/GAP/GAP_VERSION" ]; then
  echo -ne "$(cat "${QUIP_ROOT}/src/GAP/GAP_VERSION")"
  exit 0
elif [ -d "${QUIP_ROOT}/.git" ]; then
  # Sort everything as "DATE filename" list and take the last one
  GAP_VERSION=$(
    for I in $GAP_FILES
    do
      if [ -d "${QUIP_ROOT}/$(dirname "$I")" ]
      then
        cd "${QUIP_ROOT}/$(dirname "$I")"
        echo $(git_date $(basename "$I")) "$I"
        cd - >/dev/null
      fi
    done | sort -n | tail -1 )
  # filename to stderr if requested;
  # split string with parameter expansion
  if [ $VERBOSE == "true" ]; then
    echo "${GAP_VERSION##* }" >&2
  fi
  echo "${GAP_VERSION%% *}"
else
   echo "0"
   exit 0
fi
